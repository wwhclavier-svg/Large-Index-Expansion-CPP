(* test_binary_export.wl — Verify RingData binary export compatibility *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
Get["ExportBinary_IBPMatrix.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

(* Load SR212 checkpoint *)
$FamilyPath = $LIECPPPath <> "/verify/SR212/";
Print["Loading SR212 checkpoint..."];
cpFile = $FamilyPath <> "PrepareCheckpoint-SR212.wdx";
Block[{Print = If[#[[1]] === "Warning", Print[#]&, Identity]&},
  data = Import[cpFile, "WDX"];
];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
Alist = data["Family", "AList"];
ne = Length[Alist];

(* Pick a completed sector with consistent region count *)
existingRegions = data["Regions"];
targetSector = SelectFirst[SortBy[Keys[existingRegions], Total],
  True &,
  $Failed
];
Print["Target sector: ", targetSector, " (checkpoint has ", Length[existingRegions[targetSector]], " regions)"];

(* --- Export MMA version --- *)
expregMMA = Association[targetSector -> existingRegions[targetSector]];
ExportIBPMatrixBinary`ExportBinaryRingData[
  $FamilyPath <> "RingData_SR212_MMA.bin",
  expregMMA, Alist, ne, char
];

(* --- Export Singular version --- *)
(* Recompute using Singular replacement *)
ibpeqs = data["Family", "IBPEqs"];
vlist = data["Family", "VList"];
ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

bivarPrimeInfoSym = Symbol["LIECoreAlgebra`Private`bivarPrimeInfo"];
originalDownValues = DownValues[bivarPrimeInfoSym];
DownValues[bivarPrimeInfoSym] = {
  HoldPattern[bivarPrimeInfoSym[primelist_, Avar_, Bvar_, opts:OptionsPattern[]]] :>
    Module[{c = OptionValue[{opts}, Modulus], p = OptionValue[{opts}, "Parameters"],
      ls = OptionValue[{opts}, "LimitSector"], fr = OptionValue[{opts}, "FractionSymbol"]},
      SingularCoordinateRing`singularBivarPrimeInfo[primelist, Avar, Bvar,
        Modulus -> c, "Parameters" -> p, "LimitSector" -> ls, "FractionSymbol" -> fr, Verbose -> False]
    ]
};

resSingular = LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];
DownValues[bivarPrimeInfoSym] = originalDownValues;

(* Filter trivial regions to match regionsBySectors behavior *)
resSingular = Select[resSingular, Join[#["VarDep"], #["VarIndep"]] =!= {} &];
Print["Singular regions after filtering: ", Length[resSingular]];

(* Wrap Singular result in same format as existingRegions (with RecursionMatrix placeholder) *)
resSingularWrapped = Table[
  <|"CoordinateRing" -> resSingular[[i]],
    "RecursionMatrix" -> existingRegions[targetSector][[i]]["RecursionMatrix"],
    "HMatrix" -> existingRegions[targetSector][[i]]["HMatrix"],
    "HNull" -> existingRegions[targetSector][[i]]["HNull"]|>,
  {i, Length[resSingular]}
];
expregSingular = Association[targetSector -> resSingularWrapped];
ExportIBPMatrixBinary`ExportBinaryRingData[
  $FamilyPath <> "RingData_SR212_Singular.bin",
  expregSingular, Alist, ne, char
];

(* --- Compare --- *)
Print["\nComparing binary files..."];
cmd = "cmp -l " <> $FamilyPath <> "RingData_SR212_MMA.bin " <> $FamilyPath <> "RingData_SR212_Singular.bin | wc -l";
resultStr = Import["!" <> cmd, "String"];
numDiff = ToExpression[resultStr];
Print["Different bytes: ", numDiff];
If[numDiff === 0,
  Print["BINARY FILES IDENTICAL — export compatibility verified!"],
  Print["WARNING: Binary files differ by ", numDiff, " bytes"]
];
