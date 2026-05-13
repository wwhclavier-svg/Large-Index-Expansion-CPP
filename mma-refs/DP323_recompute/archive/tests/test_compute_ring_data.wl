(* test_compute_ring_data.wl — Verify ComputeRingData output consistency *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
Get["ExportBinary_IBPMatrix.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

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
ibpeqs = data["Family", "IBPEqs"];
vlist = data["Family", "VList"];
existingRegions = data["Regions"];

(* Pick the same sector as test_singular_compare *)
targetSector = SelectFirst[SortBy[Keys[existingRegions], Total],
  Length[existingRegions[#]] == 1 &,
  $Failed
];
Print["Target sector: ", targetSector];

ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

(* Run MMA *)
resMMA = LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];
Print["MMA regions: ", Length[resMMA]];

(* Run Singular *)
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
Print["Singular regions: ", Length[resSingular]];

(* Compare ComputeRingData for each corresponding region *)
allPass = True;
Do[
  If[i > Length[resSingular],
    Print["Region ", i, ": Singular missing"];
    allPass = False;
    Continue[];
  ];
  cdMMA = ExportIBPMatrixBinary`ComputeRingData[
    <|"CoordinateRing" -> resMMA[[i]]|>, Alist, ne, char];
  cdS = ExportIBPMatrixBinary`ComputeRingData[
    <|"CoordinateRing" -> resSingular[[i]]|>, Alist, ne, char];

  Print["\nRegion ", i, ":"];
  If[cdMMA["LimitSector"] =!= cdS["LimitSector"],
    Print["  FAIL LimitSector"]; allPass = False,
    Print["  PASS LimitSector"]];
  If[cdMMA["MatDim"] =!= cdS["MatDim"],
    Print["  FAIL MatDim: MMA=", cdMMA["MatDim"], ", Singular=", cdS["MatDim"]]; allPass = False,
    Print["  PASS MatDim=", cdMMA["MatDim"]]];
  If[cdMMA["FlatA"] =!= cdS["FlatA"],
    Print["  FAIL FlatA"]; allPass = False,
    Print["  PASS FlatA"]];
  If[cdMMA["FlatAinv"] =!= cdS["FlatAinv"],
    Print["  FAIL FlatAinv"]; allPass = False,
    Print["  PASS FlatAinv"]];
, {i, Length[resMMA]}];

Print["\n========================================"];
If[allPass,
  Print["ComputeRingData CHECKS PASSED"],
  Print["ComputeRingData CHECKS FAILED"]];
Print["========================================"];
