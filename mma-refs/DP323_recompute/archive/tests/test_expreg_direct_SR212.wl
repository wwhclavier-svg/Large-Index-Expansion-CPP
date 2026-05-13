(* test_expreg_direct_SR212.wl — Test modified expRegSolve2 directly *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

Print["Loading SR212 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
existingRegions = data["Regions"];

targetSector = SelectFirst[SortBy[Keys[existingRegions], Total],
  Length[existingRegions[#]] == 1 &,
  $Failed
];
If[targetSector === $Failed, Print["ERROR: No suitable sector"]; Exit[1]];

Print["Target sector: ", targetSector];

(* Call expRegSolve2 directly *)
ibpeqssub = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

Print["\nRunning expRegSolve2 with Singular GB..."];
result = LIECoreAlgebra`expRegSolve2[ibpeqssub, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];

Print["Result type: ", Head[result]];
Print["Region count: ", Length[result]];

If[!ListQ[result], Print["ERROR: result is not list"]; Exit[1]];

(* Compare with checkpoint *)
oldRegions = existingRegions[targetSector];
If[Length[result] =!= Length[oldRegions],
  Print["FAIL: Region count mismatch: checkpoint=", Length[oldRegions], ", new=", Length[result]];
  Exit[1];
];

keysToCheck = {"VarDep", "MinPoly", "VarIndep", "VarRule", "FractionRule",
  "MonomialBasis", "MonomialBasisIndex", "MonomialBasisMatrix", "VarDeg"};

allPass = True;
Do[
  regOld = oldRegions[[i]];
  regNew = result[[i]];
  Print["\nRegion ", i, ":"];
  Do[
    key = keysToCheck[[k]];
    valOld = regOld["CoordinateRing", key];
    valNew = regNew[key];
    If[Head[valOld] === List && Head[valNew] === List,
      valOld = Sort[valOld];
      valNew = Sort[valNew];
    ];
    If[Head[valOld] === Association && Head[valNew] === Association,
      valOld = KeySort[valOld];
      valNew = KeySort[valNew];
    ];
    If[valOld =!= valNew,
      Print["  FAIL: ", key];
      Print["    Old: ", Short[valOld, 50]];
      Print["    New: ", Short[valNew, 50]];
      allPass = False;
      ,
      Print["  PASS: ", key];
    ];
  , {k, Length[keysToCheck]}];
, {i, Length[result]}];

Print["\n========================================"];
If[allPass, Print["ALL CHECKS PASSED"], Print["SOME CHECKS FAILED"]];
Print["========================================"];
