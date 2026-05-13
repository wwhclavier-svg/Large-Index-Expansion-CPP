(* test_expreg_debug_SR212.wl — Debug modified expRegSolve2 *)
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

ibpeqssub = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

Print["\nRunning expRegSolve2 with Verbose..."];
result = LIECoreAlgebra`expRegSolve2[ibpeqssub, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> True];

Print["\nResult: ", Length[result], " regions"];
