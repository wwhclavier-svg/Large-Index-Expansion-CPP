(* test_buildrec_SR212.wl — Test buildRecursionMatrix on 3 regions *)
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

targetSector = {0, 1, 0, 1, 1};
ne = Length[Alist];

ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

Print["\nRunning expRegSolve2..."];
regs = LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];

Print["Got ", Length[regs], " regions"];

Do[
  Print["\nBuilding recursion matrix for region ", i, "..."];
  result = LIERegions`Private`buildRecursionMatrix[ibpeqs, regs[[i]], ne,
    "LimitSector" -> targetSector, Modulus -> char];
  Print["  Result keys: ", Keys[result]];
, {i, Length[regs]}];

Print["\nDone"];
