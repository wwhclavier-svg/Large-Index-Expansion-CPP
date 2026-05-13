(* test_inspect_3regions_SR212.wl — Inspect the 3 regions from expRegSolve2 *)
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

ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

Print["\nRunning expRegSolve2..."];
result = LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];

Print["Number of regions: ", Length[result]];
Do[
  reg = result[[i]];
  Print["\nRegion ", i, ":"];
  Print["  VarDep: ", reg["VarDep"]];
  Print["  VarIndep: ", reg["VarIndep"]];
  Print["  VarRule: ", reg["VarRule"]];
  Print["  MinPoly: ", reg["MinPoly"]];
  Print["  VarDeg: ", reg["VarDeg"]];
  Print["  FractionRule keys: ", Keys[reg["FractionRule"]]];
, {i, Length[result]}];

(* Compare with checkpoint *)
Print["\nCheckpoint region:"];
checkReg = data["Regions"][targetSector][[1, "CoordinateRing"]];
Print["  VarDep: ", checkReg["VarDep"]];
Print["  VarIndep: ", checkReg["VarIndep"]];
Print["  VarRule: ", checkReg["VarRule"]];
Print["  MinPoly: ", checkReg["MinPoly"]];
Print["  VarDeg: ", checkReg["VarDeg"]];
