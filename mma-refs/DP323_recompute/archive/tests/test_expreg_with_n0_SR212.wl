(* test_expreg_with_n0_SR212.wl — Test expRegSolve2 with n->0 like regionsBySectors does *)
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

(* Replicate regionsBySectors input preparation EXACTLY *)
ibpeqs_n0 = ibpeqs /. n -> 0;  (* Note: regionsBySectors does this due to string bug *)
ibpeqssub = LIEUtility`sectorLimitIBP[ibpeqs_n0, targetSector, vlist];

Print["Running expRegSolve2 with n->0 preprocessing..."];
result = LIECoreAlgebra`expRegSolve2[ibpeqssub, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];

Print["Result: ", Length[result], " regions"];
Do[
  Print["Region ", i, ": VarDep=", result[[i]]["VarDep"]];
, {i, Length[result]}];

(* Also test WITHOUT n->0 *)
Print["\nRunning expRegSolve2 WITHOUT n->0 preprocessing..."];
ibpeqssub2 = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];
result2 = LIECoreAlgebra`expRegSolve2[ibpeqssub2, Alist, vlist,
  "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];

Print["Result: ", Length[result2], " regions"];
Do[
  Print["Region ", i, ": VarDep=", result2[[i]]["VarDep"]];
, {i, Length[result2]}];
