(* Better minimal test of regionsBySectors *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
char = data["Config", "Modulus"];

result = LIERegions`regionsBySectors[ibpeqs, {{0, 1, 0, 1, 1}}, Alist, vlist, Modulus -> char];

Print["Head: ", Head[result]];
Print["AssociationQ: ", AssociationQ[result]];
Print["Keys: ", If[AssociationQ[result], Keys[result], "N/A"]];
Print["Normal: ", If[AssociationQ[result], Normal[result], "N/A"]];
