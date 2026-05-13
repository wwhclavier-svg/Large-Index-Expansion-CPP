(* test_regions_n0_SR212.wl — Test regionsBySectors with n->0 on one sector *)
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

(* Call regionsBySectors with n->0 pre-processing like it does internally *)
ibpeqs_n0 = ibpeqs /. n -> 0;

Print["\nCalling regionsBySectors with n->0..."];
result = LIERegions`regionsBySectors[ibpeqs_n0, {{0, 1, 0, 1, 1}}, Alist, vlist,
  Modulus -> char, "EnableFieldExtension" -> True, Verbose -> False, "Timeout" -> 300];

Print["Result type: ", Head[result]];
Print["Result: ", result];

(* Also call WITHOUT n->0 *)
Print["\nCalling regionsBySectors WITHOUT n->0..."];
result2 = LIERegions`regionsBySectors[ibpeqs, {{0, 1, 0, 1, 1}}, Alist, vlist,
  Modulus -> char, "EnableFieldExtension" -> True, Verbose -> False, "Timeout" -> 300];

Print["Result2 type: ", Head[result2]];
Print["Result2: ", result2];
