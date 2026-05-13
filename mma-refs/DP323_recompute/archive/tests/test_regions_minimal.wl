(* Minimal test of regionsBySectors *)
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

(* Call with minimal arguments *)
Print["Test 1: regionsBySectors with default options"];
result1 = LIERegions`regionsBySectors[ibpeqs, {{0, 1, 0, 1, 1}}, Alist, vlist, Modulus -> char];
Print["Head: ", Head[result1]];
Print["SameQ[result1, unevaluated]? ", SameQ[result1, Hold[LIERegions`regionsBySectors[ibpeqs, {{0, 1, 0, 1, 1}}, Alist, vlist, Modulus -> char]]]];

(* Try with explicit Timeout *)
Print["\nTest 2: with explicit Timeout"];
result2 = LIERegions`regionsBySectors[ibpeqs, {{0, 1, 0, 1, 1}}, Alist, vlist, Modulus -> char, "Timeout" -> 300];
Print["Head: ", Head[result2]];
