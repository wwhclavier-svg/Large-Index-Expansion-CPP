(* Debug regionsBySectors *)
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

(* Catch all messages *)
Check[
  result = LIERegions`regionsBySectors[ibpeqs, {{0, 1, 0, 1, 1}}, Alist, vlist, Modulus -> char];
  Print["No errors caught"];
,
  Print["Error caught: ", $MessageList]
];

Print["Result head: ", Head[result]];
Print["Result: ", Short[result, 200]];

(* Try extracting the Table expression *)
If[Head[result] === Association && Length[result] == 1 && Head[result[[1]]] === Table,
  Print["Found unevaluated Table inside Association"];
  tableExpr = result[[1]];
  Print["Table iterator: ", tableExpr[[2]]];
  Print["Evaluating Table directly..."];
  directResult = Evaluate[tableExpr];
  Print["Direct result head: ", Head[directResult]];
  Print["Direct result: ", Short[directResult, 200]];
];
