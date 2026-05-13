(* test_expreg_singular_SR212.wl — Test modified expRegSolve2 via regionsBySectors *)
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

(* Pick first completed sector with 1 region *)
targetSector = SelectFirst[SortBy[Keys[existingRegions], Total],
  Length[existingRegions[#]] == 1 &,
  $Failed
];
If[targetSector === $Failed, Print["ERROR: No suitable sector"]; Exit[1]];

Print["Target sector: ", targetSector];
Print["Checkpoint regions: ", existingRegions[targetSector]];

(* Call regionsBySectors with modified expRegSolve2 *)
Print["\nRunning regionsBySectors with Singular GB..."];
result = LIERegions`regionsBySectors[ibpeqs, {targetSector}, Alist, vlist,
  Modulus -> char, "EnableFieldExtension" -> True, Verbose -> False, "Timeout" -> 300];

If[!AssociationQ[result], Print["ERROR: result is not association"]; Exit[1]];
If[!KeyExistsQ[result, targetSector], Print["ERROR: sector not in result"]; Exit[1]];

newRegions = result[targetSector];
Print["New regions: ", Length[newRegions]];

(* Compare region count *)
If[Length[newRegions] =!= Length[existingRegions[targetSector]],
  Print["FAIL: Region count mismatch: checkpoint=", Length[existingRegions[targetSector]], ", new=", Length[newRegions]];
  Exit[1];
];

(* Compare key fields *)
keysToCheck = {"VarDep", "MinPoly", "VarIndep", "VarRule", "FractionRule",
  "MonomialBasis", "MonomialBasisIndex", "MonomialBasisMatrix", "VarDeg"};

allPass = True;
Do[
  regOld = existingRegions[targetSector][[i]];
  regNew = newRegions[[i]];
  Print["\nRegion ", i, ":"];
  Do[
    key = keysToCheck[[k]];
    valOld = regOld[key];
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
, {i, Length[newRegions]}];

Print["\n========================================"];
If[allPass, Print["ALL CHECKS PASSED"], Print["SOME CHECKS FAILED"]];
Print["========================================"];
