(* test_singular_gb_SR212.wl — Full validation of Singular GB + fixed regionsBySectors *)
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
sectorlist = data["Family", "SectorList"];
existingRegions = data["Regions"];

Print["Running regionsBySectors with Singular GB on all sectors..."];
newResult = LIERegions`regionsBySectors[ibpeqs, sectorlist, Alist, vlist,
  Modulus -> char, "EnableFieldExtension" -> True, Verbose -> False, "Timeout" -> 300];

Print["\nValidating..."];
allPass = True;

Do[
  sec = sectorlist[[i]];
  If[KeyExistsQ[existingRegions, sec],
    oldRegs = existingRegions[sec];
    If[KeyExistsQ[newResult, sec],
      newRegs = newResult[sec];
      If[Length[oldRegs] =!= Length[newRegs],
        Print["FAIL ", sec, ": region count mismatch old=", Length[oldRegs], " new=", Length[newRegs]];
        allPass = False;
        Continue[];
      ];
      Do[
        oldCR = oldRegs[[j, "CoordinateRing"]];
        newCR = newRegs[[j, "CoordinateRing"]];
        keysToCheck = {"VarDep", "MinPoly", "VarIndep", "VarRule", "FractionRule",
          "MonomialBasis", "MonomialBasisIndex", "MonomialBasisMatrix", "VarDeg"};
        Do[
          key = keysToCheck[[k]];
          valOld = oldCR[key];
          valNew = newCR[key];
          If[Head[valOld] === List && Head[valNew] === List,
            valOld = Sort[valOld];
            valNew = Sort[valNew];
          ];
          If[Head[valOld] === Association && Head[valNew] === Association,
            valOld = KeySort[valOld];
            valNew = KeySort[valNew];
          ];
          If[valOld =!= valNew,
            Print["FAIL ", sec, " region ", j, " key ", key];
            allPass = False;
          ];
        , {k, Length[keysToCheck]}];
      , {j, Length[oldRegs]}];
    ,
      Print["FAIL ", sec, ": missing in new result"];
      allPass = False;
    ]
  ]
, {i, Length[sectorlist]}];

Print["\n========================================"];
If[allPass, Print["ALL SECTORS MATCH"], Print["SOME MISMATCHES FOUND"]];
Print["========================================"];
