(* test_dp323_perf_hard.wl — Performance test on DP323 hard sector *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/DP323/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

Print["Loading DP323 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
existingRegions = data["Regions"];

(* Find a hard sector: >= 4 props, completed in checkpoint *)
hardSectors = Select[Keys[existingRegions], Total[#] >= 4 &];
Print["Hard sectors (>=4 props, completed): ", Length[hardSectors]];

If[Length[hardSectors] == 0,
  Print["No suitable hard sector found"];
  Exit[1];
];

(* Pick first hard sector *)
targetSector = hardSectors[[1]];
Print["\nTarget sector: ", targetSector, " (", Total[targetSector], " props)"];
Print["Checkpoint regions: ", Length[existingRegions[targetSector]]];

(* Run regionsBySectors on single sector with timeout *)
Print["\nRunning regionsBySectors with Singular pipeline..."];
totalTime = AbsoluteTiming[
  result = LIERegions`regionsBySectors[ibpeqs, {targetSector}, Alist, vlist,
    Modulus -> char, "EnableFieldExtension" -> True, Verbose -> True, "Timeout" -> 600];
][[1]];

Print["\n========================================"];
Print["Total time: ", totalTime, " s"];
If[AssociationQ[result] && KeyExistsQ[result, targetSector],
  Print["Regions found: ", Length[result[targetSector]]];
  Print["Result: SUCCESS"];
,
  Print["Result: FAILED or TIMEOUT"];
  Print["Result head: ", Head[result]];
];
Print["========================================"];
