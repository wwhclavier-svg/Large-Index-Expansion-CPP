(* test_dp323_perf_singular.wl — Performance test of Singular GB on DP323 hard sector *)
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
fullSectorList = data["Family", "SectorList"];
existingRegions = data["Regions"];

Print["Total sectors: ", Length[fullSectorList]];
Print["Completed sectors: ", Length[Keys[existingRegions]]];

(* Find hard sectors: completed sectors with high prop count *)
hardSectors = Select[Keys[existingRegions], Total[#] >= 4 &];
Print["Hard sectors (>=4 props): ", Length[hardSectors]];

If[Length[hardSectors] == 0,
  Print["No hard sectors found"];
  Exit[1];
];

(* Pick the first hard sector *)
targetSector = hardSectors[[1]];
Print["\nTarget sector: ", targetSector, " (", Total[targetSector], " props)"];

(* Prepare input exactly as regionsBySectors does *)
ibpeqs_n0 = ibpeqs /. n -> 0;
ibpeqssub = LIEUtility`sectorLimitIBP[ibpeqs_n0, targetSector, vlist];

Print["\nRunning expRegSolve2 with Singular GB..."];
totalTime = AbsoluteTiming[
  result = LIECoreAlgebra`expRegSolve2[ibpeqssub, Alist, vlist,
    "LimitSector" -> targetSector, Modulus -> char, Verbose -> True];
][[1]];

Print["\n========================================"];
Print["Total time: ", totalTime, " s"];
Print["Regions found: ", Length[result]];
Print["========================================"];
