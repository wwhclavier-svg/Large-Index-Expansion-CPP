(* run_scheduled.wl -- Run ScheduledRegionSolve for NP322 with SingularCoordinateRing acceleration *)
(* Usage: wolframscript -file workspace/C002-NP322-region-solve/run_scheduled.wl *)

$ProjectRoot = DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = FileNameJoin[{$ProjectRoot, "verify", "VerifyUtility"}];
$FamilyDBPath = FileNameJoin[{$ProjectRoot, "verify", "FamilyDatabase", "FamilyDatabase.wl"}];

SetDirectory[$VerifyUtilityPath];

Print["========================================"];
Print["ScheduledRegionSolve NP322 Recompute"];
Print["ProjectRoot: ", $ProjectRoot];
Print["VerifyUtility: ", $VerifyUtilityPath];
Print["========================================"];

(* Load packages *)
Print["\n[1/5] Loading packages..."];
Get["LIEWorkflow.wl"];
Get["ScheduledRegionSolver.wl"];
Get["ExportBinary_IBPMatrix.wl"];
Print["  Packages loaded."];

(* Load family config *)
Print["\n[2/5] Loading family config..."];
Get[$FamilyDBPath];
If[!KeyExistsQ[$FamilyDatabase, "NP322"],
  Print["ERROR: NP322 not found in FamilyDatabase"];
  Exit[1]
];
config = $FamilyDatabase["NP322"];
Print["  Name: ", config["Description"]];
Print["  Propagators: ", Length[config["Propagators"]]];
Print["  Modulus: ", config["Modulus"]];

(* Define family *)
Print["\n[3/5] Defining IBP family..."];
timeDefine = AbsoluteTiming[
  data = LIEDefineFamily[
    config["Propagators"], config["LoopMomenta"],
    config["ExternalMomenta"], config["KinematicRules"],
    config["TopSector"],
    "Numeric" -> config["Numeric"],
    Modulus -> config["Modulus"],
    Verbose -> True
  ];
][[1]];
Print["  Family definition time: ", timeDefine, " s"];
Print["  NE = ", data["Family", "NE"]];
Print["  NIBP = ", data["Family", "NIBP"]];
Print["  Sectors = ", Length[data["Family", "SectorList"]]];

(* Package for ScheduledRegionSolve *)
familyConfig = <|
  "Family" -> data["Family"],
  "Config" -> <|"Modulus" -> config["Modulus"]|>
|>;

(* Run ScheduledRegionSolve *)
Print["\n[4/5] Running ScheduledRegionSolve..."];
scheduleConfig = <|
  "Tiers" -> {
    <|"Name" -> "Quick",  "Timeout" -> 120,   "Order" -> "First", "MaxParallel" -> 4|>,
    <|"Name" -> "Fast",   "Timeout" -> 1200,  "Order" -> "First", "MaxParallel" -> 3|>,
    <|"Name" -> "Retry",  "Timeout" -> 12000, "Order" -> "FIFO",  "MaxParallel" -> 3|>,
    <|"Name" -> "Long",   "Timeout" -> 96400, "Order" -> "FIFO",  "MaxParallel" -> 1|>
  },
  "MaxRetriesPerTier" -> 1,
  "CacheDir" -> FileNameJoin[{$ProjectRoot, "workspace", "C002-NP322-region-solve", "cache"}],
  "OutputDir" -> $ProjectRoot,  (* Output .bin to project root for C++ tests *)
  "LogDir" -> FileNameJoin[{$ProjectRoot, "workspace", "C002-NP322-region-solve", "cache"}],
  "ScheduleMode" -> "Interleaved",
  "MaxParallelJobs" -> 4,
  "CheckpointInterval" -> 1,
  "ParallelJobs" -> 1,
  "Verbose" -> True
|>;

result = ScheduledRegionSolve["NP322", familyConfig, scheduleConfig];

(* Report results *)
Print["\n[5/5] Result summary"];
If[result === $Failed,
  Print["  FAILED: ScheduledRegionSolve returned $Failed"];
  Exit[1]
];

status = result["Status"];
failedSectors = result["FailedSectors"];

doneCount = Length[Select[status, #["State"] === "done" &]];
trivialCount = Length[Select[status, #["State"] === "trivial" &]];
timeoutCount = Length[Select[status, #["State"] === "timeout" &]];
errorCount = Length[Select[status, #["State"] === "error" &]];

Print["  Done:     ", doneCount];
Print["  Trivial:  ", trivialCount];
Print["  Timeout:  ", timeoutCount];
Print["  Error:    ", errorCount];
Print["  Failed sectors remaining: ", Length[failedSectors]];
If[Length[failedSectors] > 0,
  Print["  Failed: ", failedSectors];
];

(* Verify output files exist *)
ibpFile = FileNameJoin[{$ProjectRoot, "IBPMat_NP322.bin"}];
ringFile = FileNameJoin[{$ProjectRoot, "RingData_NP322.bin"}];
Print["\nOutput files:"];
Print["  ", ibpFile, " ", If[FileExistsQ[ibpFile], "EXISTS", "MISSING"]];
Print["  ", ringFile, " ", If[FileExistsQ[ringFile], "EXISTS", "MISSING"]];
If[FileExistsQ[ibpFile],
  Print["  IBPMat size: ", FileByteCount[ibpFile], " bytes"];
];
If[FileExistsQ[ringFile],
  Print["  RingData size: ", FileByteCount[ringFile], " bytes"];
];

Print["\n========================================"];
Print["ScheduledRegionSolve NP322 complete."];
Print["========================================"];
