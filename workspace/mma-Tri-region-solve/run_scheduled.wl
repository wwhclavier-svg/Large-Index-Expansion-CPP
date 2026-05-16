(* run_scheduled.wl -- MMA ScheduledRegionSolve for Tri (3-props, 7 sectors) *)
$ProjectRoot = DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = FileNameJoin[{$ProjectRoot, "workspace", "shared", "VerifyUtility"}];
$FamilyDBPath = FileNameJoin[{$ProjectRoot, "verify", "FamilyDatabase", "FamilyDatabase.wl"}];

SetDirectory[$VerifyUtilityPath];

Print["========================================"];
Print["ScheduledRegionSolve Tri (MMA)"];
Print["ProjectRoot: ", $ProjectRoot];
Print["========================================"];

Print["\n[1/5] Loading packages..."];
Get["LIEWorkflow.wl"];
Get["ScheduledRegionSolver.wl"];
Get["ExportBinary_IBPMatrix.wl"];
Print["  Done."];

Print["\n[2/5] Loading family config..."];
Get[$FamilyDBPath];
config = $FamilyDatabase["Tri"];
Print["  Name: ", config["Description"]];
Print["  Propagators: ", Length[config["Propagators"]]];
Print["  TopSector: ", config["TopSector"]];

Print["\n[3/5] Defining IBP family..."];
timeDefine = AbsoluteTiming[
  data = LIEDefineFamily[
    config["Propagators"], config["LoopMomenta"],
    config["ExternalMomenta"], config["KinematicRules"],
    config["TopSector"],
    "Numeric" -> config["Numeric"],
    Modulus -> config["Modulus"],
    Verbose -> False
  ];
][[1]];
Print["  Time: ", timeDefine, " s"];
Print["  NE = ", data["Family", "NE"]];
Print["  NIBP = ", data["Family", "NIBP"]];
Print["  Sectors = ", Length[data["Family", "SectorList"]]];

familyConfig = <|
  "Family" -> data["Family"],
  "Config" -> <|"Modulus" -> config["Modulus"]|>
|>;

Print["\n[4/5] Running ScheduledRegionSolve..."];
scheduleConfig = <|
  "Tiers" -> {
    <|"Name" -> "Quick",  "Timeout" -> 300,   "Order" -> "First"|>,
    <|"Name" -> "Fast",   "Timeout" -> 3600,  "Order" -> "First"|>,
    <|"Name" -> "Retry",  "Timeout" -> 86400, "Order" -> "FIFO"|>
  },
  "MaxRetriesPerTier" -> 1,
  "CacheDir" -> FileNameJoin[{$ProjectRoot, "workspace", "mma-Tri-region-solve", "cache"}],
  "OutputDir" -> FileNameJoin[{$ProjectRoot, "workspace", "mma-Tri-region-solve"}],
  "LogDir" -> FileNameJoin[{$ProjectRoot, "workspace", "mma-Tri-region-solve", "cache"}],
  "ScheduleMode" -> "Interleaved",
  "MaxParallelJobs" -> 1,
  "CheckpointInterval" -> 1,
  "ParallelJobs" -> 1,
  "Verbose" -> True
|>;

result = ScheduledRegionSolve["Tri", familyConfig, scheduleConfig];

Print["\n[5/5] Result summary"];
If[result === $Failed,
  Print["  FAILED"];
  Exit[1]
];

status = result["Status"];
doneCount = Length[Select[status, #["State"] === "done" &]];
trivialCount = Length[Select[status, #["State"] === "trivial" &]];
Print["  Done: ", doneCount];
Print["  Trivial: ", trivialCount];

ibpFile = FileNameJoin[{$ProjectRoot, "workspace", "mma-Tri-region-solve", "IBPMat_Tri.bin"}];
ringFile = FileNameJoin[{$ProjectRoot, "workspace", "mma-Tri-region-solve", "RingData_Tri.bin"}];
Print["\nOutput:"];
Print["  ", ibpFile, " ", If[FileExistsQ[ibpFile], "EXISTS", "MISSING"]];
Print["  ", ringFile, " ", If[FileExistsQ[ringFile], "EXISTS", "MISSING"]];
If[FileExistsQ[ibpFile], Print["  IBPMat size: ", FileByteCount[ibpFile], " bytes"]];
If[FileExistsQ[ringFile], Print["  RingData size: ", FileByteCount[ringFile], " bytes"]];

Print["\nDone."];
