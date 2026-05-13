(* run_scheduled.wl -- Run ScheduledRegionSolve for SR212-5m *)
(* Usage: wolframscript -file verify/SR212-5m/run_scheduled.wl *)

$ProjectRoot = DirectoryName[DirectoryName[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = FileNameJoin[{$ProjectRoot, "verify", "VerifyUtility"}];
$FamilyDBPath = FileNameJoin[{$ProjectRoot, "verify", "FamilyDatabase", "FamilyDatabase.wl"}];

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["ScheduledRegionSolver.wl"];
Get["ExportBinary_IBPMatrix.wl"];

Get[$FamilyDBPath];
config = $FamilyDatabase["SR212-5m"];

scheduleConfig = <|
  "Tiers" -> {
    <|"Name" -> "Quick", "Timeout" -> 300,  "Order" -> "FIFO"|>,
    <|"Name" -> "Fast",  "Timeout" -> 3600, "Order" -> "FIFO"|>,
    <|"Name" -> "Retry", "Timeout" -> 36000,"Order" -> "FIFO"|>
  },
  "MaxRetriesPerTier" -> 1,
  "ParallelJobs" -> 8,
  "CacheDir" -> "verify/SR212-5m/cache/",
  "OutputDir" -> "verify/SR212-5m/"
|>;

Print["Starting SR212-5m ScheduledRegionSolve (24 sectors)..."];
result = ScheduledRegionSolve["SR212-5m", config, scheduleConfig];

Print["\nExporting IBPMat and RingData..."];
ExportIBPMat["SR212-5m", config, "Output" -> "verify/SR212-5m/"];
ExportRingData["SR212-5m", config, "Output" -> "verify/SR212-5m/"];

Print["\nDone."];
