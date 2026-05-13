(* ::Package:: *)
(* run.wl -- Cascading-timeout sector recompute for DP323 *)
(* Usage: wolframscript -file run.wl *)
(* Resumable: safe to kill and restart; per-sector cache prevents data loss *)
(* Uses ScheduledRegionSolver for generic tiered timeout scheduling *)

$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$FamilyGeneratePath = $LIECPPPath <> "/verify/DP323/";
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$RunDir = $LIECPPPath <> "/verify/DP323_recompute/";

(* ---- Load packages ---- *)
SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["ScheduledRegionSolver.wl"];
SetDirectory[$RunDir];

(* ---- Load checkpoint ---- *)
Print["[Step 1] Loading checkpoint from: ", $FamilyGeneratePath];
data = Import[$FamilyGeneratePath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
If[Head[data["Regions"]] === Join,
  Print["  Note: checkpoint Regions corrupted (Join), extracting original"];
  data["Regions"] = data["Regions"][[1]];
];

(* ---- Run ScheduledRegionSolve ---- *)
result = ScheduledRegionSolver`ScheduledRegionSolve[
  "DP323",
  data,
  <|
    "Tiers" -> {
      <|"Name" -> "Fast",  "Timeout" -> 1200|>,
      <|"Name" -> "Retry", "Timeout" -> 12000|>
    },
    "CacheDir" -> $RunDir <> "cache/",
    "OutputDir" -> $RunDir <> "output/",
    "CheckpointInterval" -> 5,
    "Verbose" -> True
  |>
];

If[result === $Failed,
  Print["Another instance is already running. Exiting."];
  Exit[1]
];

Print["\nAll done. Failed sectors: ", Length[result["FailedSectors"]]];
