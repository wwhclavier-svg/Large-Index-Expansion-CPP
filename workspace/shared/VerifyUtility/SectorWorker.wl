(* SectorWorker.wl -- Standalone sector computation worker *)
(* Usage: wolframscript -file SectorWorker.wl <cacheDir> <sectorKey> <timeout> <tierName> *)
(* Used by ScheduledRegionSolver in Interleaved mode *)

If[Length[$ScriptCommandLine] < 5,
  Print["Usage: wolframscript -file SectorWorker.wl <cacheDir> <sectorKey> <timeout> <tierName>"];
  Exit[1]
];

$CacheDir          = $ScriptCommandLine[[2]];
$SectorKey         = $ScriptCommandLine[[3]];
$Timeout           = ToExpression[$ScriptCommandLine[[4]]];
$TierName          = $ScriptCommandLine[[5]];
$VerifyUtilityPath = DirectoryName[$InputFileName];

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];

(* Load shared family data *)
$SharedData     = Import[FileNameJoin[{$CacheDir, "worker_data.wdx"}], "WDX"];
$ibpeqs         = $SharedData["IBPEqs"];
$Alist          = $SharedData["AList"];
$vlist          = $SharedData["VList"];
$char           = $SharedData["Modulus"];
$sector         = ToExpression /@ Characters[$SectorKey];
$statusFile     = FileNameJoin[{$CacheDir, "status.wdx"}];
$sectorCacheDir = FileNameJoin[{$CacheDir, "SectorCache"}];
$doneFile       = FileNameJoin[{$sectorCacheDir, "sector_" <> $SectorKey <> ".done"}];

(* Skip if already completed by another worker *)
status = Import[$statusFile, "WDX"];
If[KeyExistsQ[status, $SectorKey] && MemberQ[{"done", "trivial", "cached"}, status[$SectorKey, "State"]],
  Print["Sector ", $SectorKey, " already completed (state: ", status[$SectorKey, "State"], "). Skipping."];
  Export[$doneFile, "already-done", "Text"];
  Print["[SECTOR_RESULT] key=", $SectorKey, " tier=", $TierName, " state=skipped time=0 regions=0"];
  Exit[0];
];

(* Compute strategy:
   1. Give regionsBySectors a generous Timeout (2x real timeout) so its internal
      Exit[]-based timeout never fires — TimeConstrained handles the real limit.
   2. TimeConstrained + CheckAbort gracefully returns $timeoutMarker on abort
      instead of killing the process. *)
$timeoutMarker = Unique["$Timeout"];
{elapsed, result} = AbsoluteTiming[
  CheckAbort[
    TimeConstrained[
      LIERegions`regionsBySectors[$ibpeqs, {$sector}, $Alist, $vlist,
        Modulus -> $char,
        "EnableFieldExtension" -> True,
        Verbose -> False,
        "Timeout" -> 2 * $Timeout,
        "TimeoutMarker" -> $timeoutMarker
      ],
      $Timeout,
      $timeoutMarker
    ],
    $timeoutMarker
  ]
];

If[!AssociationQ[result],
  res = <|"State" -> "error", "Time" -> elapsed, "Regions" -> 0|>,
  If[result === $timeoutMarker,
    res = <|"State" -> "timeout", "Time" -> elapsed, "Regions" -> 0|>,
    If[KeyExistsQ[result, $sector] && result[$sector] === $timeoutMarker,
      res = <|"State" -> "timeout", "Time" -> elapsed, "Regions" -> 0|>,
      If[KeyExistsQ[result, $sector],
        nRegions = Length[result[$sector]];
        res = <|"State" -> "done", "Time" -> elapsed, "Regions" -> nRegions, "Result" -> result[$sector]|>,
        res = <|"State" -> "trivial", "Time" -> elapsed, "Regions" -> 0|>
      ]
    ]
  ]
];

(* Write result data FIRST, before status update -- survives partial failures *)
If[res["State"] === "done" && res["Regions"] > 0,
  Export[FileNameJoin[{$sectorCacheDir, "sector_" <> $SectorKey <> ".wdx"}], res["Result"], "WDX"];
];

(* Lock-based atomic status update to prevent concurrent write races *)
$statusLock = $statusFile <> ".lock";
While[FileExistsQ[$statusLock], Pause[RandomReal[{0.05, 0.2}]]];
Export[$statusLock, ToString[$ProcessID], "Text"];

status = Import[$statusFile, "WDX"];
Switch[res["State"],
  "timeout",
    status[$SectorKey] = <|"State" -> $TierName <> "-timeout", "Regions" -> 0, "Time" -> res["Time"], "Tier" -> $TierName|>,
  "error",
    status[$SectorKey] = <|"State" -> $TierName <> "-error", "Regions" -> 0, "Time" -> res["Time"], "Tier" -> $TierName|>,
  "done",
    If[res["Regions"] > 0,
      status[$SectorKey] = <|"State" -> "done", "Regions" -> res["Regions"], "Time" -> res["Time"], "Tier" -> $TierName|>;,
      status[$SectorKey] = <|"State" -> "trivial", "Regions" -> 0, "Time" -> res["Time"], "Tier" -> $TierName|>;
    ];,
  "trivial",
    status[$SectorKey] = <|"State" -> "trivial", "Regions" -> 0, "Time" -> res["Time"], "Tier" -> $TierName|>
];
tmpStatus = $statusFile <> ".tmp." <> ToString[$ProcessID];
Export[tmpStatus, status, "WDX"];
RenameFile[tmpStatus, $statusFile, OverwriteTarget -> True];
DeleteFile[$statusLock];

(* Signal completion *)
Export[$doneFile, res["State"], "Text"];

Print["[SECTOR_RESULT] key=", $SectorKey, " tier=", $TierName,
      " state=", res["State"], " time=", res["Time"],
      " regions=", res["Regions"]];

(* Exit code *)
If[res["State"] === "timeout", Exit[1], Exit[0]];
