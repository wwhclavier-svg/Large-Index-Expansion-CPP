(* ::Package:: *)
(* ScheduledRegionSolver.wl -- Generic tiered timeout scheduler for region solving *)
(* Depends on: LIEWorkflow.wl, LIERegions.wl, SingularInterface.wl *)

BeginPackage["ScheduledRegionSolver`", {"LIEUtility`"}];

(* Public API *)
ScheduledRegionSolve::usage = "ScheduledRegionSolve[familyName, familyConfig, scheduleConfig] computes regions for all sectors of a family with tiered timeout scheduling.";

$DefaultScheduleConfig::usage = "Default schedule configuration for ScheduledRegionSolve.";

Begin["`Private`"];

$PackageDir = DirectoryName[$InputFileName];  (* capture at load time for worker path *)

(* ============================================================ *)
(* Default Configuration *)
(* ============================================================ *)

$DefaultScheduleConfig = <|
  "Tiers" -> {
    <|"Name" -> "Fast",    "Timeout" -> 1200,  "Order" -> "First"|>,
    <|"Name" -> "Retry",   "Timeout" -> 12000, "Order" -> "FIFO"|>
  },
  "MaxRetriesPerTier" -> 1,
  "ScheduleMode" -> "Sequential",   (* "Sequential" | "Interleaved" *)
  "MaxParallelJobs" -> 4,           (* Total concurrent workers in Interleaved mode *)
  "CacheDir" -> Automatic,
  "OutputDir" -> Automatic,
  "LogDir" -> Automatic,
  "CheckpointInterval" -> 5,
  "ParallelJobs" -> 1,
  "InstanceLockFile" -> Automatic,
  "Verbose" -> True
|>;

(* ============================================================ *)
(* P2.1 InstanceLock -- Single-instance guarantee *)
(* ============================================================ *)

AcquireInstanceLock[lockFile_String] := Module[{fh, pid, res},
  If[FileExistsQ[lockFile],
    pid = Quiet@Check[Import[lockFile, "Text"], ""];
    If[StringQ[pid] && StringLength[pid] > 0,
      res = Run["kill -0 " <> pid <> " 2>/dev/null"];
      If[res === 0,
        Print["[LOCK] Another instance is running (PID ", pid, "). Exiting."];
        Return[$Failed],
        Print["[LOCK] Stale lock found (PID ", pid, " dead). Removing."];
        DeleteFile[lockFile]
      ],
      DeleteFile[lockFile]
    ]
  ];
  Export[lockFile, ToString[$ProcessID], "Text"];
  Return[Success]
];

ReleaseInstanceLock[lockFile_String] := If[FileExistsQ[lockFile], DeleteFile[lockFile]];

(* ============================================================ *)
(* P2.2 SectorStateManager -- Atomic state persistence *)
(* ============================================================ *)

LoadStatus[file_String] := Module[{s},
  If[!FileExistsQ[file], Return[<||>]];
  s = Import[file, "WDX"];
  If[AssociationQ[s], s, <||>]
];

SaveStatusAtomic[file_String, status_Association] := Module[{tmpFile},
  tmpFile = file <> ".tmp";
  Export[tmpFile, status, "WDX"];
  RenameFile[tmpFile, file, OverwriteTarget -> True];
];

SectorToStr[sector_List] := StringJoin[ToString /@ sector];
StrToSector[str_String] := ToExpression /@ Characters[str];

(* ============================================================ *)
(* P2.4 RegionComputeEngine -- Unified compute interface *)
(* ============================================================ *)

ComputeSectorWithTimeout[sec_, ibpeqs_, Alist_, vlist_, char_, timeout_, timeoutMarker_] := Module[
  {result, elapsed, regionCount},
  {elapsed, result} = AbsoluteTiming[
    LIERegions`regionsBySectors[ibpeqs, {sec}, Alist, vlist,
      Modulus -> char,
      "EnableFieldExtension" -> True,
      Verbose -> False,
      "Timeout" -> timeout,
      "TimeoutMarker" -> timeoutMarker
    ]
  ];
  
  (* Case 1: result is not an Association -> outer failure/hang *)
  If[!AssociationQ[result],
    Return[<|"Result" -> $Failed, "Time" -> elapsed, "Regions" -> 0, "State" -> "error"|>]
  ];
  
  (* Case 2: inner timeout explicitly marked by regionsBySectors *)
  If[KeyExistsQ[result, sec] && result[sec] === timeoutMarker,
    Return[<|"Result" -> $Failed, "Time" -> elapsed, "Regions" -> 0, "State" -> "timeout"|>]
  ];
  
  (* Case 3: non-trivial success *)
  If[KeyExistsQ[result, sec],
    regionCount = Length[result[sec]];
    Return[<|"Result" -> result[sec], "Time" -> elapsed, "Regions" -> regionCount, "State" -> "done"|>]
  ];
  
  (* Case 4: trivial sector (0 regions) *)
  Return[<|"Result" -> {}, "Time" -> elapsed, "Regions" -> 0, "State" -> "done"|>]
];

(* ============================================================ *)
(* P2.3 TieredScheduler -- Generic tiered timeout scheduling *)
(* ============================================================ *)

ProcessTier[sectors_, tierConfig_, ibpeqs_, Alist_, vlist_, char_, status0_, statusFile_, sectorCacheDir_, timeoutMarker_, verboseQ_] := Module[
  {timeout = tierConfig["Timeout"], n = Length[sectors], i, sec, key, res, timeouts = {}, status = status0, sectorFile},
  
  If[verboseQ,
    Print["\n=== [Tier: ", tierConfig["Name"], "] Processing ", n, " sector(s) with timeout=", timeout, "s ==="]
  ];
  
  Do[
    sec = sectors[[i]];
    key = SectorToStr[sec];
    If[verboseQ, Print["\n[", i, "/", n, "] sector=", sec]];
    
    res = ComputeSectorWithTimeout[sec, ibpeqs, Alist, vlist, char, timeout, timeoutMarker];
    
    Switch[res["State"],
      "timeout",
      If[verboseQ, Print["  -> TIMEOUT at ", timeout, "s, will retry at higher timeout"]];
      status[key] = <|"State" -> tierConfig["Name"] <> "-timeout", "Regions" -> 0, "Time" -> res["Time"]|>;
      AppendTo[timeouts, sec];
      ,
      "error",
      If[verboseQ, Print["  -> ERROR after ", res["Time"], "s"]];
      status[key] = <|"State" -> "error", "Regions" -> 0, "Time" -> res["Time"]|>;
      AppendTo[timeouts, sec];
      ,
      "done",
      If[res["Regions"] > 0,
        If[verboseQ, Print["  -> DONE: ", res["Regions"], " regions in ", res["Time"], "s"]];
        status[key] = <|"State" -> "done", "Regions" -> res["Regions"], "Time" -> res["Time"]|>;
        (* Result saved to separate sector file immediately *)
        sectorFile = FileNameJoin[{sectorCacheDir, "sector_" <> key <> ".wdx"}];
        Export[sectorFile, res["Result"], "WDX"];
        ,
        If[verboseQ, Print["  -> TRIVIAL: 0 regions in ", res["Time"], "s"]];
        status[key] = <|"State" -> "trivial", "Regions" -> 0, "Time" -> res["Time"]|>;
      ];
    ];
    
    (* Atomic checkpoint after every sector - lightweight status only *)
    SaveStatusAtomic[statusFile, status];
    If[verboseQ && Mod[i, 5] === 0, Print["  [Status checkpoint saved]"]];
  , {i, n}];
  
   SaveStatusAtomic[statusFile, status];
  Export[FileNameJoin[{DirectoryName[statusFile], "status.json"}], status, "JSON"];
  {status, timeouts}
];

(* ============================================================ *)
(* P2.3b ProcessTiersInterleaved -- Parallel tier pipeline *)
(* ============================================================ *)

ProcessTiersInterleaved[initialQueues_, tiers_, ibpeqs_, Alist_, vlist_, char_,
  status0_, statusFile_, sectorCacheDir_, verboseQ_, maxParallel_] := Module[
  {nTiers = Length[tiers], queues, active, workerRecords, status = status0,
   workerScript, workerDataFile, sec, key, tierIdx, tierName, timeout,
   doneFile, workerResult, i, allDone},

  workerScript = FileNameJoin[{$PackageDir, "SectorWorker.wl"}];
  workerDataFile = FileNameJoin[{DirectoryName[statusFile], "worker_data.wdx"}];

  (* Export shared data for workers *)
  Export[workerDataFile, <|
    "IBPEqs" -> ibpeqs, "AList" -> Alist, "VList" -> vlist, "Modulus" -> char
  |>, "WDX"];

  (* Use pre-sorted initial queues for state-aware resume *)
  queues = initialQueues;
  active = ConstantArray[0, nTiers];
  workerRecords = ConstantArray[{}, nTiers]; (* {pid, sectorKey, startTime, timeout} per tier *)

  If[verboseQ,
    Print["\n========================================"];
    Print["Interleaved Tier Pipeline: ", StringRiffle[#["Name"]& /@ tiers, " + "]];
    Print["Max parallel workers: ", maxParallel];
    Print["========================================"];
  ];

  allDone = False;
  While[!allDone,
    (* --- 1. Reap completed workers --- *)
    totalActive = Total[active];
    Do[
      records = workerRecords[[tierIdx]];
      remaining = {};
      Do[
        rec = records[[j]];
        pid = rec[[1]]; sk = rec[[2]]; startTime = rec[[3]]; to = rec[[4]];
        (* Check if process still running *)
        If[Run["kill -0 " <> ToString[pid] <> " 2>/dev/null"] =!= 0,
          active[[tierIdx]]--;
          (* If no .done file, sector died before writing result — move to next tier *)
          If[!FileExistsQ[FileNameJoin[{sectorCacheDir, "sector_" <> sk <> ".done"}]],
            If[tierIdx < nTiers,
              AppendTo[queues[[tierIdx + 1]], StrToSector[sk]];
              If[verboseQ, Print["  [RE-QUEUE] ", sk, " -> ", tiers[[tierIdx + 1]]["Name"], " tier"];];
            ,
              status[sk] = <|"State" -> "failed"|>;
              If[verboseQ, Print["  [FAILED] ", sk, " — no more tiers"];];
            ];
          ];
        ,
        (* Process alive — check for hang: exceeded timeout + 120s grace *)
        If[SessionTime[] - startTime > to + 120,
          If[verboseQ, Print["  [HANG] ", sk, " exceeded ", to, "s+120s, killing PID=", pid]];
          Run["kill -9 " <> ToString[pid] <> " 2>/dev/null"];
          active[[tierIdx]]--;
          If[tierIdx < nTiers,
            AppendTo[queues[[tierIdx + 1]], StrToSector[sk]];
            If[verboseQ, Print["  [RE-QUEUE] ", sk, " -> ", tiers[[tierIdx + 1]]["Name"], " tier"];];
          ,
            status[sk] = <|"State" -> "failed"|>;
            If[verboseQ, Print["  [FAILED] ", sk, " — no more tiers"];];
          ];,
          AppendTo[remaining, rec];
        ];
        ];
      , {j, Length[records]}];
      workerRecords[[tierIdx]] = remaining;
    , {tierIdx, nTiers}];

    (* Scan doneFiles for completed workers — check all files in cache dir *)
    doneFiles = FileNames["sector_*.done", sectorCacheDir];
    Do[
      doneFile = doneFiles[[i]];
      sectorKey = StringReplace[FileNameTake[doneFile], {"sector_" -> "", ".done" -> ""}];
      workerResult = Import[doneFile, "Text"];
      DeleteFile[doneFile];
      (* Reload status updated by worker *)
      status = LoadStatus[statusFile];
      If[verboseQ,
        st = status[sectorKey, "State"];
        reg = status[sectorKey, "Regions"];
        tm = status[sectorKey, "Time"];
        Which[
          st === "done",       Print["  [DONE] ", sectorKey, " -> ", reg, " regions, ", Round[tm,0.1], "s"],
          StringContainsQ[st, "-timeout"], Print["  [TIMEOUT] ", sectorKey, " at ", Round[tm,0.1], "s (-> next tier)"],
          StringContainsQ[st, "-error"],   Print["  [ERROR] ", sectorKey, " at ", Round[tm,0.1], "s (-> next tier)"],
          st === "trivial",    Print["  [TRIVIAL] ", sectorKey, " 0 regions, ", Round[tm,0.1], "s"],
          st === "error",      Print["  [ERROR] ", sectorKey]
        ];
      ];
      (* If timeout or error, move to next tier *)
      tierName = StringReplace[StringReplace[status[sectorKey, "State"], "-error" -> ""], "-timeout" -> ""];
      currentTier = First[FirstPosition[tiers, _?(#["Name"] === tierName &), {Length[tiers]}]];
      If[(StringContainsQ[status[sectorKey, "State"], "-timeout"] ||
          StringContainsQ[status[sectorKey, "State"], "-error"]) && currentTier < nTiers,
        AppendTo[queues[[currentTier + 1]], StrToSector[sectorKey]];
      ];
    , {i, Length[doneFiles]}];
    (* Export readable JSON for live dashboard *)
    Export[FileNameJoin[{DirectoryName[statusFile], "status.json"}], status, "JSON"];

    (* --- 2. Escalate timeout sectors in queues --- *)
    (* Already handled by doneFile check above *)

    (* --- 3. Launch new workers (priority: low tier first) --- *)
    (* Per-tier parallelism limits: default to 1 if not specified in tier config *)
    tierMaxParallel = Table[Lookup[tiers[[i]], "MaxParallel", 1], {i, nTiers}];
    totalActive = Total[active];
    While[totalActive < maxParallel,
      launched = False;
      For[tierIdx = 1, tierIdx <= nTiers && totalActive < maxParallel && !launched, tierIdx++,
        tierSlots = tierMaxParallel[[tierIdx]] - active[[tierIdx]];
        If[Length[queues[[tierIdx]]] > 0 && tierSlots > 0 && totalActive < maxParallel,
          sec = First[queues[[tierIdx]]];
          queues[[tierIdx]] = Rest[queues[[tierIdx]]];
          key = SectorToStr[sec];
          tierName = tiers[[tierIdx]]["Name"];
          timeout = tiers[[tierIdx]]["Timeout"];

          If[verboseQ,
            Print["  [LAUNCH] ", key, " -> ", tierName, " tier (", timeout, "s)"];
          ];

          (* Launch worker in background *)
          workerLog = FileNameJoin[{DirectoryName[statusFile], "worker_" <> key <> ".log"}];
          cmd = "wolframscript -file \"" <> workerScript <> "\" \"" <>
                DirectoryName[statusFile] <> "\" \"" <> key <> "\" " <>
                ToString[timeout] <> " \"" <> tierName <> "\" >\"" <> workerLog <> "\" 2>&1 & echo $!";
          pidStr = StringTrim[RunProcess[{"bash", "-c", cmd}]["StandardOutput"]];
          If[StringQ[pidStr] && StringLength[pidStr] > 0,
            pid = ToExpression[pidStr];
            AppendTo[workerRecords[[tierIdx]], {pid, key, SessionTime[], timeout}];
            active[[tierIdx]]++;
            totalActive++;
            launched = True;
            (* Mark as running in status for live dashboard *)
            status[key] = <|"State" -> "running", "Regions" -> 0, "Time" -> 0, "Tier" -> tierName|>;
          If[verboseQ, Print["    [PID=", pid, ", active=", active, ", totalActive=", totalActive, "]"]];
          (* Remove doneFile from any previous run of this sector *)
          doneFile = FileNameJoin[{sectorCacheDir, "sector_" <> key <> ".done"}];
          If[FileExistsQ[doneFile], DeleteFile[doneFile]],
          If[verboseQ, Print["    [FAILED to launch: pidStr=", pidStr, "]"]];
        ];
      ];
    ];
    If[!launched, Break[]]; (* no tier could launch, exit *)
  ];

    (* --- 4. Check termination --- *)
    allDone = Total[Length /@ queues] === 0 && Total[active] === 0;
    (* If[verboseQ && !allDone, Print["  [LOOP] queues=", Length /@ queues, " active=", active, " allDone=", allDone]]; *)
    If[!allDone, Pause[3]];
  ];

  If[verboseQ, Print["\nInterleaved pipeline complete."]];
  Export[FileNameJoin[{DirectoryName[statusFile], "status.json"}], status, "JSON"];
  status
];

ScheduledRegionSolve[familyName_String, familyConfig_Association, scheduleConfig0_Association:<||>] := Module[
  {scheduleConfig, cacheDir, outputDir, lockFile, statusFile, 
   data, existingRegions, fullSectorList, char, ibpeqs, Alist, vlist,
   status, pending, timeouts, failedSectors, tier, i, 
   allRegions, binFile, ringFile, checkpointInterval, verboseQ,
   timeoutMarker = Unique["$Timeout"]},
  
  (* Merge schedule config with defaults *)
  scheduleConfig = Merge[{$DefaultScheduleConfig, scheduleConfig0}, Last];
  
  (* Resolve directories *)
  cacheDir = If[scheduleConfig["CacheDir"] === Automatic,
    FileNameJoin[{Directory[], "cache"}],
    scheduleConfig["CacheDir"]
  ];
  outputDir = If[scheduleConfig["OutputDir"] === Automatic,
    FileNameJoin[{Directory[], "output"}],
    scheduleConfig["OutputDir"]
  ];
  If[!FileExistsQ[cacheDir], CreateDirectory[cacheDir]];
  If[!FileExistsQ[outputDir], CreateDirectory[outputDir]];
  
  lockFile = If[scheduleConfig["InstanceLockFile"] === Automatic,
    FileNameJoin[{cacheDir, "instance.lock"}],
    scheduleConfig["InstanceLockFile"]
  ];
  statusFile = FileNameJoin[{cacheDir, "status.wdx"}];
  checkpointInterval = scheduleConfig["CheckpointInterval"];
  verboseQ = scheduleConfig["Verbose"];
  
  (* P2.1: Acquire instance lock *)
  If[AcquireInstanceLock[lockFile] === $Failed,
    Return[$Failed]
  ];
  
  result = $Failed;
  Check[
    (* Load family data from config *)
    data = familyConfig;
    existingRegions = If[KeyExistsQ[data, "Regions"], data["Regions"], <||>];
    If[Head[existingRegions] === Join, existingRegions = existingRegions[[1]]];
    fullSectorList = data["Family", "SectorList"];
    char = data["Config", "Modulus"];
    ibpeqs = data["Family", "IBPEqs"];
    Alist = data["Family", "AList"];
    vlist = data["Family", "VList"];
    
    If[verboseQ,
      Print["========================================"];
      Print["ScheduledRegionSolve: ", familyName];
      Print["Sectors: ", Length[fullSectorList]];
      Print["Modulus: ", char];
      Print["Cache:   ", cacheDir];
      Print["Output:  ", outputDir];
      Print["========================================"];
    ];
    
    (* P2.2: Load or init status *)
    status = LoadStatus[statusFile];
    If[status === <||>,
      status = <||>;
      Do[
        sec = fullSectorList[[i]];
        key = SectorToStr[sec];
        (* Mark sectors with existing checkpoint data as cached *)
        If[KeyExistsQ[existingRegions, sec] && existingRegions[sec] =!= {},
          status[key] = <|"State" -> "cached", "Regions" -> Length[existingRegions[sec]]|>,
          status[key] = <|"State" -> "pending"|>
        ]
      , {i, Length[fullSectorList]}];
      SaveStatusAtomic[statusFile, status];
    ];
    
    (* Sector cache directory *)
    sectorCacheDir = FileNameJoin[{cacheDir, "SectorCache"}];
    If[!FileExistsQ[sectorCacheDir], CreateDirectory[sectorCacheDir]];
    
    (* Identify pending sectors: skip done/trivial/cached; sort into tiers by status *)
    (* State-aware resume: sectors with "X-timeout" skip tier X and start at X+1 *)
    pendingAll = Select[fullSectorList,
      !MemberQ[{"done", "trivial", "cached"}, status[SectorToStr[#], "State"]] &&
      !(KeyExistsQ[existingRegions, #] && existingRegions[#] =!= {}) &];
    
    (* Dynamically build tier queues from config tier names *)
    tierNames = Lookup[#, "Name"] & /@ scheduleConfig["Tiers"];
    nTiers = Length[tierNames];
    pendingAllForTier = Table[{}, {nTiers}];  (* one queue per tier *)
    Do[
      sec = pendingAll[[i]];
      st = status[SectorToStr[sec], "State"];
      (* Extract source tier name from state suffix, if present *)
      matched = StringCases[st, start:__ ~~ ("-timeout" | "-error") :> start];
      If[Length[matched] > 0,
        srcTier = matched[[1]];
        srcIdx = Position[tierNames, srcTier][[1, 1]] - 1;  (* 0-indexed *)
        If[srcIdx < nTiers - 1,
          AppendTo[pendingAllForTier[[srcIdx + 2]], sec],  (* next tier *)
          AppendTo[pendingAllForTier[[nTiers]], sec]        (* last tier → will be marked failed *)
        ]
      ,
        (* No tier suffix: route by special state *)
        Which[
          st === "error",                                  AppendTo[pendingAllForTier[[nTiers]], sec],  (* bare error → last tier *)
          st === "running",                                AppendTo[pendingAllForTier[[1]], sec],       (* killed mid-run → restart from Quick *)
          st === "failed",                                 AppendTo[pendingAllForTier[[nTiers]], sec],  (* failed → last tier *)
          StringContainsQ[st, "Long-"],                    AppendTo[pendingAllForTier[[nTiers]], sec],  (* Long-* → last tier *)
          True,                                           AppendTo[pendingAllForTier[[1]], sec]        (* pending/unknown → Quick *)
        ]
      ]
    , {i, Length[pendingAll]}];
    If[verboseQ,
      Print["\nPending sectors: ", Length[pendingAll], "/", Length[fullSectorList]];
      Do[
        tierName = tierNames[[ti]];
        nPending = Length[pendingAllForTier[[ti]]];
        If[nPending > 0, Print["  ", tierName, ": ", nPending]];
      , {ti, nTiers}];
      Print["  (Skipped ", Length[fullSectorList] - Length[pendingAll], " already-computed/trivial/cached sectors)"];
    ];
    (* P2.3: Process tiers *)
    If[scheduleConfig["ScheduleMode"] === "Interleaved",
      (* Interleaved: parallel tier pipeline *)
      status = ProcessTiersInterleaved[
        pendingAllForTier, scheduleConfig["Tiers"], ibpeqs, Alist, vlist, char,
        status, statusFile, sectorCacheDir, verboseQ,
        scheduleConfig["MaxParallelJobs"]
      ];
      ,
      (* Sequential: original tier-by-tier processing *)
      seqPending = Join @@ pendingAllForTier;  (* all pending sectors flattened *)
      Do[
        tier = scheduleConfig["Tiers"][[i]];
        If[Length[seqPending] === 0, Break[]];
        {status, timeouts} = ProcessTier[
          seqPending, tier, ibpeqs, Alist, vlist, char,
          status, statusFile, sectorCacheDir, timeoutMarker, verboseQ
        ];
        seqPending = timeouts;
        If[verboseQ, Print["\n  Tier '", tier["Name"], "' complete. Success: ",
          Length[seqPending], " remaining for next tier."]];
      , {i, Length[scheduleConfig["Tiers"]]}];
    ];
    
    (* Final JSON export for live dashboard *)
    Export[FileNameJoin[{cacheDir, "status.json"}], status, "JSON"];

    (* Mark final failures: exclude cached sectors from checkpoint *)
    failedSectors = Select[fullSectorList, 
      !MemberQ[{"done", "trivial", "cached"}, status[SectorToStr[#], "State"]] &];
    If[Length[failedSectors] > 0 && verboseQ,
      Print["\nFailed sectors (", Length[failedSectors], "): ", failedSectors];
    ];
    
    (* Mark remaining non-completed sectors as "failed" in status *)
    Scan[Module[{key = SectorToStr[#]},
      If[!MemberQ[{"done", "trivial", "cached"}, status[key, "State"]],
        status[key] = <|"State" -> "failed", "Regions" -> 0,
          "Time" -> Lookup[status[key], "Time", 0], "Tier" -> Lookup[status[key], "Tier", ""]|>
      ]
    ]&, fullSectorList];
    Export[FileNameJoin[{cacheDir, "status.json"}], status, "JSON"];
    
    (* Build final regions: association keyed by sector name, each value = list of region associations *)
    allRegions = Association@Table[
      sec = fullSectorList[[i]];
      key = SectorToStr[sec];
      key -> Which[
        MemberQ[{"done", "trivial"}, status[key, "State"]],
          Module[{sf = FileNameJoin[{sectorCacheDir, "sector_" <> key <> ".wdx"}], data},
            data = If[FileExistsQ[sf], Import[sf, "WDX"], {}];
            If[AssociationQ[data], {data}, If[ListQ[data], data, {}]]
          ],
        True,
          existingRegions[sec]
      ],
      {i, Length[fullSectorList]}
    ];
    
    (* Export *)
    binFile = FileNameJoin[{outputDir, "IBPMat_" <> familyName <> ".bin"}];
    ringFile = FileNameJoin[{outputDir, "RingData_" <> familyName <> ".bin"}];
    ExportIBPMatrixBinary`ExportBinaryIBPMatrix[binFile, allRegions, char];
    ExportIBPMatrixBinary`ExportBinaryRingData[ringFile, allRegions, Alist, Length[Alist], char];
    If[verboseQ,
      Print["\nExport complete:"];
      Print["  ", binFile];
      Print["  ", ringFile];
    ];
    
    (* Save final checkpoint *)
    data["Regions"] = Association@Table[
      SectorToStr[fullSectorList[[i]]] -> allRegions[[i]],
      {i, Length[fullSectorList]}
    ];
    Export[FileNameJoin[{outputDir, "PrepareCheckpoint-" <> familyName <> ".wdx"}], data, "WDX"];
    
    Print["\nScheduledRegionSolve complete for ", familyName, "."];
    result = <|"Status" -> status, "FailedSectors" -> failedSectors|>,
    Null
  ];
  
  ReleaseInstanceLock[lockFile];
  Print["[LOCK] Released."];
  result
];

End[];
EndPackage[];
