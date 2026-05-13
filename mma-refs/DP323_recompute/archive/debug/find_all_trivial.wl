(* find_all_trivial.wl — Quickly test all 104 missing sectors with 60s timeout *)
(* Identifies trivial (0-region) sectors and updates prefill_trivial.wl *)
(* Usage: wolframscript -file find_all_trivial.wl *)

$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$FamilyGeneratePath = $LIECPPPath <> "/verify/DP323/";
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$RunDir = $LIECPPPath <> "/verify/DP323_recompute/";
$CacheDir = $RunDir <> "cache/";

SectorToStr[sector_List] := StringJoin[ToString /@ sector];

(* Load packages *)
SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$RunDir];

(* Load DP323 checkpoint *)
Print["Loading DP323 checkpoint..."];
Block[{Print = If[#[[1]] === "Warning", Print[#]&, Identity]&},
  data = Import[$FamilyGeneratePath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
];
If[Head[data["Regions"]] === Join,
  data["Regions"] = data["Regions"][[1]];
];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
fullSectors = data["Family", "SectorList"];
existingRegions = data["Regions"];
missingAll = Complement[fullSectors, Keys[existingRegions]];

Print["Missing sectors to test: ", Length[missingAll]];
Print["Modulus: ", char];

(* Categorize *)
trivialList = {};
nonTrivialList = {};
timeoutList = {};
TIMEOUT = 60; (* trivial sectors finish in <10s *)

Do[
  sec = missingAll[[i]];
  nProps = Total[sec[[1;;8]]];
  Print["\n[", i, "/", Length[missingAll], "] Testing sector: ", sec, " (", nProps, " props)"];

  {elapsed, result} = AbsoluteTiming[
    TimeConstrained[
      LIERegions`regionsBySectors[ibpeqs, {sec}, Alist, vlist,
        Modulus -> char,
        "EnableFieldExtension" -> True,
        Verbose -> False,
        "Timeout" -> TIMEOUT - 5,  (* slightly less to give margin *)
        "TimeoutMarker" -> $Timeout
      ],
      TIMEOUT
    ]
  ];

  If[result === $Aborted || !AssociationQ[result],
    Print["  -> TIMEOUT at ", Round[elapsed, 0.1], "s"];
    AppendTo[timeoutList, sec];
    ,
    If[KeyExistsQ[result, sec] && result[sec] =!= $Timeout,
      nRegions = Length[result[sec]];
      If[nRegions == 0,
        Print["  -> TRIVIAL (0 regions) in ", Round[elapsed, 0.1], "s"];
        AppendTo[trivialList, sec];
        ,
        Print["  -> NON-TRIVIAL (", nRegions, " regions) in ", Round[elapsed, 0.1], "s"];
        AppendTo[nonTrivialList, sec];
      ];
      ,
      Print["  -> TRIVIAL (no sector key) in ", Round[elapsed, 0.1], "s"];
      AppendTo[trivialList, sec];
    ];
  ];
, {i, Length[missingAll]}];

(* Report *)
Print["\n========================================"];
Print["RESULTS"];
Print["========================================"];
Print["Trivial (0 regions): ", Length[trivialList]];
Print["Non-trivial:         ", Length[nonTrivialList]];
Print["Timeout:             ", Length[timeoutList]];

If[Length[trivialList] > 0,
  Print["\nTrivial sectors:"];
  Do[Print["  ", s], {s, trivialList}];
];

If[Length[nonTrivialList] > 0,
  Print["\nNon-trivial sectors:"];
  Do[Print["  ", s, " (", Total[s[[1;;8]]], " props)"], {s, nonTrivialList}];
];

If[Length[timeoutList] > 0,
  Print["\nTimeout sectors:"];
  Do[Print["  ", s, " (", Total[s[[1;;8]]], " props)"], {s, timeoutList}];
];

(* Save results *)
Export[FileNameJoin[{$CacheDir, "trivial_scan.wdx"}],
  <|"Trivial" -> trivialList, "NonTrivial" -> nonTrivialList, "Timeout" -> timeoutList, "Timestamp" -> DateString[]|>,
  "WDX"
];

Print["\nTotal trivial: ", Length[trivialList], " / ", Length[missingAll]];
Print["Saved to: ", FileNameJoin[{$CacheDir, "trivial_scan.wdx"}]];
