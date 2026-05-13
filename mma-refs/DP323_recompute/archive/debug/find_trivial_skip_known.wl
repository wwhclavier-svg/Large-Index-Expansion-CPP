(* find_trivial_skip_known.wl — Skip known 36 trivial, test remaining 68 sectors *)
(* Each sector gets max 300s. Trivial ones finish in <10s. *)
(* Usage: wolframscript -file find_trivial_skip_known.wl *)

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

Print["Missing sectors: ", Length[missingAll]];

(* Known trivial from previous run - prefill_trivial.wl *)
knownTrivial = {
  {0,0,0,0,0,0,1,1,0,0,0}, {0,0,0,0,0,1,0,1,0,0,0}, {0,0,0,0,0,1,1,0,0,0,0},
  {0,0,0,0,0,1,1,1,0,0,0}, {0,0,0,0,1,0,0,1,0,0,0}, {0,0,0,0,1,0,1,0,0,0,0},
  {0,0,0,0,1,0,1,1,0,0,0}, {0,0,0,0,1,1,0,0,0,0,0}, {0,0,0,0,1,1,1,0,0,0,0},
  {0,0,0,1,0,0,0,1,0,0,0}, {0,0,0,1,0,0,1,0,0,0,0}, {0,0,0,1,0,0,1,1,0,0,0},
  {0,0,0,1,0,1,0,0,0,0,0}, {0,0,0,1,0,1,1,0,0,0,0}, {0,0,0,1,1,0,0,0,0,0,0},
  {0,0,0,1,1,0,0,1,0,0,0}, {0,0,0,1,1,0,1,0,0,0,0}, {0,0,0,1,1,0,1,1,0,0,0},
  {0,0,0,1,1,1,0,0,0,0,0}, {0,0,1,0,0,0,0,1,0,0,0}, {0,0,1,0,0,0,1,0,0,0,0},
  {0,0,1,0,0,0,1,1,0,0,0}, {0,0,1,0,1,0,0,0,0,0,0}, {0,0,1,0,1,0,0,1,0,0,0},
  {0,0,1,0,1,0,1,0,0,0,0}, {0,0,1,0,1,0,1,1,0,0,0}, {0,0,1,1,0,0,0,0,0,0,0},
  {0,0,1,1,1,0,0,0,0,0,0}, {0,1,0,0,0,0,0,1,0,0,0}, {0,1,0,0,0,0,1,0,0,0,0},
  {0,1,0,0,0,0,1,1,0,0,0}, {0,1,0,0,1,0,0,0,0,0,0}, {0,1,0,0,1,0,1,0,0,0,0},
  {0,1,0,1,0,0,0,0,0,0,0}, {0,1,0,1,0,0,1,0,0,0,0}, {0,1,0,1,1,0,0,0,0,0,0}
};
knownTrivialStr = SectorToStr /@ knownTrivial;

(* Filter: only test sectors NOT in known trivial list *)
toTest = Select[missingAll, !MemberQ[knownTrivialStr, SectorToStr[#]] &];
Print["Already known trivial: ", Length[knownTrivial]];
Print["Remaining to test:    ", Length[toTest]];

trivialList = {};
nonTrivialList = {};
timeoutList = {};
TIMEOUT = 300; (* each sector gets 5 min max *)

$SingularRunTimeout = TIMEOUT;

Do[
  sec = toTest[[i]];
  nProps = Total[sec[[1;;8]]];
  Print["\n[", i, "/", Length[toTest], "] Testing sector: ", sec, " (", nProps, " props) at ", DateString[]];

  {elapsed, result} = AbsoluteTiming[
    TimeConstrained[
      LIERegions`regionsBySectors[ibpeqs, {sec}, Alist, vlist,
        Modulus -> char,
        "EnableFieldExtension" -> True,
        Verbose -> False,
        "Timeout" -> TIMEOUT - 10,
        "TimeoutMarker" -> $Timeout
      ],
      TIMEOUT
    ]
  ];

  If[result === $Aborted || !AssociationQ[result],
    Print["  -> TIMEOUT/ABORT at ", Round[elapsed, 0.1], "s"];
    AppendTo[timeoutList, sec];
    ,
    If[KeyExistsQ[result, sec] && result[sec] =!= $Timeout,
      nRegions = Length[result[sec]];
      If[nRegions == 0,
        Print["  -> TRIVIAL (0 regions) in ", Round[elapsed, 0.1], "s"];
        AppendTo[trivialList, sec];
        ,
        Print["  -> NON-TRIVIAL: ", nRegions, " regions in ", Round[elapsed, 0.1], "s"];
        AppendTo[nonTrivialList, sec];
      ];
      ,
      Print["  -> TRIVIAL (no sector key) in ", Round[elapsed, 0.1], "s"];
      AppendTo[trivialList, sec];
    ];
  ];

  (* Save progress *)
  Export[FileNameJoin[{$CacheDir, "trivial_scan_progress.wdx"}],
    <|"Trivial" -> trivialList, "NonTrivial" -> nonTrivialList, "Timeout" -> timeoutList,
      "Tested" -> i, "Timestamp" -> DateString[]|>, "WDX"];
, {i, Length[toTest]}];

(* Final report *)
allTrivial = Join[knownTrivial, trivialList];
Print["\n========================================"];
Print["FINAL RESULTS"];
Print["========================================"];
Print["Known trivial:        ", Length[knownTrivial]];
Print["Newly found trivial:  ", Length[trivialList]];
Print["Total trivial:        ", Length[allTrivial]];
Print["Non-trivial:          ", Length[nonTrivialList]];
Print["Timeout:              ", Length[timeoutList]];

If[Length[trivialList] > 0,
  Print["\nNew trivial sectors:"];
  Do[Print["  ", s], {s, trivialList}];
];

If[Length[allTrivial] > 0,
  Print["\nAll trivial sectors (", Length[allTrivial], "):"];
  Do[Print["  ", SectorToStr[s]], {s, allTrivial}];
];

Export[FileNameJoin[{$CacheDir, "trivial_scan_final.wdx"}],
  <|"AllTrivial" -> allTrivial, "NonTrivial" -> nonTrivialList, "Timeout" -> timeoutList,
    "Timestamp" -> DateString[]|>, "WDX"];
Print["\nSaved to: ", FileNameJoin[{$CacheDir, "trivial_scan_final.wdx"}]];
