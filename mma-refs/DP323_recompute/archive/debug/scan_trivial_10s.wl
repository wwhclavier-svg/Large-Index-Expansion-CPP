(* scan_trivial_10s.wl — 10s cutoff: trivial finish fast, non-trivial skip *)
(* Usage: wolframscript -file scan_trivial_10s.wl *)

$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$FamilyGeneratePath = $LIECPPPath <> "/verify/DP323/";
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$RunDir = $LIECPPPath <> "/verify/DP323_recompute/";
$CacheDir = $RunDir <> "cache/";

SectorToStr[sector_List] := StringJoin[ToString /@ sector];

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$RunDir];

Print["Loading DP323 checkpoint..."];
Block[{Print = If[#[[1]] === "Warning", Print[#]&, Identity]&},
  data = Import[$FamilyGeneratePath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
fullSectors = data["Family", "SectorList"];
existingRegions = data["Regions"];
missingAll = Complement[fullSectors, Keys[existingRegions]];
Print["Missing sectors: ", Length[missingAll]];

(* Known trivial from previous run *)
knownTrivial = {
  {0,0,0,0,0,0,1,1,0,0,0},{0,0,0,0,0,1,0,1,0,0,0},{0,0,0,0,0,1,1,0,0,0,0},
  {0,0,0,0,0,1,1,1,0,0,0},{0,0,0,0,1,0,0,1,0,0,0},{0,0,0,0,1,0,1,0,0,0,0},
  {0,0,0,0,1,0,1,1,0,0,0},{0,0,0,0,1,1,0,0,0,0,0},{0,0,0,0,1,1,1,0,0,0,0},
  {0,0,0,1,0,0,0,1,0,0,0},{0,0,0,1,0,0,1,0,0,0,0},{0,0,0,1,0,0,1,1,0,0,0},
  {0,0,0,1,0,1,0,0,0,0,0},{0,0,0,1,0,1,1,0,0,0,0},{0,0,0,1,1,0,0,0,0,0,0},
  {0,0,0,1,1,0,0,1,0,0,0},{0,0,0,1,1,0,1,0,0,0,0},{0,0,0,1,1,0,1,1,0,0,0},
  {0,0,0,1,1,1,0,0,0,0,0},{0,0,1,0,0,0,0,1,0,0,0},{0,0,1,0,0,0,1,0,0,0,0},
  {0,0,1,0,0,0,1,1,0,0,0},{0,0,1,0,1,0,0,0,0,0,0},{0,0,1,0,1,0,0,1,0,0,0},
  {0,0,1,0,1,0,1,0,0,0,0},{0,0,1,0,1,0,1,1,0,0,0},{0,0,1,1,0,0,0,0,0,0,0},
  {0,0,1,1,1,0,0,0,0,0,0},{0,1,0,0,0,0,0,1,0,0,0},{0,1,0,0,0,0,1,0,0,0,0},
  {0,1,0,0,0,0,1,1,0,0,0},{0,1,0,0,1,0,0,0,0,0,0},{0,1,0,0,1,0,1,0,0,0,0},
  {0,1,0,1,0,0,0,0,0,0,0},{0,1,0,1,0,0,1,0,0,0,0},{0,1,0,1,1,0,0,0,0,0,0}
};
knownTrivialStr = SectorToStr /@ knownTrivial;

toTest = Select[missingAll, !MemberQ[knownTrivialStr, SectorToStr[#]] &];
Print["Known trivial: ", Length[knownTrivial], ", to test: ", Length[toTest]];

(* Set Singular subprocess timeout to 10s *)
$SingularTimeout = 3;

trivialList = {};
nonTrivialList = {};
TIMEOUT = 3;

Do[
  sec = toTest[[i]];
  nProps = Total[sec[[1;;8]]];
  Print["[", i, "/", Length[toTest], "] ", SectorToStr[sec], " (", nProps, " props)... "];

  {elapsed, result} = AbsoluteTiming[
    TimeConstrained[
      LIERegions`regionsBySectors[ibpeqs, {sec}, Alist, vlist,
        Modulus -> char, "EnableFieldExtension" -> True, Verbose -> False,
        "Timeout" -> TIMEOUT, "TimeoutMarker" -> $Timeout
      ],
      TIMEOUT + 3
    ]
  ];

  Which[
    result === $Aborted || !AssociationQ[result],
      Print["  skip (>", Round[elapsed, 0.1], "s)"];
      ,
    KeyExistsQ[result, sec] && result[sec] === $Timeout,
      Print["  skip (inner timeout)"];
      ,
    KeyExistsQ[result, sec] && Length[result[sec]] == 0,
      Print["  TRIVIAL (", Round[elapsed, 0.1], "s)"];
      AppendTo[trivialList, sec];
      ,
    KeyExistsQ[result, sec],
      Print["  NON-TRIVIAL: ", Length[result[sec]], " regions (", Round[elapsed,0.1], "s)"];
      AppendTo[nonTrivialList, sec];
      ,
    AssociationQ[result],
      Print["  TRIVIAL (", Round[elapsed, 0.1], "s)"];
      AppendTo[trivialList, sec];
      ,
    True,
      Print["  skip (error, ", Round[elapsed,0.1], "s)"];
  ];

  (* Save incremental *)
  Export[FileNameJoin[{$CacheDir, "scan10s_progress.wdx"}],
    <|"Trivial" -> trivialList, "NonTrivial" -> nonTrivialList,
      "Tested" -> i, "Timestamp" -> DateString[]|>, "WDX"];
, {i, Length[toTest]}];

allTrivial = Join[knownTrivial, trivialList];
Print["\n========================================"];
Print["Known trivial:  ", Length[knownTrivial]];
Print["New trivial:    ", Length[trivialList]];
Print["Total trivial:  ", Length[allTrivial]];
Print["Non-trivial:    ", Length[nonTrivialList]];
Print["Untested:       ", Length[toTest] - Length[trivialList] - Length[nonTrivialList]];

If[Length[trivialList] > 0,
  Print["\nNew trivial:"];
  Do[Print["  ", SectorToStr[s]], {s, trivialList}];
];

Export[FileNameJoin[{$CacheDir, "scan10s_result.wdx"}],
  <|"AllTrivial" -> allTrivial, "NonTrivial" -> nonTrivialList,
    "Timestamp" -> DateString[]|>, "WDX"];
Print["\nSaved."];
