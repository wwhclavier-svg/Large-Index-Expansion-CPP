(* compute_sectors_37_38.wl — 单独计算 sector 37 和 38 *)
(* Usage: wolframscript -file compute_sectors_37_38.wl *)

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

sectors = {
  {0,1,0,1,1,0,1,1,0,0,0},  (* sector 37, 5 props *)
  {0,1,0,1,1,1,1,1,0,0,0}   (* sector 38, 6 props *)
};

Do[
  sec = sectors[[i]];
  Print["\n========================================"];
  Print["Computing sector ", i + 36, ": ", sec];
  Print["Active props: ", Total[sec]];
  Print["========================================"];

  {elapsed, result} = AbsoluteTiming[
    LIERegions`regionsBySectors[ibpeqs, {sec}, Alist, vlist,
      Modulus -> char,
      "EnableFieldExtension" -> True,
      Verbose -> True,
      "TimeoutMarker" -> $Timeout
    ]
  ];

  Print["Elapsed: ", Round[elapsed, 0.1], "s"];

  (* Robust result parsing *)
  If[!AssociationQ[result],
    Print["RESULT: FAILED (non-Association result), head=", Head[result]];
    Continue[];
  ];

  If[KeyExistsQ[result, sec] && result[sec] === $Timeout,
    Print["RESULT: TIMEOUT"];
    Continue[];
  ];

  If[KeyExistsQ[result, sec],
    regions = result[sec];
    nRegions = Length[regions];
    Print["Regions found: ", nRegions];

    If[nRegions > 0,
      fname = FileNameJoin[{$CacheDir, "sector_" <> SectorToStr[sec] <> ".wdx"}];
      Export[fname, <|"Sector" -> sec, "Regions" -> regions, "Timestamp" -> DateString[]|>, "WDX"];
      Print["Saved to: ", fname];
      Do[
        coordRing = regions[[j]]["CoordinateRing"];
        Print["  Region ", j, ": ", coordRing["VarDep"], " vars, deg=", coordRing["VarDeg"]];
      , {j, nRegions}];
    ,
      Print["0 regions (trivial)"];
    ];
    Print["RESULT: ", nRegions, " regions, ", Round[elapsed, 0.1], "s"];
    ,
    Print["0 regions (trivial or empty association)"];
    Print["RESULT: 0 regions, ", Round[elapsed, 0.1], "s"];
  ];
, {i, 1, 2}];

Print["\nDone."];
