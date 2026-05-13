(* run_111110100.wl -- Single sector Long tier run for NP322 111110100 *)
(* Usage: wolframscript -file verify/NP322/run_111110100.wl *)

$ProjectRoot = "/home/ykm/Large-Index-Expansion-CPP";
$VerifyUtility = FileNameJoin[{$ProjectRoot, "verify", "VerifyUtility"}];

SetDirectory[$VerifyUtility];
Get["LIEWorkflow.wl"];

Needs["LIE`"];
Needs["SingularCoordinateRing`"];

Print["========================================"];
Print["NP322 sector 111110100 - Long tier run"];
Print["Timeout: 96400s (26.8h), single thread"];
Print["========================================"];

(* Load family *)
Get[FileNameJoin[{$ProjectRoot, "verify", "FamilyDatabase", "FamilyDatabase.wl"}]];
config = $FamilyDatabase["NP322"];
data = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"],
  config["ExternalMomenta"], config["KinematicRules"],
  config["TopSector"],
  "Numeric" -> config["Numeric"],
  Modulus -> config["Modulus"], Verbose -> False
];
Print["Family loaded. NE = ", data["Family", "NE"]];

(* Get IBP data for sector 111110100 *)
sectorList = data["Family", "SectorList"];
sector = {1,1,1,1,1,0,1,0,0};
If[!MemberQ[sectorList, sector],
  Print["ERROR: sector not in sector list"]; Exit[1];
];

ibpeqs = data["IBP", "Sectors", sector];
Alist = data["Family", "AList"];
ibpRelations = data["IBP", "Relations"];
vlist = data["Family", "VList"];
char = config["Modulus"];
ne = data["Family", "NE"];

Print["Running SectorWorker for 111110100..."];
Print["Timeout: 96400s (will be killed if exceeded)"];
Print["Start time: ", DateString[]];

cacheDir = FileNameJoin[{$ProjectRoot, "verify", "NP322", "cache"}];

result = TimeConstrained[
  regionsBySectors[{sector}, ibpeqs, Alist, vlist,
    Modulus -> char,
    "Timeout" -> 96400,
    Verbose -> True],
  96400 + 300,
  $Failed
];

If[result === $Failed || !AssociationQ[result],
  Print["\nFAILED: sector 111110100 did not complete in 96400s"];
  Exit[1];
];

If[KeyExistsQ[result, sector],
  regions = result[sector];
  nRegions = Length[regions];
  nbTotal = Total[Length[#["MonomialBasisIndex"]] & /@ regions];
  Print["\nSUCCESS: ", nRegions, " regions, total nb = ", nbTotal];
  (* Save to SectorCache *)
  Export[
    FileNameJoin[{cacheDir, "SectorCache", "sector_111110100.wdx"}],
    regions, "WDX"
  ];
  Print["Saved to SectorCache/sector_111110100.wdx"];
  (* Save status update *)
  statusFile = FileNameJoin[{cacheDir, "status.wdx"}];
  If[FileExistsQ[statusFile],
    status = Import[statusFile, "WDX"];
    status["111110100"] = <|"State" -> "done", "Regions" -> nRegions, "Time" -> 0|>;
    Export[statusFile, status, "WDX"];
    Export[FileNameJoin[{cacheDir, "status.json"}], status, "JSON"];
    Print["Status updated: 111110100 -> done"];
  ];
,
  Print["\nTrivial: no regions found for 111110100"];
  Export[
    FileNameJoin[{cacheDir, "SectorCache", "sector_111110100.wdx"}],
    {}, "WDX"
  ];
];

Print["Finish time: ", DateString[]];
Print["========================================"];
