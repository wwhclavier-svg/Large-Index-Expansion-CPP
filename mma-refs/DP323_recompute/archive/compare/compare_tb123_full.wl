(* compare_tb123_full.wl — Full TB123 A/B comparison (Singular vs MMA) *)
(* Requires: verify/TB123/PrepareCheckpoint-TB123-MMA.wdx (MMA baseline) *)

$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";
$TB123Path = $LIECPPPath <> "/verify/TB123/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["TB123"];

workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sectorList = workflow["Family", "SectorList"];

(* Load MMA baseline *)
Print["Loading MMA baseline..."];
baselineFile = $TB123Path <> "PrepareCheckpoint-TB123-MMA.wdx";
If[FileExistsQ[baselineFile],
  baseline = Import[baselineFile, "WDX"];
  If[Head[baseline["Regions"]] === Join, baseline["Regions"] = baseline["Regions"][[1]]];
  Print["MMA completed sectors: ", Length[Keys[baseline["Regions"]]]];
  ,
  Print["ERROR: MMA baseline not found at ", baselineFile];
  Print["Run generate_tb123_mma_baseline.wl first to create MMA baseline."];
  Exit[1];
];

Print["Total sectors: ", Length[sectorList]];

(* Run all sectors via Singular *)
Print["\nRunning Singular for all sectors..."];
totalTime = AbsoluteTiming[
  singularRes = regionsBySectors[ibpeqs, sectorList, Alist, vlist,
    Modulus -> char, Verbose -> False];
][[1]];
Print["Total Singular time: ", totalTime, " s"];
Print["Singular completed sectors: ", Length[Keys[singularRes]]];

(* Compare each sector *)
Print["\n=== Sector-by-sector comparison ==="];
allMatch = True;
mismatchSectors = {};

Do[
  sec = sectorList[[i]];
  If[!KeyExistsQ[baseline["Regions"], sec],
    Print["Sector ", sec, " SKIP (no MMA baseline)"];
    Continue[];
  ];
  If[!KeyExistsQ[singularRes, sec],
    Print["Sector ", sec, " MISMATCH (Singular missing)"];
    allMatch = False;
    AppendTo[mismatchSectors, sec];
    Continue[];
  ];
  
  mmaRegs = baseline["Regions"][sec];
  singRegs = singularRes[sec];
  
  If[Length[mmaRegs] =!= Length[singRegs],
    Print["Sector ", sec, " MISMATCH: region count ", Length[mmaRegs], " vs ", Length[singRegs]];
    allMatch = False;
    AppendTo[mismatchSectors, sec];
    Continue[];
  ];
  
  secMatch = True;
  Do[
    mmaRing = mmaRegs[[j]]["CoordinateRing"];
    singRing = singRegs[[j]]["CoordinateRing"];
    
    If[mmaRing["VarDeg"] =!= singRing["VarDeg"],
      Print["Sector ", sec, " Region ", j, " VarDeg MISMATCH: MMA=", mmaRing["VarDeg"], " Sing=", singRing["VarDeg"]];
      secMatch = False;
    ];
    If[mmaRing["VarIndep"] =!= singRing["VarIndep"],
      Print["Sector ", sec, " Region ", j, " VarIndep MISMATCH: MMA=", mmaRing["VarIndep"], " Sing=", singRing["VarIndep"]];
      secMatch = False;
    ];
    If[Length[mmaRing["MinPoly"]] =!= Length[singRing["MinPoly"]],
      Print["Sector ", sec, " Region ", j, " MinPoly length MISMATCH"];
      secMatch = False;
    ];
  , {j, Length[mmaRegs]}];
  
  If[secMatch,
    Print["Sector ", sec, " MATCH (", Length[mmaRegs], " regions)"];
    ,
    allMatch = False;
    AppendTo[mismatchSectors, sec];
  ];
, {i, Length[sectorList]}];

Print["\n=== SUMMARY ==="];
If[allMatch,
  Print["ALL SECTORS MATCH!"];
  ,
  Print["MISMATCH in sectors: ", mismatchSectors];
];
Print["Done."];