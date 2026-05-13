(* compare_tb123_script.wl — Compare TB123 Singular vs MMA baseline (to be run after baseline is ready) *)
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

Print["Loading MMA baseline..."];
baseline = Import[$TB123Path <> "PrepareCheckpoint-TB123-MMA.wdx", "WDX"];
If[Head[baseline["Regions"]] === Join, baseline["Regions"] = baseline["Regions"][[1]]];

mmaSectors = Keys[baseline["Regions"]];
Print["MMA baseline sectors: ", Length[mmaSectors]];

workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];

Print["Running Singular pipeline..."];
totalTime = AbsoluteTiming[
  singularRes = regionsBySectors[ibpeqs, mmaSectors, Alist, vlist,
    Modulus -> char, Verbose -> False];
][[1]];
Print["Singular total time: ", totalTime, " s"];
Print["Singular completed sectors: ", Length[Keys[singularRes]]];

Print["\n=== Comparison ==="];
allMatch = True;
mismatchList = {};

Do[
  sec = mmaSectors[[i]];
  If[!KeyExistsQ[singularRes, sec],
    Print["Sector ", sec, " MISSING in Singular"];
    allMatch = False;
    AppendTo[mismatchList, sec];
    Continue[];
  ];
  
  mmaRegs = baseline["Regions"][sec];
  singRegs = singularRes[sec];
  
  If[Length[mmaRegs] =!= Length[singRegs],
    Print["Sector ", sec, " region count mismatch: ", Length[mmaRegs], " vs ", Length[singRegs]];
    allMatch = False;
    AppendTo[mismatchList, sec];
    Continue[];
  ];
  
  secMatch = True;
  Do[
    mmaRing = mmaRegs[[j]]["CoordinateRing"];
    singRing = singRegs[[j]]["CoordinateRing"];
    
    If[mmaRing["VarDeg"] =!= singRing["VarDeg"],
      Print["Sector ", sec, " Region ", j, " VarDeg mismatch"];
      secMatch = False;
    ];
    If[mmaRing["VarIndep"] =!= singRing["VarIndep"],
      Print["Sector ", sec, " Region ", j, " VarIndep mismatch"];
      secMatch = False;
    ];
    If[Length[mmaRing["MinPoly"]] =!= Length[singRing["MinPoly"]],
      Print["Sector ", sec, " Region ", j, " MinPoly length mismatch"];
      secMatch = False;
    ];
  , {j, Length[mmaRegs]}];
  
  If[!secMatch,
    AppendTo[mismatchList, sec];
    allMatch = False;
  ];
, {i, Length[mmaSectors]}];

If[allMatch,
  Print["ALL TB123 SECTORS MATCH!"];
  ,
  Print["Mismatched sectors: ", mismatchList];
];

Print["Done."];
