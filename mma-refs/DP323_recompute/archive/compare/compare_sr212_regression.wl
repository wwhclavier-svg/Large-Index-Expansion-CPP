(* compare_sr212_regression.wl — SR212 regression test after SingularCoordinateRing fixes *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";
$SR212Path = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["SR212"];

(* Load MMA baseline *)
Print["Loading SR212 MMA baseline..."];
baseline = Import[$SR212Path <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[baseline["Regions"]] === Join, baseline["Regions"] = baseline["Regions"][[1]]];

workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sectorList = workflow["Family", "SectorList"];

Print["Running Singular pipeline..."];
totalTime = AbsoluteTiming[
  singularRes = regionsBySectors[ibpeqs, sectorList, Alist, vlist,
    Modulus -> char, Verbose -> False];
][[1]];
Print["Singular time: ", totalTime, " s"];

Print["\n=== Comparison ==="];
allMatch = True;
Do[
  sec = sectorList[[i]];
  If[!KeyExistsQ[baseline["Regions"], sec],
    Print["Sector ", sec, " SKIP (no MMA baseline)"];
    Continue[];
  ];
  If[!KeyExistsQ[singularRes, sec],
    Print["Sector ", sec, " MISSING"];
    allMatch = False;
    Continue[];
  ];
  
  mmaRegs = baseline["Regions"][sec];
  singRegs = singularRes[sec];
  
  If[Length[mmaRegs] =!= Length[singRegs],
    Print["Sector ", sec, " region count mismatch"];
    allMatch = False;
    Continue[];
  ];
  
  secMatch = True;
  Do[
    mmaRing = mmaRegs[[j]]["CoordinateRing"];
    singRing = singRegs[[j]]["CoordinateRing"];
    If[mmaRing["VarDeg"] =!= singRing["VarDeg"], secMatch = False];
    If[mmaRing["VarIndep"] =!= singRing["VarIndep"], secMatch = False];
    If[Length[mmaRing["MinPoly"]] =!= Length[singRing["MinPoly"]], secMatch = False];
  , {j, Length[mmaRegs]}];
  
  If[secMatch,
    Print["Sector ", sec, " MATCH (", Length[mmaRegs], " regions)"]
    ,
    Print["Sector ", sec, " MISMATCH"];
    allMatch = False;
  ];
, {i, Length[sectorList]}];

If[allMatch, Print["\nALL SR212 SECTORS MATCH!"]];
Print["Done."];
