(* compare_box_mma_singular.wl — Compare Box Singular vs MMA baseline *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";
$BoxPath = $LIECPPPath <> "/verify/Box/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["Box"];

workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sector = {1, 1, 1, 1};

(* Load MMA baseline *)
Print["Loading MMA baseline..."];
 baseline = Import[$BoxPath <> "PrepareCheckpoint-Box-MMA.wdx", "WDX"];
If[Head[baseline["Regions"]] === Join, baseline["Regions"] = baseline["Regions"][[1]]];

mmaRegions = baseline["Regions"][sector];
Print["MMA regions for sector ", sector, ": ", Length[mmaRegions]];

Do[
  Print["MMA Region ", i, ": VarDeg=", mmaRegions[[i]]["CoordinateRing"]["VarDeg"],
    ", VarIndep=", mmaRegions[[i]]["CoordinateRing"]["VarIndep"],
    ", MinPoly length=", Length[mmaRegions[[i]]["CoordinateRing"]["MinPoly"]]];
, {i, Length[mmaRegions]}];

(* Run Singular *)
Print["\nRunning Singular..."];
singularRes = regionsBySectors[ibpeqs, {sector}, Alist, vlist,
  Modulus -> char, Verbose -> False];
singularRegions = singularRes[sector];
Print["Singular regions: ", Length[singularRegions]];

Do[
  Print["Singular Region ", i, ": VarDeg=", singularRegions[[i]]["CoordinateRing"]["VarDeg"],
    ", VarIndep=", singularRegions[[i]]["CoordinateRing"]["VarIndep"],
    ", MinPoly length=", Length[singularRegions[[i]]["CoordinateRing"]["MinPoly"]]];
, {i, Length[singularRegions]}];

(* Compare *)
Print["\n=== Comparison ==="];
If[Length[mmaRegions] =!= Length[singularRegions],
  Print["MISMATCH: Region count ", Length[mmaRegions], " vs ", Length[singularRegions]];
];

minLen = Min[Length[mmaRegions], Length[singularRegions]];
Do[
  mmaRing = mmaRegions[[i]]["CoordinateRing"];
  singRing = singularRegions[[i]]["CoordinateRing"];
  
  match = True;
  If[mmaRing["VarDeg"] =!= singRing["VarDeg"],
    Print["Region ", i, " VarDeg MISMATCH: MMA=", mmaRing["VarDeg"], " Sing=", singRing["VarDeg"]];
    match = False;
  ];
  If[mmaRing["VarIndep"] =!= singRing["VarIndep"],
    Print["Region ", i, " VarIndep MISMATCH: MMA=", mmaRing["VarIndep"], " Sing=", singRing["VarIndep"]];
    match = False;
  ];
  If[Length[mmaRing["MinPoly"]] =!= Length[singRing["MinPoly"]],
    Print["Region ", i, " MinPoly length MISMATCH: MMA=", Length[mmaRing["MinPoly"]], " Sing=", Length[singRing["MinPoly"]]];
    match = False;
  ];
  If[match, Print["Region ", i, " MATCH"]];
, {i, minLen}];

Print["=== DONE ==="];
