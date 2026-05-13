(* compare_matrix_values.wl — Compare MonomialBasisMatrix values *)
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
baseline = Import[$BoxPath <> "PrepareCheckpoint-Box-MMA.wdx", "WDX"];
If[Head[baseline["Regions"]] === Join, baseline["Regions"] = baseline["Regions"][[1]]];

mmaRegions = baseline["Regions"][sector];

(* Run Singular *)
singularRes = regionsBySectors[ibpeqs, {sector}, Alist, vlist,
  Modulus -> char, Verbose -> False];
singularRegions = singularRes[sector];

Print["=== Comparing MonomialBasisMatrix ==="];
Do[
  mmaRing = mmaRegions[[i]]["CoordinateRing"];
  singRing = singularRegions[[i]]["CoordinateRing"];
  
  mmaMBM = mmaRing["MonomialBasisMatrix"];
  singMBM = singRing["MonomialBasisMatrix"];
  
  Print["Region ", i, ":"];
  Print["  MMA keys: ", Keys[mmaMBM]];
  Print["  Sing keys: ", Keys[singMBM]];
  
  If[Keys[mmaMBM] =!= Keys[singMBM],
    Print["  KEY MISMATCH!"];
    Continue[];
  ];
  
  keyMatch = True;
  Do[
    mmaMat = mmaMBM[key];
    singMat = singMBM[key];
    If[mmaMat =!= singMat,
      Print["  Key ", key, " MISMATCH"];
      Print["    MMA: ", Normal[mmaMat]];
      Print["    Sing: ", Normal[singMat]];
      keyMatch = False;
    ];
  , {key, Keys[mmaMBM]}];
  
  If[keyMatch, Print["  ALL KEYS MATCH"]];
, {i, Length[mmaRegions]}];

Print["=== Comparing FractionRule ==="];
Do[
  mmaFR = mmaRegions[[i]]["CoordinateRing"]["FractionRule"];
  singFR = singularRegions[[i]]["CoordinateRing"]["FractionRule"];
  
  If[mmaFR =!= singFR,
    Print["Region ", i, " FractionRule MISMATCH"];
    diffKeys = Select[Keys[mmaFR], mmaFR[#] =!= singFR[#] &];
    Print["  Diff keys: ", diffKeys];
    Do[
      Print["  ", k, ": MMA=", mmaFR[k], " Sing=", singFR[k]];
    , {k, diffKeys}];
  , 
    Print["Region ", i, " FractionRule MATCH"];
  ];
, {i, Length[mmaRegions]}];

Print["Done."];
