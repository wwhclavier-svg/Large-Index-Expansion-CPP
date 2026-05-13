(* test_box_singular_fix.wl — Validate SingularCoordinateRing Box fix *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["Box"];

(* Define family to generate ibpeqs, Alist, vlist *)
workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sector = {1, 1, 1, 1};

Print["=== Box sector ", sector, " via regionsBySectors ==="];
Print["Alist = ", Alist];
Print["vlist = ", vlist];
Print["IBPEqs length = ", Length[ibpeqs]];

time = AbsoluteTiming[
  res = regionsBySectors[ibpeqs, {sector}, Alist, vlist,
    Modulus -> char, Verbose -> True];
][[1]];
Print["Time: ", time];
Print["Return type: ", Head[res]];
Print["Keys: ", If[Head[res] === Association, Keys[res], "N/A"]];

(* regionsBySectors returns Association with sector keys *)
regions = If[Head[res] === Association, res[sector], res];
Print["Regions found for sector: ", Length[regions]];

If[Length[regions] > 0 && Head[regions] === List,
  Do[
    regionWrap = regions[[i]];
    ring = regionWrap["CoordinateRing"];
    Print["Region ", i, ":"];
    Print["  VarDeg=", ring["VarDeg"]];
    Print["  VarIndep=", ring["VarIndep"]];
    Print["  VarRule length=", Length[ring["VarRule"]]];
    If[Length[ring["VarRule"]] > 0,
      Print["  VarRule[[1]]=", ring["VarRule"][[1]]];
    ];
    Print["  MinPoly length=", Length[ring["MinPoly"]]];
    
    (* Check for corrupted Part expressions *)
    vf = FullForm[ring["VarRule"]];
    If[!FreeQ[vf, Part],
      Print["  ERROR: Part expression detected in VarRule FullForm!"];
    ];
    
    (* Check VarDeg length matches Avar count *)
    If[Length[ring["VarDeg"]] =!= Length[Alist],
      Print["  ERROR: VarDeg length ", Length[ring["VarDeg"]], 
        " != Avar count ", Length[Alist]];
    ];
  , {i, Length[regions]}];
  ,
  Print["No regions found or unexpected format!"];
];

Print["=== DONE ==="];
