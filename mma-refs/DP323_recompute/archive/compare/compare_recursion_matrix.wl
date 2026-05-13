(* compare_recursion_matrix.wl — Compare RecursionMatrix logic content *)
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

Print["=== Comparing RecursionMatrix ==="];
Do[
  mmaRec = mmaRegions[[i]]["RecursionMatrix"];
  singRec = singularRegions[[i]]["RecursionMatrix"];
  
  Print["Region ", i, ":"];
  
  opMatch = True;
  Do[
    op = opList[[j]];
    If[!KeyExistsQ[mmaRec, op] || !KeyExistsQ[singRec, op],
      If[KeyExistsQ[mmaRec, op] =!= KeyExistsQ[singRec, op],
        Print["  ", op, " key mismatch"];
        opMatch = False;
      ];
      Continue[];
    ];
    
    mmaVal = mmaRec[op];
    singVal = singRec[op];
    
    If[Head[mmaVal] === SparseArray && Head[singVal] === SparseArray,
      mmaNormal = Normal[mmaVal];
      singNormal = Normal[singVal];
      If[mmaNormal =!= singNormal,
        (* Check if diff is just mod char equivalent *)
        diff = mmaNormal - singNormal;
        diffMod = If[char =!= 0, PolynomialMod[diff, char], diff];
        If[diffMod === 0*diffMod || diffMod === ConstantArray[0, Dimensions[diffMod]],
          Print["  ", op, " mod-equivalent (display diff only)"]
          ,
          Print["  ", op, " MISMATCH"];
          Print["    MMA dims: ", Dimensions[mmaNormal]];
          Print["    Sing dims: ", Dimensions[singNormal]];
          opMatch = False;
        ]
        ,
        Print["  ", op, " MATCH"];
      ]
      ,
      If[mmaVal =!= singVal,
        Print["  ", op, " MISMATCH: MMA=", mmaVal, " Sing=", singVal];
        opMatch = False;
        ,
        Print["  ", op, " MATCH"];
      ]
    ]
  , {op, {"M1", "N1", "K1", "K2", "K1s", "K2s", "F2", "F2s", "F0", "nibp", "ne", "nb", "incre"}}];
  
  If[opMatch, Print["  ALL OPS MATCH"]];
, {i, Min[Length[mmaRegions], Length[singularRegions]]}];

Print["Done."];
