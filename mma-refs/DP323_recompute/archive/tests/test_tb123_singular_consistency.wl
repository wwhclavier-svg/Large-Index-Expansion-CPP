(* test_tb123_singular_consistency.wl — Verify TB123 Singular pipeline consistency *)
(* Runs same sectors twice via Singular and compares results *)

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
config = $FamilyDatabase["TB123"];

Print["=== Defining TB123 family ==="];
workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sectorList = workflow["Family", "SectorList"];

(* Pick 5 diverse sectors covering different complexity levels *)
testSectors = {
  {1, 0, 0, 0, 0, 0, 0},   (* simple: top sector with only l1 *)
  {0, 1, 0, 0, 0, 0, 0},   (* simple: top sector with only l2 *)
  {1, 1, 0, 0, 0, 0, 0},   (* medium: both loops active *)
  {1, 1, 1, 0, 0, 0, 0},   (* medium-hard: 3 loops *)
  {1, 1, 1, 1, 0, 0, 0}    (* harder: 4 loops *)
};

Print["Test sectors: ", testSectors];
Print["Total sectors in family: ", Length[sectorList]];

(* Filter to only those that exist in sectorList *)
testSectors = Select[testSectors, MemberQ[sectorList, #] &];
Print["Filtered test sectors: ", testSectors];

(* First run *)
Print["\n=== First Singular run ==="];
totalTime1 = AbsoluteTiming[
  run1 = regionsBySectors[ibpeqs, testSectors, Alist, vlist,
    Modulus -> char, Verbose -> True];
][[1]];
Print["Run 1 time: ", totalTime1, " s"];

(* Second run *)
Print["\n=== Second Singular run ==="];
totalTime2 = AbsoluteTiming[
  run2 = regionsBySectors[ibpeqs, testSectors, Alist, vlist,
    Modulus -> char, Verbose -> True];
][[1]];
Print["Run 2 time: ", totalTime2, " s"];

(* Compare *)
Print["\n=== Consistency check ==="];
allMatch = True;
Do[
  sec = testSectors[[i]];
  Print["Sector ", sec, ":"];
  
  If[!KeyExistsQ[run1, sec],
    Print["  RUN1: MISSING"];
    allMatch = False;
    Continue[];
  ];
  If[!KeyExistsQ[run2, sec],
    Print["  RUN2: MISSING"];
    allMatch = False;
    Continue[];
  ];
  
  regs1 = run1[sec];
  regs2 = run2[sec];
  
  If[Length[regs1] =!= Length[regs2],
    Print["  REGION COUNT: ", Length[regs1], " vs ", Length[regs2], " MISMATCH"];
    allMatch = False;
    Continue[];
  ];
  
  Print["  Region count: ", Length[regs1], " — OK"];
  
  Do[
    ring1 = regs1[[j]]["CoordinateRing"];
    ring2 = regs2[[j]]["CoordinateRing"];
    
    checks = {
      {"VarDeg", ring1["VarDeg"], ring2["VarDeg"]},
      {"VarIndep", ring1["VarIndep"], ring2["VarIndep"]},
      {"MinPoly", ring1["MinPoly"], ring2["MinPoly"]},
      {"MonomialBasisIndex", Keys[ring1["MonomialBasisMatrix"]], Keys[ring2["MonomialBasisMatrix"]]}
    };
    
    allFieldMatch = True;
    Do[
      field = checks[[k]];
      If[field[[2]] =!= field[[3]],
        Print["    Region ", j, " ", field[[1]], " MISMATCH"];
        Print["      Run1: ", field[[2]]];
        Print["      Run2: ", field[[3]]];
        allFieldMatch = False;
        allMatch = False;
      ];
    , {k, Length[checks]}];
    
    If[allFieldMatch,
      Print["    Region ", j, ": ALL FIELDS MATCH"];
    ];
  , {j, Length[regs1]}];
, {i, Length[testSectors]}];

Print["\n=== FINAL RESULT ==="];
If[allMatch,
  Print["ALL CHECKS PASSED — TB123 Singular pipeline is consistent"];
  ,
  Print["MISMATCHES FOUND — see above for details"]
];