(* Evaluate the exact Table from regionsBySectors *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
char = data["Config", "Modulus"];

sectorlist = {{0, 1, 0, 1, 1}};

(* Replicate the exact Table body from regionsBySectors *)
ne = Length[Alist];
ibpeqs_n0 = ibpeqs /. n -> 0;

(* Evaluate the Table directly *)
Print["Evaluating Table directly..."];
tableResult = Table[
  sec = sectorlist[[i]];
  Print["\n      +++  sector:  ", sec, "  (", i, "/", Length[sectorlist], ")  +++ "];
  secStart = AbsoluteTime[];
  ibpeqssub = LIEUtility`sectorLimitIBP[ibpeqs_n0, sec, vlist];
  {regTime, expregsub} = AbsoluteTiming @ TimeConstrained[
    LIECoreAlgebra`expRegSolve2[ibpeqssub, Alist, vlist, "LimitSector" -> sec, Modulus -> char, Verbose -> False],
    120,
    $Failed
  ];
  If[expregsub === $Failed,
    Print["  *** TIMEOUT after 120s, skipping sector ***"];
    Nothing
    ,
    expregsub = Select[expregsub, Join[#["VarDep"], #["VarIndep"]] =!= {} &];
    If[expregsub =!= {},
      nRegs = Length[expregsub];
      {buildTime, result} = AbsoluteTiming @ Table[
        LIERegions`Private`buildRecursionMatrix[ibpeqs, reg, ne, "LimitSector" -> sec, Modulus -> char],
        {reg, expregsub}
      ];
      Print["  >> ", nRegs, " region(s): expReg=", Round[regTime, 0.001], "s  build=", Round[buildTime, 0.001], "s"];
      sec -> result
      ,
      Print["  >> 0 regions (all trivial): expReg=", Round[regTime, 0.001], "s"];
      Nothing
    ]
  ];
  secElapsed = AbsoluteTime[] - secStart;
  Print["  [sector total: ", Round[secElapsed, 0.001], "s]"];
, {i, Length[sectorlist]}];

Print["\nTable result type: ", Head[tableResult]];
Print["Table result length: ", Length[tableResult]];
Print["Table result: ", Short[tableResult, 200]];

Print["\nCreating Association..."];
assoc = Association[tableResult];
Print["Association type: ", Head[assoc]];
Print["AssociationQ: ", AssociationQ[assoc]];
Print["Keys: ", Keys[assoc]];
