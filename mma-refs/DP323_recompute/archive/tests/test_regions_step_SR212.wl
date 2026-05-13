(* test_regions_step_SR212.wl — Step through regionsBySectors logic *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

Print["Loading SR212 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
sectorlist = {{0, 1, 0, 1, 1}}; (* Just one sector *)

ne = Length @ Alist;

(* Manually replicate regionsBySectors logic *)
Print["\nReplicating regionsBySectors logic..."];
result = Table[
  sec = sectorlist[[i]];
  Print["Sector: ", sec];
  
  ibpeqssub = LIEUtility`sectorLimitIBP[ibpeqs, sec, vlist];
  Print["  sectorLimitIBP done"];
  
  expregsub = LIECoreAlgebra`expRegSolve2[ibpeqssub, Alist, vlist, "LimitSector" -> sec, Modulus -> char, Verbose -> False];
  Print["  expRegSolve2 done, regions: ", Length[expregsub]];
  
  expregsub = Select[expregsub, Join[#["VarDep"], #["VarIndep"]] =!= {} &];
  Print["  after Select: ", Length[expregsub]];
  
  If[expregsub =!= {},
    nRegs = Length[expregsub];
    buildResult = Table[
      LIERegions`Private`buildRecursionMatrix[ibpeqs, reg, ne, "LimitSector" -> sec, Modulus -> char],
      {reg, expregsub}
    ];
    Print["  buildRecursionMatrix done"];
    sec -> buildResult
    ,
    Nothing
  ]
, {i, Length @ sectorlist}];

Print["\nResult type: ", Head[result]];
Print["Result length: ", Length[result]];

assoc = Association[result];
Print["Association type: ", Head[assoc]];
Print["Association keys: ", Keys[assoc]];
