(* blade_test_SR212-3m.m — v6: no gExp filter, reduce ALL integrals *)
(* 运行: cd build && wolframscript -file blade_test_SR212-3m.m *)

$ProjectRoot = DirectoryName[$InputFileName] <> "/..";
$BladePath = "/home/ykm/blade-workspace";
Get[FileNameJoin[{$BladePath, "Blade.wl"}]];

Get[FileNameJoin[{$ProjectRoot, "verify/FamilyDatabase/FamilyDatabase.wl"}]];

famname = "SR212-3m";
family = Symbol[StringReplace[famname, "-" -> ""]];

dimension = 4 - 2*eps;
loop = $FamilyDatabase[famname]["LoopMomenta"];
leg = $FamilyDatabase[famname]["ExternalMomenta"];
conservation = {};
replacement = $FamilyDatabase[famname]["KinematicRules"];
propagator = $FamilyDatabase[famname]["Propagators"];
topsector = $FamilyDatabase[famname]["TopSector"];
numeric = DeleteCases[$FamilyDatabase[famname]["Numeric"], "d" -> _];
dValue = "d" /. $FamilyDatabase[famname]["Numeric"];
modulus = $FamilyDatabase[famname]["Modulus"];

Print["Family: ", famname, " (Blade: ", ToString[family], ")"];
Print["TopSector: ", topsector];

BLSetReducerOptions["MaxIncrement" -> 5, "BlackBoxRank" -> 3, "BlackBoxDot" -> 1];
BLFamilyDefine[family, dimension, propagator, loop, leg,
  conservation, replacement, topsector, numeric];

(* Step 1: Top sector reduction *)
Print["\n=== Step 1: Top sector ==="];
topExps = ReplacePart[ConstantArray[0, Length[topsector]], Thread[Position[topsector, 1] -> 1]];
AbsoluteTiming[BLReduce[{BL[family, topExps]}];];
masters = BLFamilyInf[family]["Masters"];
Print["Masters: ", Length[masters]];
Do[Print["  ", m], {m, masters}];

(* Step 2: Load C++ *)
allRelFiles = FileNames["AllRelations_SR212-3m_k*.m", $ProjectRoot];
allRelFile = Last[Sort[allRelFiles]];
Print["\n=== Step 2: ", allRelFile, " ==="];
Get[allRelFile];

Do[
  Print["  (lev=", #["Lev"], ", deg=", #["Deg"], ")",
    " sols=", #["NumSolutions"],
    " |alpha|=", Length[#["Alphas"]],
    " |beta|=", Length[#["Betas"]]] & /@ $AllRelations;
];

(* Step 3: NuVerify *)
Print["\n=== Step 3: NuVerify ==="];
NE = $AllRelations[[1]]["NE"];
nuSample = {1, 2, 1, 2, 1};
Print["NE=", NE, " nu=", nuSample];
epsSol = Solve[dimension - dValue == 0, eps][[1]];

redCache = <||>;
reduceOne[e_] := Module[{t, r}, t = BL[family, e]; r = BLReduce[{t}]; If[Length[r] > 0, r[[1]], 0]];
cachedReduce[e_] := If[!KeyExistsQ[redCache, e], redCache[e] = reduceOne[e]]; redCache[e];

totalPass = 0; totalFail = 0; totalRels = 0;

Do[
  entry = $AllRelations[[i]];
  lev = entry["Lev"]; deg = entry["Deg"];
  sols = entry["NumSolutions"];
  If[sols == 0, Continue[]];

  alphas = entry["Alphas"];
  betas = entry["Betas"];
  coeffMat = entry["Coefficients"];
  nAlpha = Length[alphas];
  nBeta = Length[betas];

  Print["\n(lev=", lev, ", deg=", deg, ") sols=", sols];

  (* 收集所有唯一积分（包括 sub-sector） *)
  allGexp = {};
  Do[
    gExp = Table[nuSample[[j]] - alphas[[ai, j]], {j, 1, NE}];
    AppendTo[allGexp, gExp];
  , {ai, 1, nAlpha}];
  allGexp = DeleteDuplicates[allGexp];
  Print["  Unique integrals: ", Length[allGexp]];

  (* 预先约化所有积分 *)
  Do[cachedReduce[g], {g, allGexp}];
  Print["  Cache: ", Length[redCache], " total"];

  (* 验证 *)
  Do[
    expr = 0; idx = 1;
    Do[
      alpha = alphas[[ai]];
      Do[
        beta = betas[[bi]];
        coeff = coeffMat[[idx, si+1]];
        If[coeff != 0,
          gExp = Table[nuSample[[j]] - alpha[[j]], {j, 1, NE}];
          betaFactor = Product[nuSample[[j]]^beta[[j]], {j, 1, NE}];
          gRed = cachedReduce[gExp];
          expr = expr + coeff * betaFactor * gRed;
        ];
        idx++;
      , {bi, 1, nBeta}];
    , {ai, 1, nAlpha}];

    (* mod p check *)
    expanded = Expand[expr /. epsSol];
    passQ = True;
    Do[
      coeff = Coefficient[expanded, m];
      If[coeff =!= 0 && NumericQ[coeff] && Mod[Numerator[coeff], modulus] =!= 0,
        passQ = False;
      ];
    , {m, masters}];

    totalRels++;
    If[passQ,
      Print["  Sol ", si+1, ": PASS"];
      totalPass++,
      Print["  Sol ", si+1, ": FAIL"];
      totalFail++
    ];

  , {si, 0, sols-1}];

, {i, Length[$AllRelations]}];

Print["\n========================================"];
Print["  RESULT: ", totalPass, "/", totalRels, " PASSED"];
If[totalFail > 0, Print["  FAILED: ", totalFail]];
Print["  Cache: ", Length[redCache], " integrals"];
Print["========================================"];
