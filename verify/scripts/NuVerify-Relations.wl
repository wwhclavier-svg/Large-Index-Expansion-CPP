(* ν 采样验证 C++ 关系 *)
Get["/root/Large-Index-Expansion-MMA-Mini/KiraRuleLoader.wl"];

modulus = 179424673;
dVal = 1/3;

{allRulesSymbolic, kiraReduce} = LoadKiraRules[
  "/root/kira_tests/bub00_new/bub00",
  "d" -> dVal, "Modulus" -> modulus, "NProp" -> 2
];

(* g[ν-α] → j[ν-α]: substitute ν then g→j *)
gToJ[expr_, nu1_, nu2_] := Module[{withNu},
  withNu = expr /. {v1 -> nu1, v2 -> nu2};
  withNu /. g[a_, b_] :> j[a, b]
];

verifyRelation[relExpr_, nu1_, nu2_] := Module[{jForm, reduced},
  jForm = gToJ[relExpr, nu1, nu2];
  reduced = kiraReduce[jForm];
  reduced === 0
];

(* 加载关系 *)
Get["verify/results/bub00/Compare-CPPRelation-bub00_lev2_deg0.m"]; rel2d0 = $RelationResult;
Get["verify/results/bub00/Compare-CPPRelation-bub00_lev2_deg1.m"]; rel2d1 = $RelationResult;
Get["verify/results/bub00/Compare-CPPRelation-bub00_lev2_deg2.m"]; rel2d2 = $RelationResult;

testNuPoints = {
  {1,0},{0,1},{1,1},{1,2},{2,1},{2,0},{0,2},
  {3,3},{3,4},{4,3},{5,5},{10,10},
  {100,100},{1000,2000},{99999,77777}
};

Print["=== C++ Relation ν-Verification ==="];
Print["Modulus: ", modulus, ", d: ", dVal];
Print["Test ν: ", testNuPoints];
Print[""];

verifyConfig[label_, data_] := Module[
  {relsAssoc, nSols, passCount, failCount, sol, nu, relExpr, relStr, isPass},
  (* Relations 是 rule 列表，转为 Association *)
  relsAssoc = Association[data[["Relations"]]];
  nSols = data[["NumSolutions"]];

  Print["--- ", label, " (lev=", data[["Lev"]], ", deg=", data[["Deg"]],
    ", ", nSols, " basis) ---"];

  passCount = 0; failCount = 0;

  Do[
    (* 获取关系表达式字符串，转为表达式 *)
    relStr = relsAssoc[StringJoin["Basis", ToString[sol]]];
    relExpr = ToExpression[relStr, InputForm];

    Do[
      isPass = verifyRelation[relExpr, nu[[1]], nu[[2]]];
      If[isPass, passCount++, failCount++];
      If[!isPass, Print["  FAIL: Basis", sol, " at nu=", nu]],
      {nu, testNuPoints}
    ],
    {sol, 1, nSols}
  ];

  Print["  Passed: ", passCount, "/", nSols*Length[testNuPoints],
    If[failCount > 0, " (" <> ToString[failCount] <> " FAILED)", ""]];
  Print[""];
  {passCount, failCount}
];

r0 = verifyConfig["Config A", rel2d0];
r1 = verifyConfig["Config B", rel2d1];
r2 = verifyConfig["Config C", rel2d2];

Print["=== Overall Summary ==="];
Print["Config A (lev=2,deg=0): ", r0[[1]], "/", r0[[1]]+r0[[2]]];
Print["Config B (lev=2,deg=1): ", r1[[1]], "/", r1[[1]]+r1[[2]]];
Print["Config C (lev=2,deg=2): ", r2[[1]], "/", r2[[1]]+r2[[2]]];
totalPass = r0[[1]] + r1[[1]] + r2[[1]];
totalAll = r0[[1]]+r0[[2]] + r1[[1]]+r1[[2]] + r2[[1]]+r2[[2]];
Print["Grand Total: ", totalPass, "/", totalAll, " passed"];
