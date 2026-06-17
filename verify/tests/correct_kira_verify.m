(* correct_kira_verify.m *)
(* 正确的 Kira 验证流程（修正版）
 * 
 * 验证步骤：
 * 1. 取 LIE 关系中的符号 ν
 * 2. 在测试点 ν = (v1_test, v2_test) 处代入
 * 3. g[a__] -> j[a] 转换
 * 4. 用 Kira 规则替换所有 j[...]-化简为 0 或 master 的线性组合
 * 5. 检验结果是否为 0（即所有 master 系数是否抵消）
 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];

(* 加载 Kira 规则 *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];

(* 使用 KiraRuleLoader 上下文的 j 符号 *)
kiraJ = KiraRuleLoader`j;

Print["========================================"];
Print["正确的 Kira 验证（修正版）"];
Print["========================================"];

(* 读取 C++ 关系 *)
Get["/root/Large-Index-Expansion-CPP/verify/bub00/Relations_bub00_lev2_deg1.m"];
rel = $RelationResult;
relations = rel["Relations"];

Print["\nC++ Relations (lev=2, deg=1):"];
Print["  NumSolutions: ", rel["NumSolutions"]];

(* 提取 g 项的辅助函数 *)
extractGTerms[expr_] := Module[{terms, result},
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    result = {};
    Do[
        term = terms[[i]];
        termList = If[Head[term] === Times, List @@ term, {term}];
        gTerms = Cases[termList, g[__]];
        coef = Times @@ DeleteCases[termList, Alternatives @@ gTerms];
        Do[
            gTerm = gTerms[[j]];
            (* 提取 g[...] 的参数 *)
            args = List @@ gTerm;
            (* 转换: g[ν-α] -> j[ν-α] *)
            (* 即 g[a,b] -> j[a,b] *)
            jExpr = Apply[kiraJ, args];
            AppendTo[result, {coef, args, jExpr}]
        , {j, Length[gTerms]}]
    , {i, Length[terms]}];
    result
];

(* 验证函数 *)
verifyAtNu[relExpr_, nu1_, nu2_] := Module[
    {gTerms, substituted, jForm, reduced, masterCoefs, result},
    
    (* Step 1: 提取 g 项 *)
    gTerms = extractGTerms[relExpr];
    
    (* Step 2: 代入 ν 值 *)
    (* g[ν1, ν2-1] -> g[nu1, nu2-1] -> j[nu1, nu2-1] *)
    substituted = Map[
        Function[{term},
            {term[[1]] (* coefficient *), 
             Apply[kiraJ, term[[2]] /. {v1 -> nu1, v2 -> nu2}] (* j[nu1, nu2-nu2] = j[nu1, 0] *)
            }
        ],
        gTerms
    ];
    
    (* Step 3: 构建 j 表达式 *)
    jForm = Plus @@ (#1 * #2 & @@@ substituted);
    
    (* Step 4: Kira 约化 *)
    reduced = kiraReduce[jForm];
    
    (* Step 5: 检查是否为 0 *)
    result = (reduced === 0);
    
    {result, jForm, reduced}
];

(* 测试点 *)
testPoints = {{1, 1}, {2, 3}, {5, 7}};

Print["\n========================================"];
Print["验证非平凡关系"];
Print["========================================"];

nontrivial = Select[relations, # =!= 0 &];
Print["非平凡关系数: ", Length[nontrivial]];

totalPass = 0; totalFail = 0;

Do[
    relExpr = nontrivial[[i]];
    Print["\n关系 ", i, ": ", InputForm[relExpr]//Short];
    
    passCount = 0; failCount = 0;
    Do[
        {nu1, nu2} = testPoints[[tp]];
        {result, jForm, reduced} = verifyAtNu[relExpr, nu1, nu2];
        
        status = If[result, "PASS", "FAIL"];
        If[result, passCount++, failCount++];
        
        Print["  nu=(", nu1, ",", nu2, "): ", status];
        If[!result,
            Print["    jForm: ", InputForm[jForm]//Short];
            Print["    reduced: ", InputForm[reduced]//Short];
            (* 解释为什么不是 0 *)
            If[reduced =!= 0,
                Print["    NOTE: ", InputForm[reduced], " ≠ 0 (master 积分系数未抵消)"]
            ]
        ]
    , {tp, 1, Length[testPoints]}];
    
    Print["  --> ", passCount, "/", Length[testPoints], " PASS"];
    totalPass += passCount; totalFail += failCount;
, {i, 1, Length[nontrivial]}];

Print["\n========================================"];
Print["OVERALL: ", totalPass, " PASS, ", totalFail, " FAIL"];
Print["Result: ", If[totalFail === 0, "ALL PASS", "SOME FAIL - 需要检查关系正确性"]];
Print["========================================"];
Print["\n注意: 如果关系正确, reduced 应该全为 0"];
Print["如果 reduced ≠ 0, 说明 C++ 输出的关系可能有问题"];
Print["或者 LIE 关系需要更多 ν 点才能唯一确定"];