(* correct_kira_verify_v2.m *)
(* 正确的 Kira 验证流程（v2 修正版）
 * 
 * LIE 关系形式: sum_{alpha,beta} c_{alpha,beta} * g[nu - alpha] * nu^beta = 0
 * 
 * 验证步骤（根据 Wen-Hao Wu 的说明）:
 * 1. 取符号 ν 的 LIE 关系
 * 2. 在测试点 ν = (v1_test, v2_test) 处代入
 *    - g[v1-alpha1, v2-alpha2] -> j[nu1-alpha1, nu2-alpha2]
 *    - ν^β 因子 -> 代入数值
 * 3. g[a__] -> j[a] 转换
 * 4. 用 Kira 规则替换所有 j[...]-化简
 * 5. 检验结果是否为 0
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
Print["正确的 Kira 验证（v2 修正版）"];
Print["LIE 关系形式: sum c_{alpha,beta} * g[nu-alpha] * nu^beta = 0"];
Print["========================================"];

(* 读取 C++ 关系 *)
Get["/root/Large-Index-Expansion-CPP/verify/bub00/Relations_bub00_lev2_deg1.m"];
rel = $RelationResult;
relations = rel["Relations"];

Print["\nC++ Relations (lev=2, deg=1):"];
Print["  NumSolutions: ", rel["NumSolutions"]];
Print["\n原始关系:"];
Do[
    Print["  关系 ", i, ": ", InputForm[relations[[i]]]//Short]
, {i, 1, Length[relations]}];

(* ========================================================================
 * 解析 C++ 导出的多项式格式
 * 格式: coeff * v1^β1 * v2^β2 * g[v1-α1, v2-α2]
 * 
 * 其中 g[v1-α1, v2-α2] 表示 g[ν - α]
 *       v1^β1 * v2^β2 表示 ν^β
 * ======================================================================== *)

(* 解析一个项，返回 {系数, alpha, beta}
 * 输入: coeff * v1^β1 * v2^β2 * g[v1-α1, v2-α2]
 * 例如: 78498295*v1*v2*g[v1,v2-1]
 *)
parseTerm[term_] := Module[
    {factors, coef, v1Pow, v2Pow, gArgs, alpha1, alpha2},
    
    factors = If[Head[term] === Times, List @@ term, {term}];
    
    (* 分离系数部分和 g 部分 *)
    gPart = Cases[factors, g[__]][[1]];
    nonGPart = DeleteCases[factors, g[__]];
    
    (* 系数是所有非 g 因子的乘积 *)
    coef = Times @@ nonGPart;
    
    (* 解析 g 参数: g[v1-α1, v2-α2] *)
    gArgs = List @@ gPart;
    (* gArgs = {v1-α1, v2-α2} *)
    
    (* 解析 v 变量的指数 *)
    v1Pow = 0; v2Pow = 0;
    Do[
        f = factors[[i]];
        Which[
            Head[f] === Symbol && StringTake[SymbolName[f], 1] === "v",
                Which[
                    f === v1, v1Pow = 1,
                    f === v2, v2Pow = 1
                ],
            Head[f] === Power,
                Which[
                    f[[1]] === v1, v1Pow = f[[2]],
                    f[[1]] === v2, v2Pow = f[[2]]
                ]
        ]
    , {i, Length[factors]}];
    
    (* 提取 alpha: 从 g[v1-α1, v2-α2] 中提取 α1, α2 *)
    (* 格式是 "v1" (alpha=0), "v1-1" (alpha=1), "v1-2" (alpha=2) 等 *)
    alpha1 = 0; alpha2 = 0;
    Do[
        arg = gArgs[[k]];
        argName = If[Head[arg] === Symbol, SymbolName[arg], ""];
        If[k == 1,
            Which[
                arg === v1, alpha1 = 0,
                Head[arg] === Plus && arg[[1]] === v1, alpha1 = -arg[[2]],
                Head[arg] === Symbol, alpha1 = 0
            ],
            Which[
                arg === v2, alpha2 = 0,
                Head[arg] === Plus && arg[[1]] === v2, alpha2 = -arg[[2]],
                Head[arg] === Symbol, alpha2 = 0
            ]
        ]
    , {k, Length[gArgs]}];
    
    {coef, {alpha1, alpha2}, {v1Pow, v2Pow}}
];

(* 验证函数: 在指定 nu 点验证关系是否成立 *)
verifyRelationAtNu[relExpr_, nu1_, nu2_] := Module[
    {terms, parsed, jExpr, nuBeta, substituted, reduced, result},
    
    (* 解析关系中的各项 *)
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    
    jExpr = 0;
    Do[
        term = terms[[i]];
        {coef, alpha, beta} = parseTerm[term];
        
        (* 代入 nu 值:
         * - g[nu - alpha] -> j[nu1 - alpha1, nu2 - alpha2]
         * - nu^beta -> nu1^beta1 * nu2^beta2
         *)
        jArg1 = nu1 - alpha[[1]];
        jArg2 = nu2 - alpha[[2]];
        nuBeta = nu1^beta[[1]] * nu2^beta[[2]];
        
        termValue = coef * nuBeta * kiraJ[jArg1, jArg2];
        jExpr = jExpr + termValue;
        
        Print["    项: coeff=", coef, ", alpha=", alpha, ", beta=", beta];
        Print["      -> j[", jArg1, ",", jArg2, "] * ", nuBeta];
    , {i, Length[terms]}];
    
    Print["  代入后的 j 表达式: ", InputForm[jExpr]//Short];
    
    (* Kira 约化 *)
    reduced = kiraReduce[jExpr];
    Print["  Kira 约化后: ", InputForm[reduced]//Short];
    
    (* 检查是否为 0 *)
    result = (reduced === 0);
    
    result
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
    Print["\n========== 关系 ", i, " =========="];
    Print["原始: ", InputForm[relExpr]//Short];
    
    passCount = 0; failCount = 0;
    Do[
        {nu1, nu2} = testPoints[[tp]];
        Print["\n  测试点 nu=(", nu1, ",", nu2, "):"];
        
        result = verifyRelationAtNu[relExpr, nu1, nu2];
        
        status = If[result, "PASS", "FAIL"];
        If[result, passCount++, failCount++];
        Print["  --> ", status];
    , {tp, 1, Length[testPoints]}];
    
    Print["\n  --> ", passCount, "/", Length[testPoints], " PASS"];
    totalPass += passCount; totalFail += failCount;
, {i, 1, Length[nontrivial]}];

Print["\n========================================"];
Print["OVERALL: ", totalPass, " PASS, ", totalFail, " FAIL"];
Print["Result: ", If[totalFail === 0, "ALL PASS", "SOME FAIL"]];
Print["========================================"];
