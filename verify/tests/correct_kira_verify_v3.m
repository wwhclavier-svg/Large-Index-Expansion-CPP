(* correct_kira_verify_v3.m *)
(* 正确的 Kira 验证流程（v3 修正版）
 * 
 * LIE 关系形式: sum_{alpha,beta} c_{alpha,beta} * g[nu - alpha] * nu^beta = 0
 * 
 * C++ 导出格式: coeff * v1^β1 * v2^β2 * g[v1-α1, v2-α2]
 * 其中 "v1-α" 表示 ν1 - α，所以 g[v1-1, v2] 对应 α = (1, 0)
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

kiraJ = KiraRuleLoader`j;

Print["========================================"];
Print["正确的 Kira 验证（v3 修正版）"];
Print["========================================"];

(* 读取 C++ 关系 *)
Get["/root/Large-Index-Expansion-CPP/verify/bub00/Relations_bub00_lev2_deg1.m"];
rel = $RelationResult;
relations = rel["Relations"];

Print["\nC++ Relations (lev=2, deg=1):"];
Print["  NumSolutions: ", rel["NumSolutions"]];

(* ========================================================================
 * 正确解析 C++ 导出的多项式格式
 * 格式: coeff * v1^β1 * v2^β2 * g[v1-α1, v2-α2]
 * 
 * g 参数解析:
 *   "v1"    -> ν1 - 0 -> α1 = 0
 *   "v1-1"  -> ν1 - 1 -> α1 = 1
 *   "v1-2"  -> ν1 - 2 -> α1 = 2
 *   "v2-1"  -> ν2 - 1 -> α2 = 1
 *   "-2+v2" -> ν2 - 2 -> α2 = 2 (等价于 v2-2)
 * ======================================================================== *)

(* 从 "v1-α" 或 "-α+v1" 格式提取 α *)
extractAlphaFromTerm[term_] := Module[
    {coef, const},
    Which[
        (* 格式: v1 ± k 或 k ± v1 *)
        Head[term] === Symbol,
            Return[0],
        Head[term] === Plus,
            (* 分两种情况: v1 - k 或 -k + v1 *)
            If[term[[1]] === v1 || term[[1]] === v2,
                (* v1 - k 格式 *)
                coef = term[[1]]; const = term[[2]],
                (* -k + v1 格式 *)
                coef = term[[2]]; const = term[[1]]
            ];
            (* α = -const_part *)
            Return[-const],
        Head[term] === Power,
            Return[term[[2]]]
    ];
    0
];

(* 解析 g[...] 的参数，返回 alpha 向量 *)
parseGArgs[gargs_] := Module[
    {alpha1, alpha2},
    alpha1 = extractAlphaFromTerm[gargs[[1]]];
    alpha2 = extractAlphaFromTerm[gargs[[2]]];
    {alpha1, alpha2}
];

(* 解析一个项，返回 {系数, alpha, beta} *)
parseTerm[term_] := Module[
    {factors, coef, v1Pow, v2Pow, gPart, gArgs, alpha},
    
    factors = If[Head[term] === Times, List @@ term, {term}];
    
    (* 分离系数部分和 g 部分 *)
    gPart = Cases[factors, g[__]][[1]];
    nonGPart = DeleteCases[factors, g[__]];
    
    (* 系数是所有非 g 因子的乘积（这里包括 v1^β1 * v2^β2） *)
    coef = Times @@ nonGPart;
    
    (* 解析 v 变量的指数作为 beta *)
    v1Pow = 0; v2Pow = 0;
    Do[
        f = factors[[i]];
        Which[
            f === v1, v1Pow = v1Pow + 1,
            f === v2, v2Pow = v2Pow + 1,
            Head[f] === Power,
                Which[
                    f[[1]] === v1, v1Pow = f[[2]],
                    f[[1]] === v2, v2Pow = f[[2]]
                ]
        ]
    , {i, Length[factors]}];
    
    (* 解析 g 参数得到 alpha *)
    gArgs = List @@ gPart;
    alpha = parseGArgs[gArgs];
    
    {coef, alpha, {v1Pow, v2Pow}}
];

(* 验证函数 *)
verifyRelationAtNu[relExpr_, nu1_, nu2_] := Module[
    {terms, jExpr, result, reduced},
    
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    
    jExpr = 0;
    Do[
        term = terms[[i]];
        {coef, alpha, beta} = parseTerm[term];
        
        (* 代入: g[nu-alpha] -> j[nu-alpha] *)
        jArg1 = nu1 - alpha[[1]];
        jArg2 = nu2 - alpha[[2]];
        
        (* 代入: nu^beta -> nu1^beta1 * nu2^beta2 *)
        nuBeta = nu1^beta[[1]] * nu2^beta[[2]];
        
        termValue = coef * nuBeta * kiraJ[jArg1, jArg2];
        jExpr = jExpr + termValue;
    , {i, Length[terms]}];
    
    (* Kira 约化 *)
    reduced = kiraReduce[jExpr];
    
    (* 检查是否为 0 *)
    result = (reduced === 0);
    
    {result, jExpr, reduced}
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
        
        (* 解析并显示中间结果 *)
        terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
        Do[
            term = terms[[j]];
            {coef, alpha, beta} = parseTerm[term];
            Print["    项: coeff=", coef, ", alpha=", alpha, ", beta=", beta];
        , {j, Length[terms]}];
        
        {result, jExpr, reduced} = verifyRelationAtNu[relExpr, nu1, nu2];
        
        status = If[result, "PASS", "FAIL"];
        If[result, passCount++, failCount++];
        Print["  代入后: ", InputForm[jExpr]//Short];
        Print["  Kira后: ", InputForm[reduced]//Short];
        Print["  --> ", status];
    , {tp, 1, Length[testPoints]}];
    
    Print["  --> ", passCount, "/", Length[testPoints], " PASS"];
    totalPass += passCount; totalFail += failCount;
, {i, 1, Length[nontrivial]}];

Print["\n========================================"];
Print["OVERALL: ", totalPass, " PASS, ", totalFail, " FAIL"];
Print["========================================"];
