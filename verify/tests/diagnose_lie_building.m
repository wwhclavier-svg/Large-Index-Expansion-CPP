(* diagnose_lie_building.m *)
(* 诊断 LIE 方程构建过程 *)

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
Print["LIE 方程构建诊断"];
Print["========================================"];

(* 读取 C++ 统一格式数据 *)
Get["/root/Large-Index-Expansion-CPP/verify/bub00/Compare-CPPRelation-bub00_lev2_deg1.m"];
data = $RelationResult;

alphas = data["Alphas"];
betas = data["Betas"];
coeffs = data["Coefficients"];
relations = data["Relations"];

Print["\nAlphas: ", alphas];
Print["Betas: ", betas];
Print["\n系数矩阵维度: ", Dimensions[coeffs]];
Print["NumSolutions: ", data["NumSolutions"]];

(* ========================================================================
 * 诊断 1: 检查 C++ 输出的关系结构
 * ======================================================================== *)
Print["\n========================================"];
Print["诊断 1: C++ 输出的关系"];
Print["========================================"];

(* 从 Unified 格式重建多项式格式 *)
numAlpha = Length[alphas];
numBeta = Length[betas];
numSols = data["NumSolutions"];

Do[
    solIdx = s;
    Print["\n--- Solution ", s, " ---"];
    
    expr = 0;
    Do[
        idx = a * numBeta + b + 1; (* MMA 1-indexed *)
        coeff = coeffs[[idx, s + 1]]; (* +1 因为第一列是特解 *)
        If[coeff =!= 0,
            alpha = alphas[[a + 1]];
            beta = betas[[b + 1]];
            term = coeff;
            (* 乘以 v^beta *)
            Do[
                If[beta[[i]] > 0,
                    term = term * Symbol["v" <> ToString[i]]^beta[[i]]
                ]
            , {i, Length[beta]}];
            (* 乘以 g[nu - alpha] *)
            gArg = Table[
                If[alpha[[i]] == 0,
                    Symbol["v" <> ToString[i]],
                    Symbol["v" <> ToString[i]] - alpha[[i]]
                ]
            , {i, Length[alpha]}];
            term = term * g @@ gArg;
            expr = expr + term;
        ]
    , {a, 0, numAlpha - 1}, {b, 0, numBeta - 1}];
    
    Print["表达式: ", InputForm[expr]//Short];
    relations[[s]] = expr;
, {s, 2, numSols + 1}]; (* 从 2 开始，跳过 "Particular" *)

(* ========================================================================
 * 诊断 2: 验证时检查每项的 alpha 和 beta
 * ======================================================================== *)
Print["\n========================================"];
Print["诊断 2: 检查每项的 alpha/beta 结构"];
Print["========================================"];

(* 解析 g 项的正确方式:
 * g[v1-alpha1, v2-alpha2] 表示 g[nu - alpha]
 * 代入 nu 后: g[nu1-alpha1, nu2-alpha2] -> j[nu1-alpha1, nu2-alpha2]
 *)

testNu = {1, 1};
{nu1, nu2} = testNu;

Do[
    relExpr = relations[[i]];
    Print["\n关系 ", i, ": ", InputForm[relExpr]//Short];
    
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    
    jAfterSubst = 0;
    Do[
        term = terms[[j]];
        factors = If[Head[term] === Times, List @@ term, {term}];
        gPart = First[Cases[factors, g[__]]];
        nonGPart = DeleteCases[factors, g[__]];
        coeff = Times @@ nonGPart;
        
        gArgs = List @@ gPart;
        
        (* 提取 alpha: v1 - k -> alpha1 = k, v2 - k -> alpha2 = k *)
        alpha1 = 0; alpha2 = 0;
        Do[
            arg = gArgs[[k]];
            Which[
                k == 1 && arg === v1, alpha1 = 0,
                k == 1 && Head[arg] === Plus && arg[[1]] === v1, alpha1 = -arg[[2]],
                k == 1 && Head[arg] === Plus && arg[[2]] === v1, alpha1 = -arg[[1]],
                k == 2 && arg === v2, alpha2 = 0,
                k == 2 && Head[arg] === Plus && arg[[1]] === v2, alpha2 = -arg[[2]],
                k == 2 && Head[arg] === Plus && arg[[2]] === v2, alpha2 = -arg[[1]],
                True, If[k == 1, alpha1 = 0, alpha2 = 0]
            ]
        , {k, 2}];
        
        (* 计算 j[nu - alpha] *)
        jArg1 = nu1 - alpha1;
        jArg2 = nu2 - alpha2;
        
        jTerm = coeff * kiraJ[jArg1, jArg2];
        jAfterSubst = jAfterSubst + jTerm;
        
        Print["  项 ", j, ": coeff=", coeff, ", alpha=(", alpha1, ",", alpha2, 
              ") -> j[", jArg1, ",", jArg2, "] = ", InputForm[kiraReduce[jTerm]]//Short];
    , {j, Length[terms]}];
    
    Print["  代入后: ", InputForm[jAfterSubst]//Short];
    Print["  Kira后: ", InputForm[kiraReduce[jAfterSubst]]//Short];
, {i, 1, Length[relations]}];

(* ========================================================================
 * 诊断 3: 检查零空间自洽性（用数值代入验证）
 * ======================================================================== *)
Print["\n========================================"];
Print["诊断 3: 零空间自洽性验证（MatrixBuild 方式）"];
Print["========================================"];

(* 这个测试用另一个 nu 点来验证关系是否在数值上成立 *)
(* 如果 C++ 零空间正确，关系应该对任何 nu 都成立 *)

testNu2 = {2, 3};
{nu1, nu2} = testNu2;

Print["\n用 nu = (2, 3) 验证:"];

Do[
    relExpr = relations[[i]];
    Print["\n关系 ", i];
    
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    
    total = 0;
    Do[
        term = terms[[j]];
        factors = If[Head[term] === Times, List @@ term, {term}];
        gPart = First[Cases[factors, g[__]]];
        nonGPart = DeleteCases[factors, g[__]];
        coeff = Times @@ nonGPart;
        
        gArgs = List @@ gPart;
        
        (* 提取 alpha *)
        alpha1 = 0; alpha2 = 0;
        Do[
            arg = gArgs[[k]];
            Which[
                k == 1 && arg === v1, alpha1 = 0,
                k == 1 && Head[arg] === Plus && arg[[1]] === v1, alpha1 = -arg[[2]],
                k == 1 && Head[arg] === Plus && arg[[2]] === v1, alpha1 = -arg[[1]],
                k == 2 && arg === v2, alpha2 = 0,
                k == 2 && Head[arg] === Plus && arg[[1]] === v2, alpha2 = -arg[[2]],
                k == 2 && Head[arg] === Plus && arg[[2]] === v2, alpha2 = -arg[[1]],
                True, If[k == 1, alpha1 = 0, alpha2 = 0]
            ]
        , {k, 2}];
        
        (* j[nu - alpha] 直接代入数值，不经过 Kira *)
        jArg1 = nu1 - alpha1;
        jArg2 = nu2 - alpha2;
        
        (* 这里是关键区别：我们直接代入数值到 g，然后转成 j *)
        (* 但 j 的值需要通过 Kira 规则计算 *)
        jValue = kiraReduce[kiraJ[jArg1, jArg2]];
        
        total = total + coeff * jValue;
    , {j, Length[terms]}];
    
    Print["  直接代入后的值: ", InputForm[total]//Short];
    Print["  是否为 0: ", total === 0];
, {i, 1, Length[relations]}];
