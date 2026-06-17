(* diagnose_cpp_nullspace.m *)
(* 诊断: 为什么 C++ 零空间解不满足 Kira 规则 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];
kiraJ = KiraRuleLoader`j;

Print["============================================="];
Print["诊断: C++ 零空间 vs 物理 Kira 规则"];
Print["============================================="];

(* ==============================================================
 * 第一部分: 解析 C++ 关系
 * ============================================================== *)
Print["\n=== C++ 关系解析 (lev=2 deg=1) ==="];

(* 手动输入关系 *)
basis1 = 78498295*g[v1,v2] + 89712338*g[v1,v2-1] + 1*g[v1,v2-2];
basis2 = 78498295*v2*g[v1,v2] + 89712338*v2*g[v1,v2-1] + 1*v2*g[v1,v2-2];
basis3 = 78498295*v1*g[v1,v2] + 89712338*v1*g[v1,v2-1] + 1*v1*g[v1,v2-2];

(* 解析函数: 从 g[nu-alpha] 提取 alpha, 从系数提取 beta *)
extractTermInfo[term_] := Module[
    {factors, coef, v1Pow, v2Pow, gArgs, al1, al2},
    factors = If[Head[term] === Times, List @@ term, {term}];
    gArgs = Cases[factors, g[__]][[1]];
    gArgs = List @@ gArgs;
    
    (* 提取 alpha: g[arg1, arg2] 中的 alpha = 0 if symbol, or 0 - constant *)
    al1 = Which[Head[gArgs[[1]]] === Symbol, 0,
                Head[gArgs[[1]]] === Integer, -gArgs[[1]],
                Head[gArgs[[1]]] === Plus, -gArgs[[1,2]]];
    al2 = Which[Head[gArgs[[2]]] === Symbol, 0,
                Head[gArgs[[2]]] === Integer, -gArgs[[2]],
                Head[gArgs[[2]]] === Plus, -gArgs[[2,2]]];
    
    (* 提取非 g 因子 *)
    nonG = DeleteCases[factors, g[__]];
    coef = Times @@ nonG;
    
    (* 提取 v1, v2 指数 *)
    v1Pow = 0; v2Pow = 0;
    Do[
        f = nonG[[i]];
        Which[
            f === v1, v1Pow++,
            f === v2, v2Pow++,
            Head[f] === Power && f[[1]] === v1, v1Pow = f[[2]],
            Head[f] === Power && f[[1]] === v2, v2Pow = f[[2]]
        ]
    , {i, Length[nonG]}];
    
    {coef, {al1, al2}, {v1Pow, v2Pow}}
];

parseRelation[relExpr_] := Module[{terms, info},
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    info = extractTermInfo /@ terms;
    info
];

Print["Basis1:"];
Do[Print["  coeff=", t[[1]], " α=", t[[2]], " β=", t[[3]]], {t, parseRelation[basis1]}];

Print["Basis2:"];
Do[Print["  coeff=", t[[1]], " α=", t[[2]], " β=", t[[3]]], {t, parseRelation[basis2]}];

(* ==============================================================
 * 第二部分: Kira 残差分析
 * ============================================================== *)
Print["\n=== Kira 残差分析 ==="];

testPoints = {{1,1}, {1,2}, {2,1}, {2,3}};

evaluateAtNu[relation_, nu1_, nu2_] := Module[
    {terms, jExpr},
    terms = If[Head[relation] === Plus, List @@ relation, {relation}];
    jExpr = 0;
    Do[
        info = extractTermInfo[terms[[i]]];
        {coef, alpha, beta} = info;
        nuBeta = nu1^beta[[1]] * nu2^beta[[2]];
        jExpr = jExpr + coef * nuBeta * kiraJ[nu1 - alpha[[1]], nu2 - alpha[[2]]];
    , {i, Length[terms]}];
    kiraReduce[jExpr]
];

Print["Basis1 残差:"];
Do[
    {nu1, nu2} = tp;
    residual = evaluateAtNu[basis1, nu1, nu2];
    Print["  ν=(", nu1, ",", nu2, "): ", residual, " → ", If[residual === 0, "PASS", "FAIL"]];
, {tp, testPoints}];

Print["Basis2 残差:"];
Do[
    {nu1, nu2} = tp;
    residual = evaluateAtNu[basis2, nu1, nu2];
    Print["  ν=(", nu1, ",", nu2, "): ", residual, " → ", If[residual === 0, "PASS", "FAIL"]];
, {tp, testPoints}];

Print["Basis3 残差:"];
Do[
    {nu1, nu2} = tp;
    residual = evaluateAtNu[basis3, nu1, nu2];
    Print["  ν=(", nu1, ",", nu2, "): ", residual, " → ", If[residual === 0, "PASS", "FAIL"]];
, {tp, testPoints}];

(* ==============================================================
 * 第三部分: 检查 Kira 规则中每个 j[α,β] 的约化结果
 * ============================================================== *)
Print["\n=== Kira 规则关键数值 ==="];

(* 找出 Basis1 中出现的所有 j[ν-α] 映射到 Kira 的结果 *)
nuPoints = {{1,1}, {2,3}};
Do[
    {nu1, nu2} = nu;
    Print["\nν=(", nu1, ",", nu2, "):"];
    (* j[ν1, ν2] -> master? *)
    Print["  j[", nu1, ",", nu2, "] = ", kiraJ[nu1, nu2]];
    Print["  j[", nu1, ",", nu2-1, "] = ", kiraJ[nu1, nu2-1]];
    Print["  j[", nu1, ",", nu2-2, "] = ", kiraJ[nu1, nu2-2]];
    Print["  → j[ν-(0,0)] = ", kiraReduce[kiraJ[nu1, nu2]]];
    Print["  → j[ν-(0,1)] = ", kiraReduce[kiraJ[nu1, nu2-1]]];
    Print["  → j[ν-(0,2)] = ", kiraReduce[kiraJ[nu1, nu2-2]]];
    
    (* 这些在大指标下的值 *)
    bigNu1 = nu1; bigNu2 = nu2;
    Print["  j[", bigNu1, ",", bigNu2, "] = ", InputForm[kiraJ[bigNu1, bigNu2]]];
    Print["  j[", bigNu1, ",", bigNu2-1, "] = ", InputForm[kiraJ[bigNu1, bigNu2-1]]];
, {nu, nuPoints}];

(* ==============================================================
 * 第四部分: MMA 关系的关键项分析
 * 检查 MMA 关系中 j[1,1] 的系数和是如何消去的
 * ============================================================== *)
Print["\n=== MMA 关系关键项分析 (ν=(1,1) 处 j[1,1] 消去机制) ==="];

(* 手动展开 MMA 关系1 (来自 verify_mma_kira.m 输出):
 * -119616448*g[-1+v1,v2] + v1*g[-1+v1,v2] + 59808226*g[v1,-1+v2] + 
 * 179424670*v1*g[v1,-1+v2] + 179424671*v2*g[v1,-1+v2] + 3*g[v1,v2] + 
 * 3*v1*g[v1,v2] + 179424667*v2*g[v1,v2] 
 * 
 * 代入 ν=(1,1):
 * j[0,1] 项: -119616448 + 1 = -119616447
 * j[1,0] 项: 59808226 + 179424670 + 179424671 = 418757567 mod p
 * j[1,1] 项: 3 + 3 + 179424667 = 179424673 ≡ 0 (mod p)
 *)

Print["MMA 关系1 在 ν=(1,1):"];
Print["  j[0,1] 系数: (-119616448 + 1) = ", (-119616448 + 1)];
Print["  j[1,0] 系数: (59808226 + 179424670 + 179424671) = ", 
      59808226 + 179424670 + 179424671];
Print["  j[1,1] 系数: (3 + 3 + 179424667) = ", 3 + 3 + 179424667, 
      " ≡ 0 mod 179424673"];

Print["\nC++ Basis1 在 ν=(1,1):"];
Print["  j[1,-1] 系数: 1 (→ 0, sector 0)"];
Print["  j[1,0]  系数: 89712338 (→ 0, sector 0)"];
Print["  j[1,1]  系数: 78498295 (→ master integral, ≠ 0)"];
Print["  → 只有一个 j[1,1] 项，无法消去"];

(* ==============================================================
 * 第五部分: 验证 C++ 的 M(ν) 矩阵是否包含 MMA 关系
 * 通过检查: M(ν) * v_mma = 0 对随机 ν 成立?
 * ============================================================== *)
Print["\n=== M(ν) 矩阵自洽性检查 ==="];

(* 构建 M(ν): 根据 LIE 理论, M(ν) 的行来自 seriesCoefficient h_jk(ν-α) *)
(* 对于每个 (α, β), 列定义为 c[α,β] *)
(* 行: 系数展开的 n^kt 项 *)
(* 
   M(ν)_{row, col} = Σ_regime Σ_nimax h_jk(ν-α, nimax) 
   其中 row = (regime, nimax, j, kt)
    col = (α, β) 
*)

(* 关键问题: M(ν) 是否构建正确? *)
(* 我们通过加载 C++ 层递归系数来验证 *)

Get["/root/Large-Index-Expansion-CPP/verify/bub00/Compare-CPPResult-bup00.m"];

(* 验证展开系数本身 *)
Print["C++ 展开系数 (Kmax=4, NE=2, NB=1):"];
$ExpansionResults // InputForm // Print;

