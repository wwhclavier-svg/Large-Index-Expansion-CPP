(* diagnose_cpp_kernel.m *)
(* 直接验证 C++ 零空间向量是否满足 M(ν)·x = 0 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];
Get["LIEWorkflow.wl"];

Print["========================================"];
Print["诊断: C++ 零空间正确性验证"];
Print["========================================"];

(* 加载 Kira 规则 *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];
kiraJ = KiraRuleLoader`j;

(* ========================================================================
 * 核心问题分析
 * ========================================================================
 * C++ 输出的关系形式: sum_{alpha,beta} coeff[alpha,beta,sol] * v^beta * g[nu - alpha] = 0
 *
 * 这里的 g[nu - alpha] 是什么意思?
 * - 如果 g[ν - α] 表示 g_{ν-α}(ν)（IBP 积分），则代入 ν 后得到 j[ν-α]
 * - 但 Kira 规则是关于 j[α,β] 的，约化到 master 的线性组合
 *
 * 关键问题: C++ 方程是 M(nu)·coeff(nu) = 0
 *          Kira 方程是 j[α,β] 的约化规则
 *          这两个东西如何对应?
 *
 * 物理上: j[α,β] 是 IBP 积分，系数依赖于运动学不变量
 *         关系 M(nu)·coeff(nu) = 0 表示 coeff(nu) 是 M(nu) 的零空间向量
 *         Kira 规则告诉我们如何将 j[α,β] 约化到 master
 *
 * ======================================================================== *)

(* 读取 C++ 导出的统一格式数据 *)
cppUnifiedFile = "/root/Large-Index-Expansion-CPP/verify/bub00/Compare-CPPRelation-bub00_lev2_deg1.m";
cppData = Get[cppUnifiedFile];

alphas = cppData["Alphas"];
betas = cppData["Betas"];
coeffs = cppData["Coefficients"];  (* 18 x 4 矩阵: 18 个 (alpha,beta) 对, 4 列 (1 特解 + 3 零空间向量) *)
numAlpha = Length[alphas];
numBeta = Length[betas];
numSols = cppData["NumSolutions"];

Print["\nC++ 数据:"];
Print["  Alphas: ", alphas];
Print["  Betas: ", betas];
Print["  Coefficients 维度: ", Dimensions[coeffs]];
Print["  解数量: ", numSols];

(* ========================================================================
 * 关键问题: C++ 的零空间向量对应什么?
 *
 * Mext 结构: [col 0] = 特解, [col 1..N] = 零空间基向量
 * 每个 Mext[i][s] 是长度为 numAlpha * numBeta 的向量的第 i 个元素
 *
 * Mext 按 alpha-major 排列: i = a * numBeta + b
 *   a = alpha 索引 (0..numAlpha-1)
 *   b = beta 索引 (0..numBeta-1)
 *
 * 所以 Mext[i][s] = coeff_{a,b,s}
 * ======================================================================== *)

(* ========================================================================
 * 验证: 重新构建矩阵 M(ν)，用数值点验证 M·x = 0
 * ======================================================================== *)

(* 首先，让我理解 Kira 规则和 C++ 方程的关系 *)
Print["\n========================================"];
Print["分析: Kira 规则 vs C++ 方程"];
Print["========================================"];

(* Kira 规则告诉我们:
 * j[α,β] 可以约化到 j[1,1] 的线性组合
 *
 * C++ 方程: sum_{a,b} c_{a,b}(nu) * v1^{b1} * v2^{b2} * g[ν - α] = 0
 *
 * 关键洞察: g[ν - α] 在代入 ν 后应该变成 Kira 的 j 积分
 *
 * 但问题是: g 是 C++ 的展开系数，不是直接的 IBP 积分!
 *)

(* 检查 Kira 规则覆盖的 sector *)
Print["\nKira 规则覆盖:"];
Print["  Master: j[1,1]"];
Print["  从 j[1,1] 可约化的 sector:"];
testSectors = {
    {1,1}, {1,2}, {1,3}, {2,1}, {2,2}, {2,3},
    {3,1}, {3,2}, {3,3}, {4,1}, {4,2},
    {0,1}, {0,2}, {1,0}, {2,0}, {0,0},
    {-1,1}, {1,-1}, {0,-1}, {-1,0}
};
Do[
    sector = testSectors[[i]];
    {a, b} = sector;
    red = kiraReduce[kiraJ[a, b]];
    status = If[red === 0, "-> 0",
                MatchQ[red, _kiraJ] || Head[red] === kiraJ, "-> master",
                "-> " <> ToString[FullForm[red]]];
    Print["  j[", a, ",", b, "] ", status];
, {i, Length[testSectors]}];

(* ========================================================================
 * 直接测试: 用 Kira 规则约化 C++ 关系的每一项
 * ======================================================================== *)

Print["\n========================================"];
Print["直接验证: 代入 nu=(1,1) 并用 Kira 规则约化"];
Print["========================================"];

nuTest = {1, 1};
{nu1, nu2} = nuTest;

(* 对于每个零空间向量 (sol 1, 2, 3) *)
Do[
    solIdx = s;  (* solIdx = 1 对应 Mext 列 1 (第一个零空间向量) *)
    
    Print["\n--- 零空间向量 ", s, " ---"];
    
    totalJ = 0;
    
    (* 遍历所有 (alpha, beta) *)
    Do[
        alpha = alphas[[a + 1]];
        beta = betas[[b + 1]];
        
        (* C++ 系数: Mext[idx][solIdx], idx = a * numBeta + b *)
        idx = a * numBeta + b + 1;  (* MMA 1-indexed *)
        coeff = coeffs[[idx, solIdx + 1]];  (* +1 因为第一列是特解 *)
        
        If[coeff =!= 0,
            (* 提取 beta 多项式因子: v1^beta1 * v2^beta2 *)
            vFactor = 1;
            Do[
                If[beta[[i]] > 0,
                    vFactor = vFactor * Symbol["v" <> ToString[i]]^beta[[i]]
                ]
            , {i, Length[beta]}];
            
            (* 提取 g[ν - alpha] 参数 *)
            gArgs = Table[
                If[alpha[[i]] == 0,
                    Symbol["v" <> ToString[i]],
                    Symbol["v" <> ToString[i]] - alpha[[i]]
                ]
            , {i, Length[alpha]}];
            
            (* 代入 ν 值 *)
            substitutedV = vFactor /. {v1 -> nu1, v2 -> nu2};
            gArgsSubst = gArgs /. {v1 -> nu1, v2 -> nu2};
            
            (* 转换为 j 积分: g[ν - α] -> j[ν - α] *)
            jArg1 = nu1 - alpha[[1]];
            jArg2 = nu2 - alpha[[2]];
            
            termJ = coeff * substitutedV * kiraJ[jArg1, jArg2];
            totalJ = totalJ + termJ;
            
            redJ = kiraReduce[termJ];
            Print["  项: coeff=", coeff, 
                  ", beta=", beta,
                  ", alpha=(", alpha[[1]], ",", alpha[[2]], ")",
                  " -> j[", jArg1, ",", jArg2, "] = ", InputForm[redJ]//Short];
        ];
    , {a, 0, numAlpha - 1}, {b, 0, numBeta - 1}];
    
    Print["  代入 nu=(1,1) 后的 j 和: ", InputForm[totalJ]//Short];
    Print["  Kira 约化后: ", InputForm[kiraReduce[totalJ]]//Short];
    Print["  是否为 0: ", kiraReduce[totalJ] === 0];
    
, {s, 1, numSols}];

(* ========================================================================
 * 另一种理解: 如果 C++ 零空间正确，它应该对任何 nu 都满足 M(nu)·x = 0
 * 让我用另一个 nu 点测试
 * ======================================================================== *)

Print["\n========================================"];
Print["用 nu=(2,3) 测试"];
Print["========================================"];

nuTest2 = {2, 3};
{nu1, nu2} = nuTest2;

Do[
    solIdx = s;
    
    Print["\n--- 零空间向量 ", s, " at nu=(2,3) ---"];
    
    totalJ = 0;
    
    Do[
        alpha = alphas[[a + 1]];
        beta = betas[[b + 1]];
        
        idx = a * numBeta + b + 1;
        coeff = coeffs[[idx, solIdx + 1]];
        
        If[coeff =!= 0,
            vFactor = 1;
            Do[
                If[beta[[i]] > 0,
                    vFactor = vFactor * Symbol["v" <> ToString[i]]^beta[[i]]
                ]
            , {i, Length[beta]}];
            
            substitutedV = vFactor /. {v1 -> nu1, v2 -> nu2};
            
            jArg1 = nu1 - alpha[[1]];
            jArg2 = nu2 - alpha[[2]];
            
            termJ = coeff * substitutedV * kiraJ[jArg1, jArg2];
            totalJ = totalJ + termJ;
        ];
    , {a, 0, numAlpha - 1}, {b, 0, numBeta - 1}];
    
    Print["  代入后的 j 和: ", InputForm[totalJ]//Short];
    Print["  Kira 约化后: ", InputForm[kiraReduce[totalJ]]//Short];
    Print["  是否为 0: ", kiraReduce[totalJ] === 0];
    
, {s, 1, numSols}];
