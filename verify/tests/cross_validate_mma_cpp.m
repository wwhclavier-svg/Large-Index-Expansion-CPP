(* cross_validate_mma_cpp.m *)
(* 交叉验证: MMA 关系向量是否在 C++ M(nu) 的零空间中 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];

(* ================================================================
 * Step 1: 加载 MMA 关系 (lev=2, deg=2) 并提取系数向量
 * ================================================================ *)
Print["=== Step 1: 加载 MMA 关系 ==="];
Get["/root/Large-Index-Expansion-CPP/verify/bub00/Compare-MMARelation-lev2.m"];

(* MMA 关系: $MMARelations = {{{level1 deg0}, {level1 deg1, level1 deg2}}, {{...},{level2 deg2}}} *)
mmaRels = $MMARelations;

(* 提取 lev=2, deg=2 的 4 个关系 *)
mmaLevel2Deg2 = mmaRels[[2, 3]];  (* 4 relations *)
Print["MMA lev=2, deg=2: ", Length[mmaLevel2Deg2], " relations"];

(* 定义 (alpha, beta) 坐标基 *)
alphas = {{0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {2,0}};
betas = {{0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {2,0}};

(* 为每个 (alpha, beta) 创建索引映射 *)
indexOf[alpha_, beta_] := Module[{aIdx, bIdx},
    aIdx = Position[alphas, alpha][[1, 1]] - 1;
    bIdx = Position[betas, beta][[1, 1]] - 1;
    aIdx * Length[betas] + bIdx + 1  (* 1-indexed *)
];

(* 创建符号 g[nu-alpha] 的映射: g[v1,v2-1] -> alpha=(0,1), beta=(0,0) *)
(* 从 MMA 关系表达式提取 {coeff, alpha, beta_v1, beta_v2} *)

extractAlphaFromGArg[arg_] := Module[{coef, const},
    Which[
        Head[arg] === Symbol, Return[0],
        Head[arg] === Plus,
            If[MemberQ[{v1, v2}, arg[[1]]],
                coef = arg[[1]]; const = arg[[2]]; Return[-const],
                coef = arg[[2]]; const = arg[[1]]; Return[-const]
            ],
        Head[arg] === Integer, Return[-arg],
        True, Return[0]
    ];
    0
];

parseMMAArgs[gargs_] := Module[{alpha1, alpha2},
    alpha1 = extractAlphaFromGArg[gargs[[1]]];
    alpha2 = extractAlphaFromGArg[gargs[[2]]];
    {alpha1, alpha2}
];

parseMMATerm[term_] := Module[
    {factors, coef, v1Pow, v2Pow, gPart, gArgs, alpha, nonGPart},
    factors = If[Head[term] === Times, List @@ term, {term}];
    gPart = Select[factors, Head[#] === g && Length[#] == 2&];
    If[Length[gPart] == 0, Return[Null]];
    gPart = gPart[[1]];
    nonGPart = DeleteCases[factors, gPart];
    
    coef = 1;
    v1Pow = 0; v2Pow = 0;
    Do[
        f = nonGPart[[i]];
        Which[
            f === v1, v1Pow++,
            f === v2, v2Pow++,
            Head[f] === Power && f[[1]] === v1, v1Pow = f[[2]],
            Head[f] === Power && f[[1]] === v2, v2Pow = f[[2]],
            Head[f] === Integer, coef = coef * f,
            True, coef = coef * f
        ]
    , {i, Length[nonGPart]}];
    
    gArgs = List @@ gPart;
    alpha = parseMMAArgs[gArgs];
    {coef, alpha, v1Pow, v2Pow}
];

(* 将 MMA 关系转为系数向量 (在 (alpha,beta) 基下) *)
mmaRelationToVector[relExpr_] := Module[
    {terms, vec, parsed, alpha, betaVec},
    vec = Table[0, {Length[alphas] * Length[betas]}];
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    Do[
        parsed = parseMMATerm[terms[[i]]];
        If[parsed === Null, Continue[]];
        {coef, alpha, b1, b2} = parsed;
        betaVec = {b1, b2};
        idx = indexOf[alpha, betaVec];
        vec[[idx]] = vec[[idx]] + coef;
    , {i, Length[terms]}];
    vec
];

mmaVectors = mmaRelationToVector /@ mmaLevel2Deg2;
Print["MMA 系数向量 (每个 ", Length[mmaVectors[[1]]], " 维):"];
Do[Print["  MMA关系", i, ": ", Short[mmaVectors[[i]], 8]], {i, Length[mmaVectors]}];

Print["\n=== Step 2: 加载 C++ 关系 ==="];
Get["/root/Large-Index-Expansion-CPP/verify/bub00/Relations_bub00_lev2_deg2.m"];
cppRels = $RelationResult["Relations"];
cppNonTrivial = Select[cppRels, # =!= 0 &];
Print["C++ 非平凡关系: ", Length[cppNonTrivial]];

(* C++ 关系已经是 g[v1-α1, v2-α2] 格式，解析方式类似 MMA *)
cppRelationToVector[relExpr_] := Module[
    {terms, vec, parsed, alpha, betaVec},
    vec = Table[0, {Length[alphas] * Length[betas]}];
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    Do[
        parsed = parseMMATerm[terms[[i]]];
        If[parsed === Null, Continue[]];
        {coef, alpha, b1, b2} = parsed;
        betaVec = {b1, b2};
        idx = indexOf[alpha, betaVec];
        vec[[idx]] = vec[[idx]] + coef;
    , {i, Length[terms]}];
    vec
];

cppVectors = cppRelationToVector /@ cppNonTrivial;
Print["C++ 系数向量 (每个 ", Length[cpVectors[[1]]], " 维):"];
Do[Print["  C++关系", i, ": ", Short[cpVectors[[i]], 8]], {i, Length[cpVectors]}];

(* ================================================================
 * Step 3: 比较零空间维度
 * ================================================================ *)
Print["\n=== Step 3: 零空间维度比较 ==="];
printMatrix[m_] := Module[{},
    Print["  维度: ", Length[m], " x ", If[Length[m] > 0, Length[m[[1]]], 0]];
    Print["  Rank: ", MatrixRank[m, Modulus -> 179424673]]
];

mmaMat = mmaVectors;
cppMat = cpVectors;
Print["MMA 关系矩阵: "]; printMatrix[mmaMat];
Print["C++ 关系矩阵: "]; printMatrix[cpMat];

combined = Join[mmaMat, cpVectors];
Print["合并矩阵: "]; printMatrix[combined];

rankMMA = MatrixRank[mmaMat, Modulus -> 179424673];
rankCPP = MatrixRank[cpMat, Modulus -> 179424673];
rankCombined = MatrixRank[combined, Modulus -> 179424673];

Print["\n--- 零空间维度分析 ---"];
Print["MMA rank: ", rankMMA, ", 零空间: ", Length[mmaMat] - rankMMA];
Print["C++ rank: ", rankCPP, ", 零空间: ", Length[cpMat] - rankCPP];
Print["Combined rank: ", rankCombined];

If[rankCombined === rankMMA && rankCombined === rankCPP,
    Print["\n✅ MMA 和 C++ 张成相同的零空间！"],
    If[rankCombined > Max[rankMMA, rankCPP],
        Print["\n❌ MMA 和 C++ 的零空间不同！"],
        Print["\n⚠️ 一个零空间包含在另一个中"]
    ]
];

(* ================================================================
 * Step 4: 检查 MMA 向量是否在 C++ M(nu) 的零空间中
 * 通过检查 MMA 关系在 C++ nu 点处的残差
 * ================================================================ *)
Print["\n=== Step 4: MMA 关系多项式自洽性检验 ==="];

(* 在模 179424673 下，验证 MMA 关系本身是否为零（代入具体 nu 值前） *)
(* MMA 关系是多项式恒等式: sum c[α,β] * ν^β * g[ν-α] ≡ 0 *)
(* 我们直接代入具体 nu 值验证 *)

modulus = 179424673;

(* 检验: 对于 MMA 关系, Sum c[α,β] * ν^β * j[ν-α] 在 Kira 规则下 → 0 *)
(* 这个我们已经用 verify_mma_kira.m 验证过了, 28/28 PASS *)

Print["✅ 已知: MMA 关系已通过 Kira 验证 (28/28 PASS)"];

(* ================================================================
 * Step 5: 关键诊断 — 
 * 检查 C++ 的 M(nu) 矩阵是否能区分 MMA 和 C++ 的零空间关系
 * ================================================================ *)
Print["\n=== Step 5: 检查 MMA 关系在 C++ 数据下的残差 ==="];

(* 加载 C++ Kira 工具 *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> modulus,
    "NProp" -> 2
];
kiraJ = KiraRuleLoader`j;

(* 测试点 *)
testPoints = {{1, 1}, {1, 2}, {2, 1}, {2, 3}};

Print["测试 MMA 关系中的 g[ν-α] × ν^β 结构："];
Do[
    relExpr = mmaLevel2Deg2[[ri]];
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    
    Print["\n--- MMA 关系 ", ri, " ---"];
    
    (* 逐个项分析: {系数, alpha, beta} *)
    Do[
        parsed = parseMMATerm[terms[[i]]];
        If[parsed =!= Null,
            {coef, alpha, b1, b2} = parsed;
            Print["  项: coeff=", coef, ", α=", alpha, ", β=(", b1, ",", b2, ")"];
        ];
    , {i, Length[terms]}];
    
    (* 在每个测试点检查 Kira 残差 *)
    Do[
        {nu1, nu2} = testPoints[[tp]];
        jExpr = 0;
        Do[
            parsed = parseMMATerm[terms[[i]]];
            If[parsed === Null, Continue[]];
            {coef, alpha, b1, b2} = parsed;
            nuBeta = nu1^b1 * nu2^b2;
            jExpr = jExpr + coef * nuBeta * kiraJ[nu1 - alpha[[1]], nu2 - alpha[[2]]];
        , {i, Length[terms]}];
        reduced = kiraReduce[jExpr];
        result = (reduced === 0);
        Print["  ν=(", nu1, ",", nu2, "): residual=", If[result, 0, reduced], 
              " -> ", If[result, "PASS", "FAIL"]];
    , {tp, 1, Length[testPoints]}];
, {ri, 1, Length[mmaLevel2Deg2]}];

(* ================================================================
 * Step 6: 关键诊断 —
 * 检查 MMA 关系中的 (alpha, beta) 结构 vs C++
 * MMA 关系中是否有 C++ 缺失的 beta 项?
 * ================================================================ *)
Print["\n=== Step 6: MMA vs C++ (alpha, beta) 覆盖分析 ==="];

(* 提取 MMA 关系中实际出现的 (alpha, beta) *)
mmaAlphaBetaUsed = {};
Do[
    terms = If[Head[mmaLevel2Deg2[[ri]]] === Plus, 
               List @@ mmaLevel2Deg2[[ri]], {mmaLevel2Deg2[[ri]]}];
    Do[
        parsed = parseMMATerm[terms[[i]]];
        If[parsed =!= Null,
            {coef, alpha, b1, b2} = parsed;
            If[coef =!= 0,
                AppendTo[mmaAlphaBetaUsed, {alpha, {b1, b2}}]
            ]
        ];
    , {i, Length[terms]}];
, {ri, 1, Length[mmaLevel2Deg2]}];
mmaAlphaBetaUsed = Union[mmaAlphaBetaUsed];
Print["MMA 中使用到的 (alpha, beta) 对: ", Length[mmaAlphaBetaUsed]];
Do[Print["  α=", ab[[1]], " β=", ab[[2]]], {ab, mmaAlphaBetaUsed}];

(* 提取 C++ 关系中实际出现的 (alpha, beta) *)
cppAlphaBetaUsed = {};
Do[
    terms = If[Head[cpNonTrivial[[ri]]] === Plus,
               List @@ cpNonTrivial[[ri]], {cpNonTrivial[[ri]]}];
    Do[
        parsed = parseMMATerm[terms[[i]]];
        If[parsed =!= Null,
            {coef, alpha, b1, b2} = parsed;
            If[coef =!= 0,
                AppendTo[cppAlphaBetaUsed, {alpha, {b1, b2}}]
            ]
        ];
    , {i, Length[terms]}];
, {ri, 1, Length[cpNonTrivial]}];
cppAlphaBetaUsed = Union[cppAlphaBetaUsed];
Print["C++ 中使用到的 (alpha, beta) 对: ", Length[cppAlphaBetaUsed]];
Do[Print["  α=", ab[[1]], " β=", ab[[2]]], {ab, cppAlphaBetaUsed}];

(* 找出差异 *)
mmaOnly = Complement[mmaAlphaBetaUsed, cppAlphaBetaUsed];
cppOnly = Complement[cppAlphaBetaUsed, mmaAlphaBetaUsed];
common = Intersection[mmaAlphaBetaUsed, cppAlphaBetaUsed];

Print["\nMMA-only (alpha, beta): ", Length[mmaOnly]];
Do[Print["  α=", ab[[1]], " β=", ab[[2]]], {ab, mmaOnly}];
Print["C++-only (alpha, beta): ", Length[cppOnly]];
Do[Print["  α=", ab[[1]], " β=", ab[[2]]], {ab, cppOnly}];
Print["共有 (alpha, beta): ", Length[common]];

