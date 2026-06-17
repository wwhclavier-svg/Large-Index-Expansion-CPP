(* debug_lie_relations.m *)
(* 诊断 C++ LIE 关系与 Kira 物理规则的对应性 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];
kiraJ = KiraRuleLoader`j;

Print["========================================"];
Print["诊断: C++ LIE 零空间 vs Kira 物理规则"];
Print["========================================"];

(* ========================================
 * 背景
 * ========================================
 * C++ 输出形式: sum_{α,β} coeff[α,β,sol] * v^β * g[ν-α] = 0
 * Kira 规则形式: j[α,β] -> master j[...] 的线性组合
 *
 * 关键问题: g[ν-α] 是否对应 Kira 的 j[α,β]？
 * ======================================== *)

(* ========================================
 * 测试 1: Kira 规则本身是否正确?
 * ======================================== *)
Print["\n[Test 1] Kira 规则有效性测试"];
Print["-------------------------------------"];

(* 验证 Kira 规则的一致性 *)
(* 如果 j[2,1] = (d-3)*j[1,1]，那么 j[2,1] - (d-3)*j[1,1] = 0 应该是平凡的 *)

dVal = 1/3;
modulus = 179424673;

expr1 = kiraJ[2,1] - (dVal - 3)*kiraJ[1,1];
red1 = kiraReduce[expr1];
Print["j[2,1] - (d-3)*j[1,1] = ", red1];
Print["  应该 = 0, 结果: ", red1 === 0];

expr2 = kiraJ[2,2] - (dVal^2 - 9*dVal + 18)*kiraJ[1,1];
red2 = kiraReduce[expr2];
Print["j[2,2] - (d^2-9d+18)*j[1,1] = ", red2];
Print["  应该 = 0, 结果: ", red2 === 0];

(* ========================================
 * 测试 2: C++ 输出的零空间向量
 * ======================================== *)
Print["\n[Test 2] C++ 零空间向量"];
Print["-------------------------------------"];

cppData = Get["/root/Large-Index-Expansion-CPP/verify/bub00/Compare-CPPRelation-bub00_lev2_deg1.m"];
relations = cppData["Relations"];

Print["C++ 输出的关系:"];
ForEach[rel, relations,
    Print["  ", rel];
];

(* ========================================
 * 测试 3: 直接代入 nu 并用 Kira 规则约化
 * ======================================== *)
Print["\n[Test 3] 代入 nu=(1,1) 并用 Kira 规则约化"];
Print["-------------------------------------"];

nu1 = 1; nu2 = 1;

(* 提取每个关系的系数 *)
(* C++ 关系形式: sum coeff * v^beta * g[v1 - alpha1, v2 - alpha2] *)
(* 代入 nu: g[nu1 - alpha1, nu2 - alpha2] -> kiraJ[nu1 - alpha1, nu2 - alpha2] *)

(* Basis1: 78498295*g[v1,v2] + 89712338*g[v1,v2-1] + 1*g[v1,v2-2] *)
expr_B1 = 78498295 * kiraJ[nu1, nu2] + 89712338 * kiraJ[nu1, nu2-1] + 1 * kiraJ[nu1, nu2-2];
red_B1 = kiraReduce[expr_B1];
Print["Basis1 at nu=(1,1): ", red_B1];
Print["  = 0? ", red_B1 === 0];

(* Basis2: 78498295*v2*g[v1,v2] + 89712338*v2*g[v1,v2-1] + 1*v2*g[v1,v2-2] *)
(* 代入 v2=1 *)
expr_B2 = 78498295 * nu2 * kiraJ[nu1, nu2] + 89712338 * nu2 * kiraJ[nu1, nu2-1] + 1 * nu2 * kiraJ[nu1, nu2-2];
red_B2 = kiraReduce[expr_B2];
Print["Basis2 at nu=(1,1): ", red_B2];
Print["  = 0? ", red_B2 === 0];

(* Basis3: 78498295*v1*g[v1,v2] + 89712338*v1*g[v1,v2-1] + 1*v1*g[v1,v2-2] *)
(* 代入 v1=1 *)
expr_B3 = 78498295 * nu1 * kiraJ[nu1, nu2] + 89712338 * nu1 * kiraJ[nu1, nu2-1] + 1 * nu1 * kiraJ[nu1, nu2-2];
red_B3 = kiraReduce[expr_B3];
Print["Basis3 at nu=(1,1): ", red_B3];
Print["  = 0? ", red_B3 === 0];

(* ========================================
 * 测试 4: 用另一个 nu 点测试
 * ======================================== *)
Print["\n[Test 4] 代入 nu=(2,3) 并用 Kira 规则约化"];
Print["-------------------------------------"];

nu1 = 2; nu2 = 3;

expr_B1 = 78498295 * kiraJ[nu1, nu2] + 89712338 * kiraJ[nu1, nu2-1] + 1 * kiraJ[nu1, nu2-2];
red_B1 = kiraReduce[expr_B1];
Print["Basis1 at nu=(2,3): ", red_B1];
Print["  = 0? ", red_B1 === 0];

expr_B2 = 78498295 * nu2 * kiraJ[nu1, nu2] + 89712338 * nu2 * kiraJ[nu1, nu2-1] + 1 * nu2 * kiraJ[nu1, nu2-2];
red_B2 = kiraReduce[expr_B2];
Print["Basis2 at nu=(2,3): ", red_B2];
Print["  = 0? ", red_B2 === 0];

expr_B3 = 78498295 * nu1 * kiraJ[nu1, nu2] + 89712338 * nu1 * kiraJ[nu1, nu2-1] + 1 * nu1 * kiraJ[nu1, nu2-2];
red_B3 = kiraReduce[expr_B3];
Print["Basis3 at nu=(2,3): ", red_B3];
Print["  = 0? ", red_B3 === 0];

(* ========================================
 * 结论
 * ======================================== *)
Print["\n========================================"];
Print["结论"];
Print["========================================"];
Print["如果所有测试的 '= 0?' 都是 False，"];
Print["说明 C++ 输出的零空间向量不满足 Kira 物理规则。"];

Print["\n可能的原因:"];
Print["1. C++ 的零空间计算有 bug"];
Print["2. C++ 输出的 g[ν-α] 和 Kira 的 j[α,β] 不是同一物理量"];
Print["3. C++ 和 Kira 使用的是不同的 integral family 配置"];
