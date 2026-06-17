(* simple_kernel_test.m *)
(* 最简单的测试: 直接用 Kira 规则约化 C++ 输出的每个关系 *)

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
Print["简单测试: 直接约化 C++ 输出的关系"];
Print["========================================"];

(* C++ 输出的 Relations:
 * Basis1: 78498295*g[v1,v2] + 89712338*g[v1,v2-1] + 1*g[v1,v2-2]
 * Basis2: 78498295*v2*g[v1,v2] + 89712338*v2*g[v1,v2-1] + 1*v2*g[v1,v2-2]
 * Basis3: 78498295*v1*g[v1,v2] + 89712338*v1*g[v1,v2-1] + 1*v1*g[v1,v2-2]
 *)

(* 理解 g[ν - α] 的含义:
 * g[v1, v2]      = g[ν1, ν2]       -> j[ν1, ν2]
 * g[v1, v2-1]   = g[ν1, ν2-1]      -> j[ν1, ν2-1]
 * g[v1, v2-2]   = g[ν1, ν2-2]      -> j[ν1, ν2-2]
 *)

Print["\n用 nu=(1,1) 代入:"];

(* Basis1: 代入 nu=(1,1) *)
nu1 = 1; nu2 = 1;
expr1 = 78498295*kiraJ[nu1, nu2] + 89712338*kiraJ[nu1, nu2-1] + 1*kiraJ[nu1, nu2-2];
Print["Basis1 at nu=(1,1): ", InputForm[expr1]];
Print["  Kira 约化: ", InputForm[kiraReduce[expr1]]];
Print["  是否为 0: ", kiraReduce[expr1] === 0];

(* Basis2: 代入 nu=(1,1), 注意 v2 因子变成数值 1 *)
expr2 = 78498295*1*kiraJ[nu1, nu2] + 89712338*1*kiraJ[nu1, nu2-1] + 1*1*kiraJ[nu1, nu2-2];
Print["Basis2 at nu=(1,1): ", InputForm[expr2]];
Print["  Kira 约化: ", InputForm[kiraReduce[expr2]]];
Print["  是否为 0: ", kiraReduce[expr2] === 0];

(* Basis3: 代入 nu=(1,1), 注意 v1 因子变成数值 1 *)
expr3 = 78498295*1*kiraJ[nu1, nu2] + 89712338*1*kiraJ[nu1, nu2-1] + 1*1*kiraJ[nu1, nu2-2];
Print["Basis3 at nu=(1,1): ", InputForm[expr3]];
Print["  Kira 约化: ", InputForm[kiraReduce[expr3]]];
Print["  是否为 0: ", kiraReduce[expr3] === 0];

Print["\n用 nu=(2,3) 代入:"];

nu1 = 2; nu2 = 3;
expr1 = 78498295*kiraJ[nu1, nu2] + 89712338*kiraJ[nu1, nu2-1] + 1*kiraJ[nu1, nu2-2];
Print["Basis1 at nu=(2,3): ", InputForm[expr1]];
Print["  Kira 约化: ", InputForm[kiraReduce[expr1]]];
Print["  是否为 0: ", kiraReduce[expr1] === 0];

expr2 = 78498295*nu2*kiraJ[nu1, nu2] + 89712338*nu2*kiraJ[nu1, nu2-1] + 1*nu2*kiraJ[nu1, nu2-2];
Print["Basis2 at nu=(2,3): ", InputForm[expr2]];
Print["  Kira 约化: ", InputForm[kiraReduce[expr2]]];
Print["  是否为 0: ", kiraReduce[expr2] === 0];

expr3 = 78498295*nu1*kiraJ[nu1, nu2] + 89712338*nu1*kiraJ[nu1, nu2-1] + 1*nu1*kiraJ[nu1, nu2-2];
Print["Basis3 at nu=(2,3): ", InputForm[expr3]];
Print["  Kira 约化: ", InputForm[kiraReduce[expr3]]];
Print["  是否为 0: ", kiraReduce[expr3] === 0];

Print["\n========================================"];
Print["结论: 如果所有结果都不是 0，说明 C++ 零空间关系有问题"];
Print["========================================"];
