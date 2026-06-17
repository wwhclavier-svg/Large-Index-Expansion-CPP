(* verify_mma_kira.m *)
(* 验证 MMA LIE 关系是否通过 Kira 约化规则 *)

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

(* 加载 MMA 关系（rank=2, maxDeg=2 版本） *)
Get["/root/Large-Index-Expansion-CPP/verify/bub00/Compare-MMARelation-lev2.m"];

(* $MMARelations = {{{deg0, deg1, deg2(level1)}, {deg0, deg1, deg2(level2)}}} *)
mmaRels = $MMARelations;
Print["MMA 关系结构: ", Dimensions[mmaRels]];
Print["Level 1 deg1: ", Length[mmaRels[[1,2]]], " relations"];
Print["Level 1 deg2: ", Length[mmaRels[[1,3]]], " relations"];
Print["Level 2 deg2: ", Length[mmaRels[[2,3]]], " relations"];

(* 收集所有非平凡关系 *)
allRelations = {};
Do[
    If[Head[mmaRels[[l, d+1]]] === List,
        Do[
            rel = mmaRels[[l, d+1, r]];
            If[rel =!= 0 && rel =!= {},
                AppendTo[allRelations, {rel, l, d}]
            ],
        {r, Length[mmaRels[[l, d+1]]]}]
    ],
{l, 1, 2}, {d, 0, 2}
];
Print["\n非平凡关系总数: ", Length[allRelations]];

(* 解析函数：从 g[ν-α] 中提取 α *)
extractAlphaFromGArg[arg_] := Module[
    {coef, const},
    Which[
        Head[arg] === Symbol,
            Return[0],
        Head[arg] === Plus,
            If[MemberQ[{v1, v2}, arg[[1]]],
                (* v1 - k 格式 *)
                coef = arg[[1]]; const = arg[[2]]; Return[-const],
                (* -k + v1 格式 *)
                coef = arg[[2]]; const = arg[[1]]; Return[-const]
            ],
        Head[arg] === Integer,
            Return[-arg],
        True,
            Return[0]
    ];
    0
];

parseMMAArgs[gargs_] := Module[
    {alpha1, alpha2},
    alpha1 = extractAlphaFromGArg[gargs[[1]]];
    alpha2 = extractAlphaFromGArg[gargs[[2]]];
    {alpha1, alpha2}
];

(* 解析 MMA 关系项：提取 {coeff, alpha, beta_v1, beta_v2} *)
parseMMATerm[term_] := Module[
    {factors, coef, v1Pow, v2Pow, gPart, gArgs, alpha},
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

(* 验证：代入 ν 值后 Kira 约化 *)
verifyMMAAtNu[relExpr_, nu1_, nu2_] := Module[
    {terms, jExpr, result, reduced},
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    jExpr = 0;
    Do[
        parsed = parseMMATerm[terms[[i]]];
        If[parsed === Null, Continue[]];
        {coef, alpha, b1, b2} = parsed;
        jArg1 = nu1 - alpha[[1]];
        jArg2 = nu2 - alpha[[2]];
        nuBeta = nu1^b1 * nu2^b2;
        termValue = coef * nuBeta * kiraJ[jArg1, jArg2];
        jExpr = jExpr + termValue;
    , {i, Length[terms]}];
    reduced = kiraReduce[jExpr];
    {reduced === 0, jExpr, reduced}
];

(* 测试点 *)
testPoints = {{1, 1}, {1, 2}, {2, 1}, {2, 3}};

Print["\n================================================"];
Print["MMA LIE Relations : Kira 验证"];
Print["================================================"];

totalPass = 0; totalFail = 0;
Do[
    {relExpr, lv, dg} = allRelations[[ri]];
    
    (* 提取 α 列表 *)
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    alphasUsed = {};
    Do[
        parsed = parseMMATerm[terms[[i]]];
        If[parsed =!= Null, AppendTo[alphasUsed, parsed[[2]]]],
    {i, Length[terms]}];
    alphasUsed = Union[alphasUsed];
    
    Print["\n--- 关系 ", ri, " (level=", lv, ", deg=", dg, ") ---"];
    Print["  α 值: ", alphasUsed];
    Print["  Short: ", Short[InputForm[relExpr], 4]];
    
    passCount = 0; failCount = 0;
    Do[
        {nu1, nu2} = testPoints[[tp]];
        Print["\n  ν=(", nu1, ",", nu2, "):"];
        
        (* 显示解析的项 *)
        Do[
            parsed = parseMMATerm[terms[[i]]];
            If[parsed =!= Null,
                {coef, alpha, b1, b2} = parsed;
                Print["    项: coeff=", coef, ", α=", alpha, ", β=(", b1, ",", b2, ")"];
            ];
        , {i, Length[terms]}];
        
        {result, jExpr, reduced} = verifyMMAAtNu[relExpr, nu1, nu2];
        status = If[result, "PASS", "FAIL"];
        Print["    代入后: ", Short[InputForm[jExpr], 3]];
        Print["    Kira 后: ", Short[InputForm[reduced], 3]];
        Print["    --> ", status];
        If[result, passCount++, failCount++];
    , {tp, 1, Length[testPoints]}];
    
    Print["  --> ", passCount, "/", Length[testPoints], " PASS"];
    totalPass += passCount;
    totalFail += failCount;
, {ri, 1, Length[allRelations]}];

Print["\n================================================"];
Print["总计: ", totalPass, " PASS, ", totalFail, " FAIL"];
Print["================================================"];
