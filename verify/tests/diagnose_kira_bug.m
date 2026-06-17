(* diagnose_kira_bug.m *)
(* 独立诊断: 验证 bub00 0/10 KiraVerify FAIL 的根因 *)
(* 快速运行: 只加载 C++ 关系 + Kira 规则，测 1 个关系 *)
(* 不依赖完整的 MMA LIEReconstruct 工作流 *)

AppendTo[$Path, "/root/Large-Index-Expansion-MMA-Mini"];

Print["=== KiraVerify Bug Diagnostic ==="];

(* Step 1: Load Kira rules *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];
Print["Kira rules loaded: ", Length[kiraRules], " rules"];

(* Step 2: Load one C++ relation (lev=0, deg=0, rel=1) *)
cppRel = Get["/root/Large-Index-Expansion-CPP/verify/bub00/Relations_bub00_lev0_deg0.m"];
relation = cppRel[["Relations", 1, 1, 1]];
Print["\nC++ relation (lev=0, deg=0, rel=1):"];
Print["  InputForm: ", InputForm[relation] // Short];

(* Step 3: Extract g[...] terms *)
extractG[expr_] := Module[{terms, result},
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    result = {};
    Do[
        term = terms[[i]];
        termList = List @@ term;
        numG = Count[termList, _?((ToString[Head[#]] === "g")&)];
        Which[
            Head[term] === Times && numG > 0,
            gfound = Cases[termList, _?((ToString[Head[#]] === "g")&)][[1]];
            coef = Coefficient[term, gfound];
            alpha = gfound[[1]] - v1;
            beta = gfound[[2]] - v2;
            AppendTo[result, {coef, alpha, beta}],
            Head[term] === Times,
            Do[
                If[ToString[Head[termList[[j]]]] === "g",
                    gterm = termList[[j]];
                    coef = Times @@ DeleteCases[termList, gterm];
                    alpha = gterm[[1]] - v1;
                    beta = gterm[[2]] - v2;
                    AppendTo[result, {coef, alpha, beta}]
                ]
            , {j, Length[termList]}],
            ToString[Head[term]] === "g",
            alpha = term[[1]] - v1;
            beta = term[[2]] - v2;
            AppendTo[result, {1, alpha, beta}]
        ]
    , {i, Length[terms]}];
    result
];

gCoeffs = extractG[relation];
Print["\ng[...] terms extracted: ", Length[gCoeffs]];
Print["  Sample: ", InputForm[gCoeffs[[1]]] // Short];

(* Step 4: Verify with CORRECT method (using kiraReduce) *)
jExpr = Plus @@ (#1 * j[{#2, #3}] & @@@ gCoeffs);
Print["\nj[...] expression: ", InputForm[jExpr] // Short];
kiraReduced = kiraReduce[jExpr];
Print["kiraReduce result: ", InputForm[kiraReduced] // Short];
Print["CORRECT (using kiraReduce): ", If[kiraReduced === 0, "PASS", "FAIL"]];

(* Step 5: Verify with BUGGY method (using ReplaceAll with rule list) *)
jExprs = gCoeffs /. {c_, a_, b_} :> c * j[{a, b}];
jVectors = jExprs /. kiraRules;
result = Plus @@ jVectors;
result = PolynomialMod[result, 179424673];
Print["\nBUGGY (using ReplaceAll): ", InputForm[result] // Short];
Print["BUGGY result === 0? ", result === 0];

Print["\n=== Diagnostic Complete ==="];
