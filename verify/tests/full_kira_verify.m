(* full_kira_verify.m *)
(* 对所有非平凡 C++ 关系运行 Kira 验证 *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];

extractG[expr_] := Module[{terms, result},
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    result = {};
    Do[
        term = terms[[i]];
        termList = List @@ term;
        numG = Count[termList, _?((Head[#] === g)&)];
        Which[
            Head[term] === Times && numG > 0,
            gfound = Cases[termList, _?((Head[#] === g)&)][[1]];
            coef = Coefficient[term, gfound];
            alpha = gfound[[1]] - v1;
            beta = gfound[[2]] - v2;
            AppendTo[result, {coef, alpha, beta}],
            Head[term] === Times,
            Do[
                If[Head[termList[[j]]] === g,
                    gterm = termList[[j]];
                    coef = Times @@ DeleteCases[termList, gterm];
                    alpha = gterm[[1]] - v1;
                    beta = gterm[[2]] - v2;
                    AppendTo[result, {coef, alpha, beta}]
                ]
            , {j, Length[termList]}],
            Head[term] === g,
            alpha = term[[1]] - v1;
            beta = term[[2]] - v2;
            AppendTo[result, {1, alpha, beta}]
        ]
    , {i, Length[terms]}];
    result
];

verifyOne[relExpr_] := Module[{gCoeffs, jExpr, kr, ok},
    gCoeffs = extractG[relExpr];
    If[gCoeffs === {}, Return[True]];
    jExpr = Plus @@ (#1 * j[#2, #3] & @@@ gCoeffs);
    kr = kiraReduce[jExpr];
    ok = (kr === 0);
    If[!ok, Print["  FAIL residual: ", InputForm[kr]//Short]];
    ok
];

Print["=== Full KiraVerify for bub00 C++ Relations ===\n"];
totalPass = 0; totalFail = 0;

Do[
    file = "/root/Large-Index-Expansion-CPP/verify/bub00/Relations_bub00_lev" <> ToString[lev] <> "_deg" <> ToString[deg] <> ".m";
    If[FileExistsQ[file],
        Get[file]; rel = $RelationResult;
        rels = rel["Relations"];
        nontrivial = Select[rels, # =!= 0 &];
        If[Length[nontrivial] > 0,
            Print["Lev=", lev, " Deg=", deg, ": ", Length[nontrivial], " nontrivial relations"];
            pass = 0; fail = 0;
            Do[
                If[verifyOne[nontrivial[[i]]], pass++, fail++]
            , {i, Length[nontrivial]}];
            Print["  -> ", pass, " PASS, ", fail, " FAIL"];
            totalPass += pass; totalFail += fail;
        ]
    ]
, {lev, 0, 2}, {deg, 0, 2}];

Print["\n========================================"];
Print["OVERALL: ", totalPass, " PASS, ", totalFail, " FAIL"];
Print["Result: ", If[totalFail === 0, "ALL PASS", "SOME FAIL"]];
