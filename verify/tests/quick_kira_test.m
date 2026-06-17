(* quick_kira_test.m *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];

Get["/root/Large-Index-Expansion-CPP/verify/bub00/Relations_bub00_lev2_deg1.m"];
rel = $RelationResult;
relations = rel["Relations"];
Print["NumSolutions: ", rel["NumSolutions"]];

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

Print["\n=== KiraVerify Results ==="];
passCount = 0; failCount = 0;
Do[
    relExpr = relations[[i]];
    If[relExpr =!= 0,
        gCoeffs = extractG[relExpr];
        jExpr = Plus @@ (#1 * j[#2, #3] & @@@ gCoeffs);
        kiraReduced = kiraReduce[jExpr];
        status = If[kiraReduced === 0, "PASS", "FAIL"];
        If[status === "PASS", passCount++, failCount++];
        Print["Relation ", i, ": ", status, If[status === "FAIL", " (res: " <> ToString[InputForm[kiraReduced]]<> ")", ""]]
    ]
, {i, Length[relations]}];
Print["\nSummary: ", passCount, " PASS, ", failCount, " FAIL out of ", Length[Select[relations, # =!= 0 &]], " non-trivial relations"];
