(* Quick KiraVerify for lev=2, deg=2 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];
kiraJ = KiraRuleLoader`j;

(* Load C++ relations *)
Get["/root/Large-Index-Expansion-CPP/Relations_bub00_lev2_deg2.m"];
rel = $RelationResult;
Print["lev=2, deg=2: NumSolutions=", rel["NumSolutions"]];

relations = rel["Relations"];
nontrivial = Select[relations, # =!= 0 &];
Print["Non-trivial relations: ", Length[nontrivial]];

(* Simple term parser: coeff * v1^b1 * v2^b2 * g[...] *)
parseTerm[term_] := Module[
    {factors, gPart, coef, b1, b2, gArgs, a1, a2, f},
    factors = If[Head[term] === Times, List @@ term, {term}];
    gPart = First[Cases[factors, _g]];
    coef = Times @@ DeleteCases[factors, _g];
    
    (* beta powers *)
    b1 = 0; b2 = 0;
    If[Cases[factors, v1] =!= {}, b1 = 1];
    If[Cases[factors, v2] =!= {}, b2 = 1];
    Do[
        f = factors[[i]];
        If[Head[f] === Power,
            If[f[[1]] === v1, b1 = f[[2]]];
            If[f[[1]] === v2, b2 = f[[2]]];
        ];
    , {i, Length[factors]}];
    
    (* alpha from g args *)
    gArgs = List @@ gPart;
    a1 = 0; a2 = 0;
    If[Head[gArgs[[1]]] === Plus, a1 = -gArgs[[1, 2]]];
    If[Head[gArgs[[2]]] === Plus, a2 = -gArgs[[2, 2]]];
    If[gArgs[[1]] === v1, a1 = 0];
    If[gArgs[[2]] === v2, a2 = 0];
    If[Head[gArgs[[1]]] === Times && gArgs[[1, 1]] === -1, a1 = gArgs[[1, 2]]];
    If[Head[gArgs[[2]]] === Times && gArgs[[2, 1]] === -1, a2 = gArgs[[2, 2]]];
    
    {coef, {a1, a2}, {b1, b2}}
];

verifyAtNu[relExpr_, nu1_, nu2_] := Module[
    {terms, jVal, reduced, coef, alpha, beta, j1, j2},
    terms = If[Head[relExpr] === Plus, List @@ relExpr, {relExpr}];
    jVal = 0;
    Do[
        {coef, alpha, beta} = parseTerm[terms[[t]]];
        j1 = nu1 - alpha[[1]];
        j2 = nu2 - alpha[[2]];
        jVal += coef * nu1^beta[[1]] * nu2^beta[[2]] * kiraJ[j1, j2];
    , {t, Length[terms]}];
    reduced = kiraReduce[jVal];
    reduced === 0
];

testPoints = {{1,1}, {2,3}, {5,7}};
totalPass = 0; totalFail = 0;

Do[
    relExpr = nontrivial[[r]];
    Do[
        {nu1, nu2} = testPoints[[tp]];
        result = verifyAtNu[relExpr, nu1, nu2];
        If[result, totalPass++, totalFail++];
    , {tp, Length[testPoints]}];
, {r, Length[nontrivial]}];

Print[""];
Print["========================================"];
Print["OVERALL: ", totalPass, " PASS, ", totalFail, " FAIL"];
Print["Total checks: ", Length[nontrivial] * Length[testPoints]];
Print["========================================"];
