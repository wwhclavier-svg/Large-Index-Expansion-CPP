(* np222_ibp_char.wl — 指纹匹配 5-prop NP222 族 (2-loop 2-point, p^2=1) + 真 IBP 特征方程
   用法: wolframscript -script tools/np222_ibp_char.wl [candIndex] *)

mmaDir = "/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA";
Get[FileNameJoin[{mmaDir, "SingularInterface.wl"}]];
Get[FileNameJoin[{mmaDir, "QuotientAlgebra.wl"}]];
Get[FileNameJoin[{mmaDir, "LargeIndexExpansion.wl"}]];
Get[FileNameJoin[{mmaDir, "BlockRecursion.wl"}]];

p = 179424673;
loopmom = {l1, l2};
extmom  = {p1};
spsRep  = {p1^2 -> 1};

(* 传播子池: 6 取 5 *)
pool = {
    l1^2,               (* 1 *)
    l2^2,               (* 2 *)
    (l1 - l2)^2,        (* 3 *)
    (l1 + p1)^2,        (* 4 *)
    (l2 + p1)^2,        (* 5 *)
    (l1 - l2 + p1)^2    (* 6 *)
};

(* --- top-sector regime A 点 (relations/RelationMeta_NP222.m, NB=1) --- *)
regimePoints = {
    {179424671, 2, 179424671, 2, 179424672},
    {179424668, 179424668, 179424668, 179424668, 5},
    {112304868, 157489941, 112304868, 157489941, 134239605},
    {67119801, 21934740, 67119801, 21934740, 45185066}
};

testCandidate[dropIdx_] := Module[
    {pdlist, ibpeqs, vlist, ibpeqsN, eig, ideal, gb, vars, match, t0},

    pdlist = Delete[pool, dropIdx];
    Print["\n===== drop pool[", dropIdx, "]: ", pdlist, " ====="];

    t0 = AbsoluteTime[];
    {ibpeqs, vlist} = Quiet[LargeIndexIBP[pdlist, loopmom, extmom, spsRep]];
    If[!ListQ[ibpeqs] || Length[ibpeqs] == 0, Print["IBP gen failed"]; Return[$Failed]];
    ibpeqsN = PolynomialMod[ibpeqs /. d -> 1/13, p];
    Print["IBP gen: ", Round[AbsoluteTime[] - t0], " s, #eqs = ", Length[ibpeqsN]];

    eig = IBPEigenEquationAB[Coefficient[#, n] & /@ ibpeqsN, 5];
    ideal = Join[eig // PolynomialMod[#, p] &, Table[A[i] B[i] - 1, {i, 5}]];
    vars = Join[Array[A, 5], Array[B, 5]];

    match = Table[
        Union[PolynomialMod[
            ideal /. Thread[vars -> Join[pt, PowerMod[pt, -1, p]]], p]],
        {pt, regimePoints}
    ];
    Print["regime residuals: ", match];
    If[Union[Flatten[match]] === {0},
        Print[">>> 指纹匹配! drop pool[", dropIdx, "]"];
        t0 = AbsoluteTime[];
        gb = GroebnerBasis[ideal, vars, Modulus -> p];
        Print["GB: ", Round[AbsoluteTime[] - t0], " s, #elts = ", Length[gb]];
        Print["GB = ", gb];
        Return[<|"drop" -> dropIdx, "pdlist" -> pdlist, "ideal" -> ideal, "gb" -> gb, "vars" -> vars|>];
    ,
        Print["no match."];
        Return[$Failed];
    ];
];

$HistoryLength = 0;
idxs = If[Length[$ScriptCommandLine] >= 2,
    {ToExpression[$ScriptCommandLine[[2]]]},
    Range[6]
];
results = Quiet[testCandidate /@ idxs];
good = DeleteCases[results, $Failed];
Print["\n===== 匹配到的候选 drop: ", (#["drop"] & /@ good), " ====="];
If[Length[good] > 0,
    Put[good[[1]], "results/intermediate/NP222_trueIBP_ideal.mx"];
    Print["saved -> results/intermediate/NP222_trueIBP_ideal.mx"];
];
