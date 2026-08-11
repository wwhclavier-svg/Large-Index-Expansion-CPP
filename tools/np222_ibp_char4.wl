(* np222_ibp_char4.wl — 扩大传播子池, 代入检验 (perm × sgn × A<->B swap) 匹配 regime 点 *)

mmaDir = "/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA";
Get[FileNameJoin[{mmaDir, "SingularInterface.wl"}]];
Get[FileNameJoin[{mmaDir, "QuotientAlgebra.wl"}]];
Get[FileNameJoin[{mmaDir, "LargeIndexExpansion.wl"}]];
Get[FileNameJoin[{mmaDir, "BlockRecursion.wl"}]];

p = 179424673;
loopmom = {l1, l2};
extmom  = {p1};
spsRep  = {p1^2 -> 1};

pool = {
    l1^2, l2^2, (l1 - l2)^2, (l1 + l2)^2,
    (l1 + p1)^2, (l2 + p1)^2, (l1 - l2 + p1)^2, (l1 - l2 - p1)^2
};

(* regime b 点 (Ainv, NB=1); R1 有多重集结构用于快速筛查 *)
R1 = {89712336, 89712337, 89712336, 89712337, 179424672};
regimeB = {
    R1,
    {71769869, 71769869, 71769869, 71769869, 107654804},
    {2193473, 2193474, 2193473, 2193474, 8773895},
    {69576395, 69576396, 69576395, 69576396, 98880910}
};

testCandidate[cand_] := Module[
    {pdlist, ibpeqs, vlist, ibpeqsN, eig, ideal, vars, perms, evalAt, hits, fullHit},

    pdlist = pool[[cand]];
    {ibpeqs, vlist} = Quiet[LargeIndexIBP[pdlist, loopmom, extmom, spsRep]];
    If[!ListQ[ibpeqs] || Length[ibpeqs] == 0, Return[$Failed]];
    ibpeqsN = Quiet[PolynomialMod[ibpeqs /. d -> 1/13, p]];
    If[!ListQ[ibpeqsN], Return[$Failed]];

    eig = Quiet[IBPEigenEquationAB[Coefficient[#, n] & /@ ibpeqsN, 5]];
    If[!ListQ[eig], Return[$Failed]];
    ideal = Join[eig // PolynomialMod[#, p] &, Table[A[i] B[i] - 1, {i, 5}]];
    vars = Join[Array[A, 5], Array[B, 5]];

    (* b=target 代入: (A=1/t, B=t); swap: (A=t, B=1/t) *)
    evalAt[t_] := Union[PolynomialMod[ideal /. Thread[vars -> Join[PowerMod[t, -1, p], t]], p]];
    evalAtSwap[t_] := Union[PolynomialMod[ideal /. Thread[vars -> Join[t, PowerMod[t, -1, p]]], p]];

    (* R1 快速筛查: 30 distinct perms × sgn × swap *)
    perms = DeleteDuplicates[Permutations[R1]] /. R1[[1]] -> 0; (* placeholder 防错 *)
    perms = DeleteDuplicates[Permutations[{1, 1, 2, 2, 3}]];
    hits = {};
    Do[
        rp = R1[[perm]];
        Do[
            t = Mod[sgn rp, p];
            If[evalAt[t] === {0}, AppendTo[hits, {perm, sgn, "norm"}]; Break[]];
            If[evalAtSwap[t] === {0}, AppendTo[hits, {perm, sgn, "swap"}]; Break[]];
        , {sgn, {1, -1}}];
    , {perm, perms}];
    If[hits === {}, Return[$Failed]];

    (* 全量检验: 对命中的 (perm,sgn,mode) 验 4 个点 *)
    fullHit = SelectFirst[hits,
        Function[h,
            And @@ Table[
                rp = Mod[h[[2]] pt[[h[[1]]]], p];
                If[h[[3]] === "norm", evalAt[rp] === {0}, evalAtSwap[rp] === {0}],
            {pt, regimeB}]
        ], None];
    If[fullHit =!= None,
        Print[">>> MATCH cand=", cand, "  props=", pdlist, "  hit=", fullHit];
        Return[<|"cand" -> cand, "pdlist" -> pdlist, "hit" -> fullHit,
                 "ideal" -> ideal, "vars" -> vars|>];
    ,
        Print["partial (R1 only) cand=", cand];
        Return[$Failed];
    ];
];

$HistoryLength = 0;
cands = Subsets[Range[8], {5}];
Print["#candidates: ", Length[cands]];
results = Table[
    Print["[", i, "/", Length[cands], "] cand = ", cands[[i]]];
    testCandidate[cands[[i]]],
    {i, Length[cands]}];
good = DeleteCases[results, $Failed];
Print["\n===== 匹配候选: ", (#["cand"] & /@ good), " ====="];
If[Length[good] > 0,
    Put[good[[1]], "results/intermediate/NP222_trueIBP_ideal.mx"];
    Print["saved."];
];
