(* np222_ibp_char3.wl — 代入检验: 符号翻转 + 坐标置换下匹配 regime 点; 命中后算局部化 vdim *)

mmaDir = "/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA";
Get[FileNameJoin[{mmaDir, "SingularInterface.wl"}]];
Get[FileNameJoin[{mmaDir, "QuotientAlgebra.wl"}]];
Get[FileNameJoin[{mmaDir, "LargeIndexExpansion.wl"}]];
Get[FileNameJoin[{mmaDir, "BlockRecursion.wl"}]];

p = 179424673;
loopmom = {l1, l2};
extmom  = {p1};
spsRep  = {p1^2 -> 1};

pool = {l1^2, l2^2, (l1 - l2)^2, (l1 + p1)^2, (l2 + p1)^2, (l1 - l2 + p1)^2};

regimeB = {
    {89712336, 89712337, 89712336, 89712337, 179424672},
    {71769869, 71769869, 71769869, 71769869, 107654804},
    {2193473, 2193474, 2193473, 2193474, 8773895},
    {69576395, 69576396, 69576395, 69576396, 98880910}
};

(* NB=2 regime: A 矩阵 (5 个 2x2), 取其特征值需逐变量; 先用 b5 的纯量值 *)
nb2b5 = 139552523;

testCandidate[dropIdx_] := Module[
    {pdlist, ibpeqs, vlist, ibpeqsN, eig, ideal, vars, perms, hit, sgn, t0, evalAt},

    pdlist = Delete[pool, dropIdx];
    {ibpeqs, vlist} = Quiet[LargeIndexIBP[pdlist, loopmom, extmom, spsRep]];
    If[!ListQ[ibpeqs] || Length[ibpeqs] == 0, Print["drop ", dropIdx, ": IBP gen failed"]; Return[$Failed]];
    ibpeqsN = PolynomialMod[ibpeqs /. d -> 1/13, p];

    eig = IBPEigenEquationAB[Coefficient[#, n] & /@ ibpeqsN, 5];
    ideal = Join[eig // PolynomialMod[#, p] &, Table[A[i] B[i] - 1, {i, 5}]];
    vars = Join[Array[A, 5], Array[B, 5]];

    (* 在候选 ideal 中代入: b 点 = sgn * pt[[perm]]; A = 1/b *)
    evalAt[pt_] := Union[PolynomialMod[
        ideal /. Thread[vars -> Join[PowerMod[pt, -1, p], pt]], p]];

    perms = Permutations[Range[5]];
    hit = None;
    Do[
        hit = SelectFirst[perms,
            Function[perm, And @@ Table[evalAt[Mod[sgn pt[[perm]], p]] === {0}, {pt, regimeB}]],
            None];
        If[hit =!= None, Break[]];
    , {sgn, {1, -1}}];

    Print["\n===== drop pool[", dropIdx, "]: ", pdlist];
    If[hit =!= None,
        Print["  >>> 4 个 NB=1 regime 点全部满足 (置换 ", hit, ", sgn=", sgn, ")"];
        Return[<|"drop" -> dropIdx, "pdlist" -> pdlist, "perm" -> hit, "sgn" -> sgn,
                 "ideal" -> ideal, "vars" -> vars|>];
    ,
        Print["  no match."];
        Return[$Failed];
    ];
];

$HistoryLength = 0;
idxs = If[Length[$ScriptCommandLine] >= 2, {ToExpression[$ScriptCommandLine[[2]]]}, Range[6]];
results = testCandidate /@ idxs;
good = DeleteCases[results, $Failed];
Print["\n===== 匹配候选 drop: ", (#["drop"] & /@ good), " ====="];
If[Length[good] > 0,
    Put[good[[1]], "results/intermediate/NP222_trueIBP_ideal.mx"];
    Print["saved -> results/intermediate/NP222_trueIBP_ideal.mx"];
];
