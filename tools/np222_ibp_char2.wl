(* np222_ibp_char2.wl — 每个候选族: 解特征 variety, 在坐标置换下匹配 regime 点集 *)

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

(* regime b 点 (Ainv, NB=1) *)
regimeB = {
    {89712336, 89712337, 89712336, 89712337, 179424672},
    {71769869, 71769869, 71769869, 71769869, 107654804},
    {2193473, 2193474, 2193473, 2193474, 8773895},
    {69576395, 69576396, 69576395, 69576396, 98880910}
};
BV = Array[b, 5];

solveIdeal[ideal_, vars_] := Module[{gb, elim, sol},
    gb = GroebnerBasis[ideal, vars, Modulus -> p];
    sol = Solve[Thread[gb == 0], vars, Modulus -> p];
    {gb, vars /. sol}
];

testCandidate[dropIdx_] := Module[
    {pdlist, ibpeqs, vlist, ibpeqsN, eig, ideal, vars, gb, pts, perms, hit, t0},

    pdlist = Delete[pool, dropIdx];
    t0 = AbsoluteTime[];
    {ibpeqs, vlist} = Quiet[LargeIndexIBP[pdlist, loopmom, extmom, spsRep]];
    If[!ListQ[ibpeqs] || Length[ibpeqs] == 0, Print["drop ", dropIdx, ": IBP gen failed"]; Return[$Failed]];
    ibpeqsN = PolynomialMod[ibpeqs /. d -> 1/13, p];

    eig = IBPEigenEquationAB[Coefficient[#, n] & /@ ibpeqsN, 5];
    ideal = Join[eig // PolynomialMod[#, p] &, Table[A[i] B[i] - 1, {i, 5}]];
    vars = Join[Array[A, 5], Array[B, 5]];

    {gb, pts} = solveIdeal[ideal, vars];
    (* b-坐标 = 后 5 个分量 *)
    bpts = (#[[6 ;; 10]] & /@ pts) // DeleteDuplicates;
    Print["\n===== drop pool[", dropIdx, "]: ", pdlist];
    Print["  #rational AB-points: ", Length[pts], "  distinct b-pts: ", Length[bpts], "  (", Round[AbsoluteTime[] - t0], " s)"];
    Print["  b-points: ", bpts];

    (* 在 120 个坐标置换下匹配 regime b 点集 *)
    perms = Permutations[Range[5]];
    hit = SelectFirst[perms,
        Function[perm,
            And @@ Table[MemberQ[bpts, Mod[pt[[perm]], p]], {pt, regimeB}]
        ], None];
    If[hit =!= None,
        Print["  >>> 匹配! 置换: ", hit];
        Print["  即 bin 坐标 i = 候选坐标 ", hit];
        Return[<|"drop" -> dropIdx, "pdlist" -> pdlist, "perm" -> hit,
                 "ideal" -> ideal, "vars" -> vars, "gb" -> gb, "bpts" -> bpts|>];
    ,
        Print["  no match under permutation."];
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
    Print["saved."];
];
