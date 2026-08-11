Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LargeIndexExpansion.wl"];
char = 179424673;
pdlist = {l1^2, (l1-p1)^2, (l1-p1-p2)^2, (l2+p1+p2)^2, l2^2, (l1+l2)^2, (l2+p1)^2};
loopmom = {l1, l2}; extmom = {p1, p2};
spsRep = {p1^2 -> 2, p2^2 -> 4, p1*p2 -> -3/2};
topsector = {1,1,1,1,1,1,0};
{ne, nl, nE, nibp, Alist, vlist, ibpeqs1, sectorlist} =
  AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topsector];
Print["ne=", ne, " nibp=", nibp, " Alist=", Alist, " vlist=", vlist];
ibpeqssub = sectorLimitIBP[PolynomialMod[ibpeqs1 /. d -> 1/13, char], topsector, vlist];
aeqs0 = Coefficient[ibpeqssub, n] /. g[a__] :> (Times @@ Thread@Power[Alist, {a} - vlist]);
aeqs = Join[aeqs0 /. {1/A[a_] :> B[a]}, Table[A[i] B[i] - 1, {i, 7}]] // PolynomialMod[#, char] &;
toSing[e_] := ToString[e, InputForm] // StringReplace[#, {
    "A[" ~~ i : DigitCharacter ~~ "]" :> "A" <> i,
    "B[" ~~ i : DigitCharacter ~~ "]" :> "B" <> i}] &;
sing = StringRiffle[toSing /@ aeqs, ",\n"];
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/tools/tb123_true_ideal.txt", sing, "Text"];
Print["ideal exported, #eqs=", Length[aeqs]];
