(* sr212_true_ab.wls — SR212 真 IBP 特征系统：regime 验证 + primdec *)
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LargeIndexExpansion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BlockRecursion.wl"];

char = 179424673;
pdlist = {l1^2, (l1+p)^2, l2^2, (l2+p)^2, (l1+l2+p)^2};
loopmom = {l1, l2}; extmom = {p}; spsRep = {p^2 -> 1};

{ibpeqsfull, vlist} = LargeIndexIBP[pdlist, loopmom, extmom, spsRep];
ibpeqs = ibpeqsfull[[1 ;; 6]] /. d -> 1/13 // PolynomialMod[#, char] &;
Print["#ibpeqs=", Length[ibpeqs]];

(* 特征方程（expRegSolve2 同款构造） *)
Alist = Array[A, 5];
aeqs0 = Coefficient[ibpeqs, n] /. g[a__] :> (Times @@ Thread@Power[Alist, {a} - vlist]);
aeqs = Join[Together /@ aeqs0 // Numerator, {Times @@ Alist * f - 1}] // PolynomialMod[#, char] &;
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/tools/sr212_true_aeqs.m", aeqs, "Package"];
Print["#aeqs=", Length[aeqs]];

(* meta top-sector regime A 点代入检验 *)
regs = {
  {179424671, 2, 179424671, 2, 179424672},
  {179424668, 179424668, 179424668, 179424668, 5},
  {112304868, 157489941, 112304868, 157489941, 134239605},
  {67119801, 21934740, 67119801, 21934740, 45185066}
};
Do[
  r = regs[[ri]];
  sub = Thread[Join[Alist, {f}] -> Append[r, PowerMod[Times @@ r, -1, char]]];
  res = PolynomialMod[aeqs /. sub, char];
  Print["reg", ri, " residues=", res];
, {ri, 1, Length[regs]}];
Print["DONE"];
