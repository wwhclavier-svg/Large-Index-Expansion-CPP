(* np222_variety_points.wl — 解出 NP222 特征 variety 的 F_p 有理点, 与 regime 点比对 *)
cppRoot = "/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP";
relDir = FileNameJoin[{cppRoot, "relations_l4d3"}];
Get[FileNameJoin[{relDir, "AllRelations_NP222_k4.m"}]];
Get[FileNameJoin[{relDir, "RelationMeta_NP222.m"}]];

p = $AllRelations[[1]]["Modulus"];
ne = $AllRelations[[1]]["NE"];
sector = {1,1,1,1,1};
BV = Array[ToExpression["b" <> ToString[#]] &, ne];

buildCharAssoc[entry_] := Module[
  {alphas = entry["Alphas"], betas = entry["Betas"], coefs = entry["Coefficients"],
   nA, nB, nR, w, wmax, sel},
  nA = Length[alphas]; nB = Length[betas]; nR = entry["NumRelations"];
  w = (sector . #) & /@ betas; wmax = Max[w]; sel = Flatten@Position[w, wmax];
  Table[
    Merge[
      DeleteCases[
        Flatten@Table[
          If[coefs[[(i-1)*nB+j, r+1]] =!= 0, alphas[[i]] -> coefs[[(i-1)*nB+j, r+1]], Nothing],
          {j, sel}, {i, nA}],
        _ -> 0],
      Total],
    {r, nR}]
];
clearAssoc[assoc_] := Module[{keys = Keys[assoc], mins},
  If[Length[keys] == 0, Return[assoc]];
  mins = Min /@ Transpose[keys];
  KeyMap[# - mins &, assoc]
];
assocToPoly[assoc_] := Plus @@ KeyValueMap[#2 * Times @@ (BV^#1) &, assoc];

blocks = Select[$AllRelations, #["NumRelations"] > 0 &];
Print["blocks: ", {#["Lev"], #["Deg"], #["NumRelations"]} & /@ blocks];
assocC = clearAssoc /@ Flatten[DeleteCases[buildCharAssoc[#], <||>] & /@ blocks];
polys = DeleteCases[assocToPoly /@ assocC, 0];
Print["total char polys: ", Length[polys]];

Print["GB... ", DateString[]];
gb = GroebnerBasis[polys, BV, Modulus -> p];
Print["GB size: ", Length[gb], "  ", DateString[]];
Print["GB: ", gb];

Print["Solve... ", DateString[]];
sol = Solve[Thread[gb == 0], BV, Modulus -> p];
Print["rational points over F_p: ", Length[sol], "  ", DateString[]];

(* regime 点: NB=1 直接取 Ainv; NB=2 求特征值 *)
regimes = $RelationMeta["Regimes"];
nb1pts = Cases[regimes, r_ /; r["NB"] == 1 :> r["Ainv"][[;;,1]]];
nb2 = Cases[regimes, r_ /; r["NB"] == 2];
Print["NB=1 regime points: ", Length[nb1pts], "  NB=2 regimes: ", Length[nb2]];

(* NB=2: 每个变量的 2x2 Ainv 矩阵的特征值 (mod p) *)
If[Length[nb2] > 0,
  m1 = nb2[[1]]["Ainv"];
  Do[
    mat = Partition[m1[[i]], 2];
    cp = CharacteristicPolynomial[mat, x, Modulus -> p];
    fac = FactorList[cp, Modulus -> p];
    Print["  NB=2 var b", i, ": charpoly factors: ", fac];
    , {i, ne}];
];

(* 比对: 每个解点是否等于某个 NB=1 点 *)
pts = BV /. sol;
memberFlag = Table[
  pt = pts[[i]];
  hits = Position[nb1pts, q_ /; And @@ Table[Mod[pt[[k]] - q[[k]], p] == 0, {k, ne}], 1];
  {i, pt, If[hits === {}, "NEW", "regime-" <> ToString[hits[[1,1]]]]},
  {i, Length[pts]}];
Print[""];
Do[Print[r[[1]], "  b=", r[[2]], "  -> ", r[[3]]], {r, memberFlag}];
