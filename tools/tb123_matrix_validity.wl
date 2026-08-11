(* tb123_matrix_validity.wl — 矩阵版有效性检验
   NB>1 的 regime 中, A_i 是 nb×nb 乘法矩阵(拉平存储).
   特征多项式 p(b) 应在 B_i = Ainv_i 矩阵下给出零矩阵.
   用法: wolframscript -script tools/tb123_matrix_validity.wl TB123 4
*)

args = Rest[$ScriptCommandLine];
{family, order} = {args[[1]], ToExpression[args[[2]]]};
cppRoot = "/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP";

Get[FileNameJoin[{cppRoot, "relations", "AllRelations_" <> family <> "_k" <> ToString[order] <> ".m"}]];
Get[FileNameJoin[{cppRoot, "relations", "RelationMeta_" <> family <> ".m"}]];

p = $AllRelations[[1]]["Modulus"];
ne = $AllRelations[[1]]["NE"];
regimes = $RelationMeta["Regimes"];
sectorSums = Total[#["Sector"]] & /@ regimes;
sector = regimes[[First@Ordering[sectorSums, -1]]]["Sector"];
topRegimes = Cases[regimes, r_ /; r["Sector"] == sector];
Print["family=", family, " ne=", ne, " top regimes (all NB): ", #["NB"] & /@ topRegimes];

buildCharAssoc[entry_] := Module[
  {alphas = entry["Alphas"], betas = entry["Betas"], coefs = entry["Coefficients"],
   nA, nB, nR, w, wmax, sel},
  nA = Length[alphas]; nB = Length[betas]; nR = entry["NumRelations"];
  w = (sector . #) & /@ betas;
  wmax = Max[w];
  sel = Flatten@Position[w, wmax];
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

(* 矩阵求值: pt = B_i 矩阵列表 *)
evalAssocMat[assoc_, Bmats_, nb_] := Module[{acc = ConstantArray[0, {nb, nb}]},
  Do[
    acc += Mod[c * Fold[If[#2 === 1, #1, #1 . #2] &, IdentityMatrix[nb],
       Table[If[e[[k]] > 0, MatrixPower[Mod[Bmats[[k]], p], e[[k]]], 1], {k, ne}]], p],
    {e, Keys[assoc]}, {c, Values[assoc]}];
  Mod[acc, p]
];

blocks = Select[$AllRelations, #["NumRelations"] > 0 &];
cumul = {};
Do[
  assoc = DeleteCases[buildCharAssoc[blk], <||>];
  assocC = clearAssoc /@ assoc;
  cumul = Join[cumul, assocC];
  Print["=== (lev=", blk["Lev"], ", deg=", blk["Deg"], "): ", Length[assocC], " char polys ==="];
  Do[
    nb = rg["NB"];
    Bmats = Partition[#, nb] & /@ rg["Ainv"];
    res = Max[Max[Abs[evalAssocMat[#, Bmats, nb]]] & /@ assocC];
    Print["  regime NB=", nb, "  max |matrix residue| = ", res];
    , {rg, topRegimes}];
  , {blk, blocks}];

Print["\n=== cumulative: ", Length[cumul], " char polys ==="];
Do[
  nb = rg["NB"];
  Bmats = Partition[#, nb] & /@ rg["Ainv"];
  res = Max[Max[Abs[evalAssocMat[#, Bmats, nb]]] & /@ cumul];
  Print["  regime NB=", nb, "  max |matrix residue| = ", res];
  , {rg, topRegimes}];
Print["Done. ", DateString[]];
