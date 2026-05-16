(* Compare C++ and MMA equation matrices for bub00 *)
$ProjectDir = "/home/ykm/Large-Index-Expansion-CPP";

Get[FileNameJoin[{$ProjectDir, "families/FamilyDatabase.wl"}]];
config = $FamilyDatabase["bub00"];
mod = config["Modulus"];
ne = 2; order = 4;
vv = Table[Symbol["v" <> ToString[i]], {i, ne}];

(* Load MMA expansion data *)
mmaRaw = Get[FileNameJoin[{$ProjectDir, "workspace/playground/ExpansionMMA_bub00.m"}]];
hRaw = mmaRaw[[1, "Solutions", 1, "H"]];
hlist2 = Flatten /@ hRaw;

(* Load relation meta *)
Import[FileNameJoin[{$ProjectDir, "relations/RelationMeta_bub00.m"}]];
aV = $RelationMeta["Regimes"][[1, "A"]];
th = $RelationMeta["Regimes"][[1, "Sector"]];

(* alphas/betas for (1,1) *)
al = Flatten[Table[{a1, a2}, {a1, 0, 1}, {a2, 0, 1 - a1}], 1];
be = Flatten[Table[{b1, b2}, {b1, 0, 1}, {b2, 0, 1 - b1}], 1];
nA = Length[al]; nB = Length[be]; km = order;
nu2 = {1, 0};
dE = 1;
nb2 = 1;

(* Build MMA equation matrix *)
mm = ConstantArray[0, {nb2 * (km + 1), nA * nB}];

Do[
  al0 = al[[ai]]; be0 = be[[bi]];
  col = (ai - 1) * nB + bi;
  
  (* A^{-alpha} *)
  af = 1;
  Do[
    If[aV[[i]] == 0, If[al0[[i]] == 0, Null, af *= 0], af *= PowerMod[aV[[i]], -al0[[i]], mod]];
  , {i, ne}];
  If[af === 0, Continue[]];
  
  (* (theta+nu)^beta *)
  tb = Times @@ ((nu2 + th)^be0);
  
  (* h(nu-alpha) *)
  hv = Sum[
    (hlist2[[k+1]] /. {ToString[v1] -> (nu2[[1]] - al0[[1]]), ToString[v2] -> (nu2[[2]] - al0[[2]])}) / n^k,
    {k, 0, km}];
  
  term = tb * af * hv;
  ser = CoefficientList[term, 1/n, km + 1];
  Do[
    If[row <= Length[ser], mm[[row, col]] = PolynomialMod[ser[[row]], mod]];
  , {row, km + 1}];
, {ai, nA}, {bi, nB}];

Print["MMA eqn matrix: ", Dimensions[mm]];
Print[MatrixForm[mm]];

(* Load C++ matrix *)
Get[FileNameJoin[{$ProjectDir, "debug_eqn_matrix.m"}]];
mc = $EqnMatrix;

diff = mm - mc;
Print["Difference (MMA - C++):"];
Print[MatrixForm[diff]];
Print["Max diff: ", Max[Abs[Flatten[diff]]]];
If[Max[Abs[Flatten[diff]]] === 0, Print["*** MATRICES MATCH! ***"], Print["*** MATRICES DIFFER ***"]];
