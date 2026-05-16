(* Compare equation matrices: correctly handling MMA expansion format *)
$ProjectDir = "/home/ykm/Large-Index-Expansion-CPP";

Get[FileNameJoin[{$ProjectDir, "families/FamilyDatabase.wl"}]];
config = $FamilyDatabase["bub00"];
modulus = config["Modulus"];
ne = 2; order = 4;
vlist = Table[Symbol["v" <> ToString[i]], {i, ne}];

(* Load MMA expansion data (format: {<|"Solutions"->{<|"H"->{{h0},{h1},...}|>}|>} *)
mmaRaw = Get[FileNameJoin[{$ProjectDir, "workspace/playground/ExpansionMMA_bub00.m"}]];
hRaw = mmaRaw[[1, "Solutions", 1, "H"]];
(* hRaw[k+1] is {coeff} — unwrap *)
hlist = Flatten /@ hRaw;  
hExpr[nu_List, alpha_List] := Sum[
  (hlist[[k+1]] /. {ToString[v1] -> (nu[[1]] - alpha[[1]]), ToString[v2] -> (nu[[2]] - alpha[[2]])}) / n^k,
  {k, 0, Length[hlist]-1}];

(* Load relation meta *)
Get[FileNameJoin[{$ProjectDir, "relations/RelationMeta_bub00.m"}]];
aVals = $RelationMeta["Regimes"][[1, "A"]];
aInvVals = $RelationMeta["Regimes"][[1, "Ainv"]];
theta = $RelationMeta["Regimes"][[1, "Sector"]];

(* Generate alphas/betas for (1,1) *)
alphas = Flatten[Table[{a1, a2}, {a1, 0, 1}, {a2, 0, 1 - a1}], 1];
betas = Flatten[Table[{b1, b2}, {b1, 0, 1}, {b2, 0, 1 - b1}], 1];
nAlpha = Length[alphas]; nBeta = Length[betas]; kMax = order;
nu = {1, 0};
degEff = If[AnyTrue[theta, # =!= 0 &], 1, 0];
nb = 1;

(* Build MMA equation matrix *)
M_mma = ConstantArray[0, {nb * (kMax + 1), nAlpha * nBeta}];

Do[
  alpha = alphas[[ai]]; beta = betas[[bi]];
  col = (ai - 1) * nBeta + bi;
  
  (* A^{-alpha} — handle zero A values *)
  aFac = 1;
  Do[
    If[aVals[[i]] == 0,
      If[alpha[[i]] == 0, aFac *= 1, aFac *= 0],
      aFac *= PowerMod[aVals[[i]], -alpha[[i]], modulus]
    ];
  , {i, ne}];
  
  (* (nu + theta)^beta *)
  thetaBeta = Times @@ ((nu + theta)^beta);
  
  (* h(nu - alpha) *)
  hVal = hExpr[nu, alpha];
  
  (* Full term: (theta+nu)^beta * A^{-alpha} * h(nu-alpha) *)
  term = thetaBeta * aFac * hVal;
  
  (* Expand in 1/n *)
  series = CoefficientList[term, 1/n, kMax + 1];
  For[row = 1, row <= kMax + 1, row++,
    If[row <= Length[series],
      M_mma[[row, col]] = PolynomialMod[series[[row]], modulus];
    ];
  ];
, {ai, nAlpha}, {bi, nBeta}];

Print["MMA eqn matrix (" <> ToString[Length[M_mma]] <> " x " <> ToString[Length[M_mma[[1]]]] <> "):"];
Print[MatrixForm[M_mma]];

(* Load C++ matrix *)
Get[FileNameJoin[{$ProjectDir, "debug_eqn_matrix.m"}]];
M_cpp = $EqnMatrix;

diff = M_mma - M_cpp;
Print["\nDifference (MMA - C++):"];
Print[MatrixForm[diff]];

If[Max[Abs[Flatten[diff]]] === 0,
  Print["\n*** MATRICES MATCH! ***"],
  Print["\n*** MATRICES DIFFER (max diff = ", Max[Abs[Flatten[diff]]], ") ***"]
];
