(* Compare C++ and MMA equation matrices for bub00 (lev=1, deg=1) *)
$ProjectDir = "/home/ykm/Large-Index-Expansion-CPP";

Get[FileNameJoin[{$ProjectDir, "families/FamilyDatabase.wl"}]];
config = $FamilyDatabase["bub00"];
modulus = config["Modulus"];
ne = 2; order = 4;
vlist = Table[Symbol["v" <> ToString[i]], {i, ne}];

(* Load expansion data (MMA-generated) *)
file = FileNameJoin[{$ProjectDir, "workspace/playground/ExpansionMMA_bub00.m"}];
Get[file];
hlist = $ExpansionResults[[1, 1, "Solutions", 1, "H"]];
hExpr[nu_List, alpha_List] := Sum[
  (hlist[[k+1]] /. Thread[vlist -> (nu - alpha)]) / n^k,
  {k, 0, Length[hlist]-1}];

(* Load relation meta *)
Get[FileNameJoin[{$ProjectDir, "relations/RelationMeta_bub00.m"}]];
regime = $RelationMeta["Regimes"][[1]];
aVals = regime["A"]; aInvVals = regime["Ainv"];
theta = regime["Sector"];

(* Generate alphas and betas for (lev=1, deg=1) *)
alphas = Flatten[Table[{a1, a2}, {a1, 0, 1}, {a2, 0, 1 - a1}], 1];
betas = Flatten[Table[{b1, b2}, {b1, 0, 1}, {b2, 0, 1 - b1}], 1];
nAlpha = Length[alphas]; nBeta = Length[betas];
kMax = order;

(* Test ν point *)
nu = {1, 0};

degEff = If[AnyTrue[theta, # =!= 0 &], 1, 0];

(* Build equation matrix M: rows = nb * (k_max+1), cols = nAlpha * nBeta *)
nb = 1;
M = ConstantArray[0, {nb * (kMax + 1), nAlpha * nBeta}];

Do[
  alpha = alphas[[ai]]; beta = betas[[bi]];
  
  (* n^beta * A^{-alpha} * h(nu - alpha) *)
  aFac = Times @@ (PowerMod[aVals, -alpha, modulus]);
  hVal = hExpr[nu, alpha];  (* Σ h_k(nu-alpha)/n^k *)
  
  (* Expand in n: (theta + nu)^beta factor *)
  betaSupp = Total[Pick[beta, theta, x_ /; x =!= 0]];
  shift = degEff - betaSupp;
  
  (* Compute coefficient of n^{degEff - row} *)
  For[row = 0, row <= kMax, row++,
    targetPower = degEff - row;
    (* Need: power of 1/n in n^beta * A^{-alpha} * h(nu-alpha) = targetPower *)
    (* n^beta contributes +betaSupp, h contributes -k, so need: betaSupp - k + shift_correction = targetPower *)
    (* Where shift adjusts for (theta+nu)^beta expansion *)
    col = (ai - 1) * nBeta + bi;
    
    (* Direct computation: (theta+nu)^beta * A^{-alpha} * h(nu-alpha) *)
    term = (Times @@ ((theta + nu)^beta)) * aFac * hVal;
    
    (* Expand in 1/n and extract coefficient of 1/n^{targetPower} *)
    series = CoefficientList[term, 1/n, kMax + 1];
    If[Length[series] >= row + 1,
      val = PolynomialMod[series[[row + 1]], modulus];
      M[[row + 1, col]] = val;
    ];
  ];
, {ai, nAlpha}, {bi, nBeta}];

Print["MMA equation matrix (", Length[M], " x ", Length[M[[1]]], "):"];
Print[MatrixForm[M]];

(* Compare with C++ matrix *)
cppMatrix = Get[FileNameJoin[{$ProjectDir, "debug_eqn_matrix.m"}]];
diff = M - cppMatrix;
Print["\nDifference (MMA - C++):"];
Print[MatrixForm[diff]];
Print["\nMax difference: ", Max[Abs[Flatten[diff]]]];
If[Max[Abs[Flatten[diff]]] === 0,
  Print["\n*** MATRICES MATCH! ***"],
  Print["\n*** MATRICES DIFFER ***"]
];
