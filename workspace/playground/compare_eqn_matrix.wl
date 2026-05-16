(* Compare equation matrix: C++ vs MMA *)
(* Load C++ matrix *)
Get["/home/ykm/Large-Index-Expansion-CPP/debug_eqn_matrix.m"];
cppMatrix = $EqnMatrix;
Print["C++ matrix (", Length[cppMatrix], " x ", Length[cppMatrix[[1]]], ")"];

(* Load family config *)
$ProjectDir = "/home/ykm/Large-Index-Expansion-CPP";
Get[FileNameJoin[{$ProjectDir, "families/FamilyDatabase.wl"}]];
config = $FamilyDatabase["bub00"];
modulus = config["Modulus"];

(* Load expansion data *)
Get[FileNameJoin[{$ProjectDir, "workspace/playground/ExpansionMMA_bub00.m"}]];
ne = 2;
vlist = Table[Symbol["v" <> ToString[i]], {i, ne}];
Alist = Table[A[i], {i, ne}];

(* Extract h_k from expansion *)
hlist = $ExpansionResults[[1, 1, "Solutions", 1, "H"]];

(* Load relation meta *)
Get[FileNameJoin[{$ProjectDir, "relations/RelationMeta_bub00.m"}]];
regime = $RelationMeta["Regimes"][[1]];
aVals = regime["A"];
aInvVals = regime["Ainv"];
theta = regime["Sector"] /. {0 -> 0, 1 -> 1};

(* Build equation matrix for (lev=1, deg=1) at ν = {1,0} *)
nuSample = {FFInt[1], FFInt[0]};
(* Convert to regular integers for comparison *)
nu = {1, 0};

(* Generate alphas and betas *)
lev = 1; deg = 1; order = 4;
alphas = Flatten[Table[{a1, a2}, {a1, 0, lev}, {a2, 0, lev - a1}], 1];
betas = Flatten[Table[{b1, b2}, {b1, 0, deg}, {b2, 0, deg - b1}], 1];
Print["Alphas: ", alphas];
Print["Betas: ", betas];

nAlpha = Length[alphas]; nBeta = Length[betas];
kMax = order;

(* Build equation matrix M: rows = nb * (k_max+1), cols = nAlpha * nBeta *)
nb = 1;
M = ConstantArray[0, {nb * (kMax + 1), nAlpha * nBeta}];

Do[
  alpha = alphas[[ai]];
  beta = betas[[bi]];
  
  (* 1. g = A^(-alpha) * h_k(nu - alpha) *)
  gvec = Table[
    Module[{sum = 0},
      hl = hlist[[k + 1]];
      (*  h_k expressed as polynomial in v1, v2 *)
      hVal = hl /. {v1 -> nu[[1]] - alpha[[1]], v2 -> nu[[2]] - alpha[[2]]};
      hVal
    ],
    {k, 0, kMax}
  ];
  Print["alpha=", alpha, " beta=", beta, " g=", gvec];
  
  , {ai, nAlpha}, {bi, nBeta}
];
