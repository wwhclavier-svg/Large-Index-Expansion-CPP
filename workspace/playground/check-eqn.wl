(* Compare C++ and MMA equation matrices for bub00 (lev=1, deg=1) *)
$ProjectDir = "/home/ykm/Large-Index-Expansion-CPP";
$ProjectMMADir = "/home/ykm/Documents/Wolfram Mathematica/[Work] Program_LargeIndexExpansion_2508";

(* Load family config *)
Get[FileNameJoin[{$ProjectDir, "families/FamilyDatabase.wl"}]];
config = $FamilyDatabase["bub00"];
modulus = config["Modulus"];
ne = 2; order = 4;
vlist = Table[Symbol["v" <> ToString[i]], {i, ne}];
Alist = Table[A[i], {i, ne}];

(* Load expansion data *)
Get[FileNameJoin[{$ProjectDir, "workspace/playground/ExpansionMMA_bub00.m"}]];
hlist = $ExpansionResults[[1, 1, "Solutions", 1, "H"]];
hExpr = Sum[1/n^k * hlist[[k+1]], {k, 0, order}];

(* Load relation meta *)
Get[FileNameJoin[{$ProjectDir, "relations/RelationMeta_bub00.m"}]];
regime = $RelationMeta["Regimes"][[1]];
aVals = regime["A"]; aInvVals = regime["Ainv"];
theta = regime["Sector"] /. {0 -> 0, 1 -> 1};

(* Load C++ relation *)
Get[FileNameJoin[{$ProjectDir, "relations/AllRelations_bub00_k4.m"}]];
cppEntry = Select[$AllRelations, #Lev == 1 && #Deg == 1 &][[1]];
cppAlphas = cppEntry["Alphas"]; cppBetas = cppEntry["Betas"];
cppCoeffs = cppEntry["Coefficients"];
Print["C++ has ", cppEntry["NumSolutions"], " solutions, format: {particular, basis0, basis1}"];
Print["First coefficient: ", cppCoeffs[[1]]];

(* Convert C++ relation to expression form: Σ coeff · n^β · j[ν-α] *)
cppRel = Sum[
  Module[{coeff = cppCoeffs[[idx, 2]], alpha = cppAlphas[[ai]], beta = cppBetas[[bi]]},
    coeff * (Times @@ (vlist^beta)) * "j"[Sequence @@ alpha]
  ],
  {ai, Length[cppAlphas]}, {bi, Length[cppBetas]},
  {idx, ai * Length[cppBetas] + bi + 1}];
Print["CPP relation: ", cppRel];

(* Load MMA relation *)
Get[FileNameJoin[{$ProjectMMADir, "bub00/bub00_relation.m"}]];

(* Define seriesVerify for testing *)
seriesVerify[relationExpr_, regionInfo_, hExpr_, vlist_, Alist_, modulus_] := 
  Module[{rel, limitSector, aVals, aInvVals, vlistSym = vlist, ne = Length[vlist], 
          hlist0, secpos, alphaShift, result, terms},
    aVals = regionInfo["A"]; aInvVals = regionInfo["Ainv"];
    limitSector = regionInfo["Sector"] /. {0 -> 0, 1 -> 1};
    
    (* For each term "j"[a1,a2,...] in the relation, substitute A^(-α)·h(ν-α) *)
    terms = Cases[relationExpr, _"j", Infinity] // Union;
    result = relationExpr;
    Do[
      alpha = List @@ term;
      aFac = Times @@ (PowerMod[aVals, -alpha, modulus]);
      hShift = hExpr /. Thread[vlistSym -> vlistSym - alpha];
      sub = aFac * hShift;
      result = result /. (term :> sub);
    , {term, terms}];
    
    (* Expand in 1/n and check each order *)
    result = Collect[result, n];
    coeffs = CoefficientList[result, 1/n, order + 1];
    coeffs = PolynomialMod[#, modulus] & /@ Flatten[{coeffs}];
    coeffs
  ];

(* Test at ν = {1, 0} *)
nuTest = {1, 0};
Print["\n=== C++ relation at ν={1,0} ==="];
cppResult = seriesVerify[cppRel, regime, hExpr, vlist -> nuTest, vlist, Alist, modulus];
Print["Orders 0..4: ", cppResult];

(* Print MMA relation for comparison *)
Print["\n=== MMA relations (relFF[[1,2]] = (lev=1, deg=1)) ==="];
(* TODO: load MMA relation *)
