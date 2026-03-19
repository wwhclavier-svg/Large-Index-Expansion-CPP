(* ::Package:: *)

(* ::Title:: *)
(*ReconstructReductionRelation - Full-featured Relation Reconstruction*)
(*Based on modular 5-step structure with multi-region and algebraic extension support*)


(* ::Section:: *)
(*Package Header*)


BeginPackage["ReconstructReductionRelation`"]

(* Main exported function - compatible with reconstructRelation5 interface *)
ReconstructReductionRelation::usage = "ReconstructReductionRelation[level, degree, order, hexpnlist, areglist, topsec, vlist, opts] \
reconstructs IBP reduction relations from large-index expansions.\n\
Supports multi-region, multi-solution, and algebraic extension fields.\n\
\nParameters:\n\
- level: maximum rank level for ansatz\n\
- degree: maximum coefficient polynomial degree fro ansatz\n\
- order: expansion order to use\n\
- topsec: top sector specification\n\
- hexpnlist: {{{h00, h01, ...}, {h10, h11, ...}, ...}, ...} expansion coefficients per region per solution\n\
- areglist: list of region specifications with algebraic extension data\n\
- vlist: list of variables {v1, v2, ..., vn}\n\
- opts: algorithm options\n\
\nReturns: List of relations per level per degree"

(* Utility exports *)
GenerateAnsatz::usage = "Generate ansatz polynomials for given rank and mode."
ComputeBasisPower::usage = "Compute basis power tables for algebraic extensions."
(* Future: ReconstructLinearRelation for simplified single-region problems *)

Begin["`Private`"]


(* ::Section:: *)
(*Configuration and Global State*)


$RRRVerbose = False;
$LevelCount = 0;

SetAttributes[WithIndent, HoldAll];
WithIndent[body_] := Block[{$LevelCount = $LevelCount + 1}, body];

PrintV[args___] := If[$RRRVerbose, Print[StringRepeat["  ", $LevelCount], args]];

SetAttributes[MyTimer, HoldAll];
MyTimer[title_, body_] := Module[{result, time},
  {time, result} = AbsoluteTiming[body];
  Print[title, " ... ", time];
  result
];


(* ::Section:: *)
(*Options*)


Options[ReconstructReductionRelation] = {
  (* Computation mode *)
  Modulus -> 0,                    (* 0 = exact/Q, p = finite field F_p *)
  Verbose -> False,                (* Debug output *)
  
  (* Ansatz configuration *)
  "AnsatzMode" -> 0,               (* 0=Pyramid, 1=DotPyramid, 2=Star *)
  "CoefDeg" -> 1,                  (* Maximum degree of coefficient ansatz *)
  "LimitSector" -> Automatic,      (* Sector for n->infinity limit *)
  
  (* Algorithm strategy *)
  "Strategy" -> "Incremental",     (* "Direct", "Incremental", "RegionByRegion" *)
  "PlateauSize" -> 1,              (* Stability plateau for convergence *)
  "DegreeByDegree" -> True,        (* Add equations degree by degree *)
  
  (* Linear solver *)
  "LinearSolver" -> Automatic,     (* Automatic = reduceSolve or modSolve *)
  
  (* Numerical/Optimization *)
  "Numeric" -> {},                 (* Fixed numerical values *)
  "ExtraNTerm" -> False            (* Include extra N-dependent terms *)
};


(* ::Section:: *)
(*Utility Functions*)


(* First non-zero index in list *)
FirstNonZero[list_] := Module[{i = 1},
  While[i <= Length[list] && list[[i]] === 0, i++];
  If[i > Length[list], -1, i]
];

(* Modular reduction *)
FFReduce[expr_, 0] := expr;
FFReduce[expr_, p_] := PolynomialMod[expr, p];

(* Seed generation: all exponent vectors with total power <= maxPower *)
GenerateSeeds[maxPower_, nVars_] := If[maxPower < 0, {},
  Join @@ Table[
    SortBy[Flatten[Permutations /@ (IntegerPartitions[p + nVars, {nVars}] - 1), 1], -Reverse[#] &],
    {p, 0, maxPower}
  ]
];

(* Circled minus for subsector structure {i1,i2,...,ik} *)
CircledMinus[maxPower_, subsec_, nVars_] := 
  Normal@SparseArray[Thread[subsec -> #], nVars] & /@ GenerateSeeds[maxPower, Length[subsec]];

(* Dot-rank seed generation *)
DotRankSeeds[sector_, maxDot_, maxRank_] := Module[
  {pdSet, ispSet, nPd, nIsp, dotSeeds, rankSeeds, Seeds},
  pdSet = Flatten[Position[sector, 1]];
  ispSet = Flatten[Position[sector, 0]];
  nPd = Length[pdSet];
  nIsp = Length[ispSet];
  dotSeeds = GenerateSeeds[maxDot, nPd];
  rankSeeds = GenerateSeeds[maxRank, nIsp];
  Seeds=Join @@ Outer[
    Normal@SparseArray[
      Join[Thread[pdSet -> (1 + #1)], Thread[ispSet -> -#2]],
      nPd + nIsp
    ] &,
    dotSeeds, rankSeeds, 1
  ];
  SortBy[Seeds, SectorMO[sector,#]&]
];


(* Sector-based monomial ordering *)
Options[SectorMO]={"Reverse" -> False};
SectorMO[sector_, index_, opts:OptionsPattern[]] := Module[
  {n = Length[sector], pdSet, ispSet},
  pdSet = Flatten[Position[sector, 1]];
  ispSet = Complement[Range[n], pdSet];
  If[OptionValue["Reverse"] == True, pdSet = Reverse@pdSet];
  {
    Total[index[[pdSet]]],           (* dot *)
    -Total[index[[ispSet]]],         (* rank *)
    -index[[ispSet]],                (* isp indices *)
    index[[pdSet]]                   (* pd indices *)
  }
];

SectorMOneg[sector_, index_] := SectorMO[sector, -index];


SectorMO[index_]:=Module[{sector},
sector=(index/.{a_/;a > 0 -> 1, a_/;a <= 0 -> 0});
Join[{sector}, SectorMO[sector, index, "Reverse"->True]]];


CoefficientEquations[expr_,vars_]:=Module[{eqs},
If[expr==={},Return[{}]];
eqs=CoefficientArrays[expr,vars];
eqs=Join@@(If[Head[#]===SparseArray,Most@ArrayRules[#][[;;,2]],{#}]&/@eqs);
eqs//DeleteCases[#,0]&];


(* ::Section:: *)
(*Module 1: GenerateAnsatz - Create seed polynomial ansatz*)


GenerateAnsatz[mode_, rank_, maxDeg_, limitSector_, nVars_, vlist_] := Module[
  {seeding, fullAnsatz, levelList, secIsp, levelFilter, getPoly, monomials, seedList},
  
  (* Helper: generate monomials up to given degree *)
  seeding[n_, lev_] := Join @@ Table[
    SortBy[Flatten[Permutations /@ (IntegerPartitions[l + n, {n}] - 1), 1], -Reverse[#] &],
    {l, 0, lev}
  ];(* GenerateSeeds[] *)
  seedList = seeding[nVars, maxDeg];
  monomials[d_] := (Times @@ Thread@Power[vlist, #]) & /@ seeding[nVars, d];
  
  (* Polynomial builder: Sum["b"[rk, power] * v^power] *)
  getPoly[rk_, deg_] := Dot[("b"[rk, #] & /@ seeding[nVars, deg]), monomials[deg]];
  
  secIsp = Flatten[Position[limitSector, 0]];
  
  (* Level filter by mode *)
  levelFilter = Switch[mode,
    0, Total[vlist - List @@ #] &,      (* Pyramid: distance from vlist *)
    1, Total[List @@ # - vlist] &,      (* Dot pyramid *)
    2, Total[Abs /@ (List @@ # - vlist)] &  (* Star: L1 distance *)
  ];
  
  (* Generate full ansatz by mode *)
  fullAnsatz = Switch[mode,
    0, (* Pyramid/Extended mode *)
    Module[{rkGe, rkG, base, boundary},
      rkGe = CircledMinus[rank + 1, Range[nVars], nVars];
      rkG = Select[rkGe, Total[#] <= rank &];
      
      (* Base ansatz *)
      base = Association@Table[
        "g" @@ (vlist - rk) -> getPoly[rk, maxDeg],
        {rk, rkG}
      ];
      
      (* Boundary terms: v_i * ansatz at level-1 *)
      boundary = Merge[Table[
        rkG = Select[rkGe, #[[i]] == 0 &];
        Association@Table[
          "g" @@ (vlist + UnitVector[nVars, i] - rk) -> (
            vlist[[i]] * Dot[
              ("b"[rk - UnitVector[nVars, i], # + UnitVector[nVars, i]] & /@ seeding[nVars, maxDeg - 1]),
              monomials[maxDeg - 1]
            ]
          ),
          {rk, rkG}
        ],
        {i, secIsp}
      ], First];
      
      Join[base, boundary]
    ],
    
    1, (* Dot Pyramid *)
    Association@Table[
      "g" @@ (vlist + rk) -> getPoly[rk, maxDeg],
      {rk, CircledMinus[rank + 1, Range[nVars], nVars]}
    ],
    
    2, (* Star-like: per-sector combinations *)
    Module[{rkS},
      rkS = Join @@ Table[
        DotRankSeeds[sec, rank, rank],
        {sec, Normal@SparseArray[# -> 1 & /@ #, nVars] & /@ Subsets[Range[nVars]]}
      ];
      Association@Table["g" @@ (vlist + rk) -> getPoly[rk, maxDeg], {rk, rkS}]
    ]
  ];
  
  (* Group by level *)
  levelList = Switch[mode,
    0, If[Union[limitSector] =!= {1}, Range[0, rank], Range[rank]],
    1, Range[rank],
    2, Range@Max[levelFilter /@ Keys[fullAnsatz]]
  ];
  
  Module[{grouped},
    grouped = GroupBy[Normal[fullAnsatz], levelFilter[First[#]] &];
    Association@Table[
      level -> Association[Join @@ Table[Lookup[grouped, i, {}], {i, 0, level}]],
      {level, levelList}
    ]
  ]
];


(* ::Section:: *)
(*Module 2: Basis Power Computation with Algebraic Extensions*)


(*
  Design:
  - ComputeBasisPower: Main entry point, computes all A^(-alpha) for a region
  - Uses recursive memoization (like variablePowerTable) to avoid redundant computation
  - Handles both rational (scalar) and algebraic extension (vector) cases
*)

(* Helper: Extract monomial rules for basis conversion *)
MonomialRulesPower[expr_, var_] := If[var =!= {},
  Association[Join @@ (
    If[Head[#] === SparseArray,
      Table[Length@Cases[#[[1]], i, Infinity], {i, Length[var]}] -> #[[2]] & /@ Most@ArrayRules[#],
      {Table[0, Length[var]] -> #}
    ] & /@ CoefficientArrays[expr, var]
  )],
  <|{} -> expr|>
];

(* Build variable matrices for algebraic extension *)
BuildVariableMatrices[aList_, aRule_, aIndep_, basisMatrix_, modulus_] := Module[
  {varAMatrix, varBMatrix, nb = Length[basisMatrix]},
  
  (* A[i] matrices in basis representation *)
  varAMatrix = Association@Table[
    j -> With[{rules = MonomialRulesPower[aList[[j]] /. aRule, aIndep]},
      FFReduce[(Keys[rules] /. basisMatrix) . Values[rules], modulus]
    ],
    {j, Length[aList]}
  ];
  
  (* B[i] = A[i]^(-1) matrices *)
  varBMatrix = Association@Table[
    j -> LinearSolve[varAMatrix[j], IdentityMatrix[nb], Modulus -> modulus],
    {j, Length[aList]}
  ];
  
  {varAMatrix, varBMatrix}
];

(* Recursive power table computation (like variablePowerTable) *)
ComputePowerTable[indexList_List, varAMatrix_, varBMatrix_, nb_, modulus_:0] := Module[
  {powerTable = <||>, getPower, ne = Length[varAMatrix]},
  
  (* Base case: B^0 = identity vector *)
  powerTable[ConstantArray[0, ne]] = IdentityMatrix[nb][[1]];
  
  (* Recursive memoized getter *)
  getPower[idx_] := powerTable[idx] /; KeyExistsQ[powerTable, idx];
  
  getPower[idx_] := Module[{i, res, unit},
    If[KeyExistsQ[powerTable, idx], Return[powerTable[idx]]];
    
    (* Find first non-zero entry in idx *)
    i = FirstNonZero[idx];
    If[i < 0, Return[powerTable[ConstantArray[0, ne]]]];  (* All zero *)
    
    unit = UnitVector[ne, i];
    If[idx[[i]] > 0,
      (* idx[i] > 0: use B[i] * B^(idx - e_i) *)
      res = varBMatrix[i] . getPower[idx - unit],
      (* idx[i] < 0: use A[i] * B^(idx + e_i) *)
      res = varAMatrix[i] . getPower[idx + unit]
    ];
    
    res = FFReduce[res, modulus];
    powerTable[idx] = res;
    res
  ];
  
  (* Compute all requested powers *)
  Scan[getPower, SortBy[indexList, Total]];
  
  (* Return association: index -> power vector *)
  Association@Table[idx -> powerTable[idx], {idx, indexList}]
];

(* Main function: Compute basis power table for a region *)
ComputeBasisPower[alphaList_, areg_, modulus_:0] := Module[
  {aRule, aIndep, aDep, aList, basisMatrix, basis, nb, 
   varAMatrix, varBMatrix, indexList, powerTable, scalarPowers, signList},
  
  (* Extract data from areg *)
  aRule = areg["VarRule"];
  aIndep = areg["VarIndep"];
  aDep = aRule[[;;, 1]];
  aList = Sort@Join[aDep, aIndep];
  basisMatrix = areg["MonomialBasisMatrix"];
  basis = areg["MonomialBasis"];
  nb = Length[basis];
  
  (* Convert alpha list to index list (negative for A^(-alpha)) *)
  indexList = -# & /@ alphaList;
  
  If[nb == 1,
    (* === Rational case: scalar computation === *)
    (* Direct computation: A^(-alpha) = Product[A[i]^(-alpha[i])] *)
    scalarPowers = (aList /. aRule);
    signList = If[#>=0,1,-1]&/@indexList;
    Association@Table[
      alphaList[[k]] -> {{Times @@ Thread@Power[scalarPowers, indexList[[k]]]//PolynomialMod[#,modulus]&}},
      {k, Length[alphaList]}
    ],
    
    (* === Algebraic extension case: vector computation === *)
    (* Build variable matrices *)
    {varAMatrix, varBMatrix} = BuildVariableMatrices[
      aList, aRule, aIndep, basisMatrix, modulus
    ];
    
    (* Compute recursive power table for all indices *)
    powerTable = ComputePowerTable[
      indexList, varAMatrix, varBMatrix, nb, modulus
    ];
    
    (* Convert back: index -> alpha, only for requested alphas *)
    Association@Table[
      alphaList[[k]] -> powerTable[indexList[[k]]],
      {k, Length[alphaList]}
    ]
  ]
];


(* ::Section:: *)
(*Module 3: Expansion Table Update - Multi-region Multi-solution*)


UpdateExpansionTable[expnTable_, basePower_, hexpnList_, indexNew_, vlist_, maxOrder_, modulus_] := Module[
  {nReg = Length[hexpnList], nSol, hSeries, expnNewList, expnNew, rk},
  
  Table[
    nSol = Length[hexpnList[[j]]];
    Table[
      hSeries = hexpnList[[j, m]];
      expnNewList = PolynomialMod[#,modulus]&@Expand@Table[
        (* Compute g_alpha^(i) = p_alpha * h^(i)(nu - alpha) for all alpha in indexNew *)
        expnNew = Table[
          alpha -> Dot[
            basePower[[j]][alpha],
            hSeries[[i]] /. Dispatch@Thread[vlist -> (vlist - alpha)]
          ],
          {alpha, indexNew}
        ];
        
        (* Update or create association *)
        If[i > Length[expnTable[[j, m]]],
          Association@expnNew,
          Append[expnTable[[j, m, i]], expnNew]
        ],
        {i, Min[maxOrder + 1, Length[hSeries]]}
      ];
      expnNewList,
      {m, nSol}
    ],
    {j, nReg}
  ]
];(* table index: (ref, sol, order; alpha) alpha being key *)


(* ::Section:: *)
(*Module 4: Linear System Utilities*)


RemoveSolvedVariables[cfList0_, alphaList0_, gvList0_, levelBindep_, bindep_, deg_] := Module[
  {cfList, empty, bVars, gvList, alphaList, bVarGroup},
  
  PrintV["Eliminating solved variables..."];
  
  (* Degree truncation *)
  cfList = cfList0 /. {"b"[alpha_, beta_] /; Total[beta] > deg -> 0};
  
  (* Eliminate solutions from lower degree *)
  cfList = cfList /. {"b"[alpha_, beta_] /; Or @@ (alpha == #[[1]] && And @@ Thread[beta >= #[[2]]] & /@ levelBindep) -> 0};
  
  (* Eliminate solutions from lower seed levels *)
  cfList = cfList /. {"b"[alpha_, beta_] /; Or @@ (And @@ Thread[alpha >= #[[1]]] && And @@ Thread[beta >= #[[2]]] & /@ Join @@ bindep) -> 0};
  
  (* Remove zero entries *)
  empty = Select[Range@Length@cfList, cfList[[#]] == 0 &];
  PrintV["  Eliminated: ", Length[empty], " entries"];
  {cfList, gvList, alphaList} = Part[#, Complement[Range[Length@cfList], empty]] & /@ {cfList, gvList0, alphaList0};
  
  (* Extract b-variables *)
  bVars = Union@Cases[cfList, a_/;Head[a]=="b", Infinity];
  bVars = Reverse@SortBy[bVars, {SectorMO[-#[[1]]], Total@#[[2]], -Reverse[#[[2]]]} &];
  
  {cfList, gvList, alphaList, bVars}
];
 (* Separate numerically fixed variables *)
InitializeVariables[cfList0_, bVars_, vFixed_]:=Module[{cfList, bVarGroup=<||>, bVarReg, bSolNew},
If[vFixed =!= {},
	bVarGroup = GroupBy[bVars, Union[#[[2, vFixed]]] == {0} &];
      bVarReg = bVarGroup // If[MemberQ[Keys[#], True], #[True], {}]&;
      bSolNew = bVarGroup // Association@If[MemberQ[Keys[#], False], # -> 0 & /@ #[False], {}]&;
        cfList = cfList0 /. bSolNew
        ,
    {bVarReg, bSolNew, cfList} = {bVars, <||>, cfList0}
  ];

  PrintV["  Remaining b-vars: ", Length[bVarReg]];
  {cfList, bVarReg, bSolNew}
];

(* Setup coefficient table for b-polynomials *)
SetupCoefficientTable[cfList_, aregList_, vlist_, order_] := Module[
  {nVars = Length[vlist], nReg = Length[aregList], cfTable, limitSector, limitRule},
  
  Table[
    limitSector = aregList[[j]]["LimitSector"];
    limitRule = Table[
      vlist[[i]] -> If[limitSector[[i]] > 0, n + vlist[[i]], vlist[[i]]],
      {i, nVars}
    ];
    Transpose@Table[
      Reverse@CoefficientList[cfList[[l]] /. limitRule, n, order + 1],
      {l, Length[cfList]}
    ],
    {j, nReg}
  ]
];(* table index: (reg, order, alpha) *)



(* Linear solver dispatcher *)
DispatchLinearSolve[eqs_, vars_, modulus_, solver_] := Which[
  modulus =!= 0 && solver === Automatic,
  ReduceSolve[eqs, vars, modulus],
  solver === Automatic,
  ReduceSolve[eqs, vars, 0],
  True,
  solver[eqs, vars]
];

(* Modular/Exact reduce solve *)
ReduceSolve[eqs_, vars_, p_] := Module[
  {mat, dep, indep, rule, varRev = Reverse[vars]},
  If[Length[eqs] == 0, Return[{}]];
  mat = CoefficientArrays[eqs, varRev];
  mat = Join[mat[[2]], List /@ mat[[1]], 2];
  mat = If[p =!= 0, RowReduce[mat, Modulus -> p], RowReduce[mat]];
  mat = DeleteCases[mat, a_ /; Union[Flatten[{a}]] == {0}];
  dep = FirstNonZero /@ mat;
  If[Last[dep] > Length[vars], Return[$Failed]];
  indep = Complement[Range[Length[vars]], dep];
  rule = Thread[varRev[[dep]] -> (-mat[[;;, indep]] . varRev[[indep]] - mat[[;;, -1]])];
  rule
];



(* Back-substitute solutions *)
ResubstituteSolutions[bSol_, gvList_, cfList_, bVars_, topsec_] := Module[
  {recIbp, bVarIndep, redRel},
  
  If[bSol === <||>, Return[{{}, {}}]];
  
  (* Reconstruct relation *)
  recIbp = Sum[
    gvList[[i]] * (cfList[[i]] /. bSol),
    {i, Length[cfList]}
  ] // Collect[#, bVars] &;
  
  (* Identify independent variables *)
  bVarIndep = DeleteDuplicates@Cases[recIbp, a_ /; MemberQ[bVars, a], {0, Infinity}];
  bVarIndep = Reverse@SortBy[bVarIndep, SectorMO[topsec, -#[[1]],"Reverse"->True] &];
  
  (* Extract relation *)
  redRel = If[Length[#] >= 2, #[[2]], {}] & @ Normal@CoefficientArrays[recIbp, bVarIndep];
  
  {redRel, bVarIndep}
];


(* ::Section:: *)
(*Module 5: Equation Solver Core with State Control*)


(*
  SolveDegreeEquations: Core solver for a specific seed level and coefficient degree.
  
  This module encapsulates the state-controlled equation solving loop,
  including strategy selection (all-together vs region-by-region),
  stability detection, and plateau checking.
*)

SolveDegreeEquations[coeftable0_, expntable_, alphalist_, bVars_, 
  bVarReg0_, bSolNew0_, deg_, maxorder_, level_, vlist_, numeric_, 
  modulus_, linearsolver_, strategy_, plateausize_, degreebydegree_] := 
Module[
  {
    (* State variables *)
    state = <|"Stable" -> False, "StableLevel" -> -2, "Trivial" -> True|>,
    bSolAcc = bSolNew0,
    bSolNew = <||>,
    bVarReg = bVarReg0,
    coefTable = coeftable0,
    
    (* Loop variables *)
    k = 0,
    eqs, eqs1,
    converged = False,
    nreg = Length[expntable]
  },
  
  (* Main solving loop over orders *)
  While[k <= maxorder && !converged,
    
    (* Strategy branch: All regions together vs Region by region *)
    If[(k == 0 && deg <= 1) || strategy =!= "RegionByRegion",
      (* === Branch 1: All regions together === *)
      
      bVarReg = DeleteCases[bVarReg0, a_ /; MemberQ[Keys@bSolNew, a]];
      coefTable = Map[PolynomialMod[#, modulus]&@Expand@ReplaceAll[#, bSolNew] &, coefTable, {2}];
      {eqs, bVarReg} = SetupEquationsAll[
        coefTable, expntable, alphalist, bSolAcc, bVarReg, 
        deg, maxorder, vlist, k, level, numeric, modulus, degreebydegree
      ];
      (*PrintV["eqs  ",eqs//MatrixForm];*)
      PrintV["    + SetupEquationsAll  (",Length@eqs,"\[Times]",Length@bVarReg,") -- over"];
      
      (* [Trigger 1] Handle trivial equations *)
      If[eqs === {},
        PrintV["  (ord=", k, "): #eqs=0, skipping"];
        If[!state["Stable"] && k > 0,
          state["Stable"] = True;
          state["StableLevel"] = k - 1;
        ],
        
        (* Non-trivial: solve *)
        state["Trivial"] = False;
        bSolNew = Association@DispatchLinearSolve[eqs, bVarReg, modulus, linearsolver];
        PrintV["    + DispatchLinearSolve -- over"];
        
        (* Check for empty solution *)
        If[bSolNew === <||>,
          PrintV["  Empty solution at order ", k];
          Return[{<||>, state, state["StableLevel"]}, Module]
        ];
        
        (* Update accumulated solutions *)
        bSolAcc = Map[PolynomialMod[#, modulus]&@Expand@ReplaceAll[#, bSolNew] &, bSolAcc];
        bSolAcc = AssociateTo[bSolAcc, bSolNew];        
        bVarReg = DeleteCases[bVarReg, a_ /; MemberQ[Keys@bSolNew, a]];
        
        PrintV["  Solved: ", Length[bSolAcc], "/", Length[bVars]];
        
        (* Check stability *)
        If[Length[bSolAcc] == Length[bVars] && !state["Stable"],
          state["Stable"] = True;
          state["StableLevel"] = k;
          PrintV["  Stability achieved at order ", k];
        ];
      ],

      (* === Branch 2: Region by region === *)
      state["Trivial"] = True;
      
      Do[
        (* Skip if beyond available orders for this region *)
        If[k + 1 > Length[expntable[[j]]], Continue[]];
        
      bVarReg = DeleteCases[bVarReg0, a_ /; MemberQ[Keys@bSolNew, a]];
      coefTable = Map[FFReduce[#, modulus]&@Expand@ReplaceAll[#, bSolNew] &, coefTable, {2}];
        {eqs, bVarReg} = SetupEquationsSingle[
          coefTable[[j]], expntable[[j]], alphalist, bSolNew, bVarReg,
          deg, maxorder, vlist, k, level, numeric, modulus
        ];
        
        PrintV["  Region ", j, " (ord=", k, "): #eqs=", Length[eqs], ", #vars=", Length[bVarReg]];
        
        (*[Trigger 1] Skip empty equations *)
        If[eqs === {}, Continue[]];
        state["Trivial"] = False;
        
        (* Solve *)
        bSolNew = Association@DispatchLinearSolve[eqs, bVarReg, modulus, linearsolver];
        (*PrintV["New Sol:  ",bSolNew//MatrixForm];*)
        If[bSolNew === <||>,
          PrintV["  Empty solution at region ", j, ", order ", k];
          Return[{<||>, state, state["StableLevel"]}, Module]
        ];
        
        (* Update solutions *)
        bSolAcc = Map[FFReduce[#, modulus]@Expand@ReplaceAll[#, bSolNew] &, bSolAcc];
        bSolAcc = AssociateTo[bSolAcc, bSolNew];
        PrintV["solved:  (",Length@bSolAcc,"/",Length@bVars,")"];
        PrintV["unsolved:  (",Length@bVars-Length@bSolAcc,"/",Length@bVars,")"];
        (*bVarReg = DeleteCases[bVarReg, a_ /; MemberQ[Keys@bSolNew, a]];*)
        
        
        (* Check stability *)
        If[Length[bSolAcc] == Length[bVars] && !state["Stable"],
          state["Stable"] = True;
          state["StableLevel"] = k;
          PrintV["  Stability achieved at region ", j, ", order ", k];
        ];
        
        , {j, nreg}
      ];  (* End region loop *)
      
      (* Check if all equations were trivial *)
      If[!state["Stable"] && state["Trivial"] && k =!= 0,
        PrintV["  All trivial equations, marking stable at ", k - 1];
        state["Stable"] = True;
        state["StableLevel"] = k - 1;
      ];
    ];  (* End strategy branch *)
    
    PrintV["  Flags at order ", k, ": ", state];
    
    (* [Trigger 3] Plateau check: stable for required number of orders *)
    If[state["Stable"] && k >= state["StableLevel"] + plateausize,
      PrintV["  Plateau reached at level ", state["StableLevel"]];
      converged = True;
      Break[];
    ];
    
    k++;
  ];  (* End While *)
  (*PrintV["Sol:  ",bSolAcc//MatrixForm];*)
  (* Return solution, final state, and stable level *)
  {bSolAcc, state, state["StableLevel"]}
];


(* Setup equations for all regions together *)
SetupEquationsAll[coeftable_, expntable_, alphaList_, bSolNew_, bVarReg0_, 
  deg_, maxorder_, vlist_, k_, level_, numeric_, modulus_, degreebydegree_] :=
Module[
  {coeftableMod, relansatz, bVarReg=bVarReg0, eqs, eqs1, bmat, nreg = Length[expntable]},
  
  (* Pre-substitution: apply known solutions *)
(*    bVarReg = DeleteCases[bVarReg0, a_ /; MemberQ[Keys@bSolNew, a]];
  coeftableMod = Map[FFReduce[#, modulus]&@Expand@ReplaceAll[#, bSolNew] &, coeftable, {2}];
*)
  
  (* Build relation ansatz *)
  relansatz = PolynomialMod[#, modulus]&@Table[
    Table[
      Sum[
        coeftable[[j, kk + 1, l]] * expntable[[j, m, k - kk + 1]][alphaList[[l]]],
        {l, Length[alphaList]}, {kk, 0, Min[k, deg]}
      ],
      {m, Length[expntable[[j]]]}
    ],
    {j, nreg}
  ];
  
  PrintV["    + Convolution -- over"];
  
  (* Extract equations from coefficients *)
  eqs = Flatten[relansatz /. numeric];
  eqs = DeleteCases[eqs, 0];
  eqs = CoefficientEquations[eqs, vlist];
  
  If[eqs === {}, Return[{eqs, bVarReg}]];
  
  (* Add more equations if needed (degreeByDegree strategy) *)
(*  While[degreebydegree && Length[eqs] < Length[bVarReg] && k < maxorder,
    (* This logic is simplified; original has more complex condition *)
    Break[]
  ];*)
  (* NOTE: Dynamic equation addition (disabled but preserved for future use)
     Purpose: When the current equation system is underdetermined 
     (Length[eqs] < Length[bVarReg]), dynamically add higher-order equations.
     
     Behavior when enabled:
     - Increment order k while k <= maxorder
     - Build relansatz for higher orders using expntable
     - Extract new equations via CoefficientEquations
     - Append to eqs until sufficient or maxorder reached
     
     Current status: Disabled via If[False, ...] wrapper.
     To enable: Change If[False, ...] to If[True, ...] or remove the wrapper.
  *)
  If[False,
    While[(degreebydegree == False && (k+1 <= Min[maxorder, level+deg-1])) ||
          (degreebydegree == True && (Length@eqs < Length@bVarReg]) && (k < maxorder)),
      k = k+1; PrintV["  ->  order=", k];
      With[{relansatz1},
        MyTimer["  pre-substitution:  ",
          relansatz1 = PolynomialMod[#, modulus] & @ Table[Table[
            Sum[coeftable[[j, kk+1, l]] * expntable[[j, m, k - kk + 1]][alphaList[[l]]],
              {l, Length@alphaList}, {kk, 0, Min[k, deg]}]
          , {m, Length@expntable[[j]]}], {j, nreg}]
        ];
        MyTimer["equation extraction",
          eqs1 = Flatten[relansatz1 /. numeric];
          eqs1 = CoefficientEquations[#, vlist] & @ PolynomialMod[Expand@#, modulus] & @ eqs1;
        ];
        eqs = Join[eqs, eqs1];
      ];
    ]
  ];
(* ? *) 
  
  (* Convert to matrix form *)
  bmat = CoefficientArrays[eqs, bVarReg];
  eqs = Join[bmat[[2]], List /@ bmat[[1]], 2] . Join[bVarReg, {1}];
  
  {eqs, bVarReg}
];


(* Setup equations for single region *)
SetupEquationsSingle[coeftableReg_, expntableReg_, alphaList_, bSolNew_, bVarReg0_,
  deg_, maxorder_, vlist_, k_, level_, numeric_, modulus_] :=
Module[
  {coeftableMod, relansatz, bVarReg = bVarReg0, eqs, cvars, bmat},
  
  (* Apply known solutions *)
(*  bVarReg = DeleteCases[bVarReg, a_ /; MemberQ[Keys@bSolNew, a]];
  coeftableReg = Map[FFReduce[#, modulus]@Expand@ReplaceAll[#, bSolNew] &, coeftableReg, {1}];
*)  
  (* Build relation for region j only *)
  relansatz =  PolynomialMod[#, modulus]&@Table[
    Sum[
      coeftableReg[[kk + 1, l]] * expntableReg[[m, k - kk + 1]][alphaList[[l]]],
      {l, Length[alphaList]}, {kk, 0, Min[k, deg]}
    ],
    {m, Length[expntableReg]}
  ];
  
  (* Extract equations *)
  eqs = Flatten[relansatz /. numeric];
  eqs = CoefficientEquations[eqs, vlist];
  
  (* Extract coefficient variables if any *)
  cvars = Union@Cases[eqs, _c, Infinity] // Sort;
  If[cvars =!= {}, eqs = CoefficientArrays[eqs, cvars][[2]] . cvars];
  
  (* Update variable list *)
  (*bVarReg = DeleteCases[bVarReg, a_ /; MemberQ[Keys@bSolNew, a]]*);
  
  If[eqs === {}, Return[{eqs, bVarReg}]];
  
  bmat = CoefficientArrays[eqs, bVarReg];
  eqs = Join[bmat[[2]], List /@ bmat[[1]], 2] . Join[bVarReg, {1}];
  
  {eqs, bVarReg}
];


(* ::Section:: *)
(*Main Function: ReconstructReductionRelation*)


ReconstructReductionRelation[rankLevel_, degree_, order_, hexpnList_, aregList_, topsec_, vlistIn_, opts:OptionsPattern[]] := Module[
  {
    (* Options *)
    modulus = OptionValue[Modulus],
    verbose = OptionValue[Verbose],
    ansatzMode = OptionValue["AnsatzMode"],
    maxCoefDeg = degree,
    limitSector = OptionValue["LimitSector"] /. Automatic -> Table[1, Length[vlistIn]],
    strategy = OptionValue["Strategy"],
    plateauSize = OptionValue["PlateauSize"],
    degreeByDegree = OptionValue["DegreeByDegree"],
    linearSolver = OptionValue["LinearSolver"],
    numeric = OptionValue["Numeric"],
    
    (* Variables *)
    vlist = vlistIn,  (* Local copy to avoid context issues *)
    nVars = Length[vlistIn],
    nReg = Length[aregList],
    (* Working data *)
    ansatzData, levelList, alphaList,
    basePower, expnTable,
    vFixed,
    
    (* Per-level data *)
    ansatz0, cfList0, gvList0, rankList0, indexList0, indexNew,
    levelRelList, levelBindep, levelStableBound,
    
    (* Per-degree data *)
    cfList, gvList, bVars, bVarReg, bSolNew,
    coeffTable, bSol,
    
    (* Solving result *)
    state, stableLevel, redRel, bVarIndep,
    
    (* Output *)
    relationList = {}, bindep = {}, stableBound = {}
  },

  (* Initialization *)
  $RRRVerbose = verbose;
  $LevelCount = 0;
  
  vFixed = If[numeric =!= {},
    Flatten[Position[vlist, #] & /@ numeric[[;;, 1]]],
    {}
  ];
  
  (* ======================================================
     Step 1: Generate Ansatz
     ====================================================== *)
  PrintV["Generating ansatz..."];
  ansatzData = GenerateAnsatz[ansatzMode, rankLevel, maxCoefDeg, limitSector, nVars, vlist];
  levelList = Normal@Keys[ansatzData];
  PrintV["  Levels: ", levelList];
  
  (* ======================================================
     Step 2: Prepare Basis Power Tables (per region)
     ====================================================== *)
  PrintV["Computing basis power tables..."];
  alphaList=(vlist - List @@ #) & /@ Keys[ansatzData[[-1]]];
  basePower = Table[
    ComputeBasisPower[
      alphaList,  (* Use top level keys as reference *)
      aregList[[j]],
      modulus
    ],
    {j, nReg}
  ];
    
  
  (* ======================================================
     Step 3: Initialize Expansion Table
     ====================================================== *)
  expnTable = Table[
    ConstantArray[
      ConstantArray[<||>, order + 1],
      Length[hexpnList[[j]]]
    ],
    {j, nReg}
  ];
  bindep = {};
  

  (* ======================================================
     Main Loop: Over seed levels
     ====================================================== *)
  Do[
    WithIndent[
      PrintV["=== Seed Level: ", level, " ==="];
      
      (* Get ansatz for this level *)
      ansatz0 = ansatzData[level];
      cfList0 = Values[ansatz0];
      gvList0 = Keys[ansatz0];
      rankList0 = (vlist - List @@ #) & /@ gvList0;
      
      (* Compute new indices for this level *)
      indexList0 = rankList0;
      indexNew = Complement[
        indexList0,
        If[level > 1,
          (vlist - List @@ #) & /@ List @@ Keys[ansatzData[level - 1]],
          {}
        ]
      ];
      
      PrintV["  New indices: ", Length[indexNew]];
      
      (* Initialize per-level storage *)
      levelRelList = {};
      levelBindep = {};
      levelStableBound = {};
      
      (* ======================================================
         Step 4: Update Expansion Table
         ====================================================== *)
      expnTable = UpdateExpansionTable[
        expnTable, basePower, hexpnList, indexNew, vlist, order, modulus
      ];
      
      (* ======================================================
         Loop: Over coefficient degree
         ====================================================== *)
      Do[
        WithIndent[
          PrintV["=== Coef Degree: ", deg, " ==="];
          PrintV["+++ (lev,deg) = (",level,",",deg,") +++"];
          
          (* Remove solved variables *)
          {cfList, gvList, alphaList, bVars} = RemoveSolvedVariables[
            cfList0, rankList0, gvList0, levelBindep, bindep, deg
          ];
          PrintV["\[Alpha]List: ",alphaList,"\nb:  ",bVars];
          PrintV["    + RemoveSolvedVariables -- over"];
          {cfList, bVarReg, bSolNew} = InitializeVariables[cfList, bVars, vFixed];
          PrintV["    + InitializeVariables -- over"];
          
          If[Union[cfList] == {0},
            levelRelList = Append[levelRelList, {}];
            Continue[]
          ];
          
          (* Setup coefficient table *)
          bSol = bSolNew;
          coeffTable = SetupCoefficientTable[cfList, aregList, vlist, deg];
            (* Substitute known solutions into coefficient table *)
          coeffTable = Map[FFReduce[#, modulus]&@Expand@ReplaceAll[#, bSolNew] &, coeffTable, {3}];
          PrintV["    + SetupCoefficientTable -- over"];
          (* ======================================================
             Solving Loop: Call modular solver with state control
             ====================================================== *)
          {bSol, state, stableLevel} = MyTimer["Solve loop",
            SolveDegreeEquations[
              coeffTable, expnTable, alphaList, bVars,
              bVarReg, bSolNew, deg, order, level, vlist, 
              numeric, modulus, linearSolver,
              strategy, plateauSize, degreeByDegree
            ]
          ];
          PrintV["    + SolveDegreeEquations -- over"];
          
          (* Handle solver result *)
          If[bSol === <||>,
            PrintV["  Solver failed!"];
            Break[];
          ];
          
          (* Back-substitute *)
          {redRel, bVarIndep} = ResubstituteSolutions[
            bSol, gvList, cfList, bVars, topsec
          ];
          
          levelBindep = Join[levelBindep, bVarIndep];
          levelRelList = Append[levelRelList, redRel];
          levelStableBound = Append[levelStableBound, stableLevel];
          
          PrintV["  Relations found: ", Length[redRel]];
          PrintV["rel:  ",redRel//MatrixForm];
        ],  (* End WithIndent (degree) *)
        {deg, 0, maxCoefDeg}
      ];  (* End Do (degree) *)
      
      (* Accumulate results *)
      bindep = Append[bindep, levelBindep];
      relationList = Append[relationList, levelRelList];
      stableBound = Append[stableBound, levelStableBound];
    ],  (* End WithIndent (level) *)
    {level, levelList}
  ];  (* End Do (level) *)
  
  (* Output summary *)
  PrintV["=== Summary ==="];
  Print["Dimensions: ", Map[Length, relationList, {2}]//MatrixForm];
  Print["Stability bounds: ", stableBound//MatrixForm];
  
  relationList
];


End[]

EndPackage[]
