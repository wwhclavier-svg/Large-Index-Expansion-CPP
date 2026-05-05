(* ::Package:: *)

(* Large Index Expansion - Recursion Expansion Package (renamed) *)
(* Reconstructs series coefficients using pre-computed region data *)


BeginPackage["LIEExpand`", {"LIEUtility`"}]

(* Public Interface *)
LIEExpand::usage = "LIEExpand[order, ibpeqs, expregdata, opts] reconstructs series expansion coefficients.\n"<>
  "expregdata should be the output from LIERegions`regionsBySectors.\n"<>
  "Returns {hlist, alist} where hlist contains expansion coefficients and alist contains region information."

Begin["`Private`"]

(* ============================================================ *)
(* Dependency Loading *)
(* ============================================================ *)

(* Load utility functions *)
If[!ValueQ[reduceSolve],
  Get["LIEUtility.wl"];
];

(* ============================================================ *)
(* Utility Functions *)
(* ============================================================ *)

lastNonZero[list_]:=Module[{i=Length[list]},While[list[[i]]===0&&i>=1,i--];If[i<1,Return[0],Return[i]]]

Options[isIncompatiable]={Modulus->0};
isIncompatiable[M_,N_,OptionsPattern[]]:=Module[{char=OptionValue[Modulus]},
Join[M,N,2]//RowReduce[#,Modulus->char]&//DeleteCases[#,a_/;Union[a]=={0}]&//Select[firstNonZero/@#,#>Length[M[[1]]]&]=!={}&];

(* turn expansions with indeterminates to vectors spaning the solution space *)
LIESolutionSpace[expn_,vlist_]:=Module[{monovals,monokeys},
{monokeys,monovals}=expn//monomialRulesPower[#,vlist]&//{Reverse@Keys[#],Reverse@Values[#]}&;
monovals=If[Length[#]==2,Transpose@Join[List/@#[[1]],#[[2]],2],{#[[1]]}]&@Normal@CoefficientArrays[#,Union@Cases[#,a_/;Head[a]===c,Infinity]]&@monovals;
Return[{monokeys,monovals}];
];


(* ============================================================ *)
(* Recursion Expansion Functions *)
(* ============================================================ *)

computeFlags[M1_, N1_, nb_, ne_, char_] :=
 Module[{M1mat, N1mat, AnyNonZeroQ},
  M1mat = Flatten[M1, {{1,3},{2,4}}];
  N1mat = Flatten[N1, {{1,3},{2,4}}];
  AnyNonZeroQ[mat_] := AnyTrue[If[Head[mat]===SparseArray,mat["ExplicitValues"],Flatten[mat]], # != 0 &];

  <|
   "FullRank"->(MatrixRank[M1mat,Modulus->char]>=ne*nb),
   "Incompatiable"->isIncompatiable[M1mat,N1mat,Modulus->char],
   "N1Zero"->(!AnyNonZeroQ[N1mat])
  |>
 ]


SetAttributes[chooseIncrement, HoldFirst];
chooseIncrement[data_, flags_, incre0_, PrintL_] :=
Module[{incre=2},
PrintL["flag(M1 non-deg.):  ",flags["FullRank"]];PrintL["flag(incomp.):  ",flags["Incompatiable"]];PrintL["flag(N1==0):  ",flags["N1Zero"]];
incre=Which[incre0=!=2, incre=incre0,
flags["Incompatiable"]==True, incre=3,
flags["FullRank"]==True&&flags["N1Zero"]==True, incre=1,
True,2];
PrintL["incre = ",incre];data["incre"]= incre;
incre
]


SetAttributes[recursionLoop, HoldFirst];
recursionLoop[data_, ord_, char_, FFReduce_, PrintL_, options___] := Module[{
orderStartTime,orderDuration,vlist,
    incre = data["incre"], status, hsollist, csollist, cindepall, csol, cindep,
    ne = data["ne"], nb = data["nb"], nibp = data["nibp"],
    layerByLayer = Lookup[{options}, "LayerByLayer", True],
    zeroGen = Lookup[{options}, "ZeroGeneralSolution", False],
    linearSolver=Lookup[{options}, "linearSolver", reduceSolve],
    outputFormat=Lookup[{options}, "OutputFormat", "Vector"]
    (* "Vector" take components in cindep, "Polynomial" original form *)
  },

  Print["[DEBUG] char=", char, ", FFReduce test: ", FFReduce[1788511843969144]];
  status = "restart";
  vlist = Table["v" <> IntegerString[i], {i, ne}];
  
  
  While[status === "restart",
  incre = data["incre"];
    status = Catch[
      (* --- 1. \:521d\:59cb\:5316\:6bcf\:4e00\:8f6e\:5c1d\:8bd5\:7684\:53d8\:91cf --- *)
      ClearAll[c];
      c[i_, seed_] := Table[c[i, seed, j], {j, nb}];
      
      csol = Association@Table[c[0, Table[0, ne], j] -> If[j == 1, 1, 0], {j, nb}];
      cindepall = Table[c[0, Table[0, ne], j], {j, nb}];
      hsollist = {Table[If[j == 1, 1, 0], {j, nb}]};
      

      (* --- 2. \:9636\:6570\:4e3b\:5faa\:73af --- *)
      Do[
      orderStartTime = AbsoluteTime[];
        PrintL["====== Order: " <> ToString[i] <> " (Incre: " <> ToString[incre] <> ") ======"];
        
        (* \:751f\:6210\:5f53\:524d\:9636\:6570\:7684\:6240\:6709 Seed *)
        With[{seeds = Reverse@SortBy[seedsGenPos[incre*i - 1, ne], {Total[#], -Reverse[#]} &]},
          
          (* \:79fb\:9664\:53d7\:7eaf n \:56e0\:5b50\:5f71\:54cd\:7684\:65e7\:89e3 (cindepremove) *)
          With[{oldIndep = DeleteDuplicates[#[[2]] & /@ cindepall]},
          PrintL["remove old solution: ",Join @@ (Association@Table[c[i, #, j] -> 0, {j, nb}] & /@ oldIndep)];
            csol=AssociateTo[csol, Join @@ (Association@Table[c[i, #, j] -> 0, {j, nb}] & /@ oldIndep)];
          ];

          If[layerByLayer,
            (* \:6a21\:5f0f A: \:9010\:4e2a Seed \:6c42\:89e3 (Layer by Layer) *)
            Do[
              csol = WithIndent@solveSeedBlock[data, {s}, i, csol, incre, char, linearSolver, FFReduce, PrintL],
              {s, seeds}
            ],
            (* \:6a21\:5f0f B: \:6309\:5c42(Total Seed) \:6279\:91cf\:6c42\:89e3 *)
            Do[
              With[{layerSeeds = Select[seeds, Total[#] == SL &]},
                csol = WithIndent@solveSeedBlock[data, layerSeeds, i, csol, incre, char, linearSolver, FFReduce, PrintL]
              ],
              {SL, Reverse@Range[0, incre*i - 1]}
            ]
          ];
        ];

        (* --- 3. \:6bcf\:4e00\:9636\:7ed3\:675f\:540e\:7684\:6536\:5c3e (Coda) --- *)
        With[{res = updateCoda[i, csol, incre, ne, nb, vlist, cindepall, zeroGen, FFReduce]},
          csol = res["csol"];
          cindepall = res["cindepall"];
          AppendTo[hsollist, res["hAtOrder"]];
          orderDuration = AbsoluteTime[] - orderStartTime;
          
          PrintL["Expansion at order " <> ToString[i] <> " finished. "];
          PrintL["Pivot c: " <> ToString[cindepall]];
          Print["Time taken for Order " <> ToString[i] <> ": " <> ToString[UnitConvert[Quantity[orderDuration, "Seconds"], "s"]]];
          PrintL[Last@hsollist]
        ];

      , {i, 1, ord}];
      
      "success"
    ]; (* Catch End *)
  ]; (* While End *)
  
  If[outputFormat === "Vector",
  cindep = Select[cindepall,#[[3]]==1&&AnyTrue[#[[2]],#=!=0&]&];
  hsollist = If[cindep==={},
  {hsollist},
  CoefficientArrays[hsollist, cindep]//{#[[1]], Sequence@@Transpose[#[[2]], {2,3,1}]}&];
  ];
  
  If[status === "success", hsollist, $Failed]
]

(* --- \:6838\:5fc3\:6c42\:89e3\:5b50\:51fd\:6570 --- *)
SetAttributes[solveSeedBlock, HoldFirst];
solveSeedBlock[data_, targetSeeds_, i_, currentCsol_, incre_, char_, linearSolver_, FFReduce_, PrintL_] := Module[{
    ceqs = {}, cvars, csolnew, 
    ne = data["ne"], nb = data["nb"], nibp = data["nibp"], csolUpdated
  },
  
  Do[
    (* \:8ba1\:7b97\:590d\:5408\:975e\:9f50\:6b21\:9879 *)
    PrintL["++++++  seed:  ", seed, "  ++++++"];
    PrintL["1. Inhomogenous Part: "];
    With[{inhomog = WithIndent@calculateTotalInhomog[seed, i, currentCsol, data, incre, FFReduce, If[$Verbose == True, PrintL, Nothing]]},
      (* \:627e\:51fa\:5f53\:524d Seed \:9700\:8981\:6c42\:89e3\:7684\:53d8\:91cf\:65b9\:5411 *)
      PrintL["nh = ", inhomog];
      With[{unknowns = Range[Max[lastNonZero@seed, 1], ne]},
        PrintL["  solving block:  ", unknowns];
        AppendTo[ceqs,
          Table[FFReduce@Sum[Dot[data["M1"][[m, j]], (seed[[j]] + 1) * c[i, seed + UnitVector[ne, j]]], {j, unknowns}] + inhomog[[m]], {m, nibp}]
        ];
      ];
    ];
  , {seed, targetSeeds}];

  ceqs = DeleteCases[Flatten[ceqs], 0];
  cvars = DeleteDuplicates@Cases[ceqs, a_/;Head[a]===c, {0, Infinity}]//Reverse@SortBy[#,{-#[[1]],Total[#[[2]]],-Reverse[#[[2]]]}&]&;
  PrintL["2. Set up Eqns:  "];
  PrintL["#ceqs:  ", Length@ceqs, "  #cvars:  ", Length@cvars];
  PrintL["   vars:  ", cvars];

  (* \:8c03\:7528\:7ebf\:6027\:6c42\:89e3\:5668 *)
  csolnew = ExecuteLinearSolve[ceqs, cvars, char, linearSolver];
  
  If[csolnew === {}, 
    Print["!!! Incompatible at current incre. Restarting with incre = ", data["incre"]=incre+1];
    Pause[1];
    Throw["restart"] 
  ];
  
  PrintL["3.  Eqns Solved:  "];
  PrintL["solved:  ",csolnew];

  (* \:66f4\:65b0\:89e3\:89c4\:5219\:5217\:8868 *)
  csolnew = Map[FFReduce[# /. currentCsol] &, csolnew];
  csolUpdated = Map[FFReduce[# /. csolnew] &, currentCsol];
  AssociateTo[csolUpdated, Association@csolnew]
]

(* --- \:975e\:9f50\:6b21\:9879\:7d2f\:52a0\:5668 --- *)
calculateTotalInhomog[seed_, i_, csol_, data_, incre_, FFReduce_] := Module[{
    knowns = Range[1,#-1]&@lastNonZero@seed
  },
  (* \:8fd9\:91cc\:7684\:6bcf\:4e2a\:5b50\:51fd\:6570\:5bf9\:5e94\:4f60\:539f\:59cb\:4ee3\:7801\:4e2d\:7684\:4e00\:6bb5\:6c42\:548c\:903b\:8f91 *)  
  FFReduce[
    GetInhomogNMinus[seed, i, incre, data["ne"], data["N1"], csol, FFReduce] +
    GetInhomogNZero[seed, i, incre, data["ne"], data["F0"], data["K1"], csol, FFReduce] +
    GetInhomogNPlusMinus[seed, i, incre, data["ne"], data["F2"], csol, FFReduce] +
    GetInhomogNPlus[seed, i, incre, data, csol, FFReduce] + (* \:90a3\:4e2a\:590d\:6742\:7684\:5f20\:91cf\:6c42\:548c *)
    GetInhomogM1[seed, i, data, csol, FFReduce, knowns] +
    GetInhomogMPlus[seed, i, incre, data, csol, FFReduce]
  ]//If[#===0,ConstantArray[0, data["nibp"]],#]&
]

calculateTotalInhomog[seed_, i_, csol_, data_, incre_, FFReduce_, PrintL_] := Module[{
    knowns = Range[1,#-1]&@lastNonZero@seed, eqnkey,
    inhomogNMinus,inhomogNZero,inhomogNPluMi,inhomogNPlus,inhomogM1,inhomogMPlus
  },
  
  {inhomogNMinus,inhomogNZero,inhomogNPluMi,inhomogNPlus,inhomogM1,inhomogMPlus}={
    GetInhomogNMinus[seed, i, incre, data["ne"], data["N1"], csol, FFReduce],
    GetInhomogNZero[seed, i, incre, data["ne"], data["F0"], data["K1"], csol, FFReduce],
    GetInhomogNPlusMinus[seed, i, incre, data["ne"], data["F2"], csol, FFReduce],
    GetInhomogNPlus[seed, i, incre, data, csol, FFReduce,"ListFunction"->Table],(* \:90a3\:4e2a\:590d\:6742\:7684\:5f20\:91cf\:6c42\:548c *)
    GetInhomogM1[seed, i, data, csol, FFReduce, knowns],
    GetInhomogMPlus[seed, i, incre, data, csol, FFReduce,"ListFunction"->Table]};
  
  (* Force print for k=2, l=1, seed={1,0} comparison with C++ *)
  If[i == 2 && Total[seed] == 1 && seed == {1, 0},
    Print["[MMA DETAIL] k=", i, ", l=", Total[seed], ", seed=", seed];
    Do[
      With[{nminus = FFReduce[If[Head[inhomogNMinus] === SparseArray, Flatten[Normal[inhomogNMinus]][[m]], inhomogNMinus]],
            nzero   = FFReduce[If[Head[inhomogNZero] === SparseArray, Flatten[Normal[inhomogNZero]][[m]], inhomogNZero]],
            nplumi  = FFReduce[If[Head[inhomogNPluMi] === SparseArray, Flatten[Normal[inhomogNPluMi]][[m]], inhomogNPluMi]],
            nplus   = Table[FFReduce[If[Head[inhomogNPlus[[nl]]] === SparseArray, Flatten[Normal[inhomogNPlus[[nl]]]][[m]], inhomogNPlus[[nl]]]], {nl, Length[inhomogNPlus]}],
            m1      = FFReduce[If[Head[inhomogM1] === SparseArray, Flatten[Normal[inhomogM1]][[m]], inhomogM1]],
            mplus   = Table[FFReduce[If[Head[inhomogMPlus[[nl]]] === SparseArray, Flatten[Normal[inhomogMPlus[[nl]]]][[m]], inhomogMPlus[[nl]]]], {nl, Length[inhomogMPlus]}]},
        Print["  m=", m-1, ": NMinus=", nminus, " NZero=", nzero, " NPluMi=", nplumi,
          " NPlus=", nplus, " M1=", m1, " MPlus=", mplus,
          " Total=", FFReduce[nminus + nzero + nplumi + Total[nplus] + m1 + Total[mplus]]];
      ];
    , {m, data["nibp"]}];
  ];
  
  (*eqnkey={"N(-)","N(0)","N(\[PlusMinus])","N(n>0)","M(+1)","M(n>1)"};*)
  eqnkey = Join[{"N(-)","N(0)","N(\[PlusMinus])"}, Table["N(+" <> ToString[j] <> ")", {j, Length@inhomogNPlus}], {"M(+1)"}, Table["M(+" <> ToString[j + 1] <> ")", {j, Length@inhomogMPlus}]];
  (*PrintL[(Join@@#&/@DeleteCases[#,a_/;a[[2]]==0]&@Transpose[{List/@eqnkey,
  {inhomogNMinus,inhomogNZero,inhomogNPluMi,If[Head[inhomogNPlus]===List,Sequence@@inhomogNPlus,Nothing],
  inhomogM1,If[Head[inhomogMPlus]===List,Sequence@@inhomogMPlus,Nothing]}}])//MatrixForm];*)
  
  FFReduce@Plus[inhomogNMinus,inhomogNZero,inhomogNPluMi,Sequence@@inhomogNPlus,inhomogM1,Sequence@@inhomogMPlus]//If[#===0,ConstantArray[0, data["nibp"]],#]&
]

(* --- \:6bcf\:4e00\:9636\:7684\:6536\:5c3e\:5904\:7406 --- *)
updateCoda[i_, csol_, incre_, ne_, nb_, vlist_, cindepall_, zeroGen_, FFReduce_] := Module[{
    seedsAtI, coefs, cindep, hAtOrder, newCsol = csol
  },
  seedsAtI = seedsGenPos[incre*i, ne];
  coefs = Table[c[i, s, j], {s, seedsAtI}, {j, nb}];
  cindep = Complement[Flatten[coefs], Keys@csol];
  
  If[zeroGen,
    newCsol = AssociateTo[csol /. #, #]&@Association[# -> 0 & /@ cindep];
  ];

  hAtOrder = Total[(coefs /. newCsol) * Map[Times @@ (vlist^#[[2]]) &, coefs, {2}]];
  (*hAtOrder = Total[ReplaceRepeated[#,newCsol]*(Times@@Thread@Power[vlist,#[[1,2]]])&/@coefs]*)

  <|"csol" -> newCsol, "cindepall" -> Join[cindepall, cindep], "hAtOrder" -> hAtOrder|>
]


(*SetAttributes[LookupOrSelf, HoldAll];*)
LookupOrSelf[assoc_, keys_List] := 
  Module[{results, heldList},
    (* \:6355\:83b7\:672a\:6c42\:503c\:7684\:539f\:59cb keys *)
    heldList = Unevaluated[keys];
    
    (* Lookup \:4f1a\:6c42\:503c keys \:8fdb\:884c\:67e5\:627e\:ff0c\:4f46 kept \:539f\:59cb\:5f15\:7528 *)
    results = heldList /. {k___} :> Lookup[assoc, {k}, $uniqueMissing];
    
    (* \:5355\:904d\:5408\:5e76\:ff1a\:7528 MapThread \:6bd4\:5bf9\:7ed3\:679c\:548c\:539f\:59cb\:503c *)
    MapThread[
      If[#1 === $uniqueMissing, #2, #1] &,
      {results, heldList}
    ]
  ]


GetInhomogNMinus[seed_, i_, incre_, ne_, N1_, csolAssoc_, FFReduce_] := 
 Module[{totalSeed, validJs=Select[Range[ne], seed[[#]] > 0 &], veclist, nb=Last@Dimensions@N1},
  totalSeed = Total[seed];
  If[totalSeed - 1 > incre * (i - 1), Return[0] ];
  veclist=
  FFReduce @ Sum[
    N1[[All, j]] . (LookupOrSelf[csolAssoc, c[i - 1, seed - UnitVector[ne, j]]]),
    {j, validJs}
  ]
]

GetInhomogNZero[seed_, i_, incre_, ne_, F0_, K1_, csolAssoc_, FFReduce_] := 
 Module[{totalSeed, sumK1, nb=Last@Dimensions@F0},
  totalSeed = Total[seed];
  
  (* \:6743\:91cd\:68c0\:67e5 *)
  If[totalSeed > incre * (i - 1), Return[0]];

  (* 1. \:6784\:9020\:5f53\:524d seed \:4e0b\:7684\:7cfb\:6570\:77e9\:9635\:548c (nibp x nb x nb) *)
  (* \:8fd9\:91cc\:7684 Sum \:540c\:6837\:662f\:5f20\:91cf\:95f4\:7684\:53e0\:52a0 *)
  sumK1 = Sum[
    If[seed[[j]] >= 1, seed[[j]] * K1[[All, j]], 0], 
    {j, ne}
  ];

  (* 2. (F0 + sumK1) \:662f\:4e00\:4e2a nibp x nb x nb \:7684\:5f20\:91cf *)
  (* \:4e0e c \:5411\:91cf\:70b9\:4e58\:5f97\:5230 nibp x nb \:7684\:7ed3\:679c\:77e9\:9635 *)
  FFReduce[(F0 + sumK1) . LookupOrSelf[csolAssoc, c[i - 1, seed]]]
]

GetInhomogNPlusMinus[seed_, i_, incre_, ne_, F2_, csolAssoc_, FFReduce_] := 
 Module[{totalSeed, validKs, nb=Last@Dimensions@F2},
  totalSeed = Total[seed];
  
  If[totalSeed > incre * (i - 1), Return[0]];

  (* \:627e\:51fa\:6240\:6709\:53ef\:4ee5\:51cf 1 \:7684\:5206\:91cf k *)
  validKs = Flatten[Position[seed, _?(# > 0 &)]];
  If[validKs === {}, Return[0]];

  (* \:53cc\:91cd\:6c42\:548c\:5f20\:91cf\:5316 *)
  FFReduce @ Sum[
    (* \:7cfb\:6570: -(seed[j]+1) *)
    (* \:5f20\:91cf: F2[[All, j, k]] (\:7ef4\:5ea6 nibp x nb x nb) *)
    -(seed[[j]] + 1) * F2[[All, j, k]] . (LookupOrSelf[csolAssoc, c[i - 1, seed + UnitVector[ne, j] - UnitVector[ne, k]]]),
    {j, ne}, 
    {k, validKs}
  ]
]

GetInhomogM1[seed_, i_, data_, csolAssoc_, FFReduce_, knowns_List] := 
 Module[{M1 = data["M1"], nb=data["nb"]},
  
  (* \:53ea\:6709\:5728 knowns \:4e0d\:4e3a\:7a7a\:65f6\:8ba1\:7b97 *)
  If[knowns === {}, Return[0]];

  (* \:8ba1\:7b97 M1[[All, j]] . c[i, seed + ej] *)
  FFReduce @ Sum[
    (seed[[j]] + 1) * M1[[All, j]] . (LookupOrSelf[csolAssoc,c[i, seed + UnitVector[data["ne"], j]]] ),
    {j, knowns}
  ]
]

GetInhomogNPlus[seed_, i_, incre_, data_, csolAssoc_, FFReduce_,options___] := Module[{
    ne = data["ne"], K1 = data["K1"], F2 = data["F2"],Func=Sum,
    maxL1 = incre*(i - 1) - Total[seed],nb=data["nb"],
    posK = Flatten[Position[seed, _?(# >= 1 &)]] (* \:9884\:5148\:627e\:51fa\:975e\:96f6\:5206\:91cf\:7d22\:5f15 *),
    SumFunc = Lookup[{options}, "ListFunction", Sum]
  },
  
  If[maxL1 < 1, Return[0]];

  FFReduce @ SumFunc[
    (* --- Part A: Direct Terms (\:5408\:5e76 K1 \:548c F2 \:7684\:7ebf\:6027\:8d21\:732e) --- *)
    Sum[
      With[{
        tensor = (If[MemberQ[posK, j], K1[[All, j]] * Binomial[seed[[j]] + l1, l1 + 1], 0] + 
                  (-1)^l1 * Binomial[seed[[j]] + l1, l1] * Sum[F2[[All, j, k]] * seed[[k]], {k, posK}])
      },
        If[tensor===0, 0, 
          tensor . (LookupOrSelf[csolAssoc, c[i - 1, seed + l1 * UnitVector[ne, j]]])
        ]
      ], 
      {j, ne}
    ] +

    (* --- Part B: Splitting Terms (\:5904\:7406 l2 \:8303\:56f4) --- *)
    Sum[
      With[{vec = LookupOrSelf[csolAssoc, c[i - 1, seed + l2*UnitVector[ne, j] + (l1 - l2)*UnitVector[ne, k]]]},
        ((-1)^l2 * Binomial[seed[[j]] + l2, l2] * Binomial[seed[[k]] + l1 - l2, l1 - l2 + 1]) * (F2[[All, j, k]] . vec)
      ],
      {j, ne}, {k, posK}, 
      {l2, Join[Range[l1 - 1], {l1 + 1}]} (* \:76f4\:63a5\:5728\:8fed\:4ee3\:5668\:4e2d\:5904\:7406\:975e\:8fde\:7eed\:8303\:56f4 *)
    ]
  , {l1, 1, maxL1}]
]

GetInhomogMPlus[seed_, i_, incre_, data_, csolAssoc_, FFReduce_, options___] := Module[{
    ne = data["ne"], K1s = data["K1s"], K2s = data["K2s"], F2s = data["F2s"], nb=data["nb"],
    maxL1 = incre*i - Total[seed], SumFunc = Lookup[{options}, "ListFunction", Sum], term1},

  If[maxL1 < 2, Return[0]];

  FFReduce @ SumFunc[
    (* --- Part A: Direct Terms --- *)
    Sum[
      ((K1s[[All, j]] + (-1)^l1 * K2s[[All, j]]) * Binomial[seed[[j]] + l1, l1]) . 
      (LookupOrSelf[csolAssoc, c[i, seed + l1*UnitVector[ne, j]]]),
      {j, ne}
    ] +

    (* --- Part B: Splitting Terms --- *)
    Sum[
      With[{vec = LookupOrSelf[csolAssoc, c[i, seed + l2*UnitVector[ne, j] + (l1 - l2)*UnitVector[ne, k]]]},
        ((-1)^l2 * Binomial[seed[[j]] + l2, l2] * Binomial[seed[[k]] + l1 - l2, l1 - l2]) * (F2s[[All, j, k]] . vec)
      ],
      {j, ne}, {k, ne}, {l2, 1, l1 - 1}
    ]
  , {l1, 2, maxL1}]
]


finalize[result_, data_, outputFormat_]:=Module[{res},
If[result=!=$Failed,
res = Switch[outputFormat, "Polynomial", (# . data["basis"])&/@result, _, result],
$Failed]];


 PrintL[expr___]:=If[$Verbose==True,Print[StringRepeat["    ", $levelCount], expr]];
 SetAttributes[WithIndent, HoldAll];
WithIndent[body_] := Block[{$levelCount=$levelCount+1},body];


Options[ExpandRecursive] = {
  Modulus -> 0,
  "LinearSolver" -> reduceSolve,
  "ZeroGeneralSolution" -> True,
  "LimitSector" -> {},
  Verbose -> False,
  "Increment" -> 2,
  "LayerByLayer" -> True,
  "OutputFormat" -> "Vector"(* / "Polynomial" *)
};
ExpandRecursive[recursionMatrixData_, order_, OptionsPattern[]] :=
 Module[{FFReduce, data, flags, incre, res, char=OptionValue[Modulus], M1, N1},
 FFReduce =Function[x, If[char =!= 0, PolynomialMod[x, char], x]];
 $levelCount=0;$Verbose=OptionValue[Verbose];
 (*PrintL = makePrintL[OptionValue[Verbose], levelCount];*)

  data = recursionMatrixData;
  Print["1 ...... Got Recursion Matrices in  Coordinate Ring basis"];

  flags = computeFlags[data["M1"],data["N1"],data["nb"],data["ne"],char];
  Print["2 ...... Got Rank Flags: ",flags];

  incre = WithIndent@chooseIncrement[data, flags, OptionValue["Increment"], PrintL];
  Print["3 ...... Initial Increment Determined (",incre,")"];
    
WithIndent@recursionLoop[data, order, char, FFReduce, PrintL]
 ]


(* ============================================================ *)
(* Main Public Function: LIEExpand *)
(* ============================================================ *)

Options[LIEExpand] = {
  Verbose -> False, Modulus -> 0, 
  "Increment" -> 2, "LayerByLayer" -> True, 
  "LinearSolver" -> reduceSolve
};

LIEExpand[order_, ibpeqs_, expregdata_, OptionsPattern[]] := 
 Module[{
    hexpnAssoc, alist, hlist,
    char = OptionValue[Modulus],
    solver = OptionValue["LinearSolver"],
    verbose = OptionValue[Verbose],
    incre = OptionValue["Increment"],
    lbl = OptionValue["LayerByLayer"]
  },
  
  (* 1. \:786e\:5b9a Regions \:6570\:636e\:5df2\:4ece\:53c2\:6570\:4f20\:5165 *)
  
  (* 2. \:4f7f\:7528 Table \:66ff\:4ee3 Append\:ff0c\:63d0\:9ad8\:5185\:5b58\:5206\:914d\:6548\:7387 *)
  Print["\n    =========      Reconstruct Expansion      ========="];
  hexpnAssoc = Association @ KeyValueMap[
    Function[{sec, data},
      Print["    +++ sector: ", sec, " (", Position[Keys[expregdata], sec][[1,1]], "/", Length[expregdata], ") +++"];
      
      (* \:5bf9\:6bcf\:4e00\:4e2a region \:8fdb\:884c\:5c55\:5f00 *)
      sec -> Table[With[{basis = data[[j]]["CoordinateRing"]["MonomialBasis"]},
        Print["      reg: (", j,"/",Length[data],")    basis(",Length@basis,"):   ",basis]];
        (* \:7ed3\:6784\:5316\:53d6\:503c\:ff1a#[[4]] \:662f recursionMatrix, #[[1]] \:662f expreg *)
        ExpandRecursive[
          data[[j]]["RecursionMatrix"], (* recursionMatrix *)
          order,
          "ZeroGeneralSolution" -> False,
          "LinearSolver" -> solver,
          Verbose -> verbose,
          "Increment" -> incre,
          Modulus -> char,
          "LayerByLayer" -> lbl
        ],
        {j, Length[data]}
      ]
    ],
    expregdata
  ];
  
  (* 3. \:7ed3\:679c\:5408\:5e76\:ff1a\:5229\:7528 Flatten \:548c Values \:7b80\:5316 *)
  alist = Flatten[Map[Lookup[#,"CoordinateRing",{}]&,#,{2}]&@Values[expregdata],1];
  hlist = Flatten[Values[hexpnAssoc],1];

  {hlist, alist}
  ];


End[] (* Private *)

EndPackage[]
