(* debug_mbm.wl — Debug MonomialBasisMatrix construction *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["Box"];

workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sector = {1, 1, 1, 1};

(* Temporarily patch assembleCoordinateRing to print debug info *)
Unprotect[assembleCoordinateRing];
assembleCoordinateRing[prime_, gb_List, kb_List, mbPolysIn_List, Avar_, Bvar_, char_, limitSector_, FRsym_: "B", params_: {}] := Module[{
    allVars, ltlist, lmlist, varclass, vargen, varpar,
    ximinpoly, xdrule,
    basisIndex, basis, basisMatrix,
    frrule,
    result, mbCount, mbPolys = mbPolysIn
  },

  allVars = Join[Bvar, Avar, params];
  aOnlyGb = Select[gb, Cases[#, Alternatives @@ Bvar, Infinity] === {} &];
  ltlist = First@MonomialList[#, Avar, Lexicographic] & /@ aOnlyGb;
  lmlist = {Union@Cases[#, Alternatives @@ Avar, {0, Infinity}],
      PolynomialDeg[#, Avar]} & /@ ltlist;
  varclass = If[!MemberQ[Keys[#], True], Join[#, <|True -> {}|>], #] &@
    GroupBy[lmlist, #[[2]] > 1 &];
  {vargen, varpar} = {Sort[Join @@ (First /@ #[True])],
     Join @@ (First /@ #[False])} &@varclass;
  ximinpoly = aOnlyGb[[#]] & @ Select[Range@Length@aOnlyGb,
     lmlist[[#, 2]] > 1 && SubsetQ[vargen, lmlist[[#, 1]]] &];
  mmaGb = GroebnerBasis[aOnlyGb, Join[Bvar, Avar, params], Modulus -> char];
  xdrule = Solve[# == 0 & /@ Select[mmaGb, Cases[#, Alternatives @@ varpar, Infinity] =!= {} &], varpar, Modulus -> char];
  If[xdrule === {} || xdrule === {{}},
    Print["assembleCoordinateRing: Solve returned empty"];
    xdrule = {};
    , xdrule = xdrule[[1]];
  ];
  xdrule = Map[#[[1]] -> PolynomialReduce[#[[2]], ximinpoly, vargen,
       MonomialOrder -> Lexicographic, Modulus -> char][[2]] &, xdrule];

  (* Debug MBM *)
  basisIndex = monomialExponents[#, vargen] & /@ kb;
  mbCount = Length[kb];
  Print["DEBUG: kb=", kb, " basisIndex=", basisIndex, " mbCount=", mbCount];
  Print["DEBUG: mbPolys (raw)=", mbPolys];
  
  Module[{sortOrder},
    sortOrder = Ordering[basisIndex];
    Print["DEBUG: sortOrder=", sortOrder];
    basisIndex = basisIndex[[sortOrder]];
    kb = kb[[sortOrder]];
    mbPolys = Flatten@Table[
      mbPolys[[(sortOrder[[j]] - 1)*mbCount + sortOrder[[i]]]],
      {j, mbCount}, {i, mbCount}
    ];
    Print["DEBUG: mbPolys (sorted)=", mbPolys];
  ];
  basis = Times @@ Thread[Power[vargen, #]] & /@ basisIndex;
  
  basisMatrix = Association@Table[
    basisIndex[[j]] -> Transpose@Table[
      Module[{poly, coeffs},
        poly = mbPolys[[(j - 1)*mbCount + i]];
        Print["DEBUG: j=", j, " i=", i, " poly=", poly];
        coeffs = If[poly === 0,
          Table[0, mbCount],
          Lookup[monomialRulesPower[poly, vargen], basisIndex, 0]
        ];
        Print["DEBUG: coeffs=", coeffs];
        coeffs
      ],
      {i, mbCount}
    ],
    {j, mbCount}
  ];
  Print["DEBUG: basisMatrix=", basisMatrix];

  frrule = Join @@ Table[
    If[i =!= j,
      FRsym[i, j] -> PolynomialReduce[Avar[[i]]*Bvar[[j]], prime,
        Join[params, Bvar, Avar], Modulus -> char][[2]],
      Nothing
    ],
    {i, Length@Avar}, {j, Length@Avar}
  ];

  varDepValue = Avar /. xdrule;
  varDegValue = Select[lmlist, Cases[#[[1]], Alternatives @@ Avar, Infinity] =!= {} &][[;; , 2]];
  result = <|"VarDep" -> varDepValue,
    "MinPoly" -> ximinpoly,
    "VarIndep" -> vargen,
    "VarRule" -> xdrule,
    "FractionRule" -> frrule,
    "MonomialBasis" -> basis,
    "MonomialBasisIndex" -> basisIndex,
    "MonomialBasisMatrix" -> basisMatrix,
    "VarDeg" -> varDegValue
  |>;
  Return[result]
];

(* Run just region 2 of Box *)
singularRes = regionsBySectors[ibpeqs, {sector}, Alist, vlist,
  Modulus -> char, Verbose -> False];

Print["=== Region 2 debug ==="];
ring = singularRes[sector][[2]]["CoordinateRing"];
Print["VarIndep=", ring["VarIndep"]];
Print["MonomialBasisIndex=", ring["MonomialBasisIndex"]];
Print["MonomialBasisMatrix keys=", Keys[ring["MonomialBasisMatrix"]]];
Print["MonomialBasisMatrix values=", Values[ring["MonomialBasisMatrix"]]];
