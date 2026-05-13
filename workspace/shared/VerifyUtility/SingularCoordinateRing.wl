(* ::Package:: *)
(* SingularCoordinateRing.wl — Replace MMA CoordinateRing construction with Singular *)
(*
   Problem: bivarPrimeInfo in LIECoreAlgebra.wl spends >1200s in MMA on:
     1. quotientRingBasisMatrixPower (PolynomialReduce for monomial basis matrix)
     2. GroebnerBasis for FractionRule
   Solution: Singular computes GB, kbase, multiplication table, fraction GB.
              MMA only does fast assembly (no PolynomialReduce needed).
*)

BeginPackage["SingularCoordinateRing`", {"LIEUtility`"}]

singularBivarPrimeInfo::usage = "singularBivarPrimeInfo[primelist, Avar, Bvar, opts] uses Singular to compute CoordinateRing data. Returns same format as LIECoreAlgebra`bivarPrimeInfo.";

Begin["`Private`"]

(* Load SingularInterface *)
If[!ValueQ[SingularRun],
  Get["SingularInterface.wl"];
];

(* ============================================================ *)
(* Helper: monomial exponent vector *)
(* ============================================================ *)

monomialExponents[mon_, vars_] := Exponent[mon, #] & /@ vars;

(* ============================================================ *)
(* MMA-side assembly from Singular raw output *)
(* ============================================================ *)

assembleCoordinateRing[prime_, gb_List, kb_List, mbPolys_List, Avar_, Bvar_, char_, limitSector_,
                        FRsym_: "B", params_: {}] := Module[{
    allVars, ltlist, lmlist, varclass, vargen, varpar,
    ximinpoly, xdrule,
    basisIndex, basis, basisMatrix,
    frrule,
    result, mbCount
  },

  allVars = Join[Bvar, Avar, params];

  (* Filter to A-only polynomials, matching bivarPrimeInfo's primeA step.
     Mixed A+B polynomials confuse Solve for 2D regions and produce
     inconsistent equation sets that break the VarRule computation. *)
  aOnlyGb = Select[gb, Cases[#, Alternatives @@ Bvar, Infinity] === {} &];

  (* 1. Leading term analysis (fast pattern matching, no GB/Solve needed) *)
  ltlist = First@MonomialList[#, Avar, Lexicographic] & /@ aOnlyGb;
  lmlist = {Union@Cases[#, Alternatives @@ Avar, {0, Infinity}],
      PolynomialDeg[#, Avar]} & /@ ltlist;

  varclass = If[!MemberQ[Keys[#], True], Join[#, <|True -> {}|>], #] &@
    GroupBy[lmlist, #[[2]] > 1 &];
  {vargen, varpar} = {Sort[Join @@ (First /@ #[True])],
     Join @@ (First /@ #[False])} &@varclass;

  (* 2. MinPoly: GB generators with only independent variables and degree > 1 *)
  ximinpoly = aOnlyGb[[#]] & @ Select[Range@Length@aOnlyGb,
     lmlist[[#, 2]] > 1 && SubsetQ[vargen, lmlist[[#, 1]]] &];

  (* 3. VarRule: solve for dependent vars from GB (fast: linear equations)
     Singular lp-order GB may produce polynomials that Solve cannot handle
     directly in modular arithmetic. Recompute with MMA default order to
     obtain triangulated linear forms that Solve can process.
     Use A-only GB to avoid inconsistent mixed-A+B equations. *)
  mmaGb = GroebnerBasis[aOnlyGb, Join[Bvar, Avar, params], Modulus -> char];
  xdrule = Solve[# == 0 & /@ Select[mmaGb, Cases[#, Alternatives @@ varpar, Infinity] =!= {} &], varpar, Modulus -> char];
  (* Defensive: handle empty or inconsistent solutions gracefully *)
  If[xdrule === {} || xdrule === {{}},
    Print["assembleCoordinateRing: Solve returned empty for prime ", prime,
      " varpar=", varpar, ", mmaGb=", mmaGb];
    xdrule = {};
    ,
    xdrule = xdrule[[1]];
  ];
  xdrule = Map[#[[1]] -> PolynomialReduce[#[[2]], ximinpoly, vargen,
       MonomialOrder -> Lexicographic, Modulus -> char][[2]] &, xdrule];

  (* 4. MonomialBasis from kbase *)
  basisIndex = monomialExponents[#, vargen] & /@ kb;
  mbCount = Length[kb];

  (* Sort basis to match MMA quotientRingBasisPower lexicographic order.
     Singular kbase order may differ from MMA, causing MonomialBasisMatrix
     key order mismatches that propagate to binary export differences. *)
  Module[{sortOrder},
    sortOrder = Ordering[basisIndex];
    basisIndex = basisIndex[[sortOrder]];
    kb = kb[[sortOrder]];
    (* mbPolys is flat: (j-1)*mbCount + i for old ordering.
       Reorder both j and i dimensions to match new kb order. *)
    mbPolys = Flatten@Table[
      mbPolys[[(sortOrder[[j]] - 1)*mbCount + sortOrder[[i]]]],
      {j, mbCount}, {i, mbCount}
    ];
  ];
  basis = Times @@ Thread[Power[vargen, #]] & /@ basisIndex;

  (* 5. MonomialBasisMatrix from Singular multiplication table *)
  basisMatrix = Association@Table[
    basisIndex[[j]] -> Transpose@Table[
      Module[{poly, coeffs},
        poly = mbPolys[[(j - 1)*mbCount + i]];
        coeffs = If[poly === 0,
          Table[0, mbCount],
          Lookup[monomialRulesPower[poly, vargen], basisIndex, 0]
        ];
        coeffs
      ],
      {i, mbCount}
    ],
    {j, mbCount}
  ];

  (* 6. FractionRule — match bivarPrimeInfo exactly using PolynomialReduce *)
  frrule = Join @@ Table[
    If[i =!= j,
      FRsym[i, j] -> PolynomialReduce[Avar[[i]]*Bvar[[j]], prime,
        Join[params, Bvar, Avar], Modulus -> char][[2]],
      Nothing
    ],
    {i, Length@Avar}, {j, Length@Avar}
  ];

  (* 7. Assemble — evaluate all values explicitly before wrapping in Association
     to prevent serialization issues with Module-local variables. *)
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

(* ============================================================ *)
(* Singular template for single-prime CoordinateRing computation *)
(* ============================================================ *)

singularPrimeTemplate = StringTemplate[
 "
LIB \"primdec.lib\";
ring r = `char`, (`vars`), (`order`);
option(prot);
ideal I = `ideal`;
option(redSB);
ideal gb = std(I);
ideal kb = kbase(gb);
int vd = vdim(gb);

string gbs = string(gb[1]);
for (int j = 2; j <= size(gb); j++) { gbs = gbs + \",\" + string(gb[j]); }

string kbs = string(kb[1]);
for (int j = 2; j <= size(kb); j++) { kbs = kbs + \",\" + string(kb[j]); }

string mbs = \"\";
for (int a = 1; a <= vd; a++)
{
  for (int b = 1; b <= vd; b++)
  {
    poly p = reduce(kb[a]*kb[b], gb);
    mbs = mbs + string(p);
    if (a < vd || b < vd) { mbs = mbs + \",\"; }
  }
}
`outvars`
"
];

(* ============================================================ *)
(* Main: singularBivarPrimeInfo *)
(* ============================================================ *)

Options[singularBivarPrimeInfo] = {
  Modulus -> 0,
  "FractionSymbol" -> "B",
  "Parameters" -> {},
  "LimitSector" -> {},
  Verbose -> False,
  "MonomialOrder" -> "lp"
};

singularBivarPrimeInfo[primelist_List, Avar_List, Bvar_List, OptionsPattern[]] := Module[{
    char = OptionValue[Modulus],
    params = OptionValue["Parameters"],
    limitSector = OptionValue["LimitSector"],
    FRsym = OptionValue["FractionSymbol"],
    order, varrep, posrep,
    result = {}, i,
    prime, gbList, kbList, mbList,
    allVars
  },

  If[Length[primelist] == 0, Return[{}]];

  allVars = Join[Bvar, Avar, params];

  (* Generate ring declaration *)
  {order, varrep, posrep} = SingularRingDeclaration[
    {allVars}, {{1}},
    "PositionOverTerm" -> False,
    "MonomialOrder" -> OptionValue["MonomialOrder"]
  ];

  If[OptionValue[Verbose],
    Print["[SingularCoordinateRing] Processing ", Length[primelist], " prime component(s)..."];
  ];

  (* Process each prime individually to avoid complex list serialization *)
  Do[
    prime = primelist[[i]];
    If[char =!= 0, prime = PolynomialMod[prime, char]];

    If[OptionValue[Verbose],
      Print["  Prime ", i, "/", Length[primelist], ": ", Length[prime], " generator(s)"];
    ];

    {gbList, kbList, mbList} = SingularRun[
      singularPrimeTemplate,
      <|
        "char" -> char,
        "vars" -> StringJoin@Riffle[varrep[[All, 2]], ","],
        "order" -> order,
        "ideal" -> prime,
        "outvars" -> {"gbs", "kbs", "mbs"}
      |>,
      "Replacement" -> Join[posrep, varrep]
    ];

    (* SingularRun returns {gbs, kbs, mbs} where each is already parsed by ToExpression *)
    (* gbs, kbs, mbs are MMA lists of polynomials with original variable names restored by rep2 *)

    reg = assembleCoordinateRing[
      primelist[[i]], gbList, kbList, mbList,
      Avar, Bvar, char, limitSector, FRsym, params
    ];
    AppendTo[result, reg];
  , {i, Length[primelist]}];

  Return[result]
];

End[]
EndPackage[]
