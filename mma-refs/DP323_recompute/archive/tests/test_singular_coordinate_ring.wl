(* test_singular_coordinate_ring.wl — Validate SingularCoordinateRing against MMA bivarPrimeInfo *)
(* Usage: wolframscript -file test_singular_coordinate_ring.wl *)

$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyGeneratePath = $LIECPPPath <> "/verify/DP323/";

(* Load packages *)
SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

(* Access private bivarPrimeInfo for comparison *)
bivarPrimeInfoMMA = Symbol["LIECoreAlgebra`Private`bivarPrimeInfo"];

(* Load checkpoint *)
Print["Loading DP323 checkpoint..."];
Block[{Print = If[#[[1]] === "Warning", Print[#]&, Identity]&},
  data = Import[$FamilyGeneratePath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
fullSectorList = data["Family", "SectorList"];

(* Find a completed simple sector: 2 props with 1 region *)
existingSectors = Keys[data["Regions"]];
simpleSectors = Select[existingSectors, Total[#] == 2 &];
Print["Found ", Length[simpleSectors], " completed sectors with 2 props"];
If[Length[simpleSectors] == 0,
  Print["ERROR: No completed 2-prop sectors found"];
  Exit[1];
];

(* Helper: replicate expRegSolve2 prime extraction *)
getPrimesForSector[sec_] := Module[{
  ne = Length@Alist, aeqs0, aeqs, avar, agb, agbSymbols, agbA, agbB,
  localA, localB, localRep, agbGlobal, avarGlobal, aprimelist, outRep, resultStr
},
  aeqs0 = Coefficient[ibpeqs, "n"] /. "g"[a__] :> (Times @@ Thread@Power[Alist /. A[i_] :> "A"[i], {a} - vlist]);
  aeqs = Join[aeqs0 /. {1/"A"[a_] :> "B"[a]}, Table["A"[i] "B"[i] - 1, {i, ne}]];
  avar = {Alist /. "A"[a_] :> "B"[a], Alist /. A[i_] :> "A"[i], {}};

  agb = GroebnerBasis[aeqs, Join @@ avar, Modulus -> char];

  localA = Table["A"[i], {i, ne}];
  localB = Table["B"[i], {i, ne}];
  agb = agb /. {A[i_] :> "A"[i], B[i_] :> "B"[i]};
  agbSymbols = Union@Cases[agb, h_[i_Integer] /; StringQ[h], Infinity];
  agbA = SortBy[Select[agbSymbols, #[[0]] === "A" &], #[[1]] &];
  agbB = SortBy[Select[agbSymbols, #[[0]] === "B" &], #[[1]] &];

  localRep = Join[
    If[Length[agbB] > 0, Thread[agbB -> localB[[;; Length[agbB]]]], {}],
    If[Length[agbA] > 0, Thread[agbA -> localA[[;; Length[agbA]]]], {}]
  ];
  agbGlobal = agb /. localRep;
  avarGlobal = {localB[[;; Length[agbB]]], localA[[;; Length[agbA]]], {}};

  aprimelist = SingularMinAssPrime[agbGlobal, Join @@ avarGlobal, Modulus -> char, "MonomialOrder" -> "lp"];

  (* Fix string replacement back to MMA symbols *)
  outRep = Join[
    Table["x" <> ToString[i] -> "\"A\"[" <> ToString[i] <> "]", {i, Length[localA]}],
    Table["x" <> ToString[Length[localA] + i] -> "\"B\"[" <> ToString[i] <> "]", {i, Length[localB]}]
  ];
  resultStr = ToString[aprimelist, InputForm];
  resultStr = StringReplace[resultStr, outRep];
  aprimelist = ToExpression[resultStr];
  aprimelist = aprimelist /. (Reverse /@ localRep);

  (* Filter zero-dimensional primes *)
  aprimelist = GroebnerBasis[#, Join @@ avar, Modulus -> char] & /@ aprimelist;
  aprimelist
];

(* Helper: apply sector limit *)
applySectorLimit[sec_] := Module[{
  ne = Length@Alist, aeqs0, aeqs, avar
},
  aeqs0 = Coefficient[ibpeqs, "n"] /. "g"[a__] :> (Times @@ Thread@Power[Alist /. A[i_] :> "A"[i], {a} - vlist]);
  aeqs = Join[aeqs0 /. {1/"A"[a_] :> "B"[a]}, Table["A"[i] "B"[i] - 1, {i, ne}]];
  avar = {Alist /. "A"[a_] :> "B"[a], Alist /. A[i_] :> "A"[i], {}};
  {LIECoreAlgebra`Private`sectorLimitIBP[aeqs, sec, vlist], avar}
];

(* Test on first simple sector *)
sec = simpleSectors[[1]];
Print["\n========================================"];
Print["Testing sector: ", sec, " (", Total[sec], " props)"];
Print["========================================"];

{ibpeqssub, avar} = applySectorLimit[sec];
Print["Sector-limited IBPEqs: ", Length[ibpeqssub], " equations"];

(* Get primes directly from the limited equations *)
ne = Length@Alist;
agb = GroebnerBasis[ibpeqssub, Join @@ avar, Modulus -> char];
Print["GB size: ", Length[agb]];

localA = Table["A"[i], {i, ne}];
localB = Table["B"[i], {i, ne}];
agb2 = agb /. {A[i_] :> "A"[i], B[i_] :> "B"[i]};
agbSymbols = Union@Cases[agb2, h_[i_Integer] /; StringQ[h], Infinity];
agbA = SortBy[Select[agbSymbols, #[[0]] === "A" &], #[[1]] &];
agbB = SortBy[Select[agbSymbols, #[[0]] === "B" &], #[[1]] &];

localRep = Join[
  If[Length[agbB] > 0, Thread[agbB -> localB[[;; Length[agbB]]]], {}],
  If[Length[agbA] > 0, Thread[agbA -> localA[[;; Length[agbA]]]], {}]
];
agbGlobal = agb2 /. localRep;
avarGlobal = {localB[[;; Length[agbB]]], localA[[;; Length[agbA]]], {}};

Print["Variables: ", avarGlobal];
Print["Running SingularMinAssPrime..."];
rawResult = SingularMinAssPrime[agbGlobal, Join @@ avarGlobal, Modulus -> char, "MonomialOrder" -> "lp"];

outRep = Join[
  Table["x" <> ToString[i] -> "\"A\"[" <> ToString[i] <> "]", {i, Length[localA]}],
  Table["x" <> ToString[Length[localA] + i] -> "\"B\"[" <> ToString[i] <> "]", {i, Length[localB]}]
];
resultStr = ToString[rawResult, InputForm];
resultStr = StringReplace[resultStr, outRep];
rawResult = ToExpression[resultStr];
rawResult = rawResult /. (Reverse /@ localRep);

(* Extract zero-dimensional primes, matching expRegSolve2 logic *)
aprimelist = rawResult[[1]];
dimlist = rawResult[[2]];
zeroDimIdx = Select[Range@Length[aprimelist], dimlist[[#]] == 0 &];
aprimelist = aprimelist[[#]] & /@ zeroDimIdx;
Print["Zero-dim primes: ", Length[aprimelist], " (from ", Length[dimlist], " total)"];

aprimelist = GroebnerBasis[#, Join @@ avar, Modulus -> char] & /@ aprimelist;
Print["After GB refine: ", Length[aprimelist]];
Do[Print["  Prime ", i, ": ", Length[aprimelist[[i]]], " gens"], {i, Length[aprimelist]}];

(* Run MMA bivarPrimeInfo *)
Print["\n--- Running MMA bivarPrimeInfo ---"];
timeMMA = AbsoluteTiming[
  resMMA = bivarPrimeInfoMMA[aprimelist, avar[[2]], avar[[1]],
    "FractionSymbol" -> "B", Modulus -> char, "Parameters" -> avar[[3]], Verbose -> True];
][[1]];
Print["MMA time: ", timeMMA, "s"];

(* Run Singular bivarPrimeInfo *)
Print["\n--- Running Singular bivarPrimeInfo ---"];
timeSingular = AbsoluteTiming[
  resSingular = SingularCoordinateRing`singularBivarPrimeInfo[aprimelist, avar[[2]], avar[[1]],
    "FractionSymbol" -> "B", Modulus -> char, "Parameters" -> avar[[3]], "LimitSector" -> sec, Verbose -> True];
][[1]];
Print["Singular time: ", timeSingular, "s"];

(* Compare results *)
Print["\n========================================"];
Print["Comparison Results"];
Print["========================================"];

If[Length[resMMA] =!= Length[resSingular],
  Print["FAIL: Different number of regions: MMA=", Length[resMMA], ", Singular=", Length[resSingular]];
  Exit[1];
];

keysToCheck = {"VarDep", "MinPoly", "VarIndep", "VarRule", "FractionRule",
  "MonomialBasis", "MonomialBasisIndex", "MonomialBasisMatrix", "VarDeg"};

allPass = True;
Do[
  regMMA = resMMA[[i]];
  regS = resSingular[[i]];
  Print["\nRegion ", i, ":"];

  Do[
    key = keysToCheck[[k]];
    valMMA = regMMA[key];
    valS = regS[key];

    (* Normalize for comparison: sort rules, etc. *)
    If[Head[valMMA] === List && Head[valS] === List,
      valMMA = Sort[valMMA];
      valS = Sort[valS];
    ];
    If[Head[valMMA] === Association && Head[valS] === Association,
      valMMA = KeySort[valMMA];
      valS = KeySort[valS];
    ];

    If[valMMA =!= valS,
      Print["  FAIL: ", key];
      Print["    MMA:     ", valMMA];
      Print["    Singular:", valS];
      allPass = False;
      ,
      Print["  PASS: ", key];
    ];
  , {k, Length[keysToCheck]}];
, {i, Length[resMMA]}];

Print["\n========================================"];
If[allPass,
  Print["ALL CHECKS PASSED"];
  ,
  Print["SOME CHECKS FAILED — see above"];
];
Print["========================================"];
