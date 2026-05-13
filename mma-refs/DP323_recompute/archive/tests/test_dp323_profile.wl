(* test_dp323_profile.wl — Profile each step of expRegSolve2 on DP323 sector 37 *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

Print["Loading DP323 checkpoint..."];
$FamilyPath = $LIECPPPath <> "/verify/DP323/";
Block[{Print = If[#[[1]] === "Warning", Print[#]&, Identity]&},
  data = Import[$FamilyPath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
ne = Length[Alist];

sec = {0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0};
Print["Sector: ", sec, " (", Total[sec], " props)"];

AlistCF = Alist /. A[i_] :> "A"[i];

(* Step 1: sector limit *)
t1 = AbsoluteTiming[
  ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, sec, vlist];
][[1]];
Print["1. sectorLimitIBP: ", Round[t1, 0.01], "s"];

(* Step 2: A/B conversion *)
t2 = AbsoluteTiming[
  aeqs0 = Coefficient[ibpeqslimited, "n"] /. "g"[a__] :> (Times @@ Thread@Power[AlistCF, {a} - vlist]);
  aeqs = Join[aeqs0 /. {1/"A"[a_] :> "B"[a]}, Table["A"[i] "B"[i] - 1, {i, ne}]];
  avar = {AlistCF /. "A"[a_] :> "B"[a], AlistCF, {}};
][[1]];
Print["2. A/B conversion: ", Round[t2, 0.01], "s"];

(* Step 3: GroebnerBasis *)
t3 = AbsoluteTiming[
  agb = GroebnerBasis[aeqs, Join @@ avar, Modulus -> char];
][[1]];
Print["3. GroebnerBasis (MMA): ", Round[t3, 0.01], "s, size=", Length[agb]];

(* Step 4: prepare for Singular *)
t4 = AbsoluteTiming[
  agb2 = agb /. {A[i_] :> "A"[i], B[i_] :> "B"[i]};
  localA = Table["A"[i], {i, ne}];
  localB = Table["B"[i], {i, ne}];
  agbSymbols = Union@Cases[agb2, h_[i_Integer] /; StringQ[h], Infinity];
  agbA = SortBy[Select[agbSymbols, #[[0]] === "A" &], #[[1]] &];
  agbB = SortBy[Select[agbSymbols, #[[0]] === "B" &], #[[1]] &];
  localRep = Join[
    If[Length[agbB] > 0, Thread[agbB -> localB[[;; Length[agbB]]]], {}],
    If[Length[agbA] > 0, Thread[agbA -> localA[[;; Length[agbA]]]], {}]
  ];
  agbGlobal = agb2 /. localRep;
  avarGlobal = {localB[[;; Length[agbB]]], localA[[;; Length[agbA]]], {}};
][[1]];
Print["4. Variable prep: ", Round[t4, 0.01], "s"];

(* Step 5: SingularMinAssPrime *)
t5 = AbsoluteTiming[
  rawResult = LIECoreAlgebra`Private`SingularMinAssPrime[agbGlobal, Join @@ avarGlobal, Modulus -> char, "MonomialOrder" -> "lp"];
][[1]];
Print["5. SingularMinAssPrime: ", Round[t5, 0.01], "s"];

aprimelist = rawResult[[1]];
dimlist = rawResult[[2]];
zeroDimIdx = Select[Range@Length[aprimelist], dimlist[[#]] == 0 &];
Print["   Zero-dim primes: ", Length[zeroDimIdx], "/", Length[aprimelist]];

(* Step 6: GB refine *)
t6 = AbsoluteTiming[
  aprimelist = aprimelist[[#]] & /@ zeroDimIdx;
  aprimelist = GroebnerBasis[#, Join @@ avar, Modulus -> char] & /@ aprimelist;
][[1]];
Print["6. GB refine: ", Round[t6, 0.01], "s"];

(* Step 7: bivarPrimeInfo with Singular *)
bivarPrimeInfoSym = Symbol["LIECoreAlgebra`Private`bivarPrimeInfo"];
originalDownValues = DownValues[bivarPrimeInfoSym];
DownValues[bivarPrimeInfoSym] = {
  HoldPattern[bivarPrimeInfoSym[primelist_, Avar_, Bvar_, opts:OptionsPattern[]]] :>
    Module[{c = OptionValue[{opts}, Modulus], p = OptionValue[{opts}, "Parameters"],
      ls = OptionValue[{opts}, "LimitSector"], fr = OptionValue[{opts}, "FractionSymbol"]},
      SingularCoordinateRing`singularBivarPrimeInfo[primelist, Avar, Bvar,
        Modulus -> c, "Parameters" -> p, "LimitSector" -> ls, "FractionSymbol" -> fr, Verbose -> False]
    ]
};

t7 = AbsoluteTiming[
  resSingular = bivarPrimeInfoSym[aprimelist, avar[[2]], avar[[1]],
    "FractionSymbol" -> "B", Modulus -> char, "Parameters" -> avar[[3]], "LimitSector" -> sec];
][[1]];
Print["7. singularBivarPrimeInfo: ", Round[t7, 0.01], "s, regions=", Length[resSingular]];

DownValues[bivarPrimeInfoSym] = originalDownValues;

(* Step 8: bivarPrimeInfo with MMA (timeout 60s) *)
t8 = AbsoluteTiming[
  resMMA = TimeConstrained[
    bivarPrimeInfoSym[aprimelist, avar[[2]], avar[[1]],
      "FractionSymbol" -> "B", Modulus -> char, "Parameters" -> avar[[3]]],
    60,
    $TimedOut
  ];
][[1]];
If[resMMA === $TimedOut,
  Print["8. MMA bivarPrimeInfo: TIMEOUT after 60s"],
  Print["8. MMA bivarPrimeInfo: ", Round[t8, 0.01], "s, regions=", Length[resMMA]]
];

Print["\nTotal Singular path: ", Round[t1+t2+t3+t4+t5+t6+t7, 0.01], "s"];
