(* test_singular_SR212.wl — Validate SingularCoordinateRing on SR212 *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

bivarPrimeInfoMMA = Symbol["LIECoreAlgebra`Private`bivarPrimeInfo"];
singularMinAssPrime = Symbol["LIECoreAlgebra`Private`SingularMinAssPrime"];

Print["Loading SR212 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
fullSectorList = data["Family", "SectorList"];
existingRegions = data["Regions"];

Print["Total sectors: ", Length[fullSectorList]];
Print["Completed sectors: ", Length[Keys[existingRegions]]];

(* List completed sectors by prop count *)
Do[
  secs = Select[Keys[existingRegions], Total[#] == n &];
  If[Length[secs] > 0,
    regs = Total[Length[existingRegions[#]] & /@ secs];
    Print["  ", n, " props: ", Length[secs], " sectors, ", regs, " regions"];
  ];
, {n, 1, Max[Total /@ fullSectorList]}];

(* Pick first completed sector with fewest props and 1 region *)
targetSector = SelectFirst[SortBy[Keys[existingRegions], Total],
  Length[existingRegions[#]] == 1 &,
  $Failed
];
If[targetSector === $Failed,
  Print["ERROR: No suitable 1-region sector found"];
  Exit[1];
];

Print["\nTarget sector: ", targetSector, " (", Total[targetSector], " props, ", Length[existingRegions[targetSector]], " region)"];

(* Replicate expRegSolve2 core logic for this sector *)
ne = Length@Alist;
AlistCF = Alist /. A[i_] :> "A"[i];

(* sectorLimitIBP is applied to ORIGINAL ibpeqs, then expRegSolve2 does A/B conversion *)
ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

aeqs0 = Coefficient[ibpeqslimited, "n"] /. "g"[a__] :> (Times @@ Thread@Power[AlistCF, {a} - vlist]);
ibpeqssub = Join[aeqs0 /. {1/"A"[a_] :> "B"[a]}, Table["A"[i] "B"[i] - 1, {i, ne}]];
avar = {AlistCF /. "A"[a_] :> "B"[a], AlistCF, {}};

Print["Sector-limited equations: ", Length[ibpeqssub]];

(* Compute GB and primes exactly as expRegSolve2 does *)
agb = GroebnerBasis[ibpeqssub, Join @@ avar, Modulus -> char];
Print["Global GB size: ", Length[agb]];

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

Print["Variables for Singular: ", avarGlobal];
Print["Running SingularMinAssPrime..."];
rawResult = singularMinAssPrime[agbGlobal, Join @@ avarGlobal, Modulus -> char, "MonomialOrder" -> "lp"];
Print["Raw result head: ", Head[rawResult], ", length: ", Length[rawResult]];
Print["Raw result: ", Short[rawResult, 100]];

aprimelist = rawResult[[1]];
dimlist = rawResult[[2]];
Print["Primes: ", Length[aprimelist], ", dims: ", dimlist];

zeroDimIdx = Select[Range@Length[aprimelist], dimlist[[#]] == 0 &];
Print["Zero-dim primes: ", Length[zeroDimIdx]];

If[Length[zeroDimIdx] == 0,
  Print["ERROR: No zero-dimensional primes found"];
  Exit[1];
];

aprimelist = aprimelist[[#]] & /@ zeroDimIdx;
aprimelist = GroebnerBasis[#, Join @@ avar, Modulus -> char] & /@ aprimelist;

Print["\n=== MMA bivarPrimeInfo ==="];
timeMMA = AbsoluteTiming[
  resMMA = bivarPrimeInfoMMA[aprimelist, avar[[2]], avar[[1]],
    "FractionSymbol" -> "B", Modulus -> char, "Parameters" -> avar[[3]], Verbose -> True];
][[1]];
Print["MMA result: ", Length[resMMA], " region(s), time=", timeMMA, "s"];

Print["\n=== Singular bivarPrimeInfo ==="];
timeSingular = AbsoluteTiming[
  resSingular = SingularCoordinateRing`singularBivarPrimeInfo[aprimelist, avar[[2]], avar[[1]],
    "FractionSymbol" -> "B", Modulus -> char, "Parameters" -> avar[[3]], "LimitSector" -> targetSector, Verbose -> True];
][[1]];
Print["Singular result: ", Length[resSingular], " region(s), time=", timeSingular, "s"];

(* Compare *)
Print["\n========================================"];
Print["Comparison"];
Print["========================================"];

If[Length[resMMA] =!= Length[resSingular],
  Print["FAIL: Different region counts: MMA=", Length[resMMA], ", Singular=", Length[resSingular]];
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
      Print["    MMA:     ", Short[valMMA, 50]];
      Print["    Singular:", Short[valS, 50]];
      allPass = False;
      ,
      Print["  PASS: ", key];
    ];
  , {k, Length[keysToCheck]}];
, {i, Length[resMMA]}];

Print["\n========================================"];
If[allPass, Print["ALL CHECKS PASSED"], Print["SOME CHECKS FAILED"]];
Print["========================================"];
