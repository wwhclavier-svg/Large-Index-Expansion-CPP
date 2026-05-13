(* test_agb_compare_SR212.wl — Compare agbGlobal between MMA and Singular paths *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

sgb = Symbol["LIECoreAlgebra`Private`SingularGroebnerBasisDefault"];
sma = Symbol["LIECoreAlgebra`Private`SingularMinAssPrime"];

Print["Loading SR212 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];

targetSector = {0, 1, 0, 1, 1};

ne = Length@Alist;
AlistCF = Alist /. A[i_] :> "A"[i];
ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];
aeqs0 = Coefficient[ibpeqslimited, "n"] /. "g"[a__] :> (Times @@ Thread@Power[AlistCF, {a} - vlist]);
ibpeqssub = Join[aeqs0 /. {1/"A"[a_] :> "B"[a]}, Table["A"[i] "B"[i] - 1, {i, ne}]];
avar = {AlistCF /. "A"[a_] :> "B"[a], AlistCF, {}};

(* --- MMA path --- *)
agbMMA = GroebnerBasis[ibpeqssub, Join @@ avar, Modulus -> char];
agbMMA = agbMMA /. {A[i_] :> "A"[i], B[i_] :> "B"[i]};
agbSymbolsMMA = Union@Cases[agbMMA, h_[i_Integer] /; StringQ[h], Infinity];
agbAMMA = SortBy[Select[agbSymbolsMMA, #[[0]] === "A" &], #[[1]] &];
agbBMMA = SortBy[Select[agbSymbolsMMA, #[[0]] === "B" &], #[[1]] &];
localAMMA = Table["A"[i], {i, ne}];
localBMMA = Table["B"[i], {i, ne}];
localRepMMA = Join[
  If[Length[agbBMMA] > 0, Thread[agbBMMA -> localBMMA[[;; Length[agbBMMA]]]], {}],
  If[Length[agbAMMA] > 0, Thread[agbAMMA -> localAMMA[[;; Length[agbAMMA]]]], {}]
];
agbGlobalMMA = agbMMA /. localRepMMA;
avarGlobalMMA = {localBMMA[[;; Length[agbBMMA]]], localAMMA[[;; Length[agbAMMA]]], {}};

(* --- Singular path --- *)
agbS = sgb[ibpeqssub, Join @@ avar, Modulus -> char, "MonomialOrder" -> "lp"];
agbS = agbS /. {A[i_] :> "A"[i], B[i_] :> "B"[i]};
agbSymbolsS = Union@Cases[agbS, h_[i_Integer] /; StringQ[h], Infinity];
agbAS = SortBy[Select[agbSymbolsS, #[[0]] === "A" &], #[[1]] &];
agbBS = SortBy[Select[agbSymbolsS, #[[0]] === "B" &], #[[1]] &];
localAS = Table["A"[i], {i, ne}];
localBS = Table["B"[i], {i, ne}];
localRepS = Join[
  If[Length[agbBS] > 0, Thread[agbBS -> localBS[[;; Length[agbBS]]]], {}],
  If[Length[agbAS] > 0, Thread[agbAS -> localAS[[;; Length[agbAS]]]], {}]
];
agbGlobalS = agbS /. localRepS;
avarGlobalS = {localBS[[;; Length[agbBS]]], localAS[[;; Length[agbAS]]], {}};

Print["MMA agbGlobal length: ", Length[agbGlobalMMA]];
Print["Singular agbGlobal length: ", Length[agbGlobalS]];

(* Check if ideals are equal *)
combined = GroebnerBasis[Join[agbGlobalMMA, agbGlobalS], Join @@ avarGlobalMMA, Modulus -> char, MonomialOrder -> Lexicographic];
newMMA = Select[combined, PolynomialReduce[#, agbGlobalMMA, Join @@ avarGlobalMMA, Modulus -> char, MonomialOrder -> Lexicographic][[2]] =!= 0 &];
newS = Select[combined, PolynomialReduce[#, agbGlobalS, Join @@ avarGlobalS, Modulus -> char, MonomialOrder -> Lexicographic][[2]] =!= 0 &];
Print["Combined GB new wrt MMA: ", Length[newMMA]];
Print["Combined GB new wrt Singular: ", Length[newS]];

(* Compare minAssGTZ results *)
Print["\nMMA minAssGTZ..."];
minMMA = sma[agbGlobalMMA, Join @@ avarGlobalMMA, Modulus -> char, "MonomialOrder" -> "lp"];
Print["Type: ", Head[minMMA], ", Length: ", Length[minMMA]];
If[ListQ[minMMA] && Length[minMMA] >= 2,
  Print["Num primes: ", Length[minMMA[[1]]]];
  Print["Dimensions: ", minMMA[[2]]];
  Print["Zero-dim count: ", Count[minMMA[[2]], 0]];
];

Print["\nSingular minAssGTZ..."];
minS = sma[agbGlobalS, Join @@ avarGlobalS, Modulus -> char, "MonomialOrder" -> "lp"];
Print["Type: ", Head[minS], ", Length: ", Length[minS]];
If[ListQ[minS] && Length[minS] >= 2,
  Print["Num primes: ", Length[minS[[1]]]];
  Print["Dimensions: ", minS[[2]]];
  Print["Zero-dim count: ", Count[minS[[2]], 0]];
];

(* Now run FULL expRegSolve2 with MMA GB to see what it returns *)
Print["\n--- Full MMA path ---"];
resMMA = LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist, "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];
(* Wait, expRegSolve2 now uses Singular GB! We need to test with backup. *)
Print["Skipping - expRegSolve2 already modified"];
