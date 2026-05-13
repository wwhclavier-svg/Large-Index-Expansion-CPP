(* test_gb_equivalence_SR212.wl — Verify SingularGroebnerBasisDefault equals MMA GroebnerBasis *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

(* SingularGroebnerBasisDefault is loaded into LIECoreAlgebra`Private` context *)
sgb = Symbol["LIECoreAlgebra`Private`SingularGroebnerBasisDefault"];

Print["Loading SR212 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
existingRegions = data["Regions"];

(* Pick first completed sector with 1 region *)
targetSector = SelectFirst[SortBy[Keys[existingRegions], Total],
  Length[existingRegions[#]] == 1 &,
  $Failed
];
If[targetSector === $Failed, Print["ERROR: No suitable sector"]; Exit[1]];

Print["Target sector: ", targetSector];

ne = Length@Alist;
AlistCF = Alist /. A[i_] :> "A"[i];
ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];
aeqs0 = Coefficient[ibpeqslimited, "n"] /. "g"[a__] :> (Times @@ Thread@Power[AlistCF, {a} - vlist]);
ibpeqssub = Join[aeqs0 /. {1/"A"[a_] :> "B"[a]}, Table["A"[i] "B"[i] - 1, {i, ne}]];
avar = {AlistCF /. "A"[a_] :> "B"[a], AlistCF, {}};

Print["Computing MMA GB..."];
timeMMA = AbsoluteTiming[gbMMA = GroebnerBasis[ibpeqssub, Join @@ avar, Modulus -> char]][[1]];
Print["MMA GB: ", Length[gbMMA], " polynomials, time=", timeMMA, "s"];

Print["Computing Singular GB..."];
timeS = AbsoluteTiming[gbS = sgb[ibpeqssub, Join @@ avar, Modulus -> char, "MonomialOrder" -> "lp"]][[1]];
Print["Singular GB: ", Length[gbS], " polynomials, time=", timeS, "s"];

(* Test equivalence via combined GB *)
Print["\nTesting ideal equality via combined GB..."];
gbCombined = GroebnerBasis[Join[gbMMA, gbS], Join @@ avar, Modulus -> char];
Print["Combined GB: ", Length[gbCombined], " polynomials"];
gbMMA_fromCombined = Select[gbCombined, PolynomialReduce[#, gbMMA, Join @@ avar, Modulus -> char][[2]] =!= 0 &];
gbS_fromCombined = Select[gbCombined, PolynomialReduce[#, gbS, Join @@ avar, Modulus -> char][[2]] =!= 0 &];
Print["Combined has new wrt MMA: ", Length[gbMMA_fromCombined]];
Print["Combined has new wrt Singular: ", Length[gbS_fromCombined]];

If[Length[gbMMA_fromCombined] == 0 && Length[gbS_fromCombined] == 0,
  Print["PASS: GBs are equivalent (same ideal)"];
,
  Print["FAIL: GBs differ"];
  Exit[1];
];
