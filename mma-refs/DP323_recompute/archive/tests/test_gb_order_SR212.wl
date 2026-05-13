(* test_gb_order_SR212.wl — Investigate monomial order differences *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

sgb = Symbol["LIECoreAlgebra`Private`SingularGroebnerBasisDefault"];

Print["Loading SR212 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
existingRegions = data["Regions"];

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
vars = Join @@ avar;

(* Test different monomial orders in MMA *)
orders = {Lexicographic, DegreeLexicographic, DegreeReverseLexicographic};
Do[
  Print["\n=== MMA ", orders[[i]], " ==="];
  gb = GroebnerBasis[ibpeqssub, vars, Modulus -> char, MonomialOrder -> orders[[i]]];
  Print["Length: ", Length[gb]];
  Print["LTs: ", MonomialList[#, vars, orders[[i]]][[1]] & /@ gb];
, {i, Length[orders]}];

(* Test Singular with lp *)
Print["\n=== Singular lp ==="];
gbS = sgb[ibpeqssub, vars, Modulus -> char, "MonomialOrder" -> "lp"];
Print["Length: ", Length[gbS]];
Print["LTs (MMA Lex): ", MonomialList[#, vars, Lexicographic][[1]] & /@ gbS];

(* Compare ideal equality for each *)
Print["\n=== Ideal equality tests ==="];
gbLex = GroebnerBasis[ibpeqssub, vars, Modulus -> char, MonomialOrder -> Lexicographic];
gbDp = GroebnerBasis[ibpeqssub, vars, Modulus -> char, MonomialOrder -> DegreeReverseLexicographic];

testEq[gb1_, gb2_, name_] := Module[{combined, new1, new2},
  combined = GroebnerBasis[Join[gb1, gb2], vars, Modulus -> char, MonomialOrder -> Lexicographic];
  new1 = Select[combined, PolynomialReduce[#, gb1, vars, Modulus -> char, MonomialOrder -> Lexicographic][[2]] =!= 0 &];
  new2 = Select[combined, PolynomialReduce[#, gb2, vars, Modulus -> char, MonomialOrder -> Lexicographic][[2]] =!= 0 &];
  Print[name, ": new wrt 1=", Length[new1], ", new wrt 2=", Length[new2], ", equal=", Length[new1] == 0 && Length[new2] == 0];
];

testEq[gbLex, gbS, "MMA-Lex vs Singular-lp"];
testEq[gbLex, gbDp, "MMA-Lex vs MMA-Dp"];
testEq[gbS, gbDp, "Singular-lp vs MMA-Dp"];
