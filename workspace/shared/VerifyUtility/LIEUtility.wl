(* ::Package:: *)

(* Large Index Expansion - Utility Package *)
(* Shared utility functions for LIE packages *)


BeginPackage["LIEUtility`"]

(* Exported utility functions *)
firstNonZero::usage = "firstNonZero[list] returns the index of the first non-zero element."

monomialRulesPower::usage = "monomialRulesPower[expr, varlist] extracts monomial rules with power representation."

seedsGenPos::usage = "seedsGenPos[n, nprop] generates seed positions up to level n for nprop propagators."

sectorLimitIBP::usage = "sectorLimitIBP[ibpeqs, sector, vlist] applies sector limit to IBP equations."

reduceSolve::usage = "reduceSolve[eqs, vars, opts] solves linear equations using row reduction."
reduceSolveFP::usage = "reduceSolveFP[eqs, vars, opts] solves with floating point precision."

ExecuteLinearSolve::usage = "ExecuteLinearSolve[eqs, vars, char, solver] executes appropriate linear solver."

PolynomialDeg::usage = "PolynomialDeg[poly, var] returns the total degree of polynomial poly with respect to variables var."

Begin["`Private`"]

(* ============================================================ *)
(* First Non-Zero Index *)
(* ============================================================ *)

firstNonZero[list_]:=Module[{i=1},
  While[list[[i]]===0 && i <= Length@list, i++];
  If[i>Length@list, Print["zero list"]; Return[-1], Return[i]]
];

(* ============================================================ *)
(* Monomial Rules Extraction *)
(* ============================================================ *)

monomialRules[expr_,var_]:=If[var=!={},Association[Join@@
(If[Head[#]===SparseArray,Product[var[[i]]^Length@Cases[#[[1]],i,Infinity],{i,Length@var}]->#[[2]]&/@Most@ArrayRules[#],{1->#}]&/@CoefficientArrays[expr,var])],<|1->expr|>];

monomialRules[expr_List,var_]:=Module[{rulelist,monlist},
  rulelist=monomialRules[#,var]&/@expr;
  monlist=Union@@(Keys/@rulelist);
  Association@Table[mon->Table[If[MemberQ[Keys@rule,mon],rule[mon],0],{rule,rulelist}],{mon,monlist}]
];

monomialRulesPower[expr_,var_]:=If[var=!={},Association[Join@@
(If[Head[#]===SparseArray,Table[Length@Cases[#[[1]],i,Infinity],{i,Length@var}]->#[[2]]&/@Most@ArrayRules[#],{Table[0,Length@var]->#}]&/@CoefficientArrays[expr,var])],<|{}->expr|>];

monomialRulesPower[expr_List,var_]:=Module[{rulelist,monlist},
  rulelist=monomialRulesPower[#,var]&/@expr;
  monlist=Union@@(Keys/@rulelist);
  Association@Table[mon->Table[If[MemberQ[Keys@rule,mon],rule[mon],0],{rule,rulelist}],{mon,monlist}]
];

(* ============================================================ *)
(* Linear Solvers *)
(* ============================================================ *)

Options[reduceSolve]={Modulus->0};
reduceSolve[eqs0_,vars_,OptionsPattern[]]:=Module[
  {char=OptionValue[Modulus],mat,dep,indep,rule,varrev=Reverse@vars,eqs=DeleteCases[eqs0,True|False]},
  If[eqs==={},Return[{}]];
  mat=CoefficientArrays[eqs,varrev]//Join[#[[2]],List/@#[[1]],2]&;
  mat=If[char=!=0,Map[PolynomialMod[#,char]&,mat,{2}],mat];
  mat=RowReduce[mat,Modulus->char]//DeleteCases[#,a_/;Union[a]=={0}]&;
  dep=firstNonZero/@mat;
  If[Last@dep>Length@vars,Print["Incompatiable!"];Return[{}]];
  indep=Complement[Range[Length@vars],dep];
  rule=Thread[varrev[[dep]]->(-mat[[;;,indep]] . varrev[[indep]]-mat[[;;,-1]])];
  Return[rule];
];

Options[partialReduceSolve]={Modulus->0};
partialReduceSolve[eqs_,vars_,params_,OptionsPattern[]]:=Module[
  {char=OptionValue[Modulus],mat,matext,transf,dep,indep,rule,
   nv=Length@vars,np=Length@params,varrev=Reverse@vars,parrev=Reverse@params},
  matext=Join[CoefficientArrays[eqs,varrev][[2]],IdentityMatrix[Length@eqs],2];
  matext=If[char=!=0,Map[PolynomialMod[#,char]&,matext,{2}],matext];
  matext=RowReduce[matext,Modulus->char];
  transf=Simplify@matext[[;;,nv+1;;-1]];
  mat=CoefficientArrays[eqs,Join[varrev,parrev]]//Join[#[[2]],List/@#[[1]],2]&;
  mat=PolynomialMod[transf . mat,char]//DeleteCases[#,a_/;Union[a[[1;;nv]]]=={0}]&;
  dep=Select[firstNonZero/@mat,#<=nv&];
  indep=Complement[Range[nv+np],dep];
  rule=Thread[varrev[[dep]]->(-mat[[;;,indep]] . Join[varrev,parrev][[indep]]-mat[[;;,-1]])];
  Return[rule];
];

Options[reduceSolveFP]={Precision->30};
reduceSolveFP[eqs_,vars_,OptionsPattern[]]:=Module[
  {mat,dep,indep,rule,varrev=Reverse@vars,prec=OptionValue[Precision]},
  mat=CoefficientArrays[eqs,varrev]//Join[#[[2]],List/@#[[1]],2]&;
  mat=Chop[#,10^(-prec)]&@RowReduce[mat]//DeleteCases[#,a_/;Union[a]=={0}]&;
  dep=firstNonZero/@mat;
  If[Last@dep>Length@vars,Print["Incompatiable!"];Return[{}]];
  indep=Complement[Range[Length@vars],dep];
  rule=Thread[varrev[[dep]]->(-mat[[;;,indep]] . varrev[[indep]]-mat[[;;,-1]])];
  Return[rule];
];

ExecuteLinearSolve[eqs_, vars_, char_, solver_] := Which[
  char =!= 0, reduceSolve[# == 0 & /@ eqs, vars, Modulus -> char],
  solver === FFSparseSolve, FFSparseSolve[# == 0 & /@ eqs, Reverse@vars],
  True, solver[# == 0 & /@ eqs, vars]
];

(* ============================================================ *)
(* Seed Generation *)
(* ============================================================ *)

(* generating seeds at levels in List n_List or below level n_Integer of nprop_Integer propagators *)
seedsGenPos[n_Integer,nprop_]:=Join@@Table[
  DeleteDuplicates@(Join@@(Permutations/@(IntegerPartitions[i+nprop,{nprop}]-1))),
  {i,0,n}
];

seedsGenPos[n_List,nprop_]:=Join@@Table[
  DeleteDuplicates@(Join@@(Permutations/@(IntegerPartitions[n[[i]]+nprop,{nprop}]-1))),
  {i,Length@n}
];

(* ============================================================ *)
(* Sector Limit IBP *)
(* ============================================================ *)

coefficientReplace[expr_,varhead_,rule_]:=Module[
  {varlist,res},
  varlist=Union@Cases[expr,a_/;Head[a]===varhead,{0,Infinity}];
  res=monomialRules[expr,varlist]//Total[Thread@Times[Keys[#],(Values[#]/.rule)]]&
];

sectorLimitIBP[ibpeqs_,sector_,vlist_]:=coefficientReplace[#,"g",Thread[vlist->sector*"n"+vlist]]&/@ibpeqs;


(* ============================================================ *)
(* Polynomial Degree *)
(* ============================================================ *)

PolynomialDeg[poly_,var_]:=If[poly===0,-\[Infinity],Total[CoefficientRules[poly,var,DegreeReverseLexicographic][[1,1]]]];


End[] (* Private *)

EndPackage[]
