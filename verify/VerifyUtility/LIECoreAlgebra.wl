(* ::Package:: *)

(* Large Index Expansion - Core Algebra Package *)
(* Contains expRegSolve2 and related algebraic geometry functions *)
(* Depends on Singular CAS for primary decomposition *)


BeginPackage["LIECoreAlgebra`", {"LIEUtility`"}]

(* Public Interface *)
expRegSolve2::usage = "expRegSolve2[ibpeqs, Alist, vlist, opts] solves for expansion regions using Groebner basis and primary decomposition.
"<>
  "Returns list of region associations with VarRule, VarIndep, VarDep, MinPoly, MonomialBasis, etc."

expRegSolve::usage = "expRegSolve[ibpeqs, Alist, vlist, opts] alternative region solver using product symbol."

bivarPrimeInfo::usage = "bivarPrimeInfo[primelist, Avar, Bvar, opts] extracts info from bivariate prime ideals."

primaryIdealInfo::usage = "primaryIdealInfo[pridlist, xvar, auxvar, opts] extracts info from primary ideals."

quotientRingBasisPower::usage = "quotientRingBasisPower[gblist, varlist] computes monomial basis for quotient ring."

quotientRingBasisMatrixPower::usage = "quotientRingBasisMatrixPower[gblist, varlist, monomialBasisPower] computes basis matrix."

Begin["`Private`"]

(* ============================================================ *)
(* Load Dependencies *)
(* ============================================================ *)

If[!ValueQ[firstNonZero],
  Get["LIEUtility.wl"];
];

(* Load SingularInterface for SingularMinAssPrime *)
If[!FileExistsQ["SingularInterface.wl"],
  Print["Warning: SingularInterface.wl not found. expRegSolve2 may not work."];
];
If[!ValueQ[SingularMinAssPrime],
  Get["SingularInterface.wl"];
];

(* ============================================================ *)
(* Quotient Ring Basis Functions *)
(* ============================================================ *)

(* Monomial Order Conversions from QuotientAlgebra.wl *)
MOShort={"Lexicographic"->"lp","NegativeLexicographic"->"ls","DegreeLexicographic"->"Dp",
"DegreeReverseLexicographic"->"dp","NegativeDegreeLexicographic"->"Ds","NegativeDegreeReverseLexicographic"->"ds"};
MOLong=Reverse/@MOShort;

(* Leading Term and Leading Monomial from QuotientAlgebra.wl *)
LT[exp_,var_,order_]:=MonomialList[exp,var,order/.MOLong][[1]];
LM[exp_,var_,order_]:=LT[exp,var,order]//If[Head[#]===Times,Select[#,Cases[#,a_/;MemberQ[var,a],{0,1}]=!={}&],If[Head[#]===Integer,1,#]]&;

Options[quotientRingBasisPower]={};
quotientRingBasisPower[gblist_,varlist_]:=Module[
  {n=Length@varlist,lm,basis,vardeg},
  lm=(LM[#,varlist,"lp"]&/@gblist);
  vardeg=Table[Max[Exponent[#,var]&/@Select[lm,ContainsAll[{var},Cases[#,a_/;MemberQ[varlist,a],Infinity]]&]],{var,varlist}];
  basis=If[vardeg=!={},vardeg//Flatten[Outer[List,Sequence@@Table[Range[0,#[[j]]-1],{j,Length[#]}]],Length[#]-1]&,{{}}];
  Return[basis]
];

Options[quotientRingBasisMatrixPower]={};
quotientRingBasisMatrixPower[gblist_,varlist_,monomialBasisPower_]:=Module[
  {res,mb=monomialBasisPower},
  res=Association@Table[
    mb[[j]]->Transpose@Table[
      Table[If[MemberQ[Keys[#],bas],#[bas],0],{bas,mb}]&@
      monomialRulesPower[#,varlist]&@
      PolynomialReduce[#,gblist,varlist,MonomialOrder->Lexicographic][[2]]&@
      (Times@@Thread@Power[varlist,mb[[i]]+mb[[j]]]),
      {i,Length@mb}
    ],
    {j,Length@mb}
  ];
  Return[res]
];

(* ============================================================ *)
(* Primary Ideal Info - Original Implementation *)
(* ============================================================ *)

Options[primaryIdealInfo]={Modulus->0,"FractionSymbol"->fr};
primaryIdealInfo[pridlist_,xvar_,auxvar_,OptionsPattern[]]:=Module[
  {char=OptionValue[Modulus],prid,ltlist,lmlist,varclass,vargen,varpar,ximinpoly,xdrule,frlist,frgb,frrule,FR=OptionValue["FractionSymbol"],basisIndex,basis,basisMatrix,res={}},
  Table[
    prid=pridlist[[ii]]//If[auxvar=!={},Select[#,Cases[#,Alternatives@@auxvar,Infinity]=={}&],#]&;
    Print["  ====== primary component (",ii,") ======"];
    (* determine dependent variables, independent variables, and minimal polynomials *)
    ltlist=First@MonomialList[#,xvar,Lexicographic]&/@prid;
    lmlist={Union@Cases[#,Alternatives@@xvar,{0,Infinity}],PolynomialDeg[#,xvar]}&/@ltlist;
    Print["leading degree: ",lmlist//MatrixForm];
    varclass=If[!MemberQ[Keys[#],True],Join[#,<|True->{}|>],#]&@GroupBy[#,#[[2]]>1&]&@lmlist;
    {vargen,varpar}={Sort[Join@@(First/@#[True])],Join@@(First/@#[False])}&@varclass;
    Print["generator: ",vargen,"  parametrized: ",varpar];
    ximinpoly=prid[[#]]&@Select[Range@Length@prid,lmlist[[#,2]]>1&];
    xdrule=Solve[#==0&/@Select[prid,Cases[#,Alternatives@@varpar,Infinity]=!={}&],varpar,Modulus->char][[1]];
    xdrule=xdrule//Map[#[[1]]->PolynomialReduce[#[[2]],ximinpoly,vargen,MonomialOrder->Lexicographic,Modulus->char][[2]]&,#]&;
    Print["solution:  ",xvar/.xdrule,"  with Min.Poly. ",ximinpoly];
    (* fractions *)
    frlist=Join@@Table[FR[i,j],{i,Length@xvar},{j,Length@xvar}]//DeleteCases[#,FR[a_,a_]]&;
    frgb=GroebnerBasis[Join[prid,Join@@Table[If[i=!=j,xvar[[j]]*FR[i,j]-xvar[[i]],Nothing],{i,Length@xvar},{j,Length@xvar}]],Join[auxvar,frlist,xvar],Modulus->char];
    frrule=Solve[#==0&/@Select[frgb,Cases[#,a_/;Head[a]===FR,Infinity]=!={}&],frlist,Modulus->char][[1]];
    frrule=frrule//Map[#[[1]]->PolynomialReduce[#[[2]],ximinpoly,vargen,MonomialOrder->Lexicographic,Modulus->char][[2]]&,#]&;
    (* coordinate ring as algebra *)
    basisIndex=quotientRingBasisPower[ximinpoly,vargen];
    basis=Times@@Thread[Power[vargen,#]]&/@basisIndex;
    basisMatrix=quotientRingBasisMatrixPower[ximinpoly,vargen,basisIndex];
    res=Join[res,{<|"VarDep"->varpar,"MinPoly"->ximinpoly,"VarIndep"->vargen,"VarRule"->xdrule,"FractionRule"->frrule,"MonomialBasis"->basis,"MonomialBasisIndex"->basisIndex,"MonomialBasisMatrix"->basisMatrix,"VarDeg"->lmlist[[;;,2]]|>}],
    {ii,Length@pridlist}
  ];
  Return[res]
];

(* ============================================================ *)
(* Bivariate Prime Info - Original Implementation *)
(* ============================================================ *)

Options[bivarPrimeInfo]={Modulus->0,"FractionSymbol"->fr,"Parameters"->{},Verbose->False};
bivarPrimeInfo[primelist_,Avar_,Bvar_,OptionsPattern[]]:=Module[
  {char=OptionValue[Modulus],prime,primeA,ltlist,lmlist,varclass,vargen,varpar,ximinpoly,xdrule,frlist,frgb,fractionrule,FR=OptionValue["FractionSymbol"],basisIndex,basis,basisMatrix,res={},auxvar=OptionValue["Parameters"],PrintF},
  PrintF[expr___]:=If[OptionValue[Verbose]===True,Print[expr]];
  Table[
    prime=primelist[[ii]];
    primeA=prime//Select[#,Cases[#,Alternatives@@Bvar,Infinity]=={}&]&;
    PrintF["  ====== primary component (",ii,") ======"];
    ltlist=First@MonomialList[#,Avar,Lexicographic]&/@primeA;
    lmlist={Union@Cases[#,Alternatives@@Avar,{0,Infinity}],PolynomialDeg[#,Avar]}&/@ltlist;
    PrintF["leading degree: ",lmlist//MatrixForm];
    varclass=If[!MemberQ[Keys[#],True],Join[#,<|True->{}|>],#]&@GroupBy[#,#[[2]]>1&]&@lmlist;
    {vargen,varpar}={Sort[Join@@(First/@#[True])],Join@@(First/@#[False])}&@varclass;
    PrintF["generator: ",vargen,"  parametrized: ",varpar];
    ximinpoly=prime[[#]]&@Select[Range@Length@primeA,lmlist[[#,2]]>1&];
    xdrule=Solve[#==0&/@Select[primeA,Cases[#,Alternatives@@varpar,Infinity]=!={}&],varpar,Modulus->char][[1]];
    xdrule=xdrule//Map[#[[1]]->PolynomialReduce[#[[2]],ximinpoly,vargen,MonomialOrder->Lexicographic,Modulus->char][[2]]&,#]&;
    PrintF["solution:  ",Avar/.xdrule,"  with Min.Poly. ",ximinpoly];
    (* fractions *)
    fractionrule=Join@@Table[If[i=!=j,FR[i,j]->PolynomialReduce[Avar[[i]]*Bvar[[j]],prime,Join[auxvar,Bvar,Avar],Modulus->char][[2]],Nothing],{i,Length@Avar},{j,Length@Avar}];
    (* coordinate ring as algebra *)
    basisIndex=quotientRingBasisPower[ximinpoly,vargen];
    basis=Times@@Thread[Power[vargen,#]]&/@basisIndex;
    basisMatrix=quotientRingBasisMatrixPower[ximinpoly,vargen,basisIndex];
    res=Join[res,{<|"VarDep"->Avar/.xdrule,"MinPoly"->ximinpoly,"VarIndep"->vargen,"VarRule"->xdrule,"FractionRule"->fractionrule,"MonomialBasis"->basis,"MonomialBasisIndex"->basisIndex,"MonomialBasisMatrix"->basisMatrix,"VarDeg"->lmlist[[;;,2]]|>}],
    {ii,Length@primelist}
  ];
  Return[res]
];

(* ============================================================ *)
(* expRegSolve2 - Main Region Solver *)
(* ============================================================ *)

Options[expRegSolve2]={Modulus->0,"LimitSector"->{},"Parameters"->{},"ZeroDim"->True,Verbose->False};
expRegSolve2[ibpeqs_,Alist_,vlist_,OptionsPattern[]]:=Module[
  {time,ne=Length@Alist,aeqs0,aeqs,agb,aprimelist,avar,reglist,reg,areg,aminpoly,apar,
   Blist,fagb,fracpar,ltlist,lmlist,vargen,varpar,varclass,char=OptionValue[Modulus],
   res={},basis,basisIndex,basisMatrix,limitSector,param=OptionValue["Parameters"],PrintF,
   localA,localB,localRep,agbGlobal,avarGlobal,agbSymbols,agbA,agbB,AlistCF},
  PrintF[expr___]:=If[OptionValue[Verbose]===True,Print[expr]];
  
  (* Convert Alist to context-free format "A"[i] *)
  AlistCF=Alist /. A[i_]:>"A"[i];
  
  If[OptionValue["LimitSector"]=!={},limitSector=OptionValue["LimitSector"],limitSector=Table[1,ne]];
  aeqs0=Coefficient[ibpeqs,"n"]/."g"[a__]:>(Times@@Thread@Power[AlistCF,{a}-vlist]);
  aeqs=Join[aeqs0/.{1/"A"[a_]:>"B"[a]},Table["A"[i]"B"[i]-1,{i,ne}]];
  avar={AlistCF/."A"[a_]:>"B"[a],AlistCF,param};

  time=AbsoluteTiming[agb=GroebnerBasis[aeqs,Join@@avar,Modulus->char]][[1]];
  PrintF["solving for GB in ",time," s"];
  
  (* Extract context-free symbols from agb for SingularInterface *)
  localA=Table["A"[i],{i,ne}];
  localB=Table["B"[i],{i,ne}];
  

  agb=agb/.{A[i_]:>"A"[i],B[i_]:>"B"[i]};
  

  agbSymbols=Union@Cases[agb,h_[i_Integer]/;StringQ[h],Infinity];
  agbA=SortBy[Select[agbSymbols,#[[0]]==="A"&],#[[1]]&];
  agbB=SortBy[Select[agbSymbols,#[[0]]==="B"&],#[[1]]&];
  

  localRep=Join[
    If[Length[agbB]>0,Thread[agbB->localB[[;;Length[agbB]]]],{}],
    If[Length[agbA]>0,Thread[agbA->localA[[;;Length[agbA]]]],{}],
    If[param=!={},Thread[param->("param"[#]&/@Range@Length@param)],{}]
  ];
  agbGlobal=agb/.localRep;
  avarGlobal={localB[[;;Length[agbB]]],localA[[;;Length[agbA]]],param/.localRep};
  time=AbsoluteTiming[aprimelist=SingularMinAssPrime[agbGlobal,Join@@avarGlobal,Modulus->char,"MonomialOrder"->"lp"]][[1]];

  PrintF["find minimal primes in ",time," s"];
  
  (* Fix SingularRun output: convert x1, x2, ... back to "A"[1], "B"[1], ... *)
  Module[{outRep, resultStr},
    (* Create string replacement rules: x1 -> "A"[1], etc. *)
    outRep=Join[
      Table["x"<>ToString[i]->"\"A\"["<>ToString[i]<>"]",{i,Length[localA]}],
      Table["x"<>ToString[Length[localA]+i]->"\"B\"["<>ToString[i]<>"]",{i,Length[localB]}]
    ];
    

    resultStr=ToString[aprimelist,InputForm];
    

    resultStr=StringReplace[resultStr,outRep];
    

    aprimelist=ToExpression[resultStr];
  ];
  
  (* Convert back to original context for downstream processing *)
  (* This ensures bivarPrimeInfo receives symbols in the expected context *)
  aprimelist=aprimelist/.(Reverse/@localRep);
  
  If[OptionValue["ZeroDim"]==True,
    PrintF["pos-dim:  ",aprimelist[[2]]];
    aprimelist=aprimelist[[1,#]]&/@Select[Range@Length[aprimelist[[1]]],aprimelist[[2,#]]==0&]
  ];
  aprimelist=GroebnerBasis[#,Join@@avar,Modulus->char]&/@aprimelist;
  res=bivarPrimeInfo[aprimelist,avar[[2]],avar[[1]],"FractionSymbol"->"B",Modulus->char,"Parameters"->avar[[3]],Verbose->OptionValue[Verbose]];
  res=Table[Append[res[[i]],"LimitSector"->limitSector],{i,Length@res}];
  PrintF["sol: ",Table[res[[i]]//Times@@#["VarDeg"]&,{i,Length@res}]//{#,Total[#]}&];
  Return[res];
];

(* ============================================================ *)
(* Alternative expRegSolve *)
(* ============================================================ *)

Options[expRegSolve]={"ProductSymbol"->False,Modulus->0,Verbose->False};
expRegSolve[ibpeqs_,Alist_,vlist_,OptionsPattern[]]:=Module[
  {ne=Length@Alist,char=OptionValue[Modulus],aeqs0,aeqs,agb,avar,reglist,reg,areg,aminpoly,apar,
   Blist,fagb,fracpar,ltlist,lmlist,vargen,varpar,varclass,res={},basis,basisIndex,basisMatrix,PrintF,
   localA,localf,localRep,aeqsGlobal,avarGlobal,aeqsSymbols,aeqsA,aeqsF,AlistCF},
  PrintF[a__]:=If[OptionValue[Verbose]===True,Print[a]];
  
  (* Convert Alist to context-free format "A"[i] *)
  AlistCF=Alist /. A[i_]:>"A"[i];
  
  aeqs0=Coefficient[ibpeqs,"n"]/."g"[a__]:>(Times@@Thread@Power[AlistCF,{a}-vlist]);
  avar=If[OptionValue["ProductSymbol"]==True,Join[AlistCF,{"f"[]}],Join[{"f"[]},AlistCF]];
  aeqs=Join[aeqs0//Together//Numerator,{Times@@AlistCF*"f"[]-1}];
  
  (* Extract context-free symbols from aeqs for SingularInterface *)
  localA=Table["A"[i],{i,ne}];
  localf="f"[];
  

  aeqs=aeqs/.A[i_]:>"A"[i];
  

  aeqsSymbols=Union@Cases[aeqs,h_[i_Integer]/;StringQ[h],Infinity];
  aeqsA=SortBy[Select[aeqsSymbols,#[[0]]==="A"&],#[[1]]&];
  aeqsF=Select[aeqsSymbols,#[[0]]==="f"&];
  

  localRep=Join[
    If[Length[aeqsA]>0,Thread[aeqsA->localA[[;;Length[aeqsA]]]],{}],
    If[Length[aeqsF]>0,Thread[aeqsF->{localf}],{}]
  ];
  aeqsGlobal=aeqs/.localRep;
  avarGlobal=avar/.localRep;
  reglist=SingularPrimaryDecompLex[aeqsGlobal,avarGlobal,Modulus->char];
  
  (* Fix SingularRun output: convert x1, x2, ... back to "A"[1], "A"[2], ... *)
  Module[{outRep, resultStr},

    outRep=Table["x"<>ToString[i]->"\"A\"["<>ToString[i]<>"]",{i,Length[localA]}];
    If[Length[aeqsF]>0,
      outRep=Append[outRep,"x"<>ToString[Length[localA]+1]->"\"f\"[]"];
    ];
    

    resultStr=ToString[reglist,InputForm];
    resultStr=StringReplace[resultStr,outRep];
    reglist=ToExpression[resultStr];
  ];
  
  (* Convert back to original context *)
  reglist=reglist/.(Reverse/@localRep);
  Print["Primary Decomposition:\n",Join[{Table["#"<>ToString[i]<>" (dim="<>ToString[#[[i]]]<>")",{i,Length[#]}]&@reglist[[2]]},Transpose[reglist[[1]]]]//MatrixForm];
  reglist={reglist[[1,#]],reglist[[2,#]]}&@Select[Range@Length@reglist[[2]],reglist[[2,#]]<1&];
  res=primaryIdealInfo[reglist[[1]],Alist,{"f"},"FractionSymbol"->"B",Modulus->char];
  Print["sol: ",Table[res[[i]]//Times@@#["VarDeg"]&,{i,Length@res}]//{#,Total[#]}&];
  Return[res];
];

End[] (* Private *)

EndPackage[]
