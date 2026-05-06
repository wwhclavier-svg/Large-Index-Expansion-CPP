(* ::Package:: *)

(* Large Index Expansion - Region Solving Package *)
(* Computes expansion regions for Feynman integrals *)


BeginPackage["LIERegions`", {"LIEUtility`", "LIECoreAlgebra`"}]

(* Public Interface *)
regionsBySectors::usage = "regionsBySectors[ibpeqslim, sectorlist, Alist, vlist, opts] computes expansion regions for given sectors.\n"<>
  "Returns an Association with sector keys, each containing a list of region data with CoordinateRing and RecursionMatrix."

Begin["`Private`"]

(* ============================================================ *)
(* Dependency Loading *)
(* ============================================================ *)

(* Dependencies are auto-loaded via BeginPackage *)

(* ============================================================ *)
(* Utility Functions *)
(* ============================================================ *)

lastNonZero[list_]:=Module[{i=Length[list]},While[list[[i]]===0&&i>=1,i--];If[i<1,Return[0],Return[i]]]

(* ============================================================ *)
(* Coefficient Matrix Functions *)
(* ============================================================ *)

IBPEigenEquationAB[ibpeqs_,ne_]:=ibpeqs/."g"[a__]:>(Product[If[#>=0,"A"[i]^#,"B"[i]^(-#)]&@({a}[[i]]-("v"<>ToString[i])),{i,ne}]);

Options[IBPCoefficientMatrix]={Modulus->0};
IBPCoefficientMatrix[ibpeqs_,ne_,OptionsPattern[]]:=Module[{char=OptionValue[Modulus],vlist,vSym,convRule,eigeneqs,nf=Length@ibpeqs,
coefruleN,N0coef,N0coefSym,N1coef,N2coef,F0Table,F1Table,F2Table,f0Table,f2Table,F2DTable,F2TTable,TakeValue,IBPFrules},
(* 1.preparations *)
vlist=Table["v"<>ToString[i],{i,ne}];
vSym=Table[Unique["v"],{i,ne}];
convRule=Thread[vlist->vSym];
TakeValue[dict_,key_]:=If[MemberQ[Keys[dict],key],dict[key],Table[0,nf]];
(* 2a.extract coefficients *)
eigeneqs=IBPEigenEquationAB[ibpeqs,ne];
coefruleN=monomialRulesPower[eigeneqs, Join[Array["B",ne],Array["A",ne]]];(* ~ N-part *)
(* 2b.data tables 1 *)
N0coef=(Table[0,2ne]//TakeValue[coefruleN,#]&);
N1coef=Transpose@Table[UnitVector[2ne,ne+i]//TakeValue[coefruleN,#]&,{i,ne}];
N2coef=Transpose[#,{2,3,1}]&@Table[(UnitVector[2ne,j]+UnitVector[2ne,i+ne])//TakeValue[coefruleN,#]&,{j,ne},{i,ne}];
(* 2c.data tables 2 - 使用符号变量计算系数 *)
N0coefSym=N0coef/.convRule;
F0Table=N0coef/.(#->0&/@vlist);
F2DTable=Coefficient[#,vSym]&/@N0coefSym;
f0Table=Total/@F2DTable;
F1Table=N1coef/.(#->1&/@vlist);
f2Table=N2coef/.(#->1&/@vlist);
F2Table=f2Table +(DiagonalMatrix/@F2DTable);
(* 3.collect data *)
IBPFrules=<|"F0"->F0Table,"F1"->F1Table,"F2"->F2Table,"f2"->f2Table,"F2D"->F2DTable,"f0"->f0Table|>;
If[char=!=0,IBPFrules=Association@Thread[Keys@IBPFrules->(Values@IBPFrules//PolynomialMod[#,char]&)]];
Return[IBPFrules];
];
(* Keys: {F0, F1, F2, f2=F2O, F2D, f0=F2T}, F2O:Off-diagonal, F2D:Diagonal, F2T:Trace *)
(* Old Notation substitution: Fd->F0, Fc->F2D, F2->f2, F0->f0 *)


Options[RecursionCoefficientMatrix]={Modulus->0,"LimitSector"->{}};
RecursionCoefficientMatrix[ibpcoefmat_,ABrule_,ne_,nf_,OptionsPattern[]]:=Module[{char=OptionValue[Modulus],secpos=OptionValue["LimitSector"],F0,F1,f0,f2,F2D,F2w,F1w,K1Table,K2Table,F2wsub,K1subTable,K2subTable,M1Table,M2Table,N0Table,NpTable,N1Table,NpmTable,RMrules},
(* 1.read data *)
If[secpos=={},secpos=Range[ne],secpos=Flatten@Position[secpos,1]];
{F0,F1,f2,F2D,f0}=ibpcoefmat/@{"F0","F1","f2","F2D","f0"};
(* 2a.eigenvalue-weighted ibp matrix *)
F2w=Table[f2[[m,i,j]]*("B"[j,i]/.ABrule),{m,nf},{i,ne},{j,ne}]//Expand;
F2wsub=Table[If[MemberQ[secpos,j],F2w[[m,i,j]],0],{m,nf},{i,ne},{j,ne}]//Expand;
F1w=Table[F1[[m,i]]*("A"[i]/.ABrule),{m,nf},{i,ne}]//Expand;
(* 2b.auxilliary matrix: K1,K2 *)
K1Table=Table[Table[Sum[F2w[[m,j,i]],{j,ne}]+F1w[[m,i]],{i,ne}],{m,nf}];
K1subTable=Table[If[MemberQ[secpos,i],K1Table[[m,i]],0],{m,nf},{i,ne}];
K2Table=Table[Table[Sum[F2w[[m,i,j]],{j,ne}],{i,ne}],{m,nf}];
K2subTable=Table[Table[Sum[If[MemberQ[secpos,j],F2w[[m,i,j]],0],{j,ne}],{i,ne}],{m,nf}];
(* 2c.matrix for specific order: M1,N1 *)
M1Table=Table[Table[K1subTable[[m,i]]-K2subTable[[m,i]],{i,ne}],{m,nf}];
N1Table=Table[Table[(F2D[[m,i]]+K1Table[[m,i]]),{i,ne}],{m,nf}];
(* 3.collect data *)
RMrules=<|"M1"->M1Table,"N1"->N1Table,"K1"->K1Table,"K2"->K2Table,"K1s"->K1subTable,"K2s"->K2subTable,"F2"->F2w,"F2s"->F2wsub,"F0"->F0|>;
If[char=!=0,RMrules=Association@Thread[Keys@RMrules->(Values@RMrules//PolynomialMod[#,char]&)]];
Return[RMrules];
]
(* Keys: {N1,F0,K1,F2,M1,K1s,K2s,F2s} *)


(* ============================================================ *)
(* Recursion Matrix Construction *)
(* ============================================================ *)

recursionMatrixCompanion[recursionMatrix_, expreg_, char_] := Module[{
    FFReduce, ConvertA2M, nb, M1, ne, nibp
  },
  (* 1.\:5c40\:90e8\:5efa\:7acb\:4e34\:65f6\:7684\:8f6c\:6362\:73af\:5883 *)
  FFReduce[x_] := If[char =!= 0, PolynomialMod[x, char], x];
  nb = Length[expreg["MonomialBasis"]];
  
  (* \:501f\:7528\:539f\:6709\:7684\:8f6c\:6362\:903b\:8f91\:ff0c\:4f46\:4ec5\:5728\:6b64\:8fd0\:884c\:4e00\:6b21 *)
  With[{basisMatrix = Association[Map[FFReduce, expreg["MonomialBasisMatrix"]]], 
        aindep = expreg["VarIndep"]},
    ConvertA2M[expr_] := If[expr === 0, 
      SparseArray[{}, {nb, nb}], 
      SparseArray@FFReduce[
        (Sum[m[Keys[#][[i]]] * Values[#][[i]], {i, Length[#]}] & @ 
         monomialRulesPower[expr, aindep]) /. m[a_] :> basisMatrix[a]
      ]
    ]
  ];

  (* 2.\:6267\:884c\:5168\:91cf\:8f6c\:6362\:5e76\:7a00\:758f\:5316 *)
  M1 = Map[ConvertA2M, recursionMatrix["M1"], {2}];
  {nibp, ne} = Dimensions[M1][[1;;2]];

  <|"M1" -> SparseArray@M1,
      "N1" -> SparseArray@Map[ConvertA2M, recursionMatrix["N1"], {2}],
      "K1" -> SparseArray@Map[ConvertA2M, recursionMatrix["K1"], {2}],
      "K1s" -> SparseArray@Map[ConvertA2M, recursionMatrix["K1s"], {2}],
      "K2s" -> SparseArray@Map[ConvertA2M, recursionMatrix["K2s"], {2}],
      "F2" -> SparseArray@Map[ConvertA2M, recursionMatrix["F2"], {3}],
      "F2s" -> SparseArray@Map[ConvertA2M, recursionMatrix["F2s"], {3}],
      "F0" -> SparseArray@Map[ConvertA2M, recursionMatrix["F0"], {1}],
      "nibp" -> nibp, "ne" -> ne, "nb" -> nb,  "incre" -> 2|>
]


recursionRankCheck[expReg_,recMat_,char_]:=Module[
{ne,nibp,aindep,basis,nb,basisIndex,basisMatrix,a2vrule,ConvertA2M,ConvertA2V,
	M1,N1, rankm,hmat,hmat2,hnull,rankh,Xsol},
aindep=expReg["VarIndep"];
basis=expReg["MonomialBasis"];
nb=Length@basis;
basisIndex=expReg["MonomialBasisIndex"];
basisMatrix=expReg["MonomialBasisMatrix"];
a2vrule=Table[m[basisIndex[[i]]]->UnitVector[nb,i],{i,nb}];
ConvertA2M[expr_]:=If[expr=!=0,(Sum[m[Keys[#][[i]]]*Values[#][[i]],{i,Length[#]}]&@monomialRulesPower[expr,aindep])/.m[a_]:>basisMatrix[a],ConstantArray[0,{nb,nb}]];
ConvertA2V[expr_]:=If[expr=!=0,(Sum[m[Keys[#][[i]]]*Values[#][[i]],{i,Length[#]}]&@monomialRulesPower[expr,aindep])/.a2vrule,ConstantArray[0,{nb}]];

(* 3.determine basis behavior from (M1,N1) *)
M1=Map[ConvertA2M,#,{2}]&@recMat["M1"]//Flatten[#,{{1,3},{2,4}}]&;
N1=Map[ConvertA2V,#,{2}]&@recMat["N1"]//Flatten[#,{{1,3},{2}}]&;
{nibp,ne}=Dimensions[recMat["M1"]][[1;;2]];
rankm=(MatrixRank[M1]/nb);
hmat=Table[PolynomialMod[#,char]&@Sum[ConvertA2M[recMat["M1"][[i,k]]] . ConvertA2V[recMat["N1"][[j,k]]],{k,ne}] . basis,{i,nibp},{j,nibp}];
hmat2=hmat//Map[ConvertA2M,#,{2}]&//Flatten[#,{{1,3},{2,4}}]&;
hnull=RowReduce@(NullSpace[hmat2] . M1)//If[char=!=0,RowReduce[#,Modulus->char],LatticeReduce[#]//#*LCM@@Flatten@Denominator[#]&]&//Map[ArrayReshape[#,{ne,nb}] . basis&,#]&;
rankh=MatrixRank[hmat//Map[ConvertA2M,#,{2}]&//Flatten[#,{{1,3},{2,4}}]&,Modulus->char]/nb;
(* solve linear system *)
Xsol=LinearSolve[M1,-N1,Modulus->char];
Xsol=Transpose[(ArrayReshape[#,{ne,nb}] . basis)&/@Transpose[Xsol]];
Print["M(+1):  ",recMat["M1"]//MatrixPlot,"    N(-1):  ",recMat["N1"]//MatrixPlot,"    H=M.N: ",hmat//MatrixPlot,"    sol. for M1.X=N1:  ",Xsol//MatrixPlot,"  (non-empty if compatiable)"];
Print["    corank=",ne-rankm,"    (rank(M),rank(H))=",{rankm,rankh}];
(* one-loop case (M1 is square matrix): if M1 is full rank, then N1 is always compatiable i.e. lying in the column space of M1 *)
{hmat,hnull}
];


Options[buildRecursionMatrix]={"LimitSector"->{},Modulus->0,"RankCheck"->False};
buildRecursionMatrix[ibpeqs_,expreg_,ne_,OptionsPattern[]]:=Module[{ConvertA2M,ConvertA2V,a2vrule,nb,basis,basisIndex,basisMatrix,aindep,M1,N1,Xsol,ibpMat,recMat,Alist,recurmat,recMatComp,nibp=Length@ibpeqs,sector=OptionValue["LimitSector"],hmat={},hmat2,hnull={},rankm,rankh,char=OptionValue[Modulus]},

(* 1.extract coefficient matrices *)
ibpMat=IBPCoefficientMatrix[ibpeqs,ne,Modulus->char];
(*Print["char. eqns:  ",IBPCharEqn[ibpcoefmat,Modulus->char]//Normal@CoefficientArrays[#,Alist][[2]]&//MatrixForm];*)
recMat=RecursionCoefficientMatrix[ibpMat,expreg//Join[#["VarRule"],#["FractionRule"]]&,ne,nibp,"LimitSector"->sector,Modulus->char];
recMatComp=recursionMatrixCompanion[recMat,expreg,char];
If[OptionValue["RankCheck"],{hmat,hnull}=recursionRankCheck[expreg,recMat,char]];

<|"CoordinateRing"->expreg,"RecursionMatrix"->recMatComp,"HMatrix"->hmat,"HNull"->hnull|>
]


(* ============================================================ *)
(* Main Public Function: regionsBySectors *)
(* ============================================================ *)

Options[regionsBySectors] = {
  "EnableFieldExtension" -> True, 
  Modulus -> 0,
  Verbose -> False,
  "Timeout" -> 120
}; 

regionsBySectors[ibpeqslim_, sectorlist_, Alist_, vlist_, OptionsPattern[]] := 
 Module[{i, sec, ne, ibpeqs = ibpeqslim /. "n" -> 0, ibpeqssub, 
   expregsub, expregdata, ibpmatsub, recmatsub, hmatsub, hnulsub, 
   char = OptionValue[Modulus], regAssoc},
  
  ne = Length @ Alist;
  
  expregdata = Association @ Table[
      sec = sectorlist[[i]];
      Print["\n      +++  sector:  ", sec, "  (", i, "/", Length @ sectorlist, ")  +++ "];

      ibpeqssub = sectorLimitIBP[ibpeqs, sec, vlist];

      (* \:83b7\:53d6\:4ee3\:6570\:533a\:57df\:4fe1\:606f — \:5e26\:8ba1\:65f6 *)
      {regTime, expregsub} = AbsoluteTiming @ TimeConstrained[
        expRegSolve2[ibpeqssub, Alist, vlist, "LimitSector" -> sec, Modulus -> char,Verbose -> OptionValue[Verbose]],
        OptionValue["Timeout"],
        $Failed
      ];
      If[expregsub === $Failed,
        Print["  *** TIMEOUT after ", OptionValue["Timeout"], "s, skipping sector ***"];
        Nothing
        ,
        (* \:8fc7\:6ee4\:5e73\:51e1\:533a\:57df *)
        expregsub = Select[expregsub, Join[#["VarDep"], #["VarIndep"]] =!= {} &];

        If[OptionValue["EnableFieldExtension"] == False,
          expregsub = Select[expregsub, #["VarIndep"] == {} &]
        ];

        If[expregsub =!= {},
          nRegs = Length[expregsub];
          Print["  >> ", nRegs, " region(s) in ", Round[regTime, 0.001], " s  (", If[nRegs>0, Round[regTime/nRegs, 0.001], 0], " s/region)"];
          (* \:83b7\:53d6\:77e9\:9635\:5316\:9012\:5f52\:5173\:7cfb *)
          sec -> Table[buildRecursionMatrix[ibpeqs, reg, ne, "LimitSector" -> sec, Modulus -> char], {reg, expregsub}]
          ,
          Print["  >> 0 regions (all trivial) in ", Round[regTime, 0.001], " s"];
          Nothing
        ]
      ]
    , {i, Length @ sectorlist}];

  Print["#Non-triv Sectors = ", Length @ Keys @ expregdata, "\n"];
  expregdata
];


End[] (* Private *)

EndPackage[]
