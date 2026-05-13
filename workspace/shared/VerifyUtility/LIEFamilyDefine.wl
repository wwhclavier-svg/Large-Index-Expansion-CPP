(* ::Package:: *)

(* Large Index Expansion - Family Definition Package *)
(* Defines IBP families and generates IBP relations *)


BeginPackage["LIEFamilyDefine`", {"LIEUtility`"}]

(* Public Interface *)
AsyIBPFamilyDefine::usage = "AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topsector] defines an IBP family.\n"<>
  "Returns {ne, nl, nE, nibp, Alist, vlist, ibpeqsfull, sectorlist}."

LargeIndexIBP::usage = "LargeIndexIBP[pdlist, loopmom, extmom, spsRep, opts] generates IBP relations in large index limit.\n"<>
  "Returns {ibpeqs, vlist}."

relationGeneral::usage = "relationGeneral[famname, pdlist, loopmom, extmom, spsRep] generates general IBP relations."

Begin["`Private`"]

(* ============================================================ *)
(* Load Dependencies *)
(* ============================================================ *)

If[!ValueQ[firstNonZero],
  Get["LIEUtility.wl"];
];

(* ============================================================ *)
(* Basic Definitions *)
(* ============================================================ *)



mon2F[poly_,zlist_,fsymb_]/;Head[poly]=!=List:=Total[(fsymb@@Exponent[#,zlist])*(#/.(#->1&/@zlist))&/@If[Head[poly]===Plus,List@@poly,List@poly]];
mon2F[vec_List,zlist_,fsymb_]:=mon2F[#,zlist,fsymb]&/@vec;
(* Derivative with respect to loop momenta *)
derivL[pdlist_,loopmom_,extmom_,spsRep_,sprule_]:=(Expand[Expand[Outer[Times,D[#,{loopmom}],Join[loopmom,extmom]]]/.spsRep/.sprule])&/@pdlist
(* Derivative with respect to external momenta *)
derivP[pdlist_,loopmom_,extmom_,spsRep_,sprule_]:=(Expand[Expand[Outer[Times,If[extmom=!={},D[#,{extmom}],{}],extmom]]/.spsRep/.sprule])&/@pdlist;
(*derivative of propagator denominators in terms of loop momenta or external momenta: q_j dD_i/dl_k, p_j dD_i/dp_k*)


(* Generate IBP relations *)
genIBP[derivl_]:=Module[{genRel},
genRel[j_Integer,k_Integer,ind_,derivmat_]:=Module[{eqmon,eqf,npd=Length@derivmat},
eqmon=If[j==k,"d",0]+Sum[-ind[[i]]*derivmat[[i,j,k]]/"z"[i],{i,npd}];
eqmon=Expand[eqmon/(Times@@Thread@Power["z"/@Range[npd],ind])];
eqf=mon2F[eqmon,"z"/@Range[npd],"f"]/.{"f"[vec__]:>"f"@@(-{vec})};
Return[eqf//Collect[#,_f]&];
];
genRel[#[[1]],#[[2]],Array["n"<>ToString[#]&,Length@derivl],derivl]&/@(Join@@Array[{#1,#2}&,Dimensions[derivl][[2;;3]]])
];

(* Generate Lorentz Invariance relations *)
genLI[derivp_,extmom_,spsRep_]:=Module[{genRel,ne=Length@extmom},
genRel[j_Integer,k_Integer,ind_,derivmat_]:=Module[{eqLI,eqLIf,npd=Length@derivmat},eqLI=Sum[-ind[[i1]]*Sum[(Times@@extmom[[{j,i2}]]/.spsRep)*derivmat[[i1,i2,k]]-(Times@@extmom[[{k,i2}]]/.spsRep)*derivmat[[i1,i2,j]],{i2,ne}]/"z"[i1],{i1,npd}];
eqLI=Expand[eqLI/(Times@@Thread@Power["z"/@Range[npd],ind])];
eqLIf=mon2F[eqLI,"z"/@Range[npd],"f"]/.{"f"[vec__]:>"f"@@(-{vec})};
Return[eqLIf];
];
genRel[#[[1]],#[[2]],Array["n"<>ToString[#]&,Length@derivp],derivp]&/@(Join@@Table[Table[{i,j},{j,i+1,ne}],{i,ne}])
];

(* Convert scalar products to propagator denominators *)
SP2PD[pdlist_,loopmom_,extmom_,spsRep_]:=Module[{scalarproducts,pd2sp,sp2pd,sp2pdrule,npd=Length@pdlist},
scalarproducts=Join[Union@Flatten@Outer[Times,loopmom,loopmom],Flatten@Outer[Times,loopmom,extmom]];
pd2sp=Expand[pdlist]/.Thread[Rule[scalarproducts,Array["sp",npd]]]/.spsRep;(*express pdlist in terms of sp[i]*)
sp2pd=LinearSolve[#[[2]],-#[[1]]+Array["z",npd]]&@CoefficientArrays[Expand[pdlist]/.Thread[Rule[scalarproducts,Array["sp",npd]]]/.spsRep,Array["sp",npd]];(*express scalarproduces in terms of z[i]*)
sp2pdrule=Thread[Rule[scalarproducts,sp2pd]];
Return[sp2pdrule];
];



(* ============================================================ *)
(* Core IBP Functions *)
(* ============================================================ *)

Options[relationGeneral] = {};
relationGeneral[famname_,pdlist_,loopmom_,extmom_,spsRep_,OptionsPattern[]]:=Module[
  {sprule,derivl,derivp,relationgeneral,relationspecific,npd=Length@pdlist},
  sprule=SP2PD[pdlist,loopmom,extmom,spsRep];
  derivl=derivL[pdlist,loopmom,extmom,spsRep,sprule];
  derivp=derivP[pdlist,loopmom,extmom,spsRep,sprule];
  relationgeneral=genIBP[derivl]~Join~genLI[derivp,extmom,spsRep];
  Return[relationgeneral/.{"f"[ind__]:>j@@Join[{famname},{ind}]}];
];


Options[LargeIndexIBP] = {"LimitSector"->{}};
LargeIndexIBP[pdlist_,loopmom_,extmom_,spsRep_,OptionsPattern[]]:=Module[
  {ne=Length@pdlist,ibpeqs,nlist,vlist,Alist,limitSector,limitRule},
  If[OptionValue["LimitSector"]=!={},limitSector=OptionValue["LimitSector"],limitSector=Table[1,ne]];
  nlist=Table["n"<>ToString[i],{i,ne}];
  vlist=Table["v"<>ToString[i],{i,ne}];
  limitRule=Table[nlist[[i]]->If[limitSector[[i]]>0,"n"+vlist[[i]],vlist[[i]]],{i,ne}];
  ibpeqs=relationGeneral["f",pdlist,loopmom,extmom,spsRep];
  ibpeqs=ibpeqs/.{j["f",a__]:>"g"@@({a}-nlist+vlist)}/.limitRule;
  Print["#IBP&LI: ",Length@ibpeqs];
  Return[{ibpeqs,vlist}];
];


Options[AsyIBPFamilyDefine] = {};
AsyIBPFamilyDefine[pdlist_,loopmom_,extmom_,spsRep_,topsector_]:=Module[
  {ne,nl,nE,nibp,Alist,ibpeqsfull,vlist,sectorlist},
  ne=Length@pdlist;nl=Length@loopmom;nE=Length@extmom;nibp=nl(nl+nE);
  Alist=Table["A"[i],{i,ne}];
  (* full IBP equations *)
  {ibpeqsfull,vlist}=LargeIndexIBP[pdlist,loopmom,extmom,spsRep]//{#[[1,1;;nibp]],#[[2]]}&;
  (* subsectors *)
  sectorlist=Normal@SparseArray[#->1&/@#,ne]&/@Subsets[Flatten@Position[topsector,1],{nl,ne}];
  Return[{ne,nl,nE,nibp,Alist,vlist,ibpeqsfull,sectorlist}];
];


End[] (* Private *)

EndPackage[]
