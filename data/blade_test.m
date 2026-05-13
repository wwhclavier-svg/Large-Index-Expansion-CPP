(* ::Package:: *)

$BladePath="/home/ykm/\:4e0b\:8f7d/openclaw-docker/workspace/blade";
Get[FileNameJoin[{$BladePath,"Blade.wl"}]];
Get[FileNameJoin[{NotebookDirectory[],"../verify/FamilyDatabase/FamilyDatabase.wl"}]];


famname="SR212-3m";
family = SR212m3;
dimension = 4-2*eps;
loop = $FamilyDatabase[famname]["LoopMomenta"];
leg = $FamilyDatabase[famname]["ExternalMomenta"];
conservation = {};
replacement = $FamilyDatabase[famname]["KinematicRules"](*{p1^2\[Rule]s,p2^2\[Rule]t,p3^2\[Rule]u, p1 p2\[Rule]-(s+t)/2, p2 p3\[Rule]-(t+u)/2, p1 p3\[Rule](s+2t+u)/2}*);
propagator = $FamilyDatabase[famname]["Propagators"];
topsector = $FamilyDatabase[famname]["TopSector"];
numeric =$FamilyDatabase[famname]["Numeric"]; 
modulus=$FamilyDatabase[famname]["Modulus"];
Print["numeric:  ",numeric,"  modulus:  ",modulus];
(*Save configuration*)
BLFamilyDefine[family,dimension, propagator,loop,leg,conservation, replacement,topsector,numeric]


Get[FileNameJoin[{NotebookDirectory[],"../AllRelations_"<>famname<>"_k5.m"}]];
Print[levdeg={#["Lev"],#["Deg"]}&/@$AllRelations];
Print[Length/@{#["Alphas"],#["Betas"]}&/@$AllRelations];
Print[matdim=Dimensions@#["Coefficients"]&/@$AllRelations];
Print[MatrixForm@SparseArray[Thread[Map[#+{0,1}&,levdeg]->matdim[[;;,2]]-1]//DeleteCases[#,a_/;a[[1,1]]==0]&]];


$AllRelations[[1]]["NE"]


Print[Table[$AllRelations[[i]]["Coefficients"]//Dimensions,{i,Length@$AllRelations}]];
vlist=Array["v"<>ToString[#]&,$AllRelations[[1]]["NE"]]
relcpp=Association@Table[
Block[{lev,deg,alphas,betas,rel},
{lev,deg}=$AllRelations[[i]]//{#["Lev"],#["Deg"]}&;
{alphas,betas}=$AllRelations[[i]]//{#["Alphas"],#["Betas"]}&;
rel=PolynomialMod[#,modulus]&@Sum[Sum[
$AllRelations[[i]]["Coefficients"][[(j-1)*Length@betas+k]]* 
Times@@Thread@Power[vlist,betas[[k]]],{k,Length@betas}]
"j"@@(vlist-alphas[[j]]),
{j,Length@alphas}];
{lev,deg}->Rest@rel
],{i,Length@$AllRelations}]





nuSample=RandomInteger[{1,0},5]


nuSample={1,0,1,0,1};


Print["j-target:  "];
jvar=Cases[relcpp,a_/;Head[a]==="j",Infinity]//DeleteDuplicates;
jTarget=jvar/.Thread[vlist->nuSample]


BLSetReducerOptions["MaxIncrement"->2,"IntegralOrdering"->3]


BLFamilyInf[SR212m3]


Print["BL-Reduction:  "];
res=BLReduce[jTarget/."j"[a__]:>BL[family,{a}]];


res


bladeTable=Thread[jTarget-> (res/.BL[family,a__]:>"j"@@a)];


Print["BL-Reduction Rule:  "];
bladeRule=(bladeTable/.Solve[dimension-("d"/.$FamilyDatabase[famname]["Numeric"])==0,eps][[1]])//Map[(#[[1]]->PolynomialMod[#[[2]],modulus])&,#]&


Print["BL-Verification:  "];
relveri=(Values@relcpp/.Thread[vlist->nuSample]/.bladeRule)//PolynomialMod[#,modulus]&;
relveri//MatrixForm


(* ::Subsubsection:: *)
(*Check With Series*)


Get[FileNameJoin[{NotebookDirectory[],"../ExpansionMMA_bub00.m"}]]


order=4;
hlist=$ExpansionResults[[1,1]]["Solutions"][[1]]["H"]//Sum[1/n^k*(#[[k+1]]/.{v1->"v1",v2->"v2"}),{k,0,Length[#]-1}]&;


secpos={1,1};vlist={"v1","v2"};Alist=Array[A,2];varRule={A[2]->59808223,A[1]->59808223};


seriesVerify[relations_,secpos_,degree_,order_]:=Module[{rel},
rel=relations/.{"j"[a__]:>jshift@@(vlist-{a})};
rel=(rel/n^degree)/.Thread[vlist->vlist+secpos*n];
rel=rel/.{jshift[a__]:>(PolynomialMod[#,modulus]&@Times@@Thread@Power[Alist/.varRule, -{a}]*(hlist/.Thread[vlist -> vlist-{a}]))};
rel=rel//CoefficientList[#,1/n,order+1]&//PolynomialMod[#,modulus]&]


seriesVerify[relmma,secpos,1,order]


(* ::Subsubsection:: *)
(*Compare With MMA*)


\!\(\*
TagBox[
RowBox[{
RowBox[{"relmma", "=", 
TagBox[
TagBox[
RowBox[{"(", "", 
TagBox[GridBox[{
{
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"59808225", "+", "v1"}], ")"}], " ", 
RowBox[{"g", "[", 
RowBox[{
RowBox[{
RowBox[{"-", "1"}], "+", "v1"}], ",", "v2"}], "]"}]}], "+", 
RowBox[{
RowBox[{"(", 
RowBox[{"59808226", "+", 
RowBox[{"179424670", " ", "v1"}], "+", 
RowBox[{"179424671", " ", "v2"}]}], ")"}], " ", 
RowBox[{"g", "[", 
RowBox[{"v1", ",", 
RowBox[{
RowBox[{"-", "1"}], "+", "v2"}]}], "]"}]}], "+", 
RowBox[{
RowBox[{"(", 
RowBox[{"3", "+", 
RowBox[{"3", " ", "v1"}], "+", 
RowBox[{"179424667", " ", "v2"}]}], ")"}], " ", 
RowBox[{"g", "[", 
RowBox[{"v1", ",", "v2"}], "]"}]}]}]},
{
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"179424672", "+", "v2"}], ")"}], " ", 
RowBox[{"g", "[", 
RowBox[{
RowBox[{
RowBox[{"-", "1"}], "+", "v1"}], ",", "v2"}], "]"}]}], "+", 
RowBox[{
RowBox[{"(", 
RowBox[{"59808223", "+", 
RowBox[{"2", " ", "v1"}], "+", "v2"}], ")"}], " ", 
RowBox[{"g", "[", 
RowBox[{"v1", ",", 
RowBox[{
RowBox[{"-", "1"}], "+", "v2"}]}], "]"}]}], "+", 
RowBox[{
RowBox[{"(", 
RowBox[{"179424670", "+", 
RowBox[{"3", " ", "v2"}]}], ")"}], " ", 
RowBox[{"g", "[", 
RowBox[{"v1", ",", "v2"}], "]"}]}]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
Column], "", ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]}], ";"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)
relmma=relmma/.{g->"j"}/.{v1->"v1",v2->"v2"}


relcpp[{1,1}]/.Thread[{"v1","v2"}->nuSample]//CoefficientArrays[#,jTarget][[2]]&//RowReduce[#,Modulus->modulus]&//MatrixForm


relmma/.Thread[{"v1","v2"}->nuSample]//CoefficientArrays[#,jTarget][[2]]&//RowReduce[#,Modulus->modulus]&//MatrixForm


(relmma/.Thread[{"v1","v2"}->{2,3}]/.bladeRule)//PolynomialMod[#,modulus]&
