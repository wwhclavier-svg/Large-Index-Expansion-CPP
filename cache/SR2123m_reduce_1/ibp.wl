AppendTo[$Path,"/home/ykm/blade-workspace/BLAddOns/BladeIBP"];
AppendTo[$Path,"/usr/local/src/fflow/finiteflow-master/mathlink"];
AppendTo[$LibraryPath,"/usr/local/src/fflow/finiteflow-master"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
<<BladeIBP.m;
<<FiniteFlow.m
<<topology.wl;
IntegralOrdering = 1;
CreateDirectory["results"];
Put[ZeroSectors[SR2123m],"results/zerosectors"];
If[NonZeroSectors[SR2123m]==={}, Print["this family is trivial."];Quit[]];

Module[{rank,maxdots,GetSectors,GetSeeds,cutposi,learn,maplearn,filters,prefix,famcount,selectedfiles,lirankmin,rawmasters},

rank=3;
maxdots=1;
GetSectors=Automatic;
cutposi=Position[{0, 0, 0, 0, 0},1]//Flatten;
learn={3 -> {1, 3}, 2 -> {1, 3}, 1 -> {1, 3}, _ -> {1, 3}};
maplearn=Map[{1,0}+#&,learn//Association]//Normal;
filters[learn_]:=Join[If[learn=!={},{If[-Plus@@Select[#[[2]],#<0&]>(Length@Select[#[[2]],#1>0&]/.learn)[[2]],False,True]&,
If[(Plus@@Select[#[[2]],#>0&]-Length@Select[#[[2]],#>0&])>(Length@Select[#[[2]],#1>0&]/.learn)[[1]],False,True]&},{}],
If[False,{If[(#[[2]][[cutposi]]=!=ConstantArray[1,Length[cutposi]]),False,True]&},{}]];
lirankmin[_]:=0;
GetSeeds["IBP"][sector_]:=GenSeeds[sector,{0,maxdots},0,rank,filters[learn]];
GetSeeds["LI"][sector_]:=GenSeeds[sector,{0,maxdots},lirankmin[sector],rank,filters[learn]];
GetSeeds["SR"][sector_]:=GenSeeds[sector,{0,0},0,rank,filters[learn]];
GetSeeds["Map"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank,filters[maplearn]];
rawmasters=<|BL[SR2123m, {1, 0, 1, 0, 1}] -> {BL[SR2123m, {1, -1, 1, -1, 1}], BL[SR2123m, {1, 0, 1, 0, 1}]}, BL[SR2123m, {0, 0, 1, 0, 1}] -> {BL[SR2123m, {0, 0, 1, 0, 1}]}, BL[SR2123m, {1, 0, 0, 0, 1}] -> {BL[SR2123m, {1, 0, 0, 0, 1}]}, BL[SR2123m, {1, 0, 1, 0, 0}] -> {BL[SR2123m, {1, 0, 1, 0, 0}]}|>;
GetSeeds["SubSym"][sector_]:=Lookup[rawmasters,sector,{}];
GetSeeds["ExtMap"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank];
FastGenIds[SR2123m,GetSeeds,"Directory"->"ibps","LaunchKernels"->8,"GetSectors"->GetSectors];
];

Module[{exints,usints,usZeroIntegrals},

exints = {};
Export["ibps/exints_def.m",exints];
Export["ibps/ids_SR2123m_exints.mx",(-#[[1]]+#[[2]])&/@exints,"MX"];
Export["ibps/ints_SR2123m_exints.mx",IntegralsIn[exints],"MX"];
usints = {};
Export["ibps/usints_def.m",usints];
Export["ibps/ids_SR2123m_usints.mx",(-#[[1]]+#[[2]])&/@usints,"MX"];
Export["ibps/ints_SR2123m_usints.mx",IntegralsIn[usints],"MX"];


If[!({eps}==={}),SerializeFastIds[FileNames["ibps/ids_SR2123m_*.mx"], {s -> 3/2, msq -> 1}, {eps}, "UserDefinedInts"->(First/@usints),"ExtraInts"->(First/@exints)]];
];
CloseKernels[];Pause[0.1];
Quit[];