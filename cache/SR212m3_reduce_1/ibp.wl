AppendTo[$Path,"/home/ykm/\:4e0b\:8f7d/openclaw-docker/workspace/blade/BLAddOns/BladeIBP"];
AppendTo[$Path,"/usr/local/src/fflow/finiteflow-master/mathlink"];
AppendTo[$LibraryPath,"/usr/local/src/fflow/finiteflow-master"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
<<BladeIBP.m;
<<FiniteFlow.m
<<topology.wl;
IntegralOrdering = 3;
CreateDirectory["results"];
Put[ZeroSectors[SR212m3],"results/zerosectors"];
If[NonZeroSectors[SR212m3]==={}, Print["this family is trivial."];Quit[]];

Module[{rank,maxdots,GetSectors,GetSeeds,cutposi,learn,maplearn,filters,prefix,famcount,selectedfiles,lirankmin,rawmasters},

rank=3;
maxdots=1;
GetSectors=Automatic;
cutposi=Position[{0, 0, 0, 0, 0},1]//Flatten;
learn={_ -> {1, 2}};
maplearn=Map[{1,0}+#&,learn//Association]//Normal;
filters[learn_]:=Join[If[learn=!={},{If[-Plus@@Select[#[[2]],#<0&]>(Length@Select[#[[2]],#1>0&]/.learn)[[2]],False,True]&,
If[(Plus@@Select[#[[2]],#>0&]-Length@Select[#[[2]],#>0&])>(Length@Select[#[[2]],#1>0&]/.learn)[[1]],False,True]&},{}],
If[False,{If[(#[[2]][[cutposi]]=!=ConstantArray[1,Length[cutposi]]),False,True]&},{}]];
lirankmin[_]:=0;
GetSeeds["IBP"][sector_]:=GenSeeds[sector,{0,maxdots},0,rank,filters[learn]];
GetSeeds["LI"][sector_]:=GenSeeds[sector,{0,maxdots},lirankmin[sector],rank,filters[learn]];
GetSeeds["SR"][sector_]:=GenSeeds[sector,{0,0},0,rank,filters[learn]];
GetSeeds["Map"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank,filters[maplearn]];
rawmasters=<|BL[SR212m3, {1, 0, 1, 0, 1}] -> {BL[SR212m3, {2, 0, 1, 0, 1}], BL[SR212m3, {1, 0, 1, 0, 1}]}, BL[SR212m3, {0, 0, 1, 0, 1}] -> {BL[SR212m3, {0, 0, 1, 0, 1}]}, BL[SR212m3, {1, 0, 0, 0, 1}] -> {BL[SR212m3, {1, 0, 0, 0, 1}]}, BL[SR212m3, {1, 0, 1, 0, 0}] -> {BL[SR212m3, {1, 0, 1, 0, 0}]}|>;
GetSeeds["SubSym"][sector_]:=Lookup[rawmasters,sector,{}];
GetSeeds["ExtMap"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank];
FastGenIds[SR212m3,GetSeeds,"Directory"->"ibps","LaunchKernels"->8,"GetSectors"->GetSectors];
];

Module[{exints,usints,usZeroIntegrals},

exints = {};
Export["ibps/exints_def.m",exints];
Export["ibps/ids_SR212m3_exints.mx",(-#[[1]]+#[[2]])&/@exints,"MX"];
Export["ibps/ints_SR212m3_exints.mx",IntegralsIn[exints],"MX"];
usints = {};
Export["ibps/usints_def.m",usints];
Export["ibps/ids_SR212m3_usints.mx",(-#[[1]]+#[[2]])&/@usints,"MX"];
Export["ibps/ints_SR212m3_usints.mx",IntegralsIn[usints],"MX"];


If[!({eps}==={}),SerializeFastIds[FileNames["ibps/ids_SR212m3_*.mx"], {s -> 3/2, msq -> 1, "d" -> 1/7}, {eps}, "UserDefinedInts"->(First/@usints),"ExtraInts"->(First/@exints)]];
];
CloseKernels[];Pause[0.1];
Quit[];