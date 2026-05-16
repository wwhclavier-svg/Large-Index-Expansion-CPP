AppendTo[$Path,"/home/ykm/blade-clean/BLAddOns/BladeIBP"];
AppendTo[$Path,"/usr/local/src/fflow/finiteflow-master/mathlink"];
AppendTo[$LibraryPath,"/usr/local/src/fflow/finiteflow-master"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
<<BladeIBP.m;
<<FiniteFlow.m
<<topology.wl;
IntegralOrdering = 1;
CreateDirectory["results"];
Put[ZeroSectors[bub00],"results/zerosectors"];
If[NonZeroSectors[bub00]==={}, Print["this family is trivial."];Quit[]];

Module[{rank,maxdots,GetSectors,GetSeeds,cutposi,learn,maplearn,filters,prefix,famcount,selectedfiles,lirankmin,rawmasters},

rank=2;
maxdots=3;
GetSectors=Automatic;
cutposi=Position[{0, 0},1]//Flatten;
learn={2 -> {3, 2}, 1 -> {3, 2}, 0 -> {3, 2}, _ -> {3, 2}};
maplearn=Map[{1,0}+#&,learn//Association]//Normal;
filters[learn_]:=Join[If[learn=!={},{If[-Plus@@Select[#[[2]],#<0&]>(Length@Select[#[[2]],#1>0&]/.learn)[[2]],False,True]&,
If[(Plus@@Select[#[[2]],#>0&]-Length@Select[#[[2]],#>0&])>(Length@Select[#[[2]],#1>0&]/.learn)[[1]],False,True]&},{}],
If[False,{If[(#[[2]][[cutposi]]=!=ConstantArray[1,Length[cutposi]]),False,True]&},{}]];
lirankmin[_]:=0;
GetSeeds["IBP"][sector_]:=GenSeeds[sector,{0,maxdots},0,rank,filters[learn]];
GetSeeds["LI"][sector_]:=GenSeeds[sector,{0,maxdots},lirankmin[sector],rank,filters[learn]];
GetSeeds["SR"][sector_]:=GenSeeds[sector,{0,0},0,rank,filters[learn]];
GetSeeds["Map"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank,filters[maplearn]];
rawmasters=<|BL[bub00, {1, 1}] -> {BL[bub00, {1, 1}]}|>;
GetSeeds["SubSym"][sector_]:=Lookup[rawmasters,sector,{}];
GetSeeds["ExtMap"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank];
FastGenIds[bub00,GetSeeds,"Directory"->"ibps","LaunchKernels"->8,"GetSectors"->GetSectors];
];

Module[{exints,usints,usZeroIntegrals},

exints = {};
Export["ibps/exints_def.m",exints];
Export["ibps/ids_bub00_exints.mx",(-#[[1]]+#[[2]])&/@exints,"MX"];
Export["ibps/ints_bub00_exints.mx",IntegralsIn[exints],"MX"];
usints = {};
Export["ibps/usints_def.m",usints];
Export["ibps/ids_bub00_usints.mx",(-#[[1]]+#[[2]])&/@usints,"MX"];
Export["ibps/ints_bub00_usints.mx",IntegralsIn[usints],"MX"];


If[!({eps}==={}),SerializeFastIds[FileNames["ibps/ids_bub00_*.mx"], {s -> 3, msq -> 0, "d" -> 1/3}, {eps}, "UserDefinedInts"->(First/@usints),"ExtraInts"->(First/@exints)]];
];
CloseKernels[];Pause[0.1];
Quit[];