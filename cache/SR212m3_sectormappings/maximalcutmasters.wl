AppendTo[$Path,"/home/ykm/\:4e0b\:8f7d/openclaw-docker/workspace/blade/BLAddOns/BladeIBP"];
AppendTo[$Path,"/usr/local/src/fflow/finiteflow-master/mathlink"];
AppendTo[$LibraryPath,"/usr/local/src/fflow/finiteflow-master"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
<<BladeIBP.m;
<<FiniteFlow.m;
<<topology.wl;


IntegralOrdering=3;
family=SR212m3;


FFLinearSolveLearn[eqs_,paras_,allvars_,neededvars_]:=Module[{graph,in,sys},
(*initialize fflow*)
FiniteFlow`Private`FFGraphId = Association[{}];
FiniteFlow`Private`FFAlgId = Association[{}];
FiniteFlow`Private`FFGraphInputs = Association[{}];
FFNThreads = Automatic;
(*solve*)
FFSparseSolve[eqs,allvars,"Parameters"->paras,"NeededVars"->neededvars,"IndepVarsOnly"->True]
];


SetOptions[AutoDetermine,"CutNoDot" -> False, "MaxIncrement" -> 2, "StartingPower" -> 2, "HighestPower" -> 5, "CheckSymmetry" -> True, "Nthreads" -> 8];


SetOptions[AutoDetermine,"GenCutIds" -> GenCutIds, "LinearReduce" -> FFLinearSolveLearn];


SetOptions[AutoDetermine,"Numeric"->{s -> 3/2, msq -> 1, "d" -> 1/7}];


Block[{mapdir,flag},
mapdir = "results";
If[!DirectoryQ@mapdir, CreateDirectory[mapdir]];
flag=AutoDetermine[family];
If[!FreeQ[flag,$Failed],Print["error: AutoDetermine failed"];Abort[]];
Put[MappedSectors[family],FileNameJoin[{mapdir,"mappedsectors"}]];
Put[UniqueSectors[family],FileNameJoin[{mapdir,"uniquesectors"}]];
Put[ZeroSectors[family],FileNameJoin[{mapdir,"zerosectors"}]];
Put[FastMIs[#]&/@NonZeroSectors[family]//Flatten//SortIntegrals,FileNameJoin[{mapdir,"rawmastersnomap"}]];
Put[FastMIs[#]&/@UniqueSectors[family]//Flatten//SortIntegrals,FileNameJoin[{mapdir,"rawmasters"}]];
Print["# rawmasters(nomap): ", Length@(Get@FileNameJoin[{mapdir,"rawmastersnomap"}]),", # rawmasters: ", Length@(Get@FileNameJoin[{mapdir,"rawmasters"}]) ];
];


CloseKernels[];Pause[0.1];


Quit[];