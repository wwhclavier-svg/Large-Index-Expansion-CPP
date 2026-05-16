AppendTo[$Path,"/home/ykm/blade-clean/BLAddOns/BladeIBP"];
AppendTo[$Path,"/usr/local/src/fflow/finiteflow-master/mathlink"];
AppendTo[$LibraryPath,"/usr/local/src/fflow/finiteflow-master"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
If[!DirectoryQ["ibps"], Put[{{}, {}}, "results/table"]; Abort[]];
<<BladeIBP.m;
<<FiniteFlow.m
<<topology.wl;
IntegralOrdering = 1;

Module[
{resdir,usints,exints,all,target,eqsfiles,vars,eqs,nids,position,ibplearn,mis,nonmis,indep,graph,in,ibps,depvarsorder,
time,pid,ps,eva,nonzero,sparseinfo, alleqs, fflearn, sol, rule,mem1,mem2},

resdir = "results";
If[!DirectoryQ@resdir,CreateDirectory[resdir]];
usints = Import["ibps/usints_def.m"];
exints = Import["ibps/exints_def.m"];
all = Join[First/@usints,GetAllInts[FileNames["ibps/ints_bub00_*.mx"]], First/@exints];
target = Get["target"];
If[target === 0,
target = Select[all, IntRank[#]<=Blade`SolveNoBL`IBPRank && IntDots[#]<=Blade`SolveNoBL`IBPDot &],
target = Select[EliminateZeroSectors[target], #=!=0&]];
If[target === {}, Print["target is empty."];Put[{{},{}}, "results/table"];Abort[]];
eqsfiles = FileNames["ibps/sids_*.json"];
vars = {eps};
indep = Get@FileNameJoin@{resdir, "indep"};
FFNewGraph[graph,in,vars];
If[!(vars==={}),
FFAlgJSONSparseSolver[graph,ibps,{in},"systemC.json"],
eqsfiles = FileNames["ibps/ids_*.mx"];
eqs=Join@@((Import[#]/.{s -> 3, msq -> 0, "d" -> 1/3})&/@eqsfiles);
eqs=Thread[eqs==0];
eqs=DeleteCases[eqs,True];
eqs=eqs[[indep]];
target=Get["target_rec"];
FFAlgSparseSolver[graph,ibps,{in},{},eqs,all,"NeededVars"->target]];
FFSolverOnlyHomogeneous[graph,ibps];
FFGraphOutput[graph,ibps];
ibplearn = FFSparseSolverLearn[graph,all];
{nonmis, mis} = {"DepVars", "IndepVars"}/.ibplearn;
Print["number of integrals" -> Length/@{nonmis,mis}];
sol = Switch[Length[vars],
0, FFReconstructNumeric[graph,vars,"MaxPrimes"->200,"MaxDegree"->1000,"PrintDebugInfo"->Automatic],
1, FFParallelReconstructUnivariate[graph,vars,"NThreads"->8,"MaxPrimes"->200,"MaxDegree"->1000,"PrintDebugInfo"->Automatic],
_, FFReconstructFunction[graph,vars,"NThreads"->8,"MaxPrimes"->200,"MaxDegree"->1000,"PrintDebugInfo"->Automatic]];
sol = Partition[sol,Length@mis];
sol = Join[sol, IdentityMatrix[Length@mis]];
rule = Thread[Join[nonmis,mis]->sol];
Put[{mis,rule},FileNameJoin@{resdir,"table"}];
FFDeleteGraph[graph];];
Quit[];