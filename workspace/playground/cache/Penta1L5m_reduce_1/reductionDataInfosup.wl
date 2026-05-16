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
all = Join[First/@usints,GetAllInts[FileNames["ibps/ints_Penta1L5m_*.mx"]], First/@exints];
target = Get["target"];
If[target === 0,
target = Select[all, IntRank[#]<=Blade`Trimming`IBPRank && IntDots[#]<=Blade`Trimming`IBPDot &],
target = Select[EliminateZeroSectors[target], #=!=0&]];
If[target === {}, Print["target is empty."];Put[{{},{}}, "results/table"];Abort[]];

eqsfiles = FileNames["ibps/sids_*.json"];
vars = {eps};

If[!FileExistsQ[FileNameJoin@{resdir, "indep"}],Print["file not found -> indep"];Abort[]];
indep = Get@FileNameJoin@{resdir, "indep"};
alleqs = Join@@(Import/@eqsfiles)[[All,2]];
alleqs = alleqs[[indep]];
Export["ibps/all.json",{Length[alleqs],alleqs},"RawJSON","Compact"->True];
FFSparseSystemToJSON["systemC.json",Length[alleqs],all,vars,{"ibps/all.json"},"NeededVars"->target];

FFDeleteGraph[graph];];
Quit[];