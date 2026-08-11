(* 通用 C++ 驱动 e2e: regionsBySectors -> cppExpandAndReconstruct -> cone -> reduce
   用法: wolframscript -script e2e_generic.wls <family> [order] [rMax] [sMax]
   输出: results/lie_e2e/<family>_rule.m + 打印关键结果 *)
fam = $ScriptCommandLine[[2]];
order = If[Length[$ScriptCommandLine] >= 3, ToExpression[$ScriptCommandLine[[3]]], 6];
rMax = If[Length[$ScriptCommandLine] >= 4, ToExpression[$ScriptCommandLine[[4]]], 2];
sMax = If[Length[$ScriptCommandLine] >= 5, ToExpression[$ScriptCommandLine[[5]]], 1];

Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SingularInterface.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/QuotientAlgebra.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LargeIndexExpansion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BlockRecursion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SymbolicBTForm.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LIE-CPP-Interface.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BaseFamilyDatabase.wl"];

{pdlist, loopmom, extmom, spsRep, topSector, numeric, char} =
  familySetting[$SpecifiedFamilyDatabase, fam];
numericD = {d -> (numeric["d"] /. d -> d)};
{ne, nl, nE, nibp, Alist, vlist, ibpeqs1, sectorlist} =
  AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topSector];
ibpeqsNum = PolynomialMod[ibpeqs1 /. numeric /. numericD, char];
Print["fam=", fam, " ne=", ne, " topSector=", topSector, " char=", char];

{regionData, sectorList} = regionsBySectors[ibpeqsNum, sectorlist, Alist, vlist, Modulus -> char];
Print["master sectors: ", Length[Keys[regionData]]];

rel = cppExpandAndReconstruct[fam, regionData, order, 3, 2, "Verbose" -> False];
Print["rel dims: ", Dimensions[rel]];

relForCone = Flatten /@ rel[[1 ;; Min[2, Length[rel]], 1 ;; Min[2, Length[rel[[1]]]]]];
blockDataList = setupGeneratingCone[relForCone, {topSector}, vlist, Modulus -> char, Verbose -> False];
Print["cone sectors: ", If[AssociationQ[blockDataList], Length[Keys[blockDataList]], "non-assoc"]];

targIntgList = seedsGenerate[topSector, rMax, sMax];
reductionRule = concreteFamilyReduce[blockDataList, targIntgList, vlist, "FullReduce" -> True];
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/lie_e2e/" <> fam <> "_rule.m", reductionRule];
Print["rules exported, #sectors: ", Length[reductionRule]];
