(* bub00 正确管线验证: regionsBySectors -> cppExpandAndReconstruct -> cone -> reduce *)
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SingularInterface.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/QuotientAlgebra.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LargeIndexExpansion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BlockRecursion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SymbolicBTForm.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LIE-CPP-Interface.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BaseFamilyDatabase.wl"];
fam = "bub00"; order = 6;
{pdlist, loopmom, extmom, spsRep, topSector, numeric, char} =
  familySetting[$SpecifiedFamilyDatabase, fam];
numericD = {d -> (numeric["d"] /. d -> d)};
{ne, nl, nE, nibp, Alist, vlist, ibpeqs1, sectorlist} =
  AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topSector];
ibpeqsNum = PolynomialMod[ibpeqs1 /. numeric /. numericD, char];
Print["ne=", ne, " topSector=", topSector, " char=", char];

(* 1. region 确定 (MMA) *)
{regionData, sectorList} = regionsBySectors[ibpeqsNum, sectorlist, Alist, vlist, Modulus -> char];
Print["master sectors: ", Length[Keys[regionData]]];

(* 2. C++ 展开 + 关系重构 *)
rel = cppExpandAndReconstruct[fam, regionData, order, 3, 2, "Verbose" -> False];
Print["rel dims: ", Dimensions[rel]];
Print["nonzero blocks: ", Table[{l, dd, Length[DeleteCases[rel[[l,dd]], 0] /. {} -> Nothing]}, {l, Length[rel]}, {dd, Length[rel[[1]]]}]];

(* 3. cone *)
relForCone = Flatten /@ rel[[1 ;; 2, 1 ;; 3]];
blockDataList = setupGeneratingCone[relForCone, {topSector}, vlist, Modulus -> char, Verbose -> False];
Print["cone sectors: ", If[AssociationQ[blockDataList], Length[Keys[blockDataList]], "non-assoc"]];

(* 4. reduce *)
targIntgList = seedsGenerate[topSector, 2, 1];
reductionRule = concreteFamilyReduce[blockDataList, targIntgList, vlist, "FullReduce" -> True];
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/intermediate/bub00_rule.m", reductionRule];
