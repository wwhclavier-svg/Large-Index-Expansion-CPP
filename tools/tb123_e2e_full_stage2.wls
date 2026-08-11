(* TB123 完整 e2e 第二段: ImportRegions -> expand -> relations -> cone -> reduce *)
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SingularInterface.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/QuotientAlgebra.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LargeIndexExpansion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BlockRecursion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SymbolicBTForm.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LIE-CPP-Interface.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BaseFamilyDatabase.wl"];
char = 179424673; order = 4; rMax = 2; sMax = 1;
{pdlist, loopmom, extmom, spsRep, topSector, numeric, char} =
  familySetting[$SpecifiedFamilyDatabase, "TB123"];
numericD = {d -> (numeric["d"] /. d -> d)};
{ne, nl, nE, nibp, Alist, vlist, ibpeqs1, sectorlist} =
  AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topSector];
ibpeqsNum = PolynomialMod[ibpeqs1 /. numeric /. numericD, char];
regions = First[Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/intermediate/TB123_regions.m"]];
Print["imported sectors: ", Length[Keys[regions]], "  ", DateString[]];

t1 = AbsoluteTime[];
{hlist, alist} = reconstructExpansionR[ibpeqsNum, order, sectorlist, Alist, vlist,
  Modulus -> char, "LayerByLayer" -> True, "ImportRegions" -> regions];
tExpand = AbsoluteTime[] - t1;
Print["#regions: ", Length[alist], " expand: ", tExpand, " s  ", DateString[]];
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/intermediate/TB123_ha.m", {hlist, alist}];

t2 = AbsoluteTime[];
rel = reconstructRelation5[3, 2, topSector, hlist, alist, vlist,
  "coefdeg" -> 3, Modulus -> 0];
tRel = AbsoluteTime[] - t2;
Print["Relations dims: ", Dimensions[rel], "  Time: ", tRel, " s"];
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/intermediate/TB123_rel.m", rel];

t3 = AbsoluteTime[];
nLev = Min[2, Length[rel]];
nDeg = Min[2, Length[rel[[1]]]];
relForCone = Flatten /@ rel[[1 ;; nLev, 1 ;; nDeg]];
blockDataList = setupGeneratingCone[relForCone, {topSector}, vlist,
  Modulus -> 0, Verbose -> False];
tCone = AbsoluteTime[] - t3;
If[AssociationQ[blockDataList],
  Print["cone: ", Length[Keys[blockDataList]], " sectors. Time: ", tCone, " s"],
  Print["cone: non-Association. Time: ", tCone, " s"]];

targIntgList = seedsGenerate[topSector, rMax, sMax];
Print["#targets: ", Length[targIntgList]];
t5 = AbsoluteTime[];
reductionRule = concreteFamilyReduce[blockDataList, targIntgList, vlist, "FullReduce" -> True];
tReduce = AbsoluteTime[] - t5;
Print["reduce: ", tReduce, " s"];

Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/lie_e2e/TB123.m", {"reductionRule" -> reductionRule}];
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/lie_e2e/TB123_timing.json", <|
  "family" -> "TB123", "expand" -> tExpand, "reconstruct" -> tRel,
  "cone" -> tCone, "reduce" -> tReduce, "n_targets" -> Length[targIntgList],
  "total" -> tExpand + tRel + tCone + tReduce|>];
Print["ALL DONE. ", DateString[]];
