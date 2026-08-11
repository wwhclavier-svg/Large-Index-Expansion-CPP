(* TB123 e2e 第二段: 从存档 regions 接续 *)
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SingularInterface.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/QuotientAlgebra.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LargeIndexExpansion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/BlockRecursion.wl"];
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/SymbolicBTForm.wl"];

char = 179424673;
order = 4; r = 2; s = 1;
pdlist = {l1^2, (l1-p1)^2, (l1-p1-p2)^2, (l2+p1+p2)^2, l2^2, (l1+l2)^2, (l2+p1)^2};
topSector = {1,1,1,1,1,1,0}; sectors = {topSector};
{ne, nl, nE, nibp, Alist, vlist, ibpeqs1, sectorlist} =
  AsyIBPFamilyDefine[pdlist, {l1,l2}, {p1,p2}, {p1^2->2, p2^2->4, p1*p2->-3/2}, topSector];
ibpeqsNum = PolynomialMod[ibpeqs1 /. d -> 1/13, char];
regions = Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/intermediate/TB123_regions.m"];
Print["regions imported: ", DateString[]];

t1 = AbsoluteTime[];
{hlist, alist} = reconstructExpansionR[ibpeqsNum, order, sectorlist, Alist, vlist,
  Modulus -> char, "LayerByLayer" -> True, "ImportRegions" -> regions];
tExpand = AbsoluteTime[] - t1;
Print["  #regions: ", Length[alist], " expansion time: ", tExpand, " s  ", DateString[]];

t2 = AbsoluteTime[];
rel = reconstructRelation5[3, 2, topSector, hlist, alist, vlist,
  "coefdeg" -> 3, Modulus -> 0];
tRel = AbsoluteTime[] - t2;
Print["  Relations: ", Dimensions[rel], " Time: ", tRel, " s  ", DateString[]];

t3 = AbsoluteTime[];
nLev = Min[2, Length[rel]];
nDeg = Min[2, Length[rel[[1]]]];
relForCone = Flatten /@ rel[[1 ;; nLev, 1 ;; nDeg]];
blockDataList = setupGeneratingCone[relForCone, sectors, vlist,
  Modulus -> 0, Verbose -> False];
tCone = AbsoluteTime[] - t3;
If[AssociationQ[blockDataList],
  Print["  setupGeneratingCone: ", Length[Keys[blockDataList]], " sectors. Time: ", tCone, " s"],
  Print["  setupGeneratingCone: returned non-Association. Time: ", tCone, " s"]];

targIntgList = seedsGenerate[topSector, r, s];
Print["  #targets: ", Length[targIntgList]];

t5 = AbsoluteTime[];
reductionRule = concreteFamilyReduce[blockDataList, targIntgList, vlist,
  "FullReduce" -> True];
tReduce = AbsoluteTime[] - t5;
Print["  concreteFamilyReduce Time: ", tReduce, " s"];

Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/lie_e2e/TB123.m", {"reductionRule" -> reductionRule}];
Export["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP/results/lie_e2e/TB123_timing.json", <|
  "family" -> "TB123", "expand" -> tExpand, "reconstruct" -> tRel,
  "cone" -> tCone, "reduce" -> tReduce, "n_targets" -> Length[targIntgList]|>];
Print["ALL DONE. ", DateString[]];
