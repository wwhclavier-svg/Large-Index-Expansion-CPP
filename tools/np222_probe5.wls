(* np222_probe5.wls — 假设检验：5-prop NP222 = 7-prop 族的前 5 个传播子 *)
(* 指纹：ne=5, nibp=6, 9 个 master sector 与 bin 一致 *)
Get["/home/ykm/work-knowledge/Projects/Large-Index-Expansion-MMA/LargeIndexExpansion.wl"];

pdlist = {(l2 + p1)^2, (l1 - l2 - p1 - p2)^2, l2^2, (l1 - p2)^2, (l1 - l2)^2};
loopmom = {l1, l2};
extmom = {p1, p2};
numeric = {s -> 1, s2 -> 0, d -> 1/13};
spsRep = {p1^2 -> 0, p2^2 -> s2, p1*p2 -> (s - s2)/2} /. numeric;
topsector = {1, 1, 1, 1, 1};

Print["=== AsyIBPFamilyDefine: first-5-props candidate ==="];
{ne, nl, nE, nibp, Alist, vlist, ibpeqs1, sectorlist} =
  AsyIBPFamilyDefine[pdlist /. numeric, loopmom, extmom, spsRep, topsector];

Print["ne=", ne, "  nl=", nl, "  nE=", nE, "  nibp=", nibp];
Print["#masters=", Length[sectorlist]];
Print["masters=", sectorlist];

expected = {{1,1,0,0,1},{1,0,1,0,1},{0,1,0,1,1},{0,0,1,1,1},
            {1,1,1,0,1},{1,1,0,1,1},{1,0,1,1,1},{0,1,1,1,1},{1,1,1,1,1}};
Print["MATCH_SECTORS=", Sort[sectorlist] === Sort[expected]];

Export["/home/ykm/Large-Index-Expansion-NOTES/Code/Large-Index-Expansion-CPP/tools/np222_probe5_result.m",
  <|"ne"->ne, "nibp"->nibp, "masters"->sectorlist,
    "match"->(Sort[sectorlist]===Sort[expected]),
    "ibpeqs"->(ibpeqs1 /. d -> 1/13)|>, "Package"];
Print["DONE"];
