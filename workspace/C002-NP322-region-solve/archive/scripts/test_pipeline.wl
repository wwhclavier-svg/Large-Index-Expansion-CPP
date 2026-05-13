(* test_pipeline.wl -- Test symbolic IBP generation + region solving on TB123, SR212, NP222, NP322 *)
(* v3: Uses proper symbolic expressions (not strings) matching FamilyDatabase format *)

$LIECPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPath <> "/verify/VerifyUtility/";
SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];

testFamily[name_, pdlist_, loopmom_, extmom_, spsRep_, topsector_, numeric_] := Module[
  {wf, ibpeqs, sectorlist, Alist, vlist, top, result, nNonTrivial, nRegions},
  Print["\n========== ", name, " =========="];
  wf = LIEDefineFamily[pdlist, loopmom, extmom, spsRep, topsector,
    "Numeric" -> numeric, Modulus -> 179424673, Verbose -> False];

  {ibpeqs, sectorlist, Alist, vlist} = wf["Family"] /@ {"IBPEqs","SectorList","AList","VList"};
  top = topsector;
  Print["Total sectors: ", Length[sectorlist], " | nProp=", Length[pdlist], " nLoop=", Length[loopmom], " nExt=", Length[extmom]];

  result = regionsBySectors[ibpeqs, {top}, Alist, vlist,
    Modulus -> 179424673, "EnableFieldExtension" -> True, Verbose -> False, "Timeout" -> 120];

  nNonTrivial = Length[Keys[result]];
  nRegions = Total[Length /@ Values[result]];
  Print["Top sector regions - non-trivial: ", nNonTrivial, ", total regions: ", nRegions];
  If[nNonTrivial > 0,
    Print["  Region counts: ", Length /@ Values[result]];
    Print["  Sector key: ", Keys[result][[1]]];
  ];
  nRegions
];

(* TB123: 7 propagators, 2L2P, top sector *)
r1 = testFamily["TB123",
  {-l1^2, -(l1-p1)^2, -(l1-p1-p2)^2, -(l2+p1+p2)^2, -l2^2, -(l1+l2)^2, -(l2+p1)^2},
  {l1, l2}, {p1, p2},
  {p1^2 -> s1, p2^2 -> s2, p1*p2 -> (s - s1 - s2)/2},
  {1,1,1,1,1,1,1}, {s -> 3, s1 -> 2, s2 -> 4}];

(* SR212: 5 propagators, 2L1P, top sector *)
r2 = testFamily["SR212",
  {-l1^2, -(l1+p)^2, -l2^2, -(l2+p)^2, -(l1+l2+p)^2},
  {l1, l2}, {p},
  {p^2 -> s},
  {1,1,1,1,1}, {s -> 1, "d" -> 1/13}];

(* NP222: 7 propagators, 2L2P, top sector *)
r3 = testFamily["NP222",
  {-(l2+p1)^2+msq, -(l1-l2-p1-p2)^2, -l2^2+msq, -(l1-p2)^2+msq, -(l1-l2)^2, -l1^2+msq, -(l1+p1)^2},
  {l1, l2}, {p1, p2},
  {p1^2 -> 0, p2^2 -> s2, p1*p2 -> (s - s2)/2},
  {1,1,1,1,1,1,1}, {s -> 1, s2 -> 0, msq -> 0, "d" -> 1/13}];

(* NP322: 9 propagators, 2L3P, top sector *)
r4 = testFamily["NP322",
  {-(k2+p1)^2+msq, -(k1-k2-p1-p2)^2, -k2^2+msq, -(k1-p2)^2+msq, -(k1-k2)^2, -k1^2+msq, -(k1-k2+p4)^2, -(k1+p1)^2, -(k1+p4)^2},
  {k1, k2}, {p1, p2, p4},
  {p1^2 -> m1, p2^2 -> m2, p4^2 -> 0, p1*p2 -> (s - m1 - m2)/2, p1*p4 -> (t - m1)/2, p2*p4 -> -(s + t - m1)/2},
  {1,1,1,1,1,1,1,1,1}, {s -> 1, t -> 5, m1 -> 0, m2 -> 0, msq -> 0, "d" -> 1/7}];

Print["\n========== SUMMARY =========="];
Print["TB123: ", r1, " regions (top sector)"];
Print["SR212: ", r2, " regions (top sector)"];
Print["NP222: ", r3, " regions (top sector)"];
Print["NP322: ", r4, " regions (top sector)"];
