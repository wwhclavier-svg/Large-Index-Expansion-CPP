(* test_regions_direct_SR212.wl — Run regionsBySectors exactly as VerifyExpand-Prepare does *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/verify/SR212/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

Print["Loading SR212 checkpoint..."];
data = Import[$FamilyPath <> "PrepareCheckpoint-SR212.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
sectorlist = data["Family", "SectorList"];

Print["Total sectors: ", Length[sectorlist]];

(* Run exactly as LIESolveRegions does *)
result = LIERegions`regionsBySectors[ibpeqs, sectorlist, Alist, vlist,
  Modulus -> char, "EnableFieldExtension" -> True, Verbose -> False];

Print["\nResults:"];
Do[
  sec = sectorlist[[i]];
  If[KeyExistsQ[result, sec],
    n = Length[result[sec]];
    Print[sec, " -> ", n, " region(s)"];
  ,
    Print[sec, " -> NOT FOUND"];
  ];
, {i, Length[sectorlist]}];

Print["\nDone"];
