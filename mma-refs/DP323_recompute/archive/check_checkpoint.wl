(* check_checkpoint.wl — 检查 DP323 checkpoint 缺了哪些 sector *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$FamilyGeneratePath = $LIECPPPath <> "/verify/DP323/";
$RunDir = $LIECPPPath <> "/verify/DP323_recompute/";

SectorToStr[sector_List] := StringJoin[ToString /@ sector];
SectorActiveProps[sec_] := Total[sec];

(* 加载现有数据 *)
cpFile = $FamilyGeneratePath <> "PrepareCheckpoint-DP323.wdx";
data = Import[cpFile, "WDX"];

If[Head[data["Regions"]] === Join,
  data["Regions"] = data["Regions"][[1]];
];

existingRegions = data["Regions"];
fullSectorList = data["Family", "SectorList"];
existingSectors = Keys[existingRegions];
missingSectors = Complement[fullSectorList, existingSectors];

Print["========================================"];
Print["DP323 Checkpoint Status"];
Print["========================================"];
Print["Total sectors defined: ", Length[fullSectorList]];
Print["Sectors in checkpoint: ", Length[existingSectors]];
Print["Sectors MISSING:       ", Length[missingSectors]];
Print["Total regions in CP:   ", Total[Length /@ Values[existingRegions]]];
Print[""];

Print["Missing sectors by active prop count:"];
missingByCount = KeySort[Counts[SectorActiveProps /@ missingSectors]];
Do[
  nRegAll = Count[SectorActiveProps /@ fullSectorList, n];
  nRegPres = Count[SectorActiveProps /@ existingSectors, n];
  Print["  ", n, " props: ", If[KeyExistsQ[missingByCount, n], missingByCount[n], 0],
        " missing (", nRegPres, "/", nRegAll, " present)"];
, {n, 1, 11}];

Print["\nMissing sectors breakdown:"];
Do[
  nProps = SectorActiveProps[s];
  Print["  ", SectorToStr[s], " (", nProps, " props)"];
, {s, missingSectors}];
