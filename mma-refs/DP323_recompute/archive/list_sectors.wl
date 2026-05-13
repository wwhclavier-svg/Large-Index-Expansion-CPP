(* list_sectors.wl — List sectors in DP323 checkpoint *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$FamilyGeneratePath = $LIECPPPath <> "/verify/DP323/";

data = Import[$FamilyGeneratePath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

fullSectorList = data["Family", "SectorList"];
existingRegions = data["Regions"];
existingSectors = Keys[existingRegions];

Print["Total sectors: ", Length[fullSectorList]];
Print["Completed sectors: ", Length[existingSectors]];
Print["\nCompleted sectors by prop count:"];
Do[
  secs = Select[existingSectors, Total[#] == n &];
  regs = Total[Length[existingRegions[#]] & /@ secs];
  If[Length[secs] > 0,
    Print["  ", n, " props: ", Length[secs], " sectors, ", regs, " regions"];
  ];
, {n, 1, 11}];

Print["\nSample sectors with 2 props (completed):"];
sample2 = Select[existingSectors, Total[#] == 2 &];
Do[Print["  ", s, " -> ", Length[existingRegions[s]], " region(s)"], {s, Take[sample2, UpTo[5]]}];
