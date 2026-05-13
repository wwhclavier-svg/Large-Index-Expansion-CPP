(* list_bub00_sectors.wl *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$data = Import[$LIECPPPath <> "/verify/bub00/PrepareCheckpoint-bub00.wdx", "WDX"];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

existingRegions = data["Regions"];
Print["Completed sectors: ", Length[Keys[existingRegions]]];
Do[
  secs = Select[Keys[existingRegions], Total[#] == n &];
  If[Length[secs] > 0,
    regs = Total[Length[existingRegions[#]] & /@ secs];
    Print["  ", n, " props: ", Length[secs], " sectors, ", regs, " regions"];
    Do[Print["    ", s, " -> ", Length[existingRegions[s]], " region(s)"], {s, Take[secs, UpTo[3]]}];
  ];
, {n, 1, 4}];
