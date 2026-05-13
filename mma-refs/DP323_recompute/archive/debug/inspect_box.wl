(* inspect_box.wl *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
data = Import[$LIECPPPath <> "/verify/Box/PrepareCheckpoint-Box.wdx", "WDX"];
Print["Keys: ", Keys[data]];
Print["Regions head: ", Head[data["Regions"]]];
If[AssociationQ[data["Regions"]],
  Print["Region keys (first 3): ", Take[Keys[data["Regions"]], 3]];
  Print["Region keys heads: ", Head /@ Take[Keys[data["Regions"]], 3]];
];
