(* inspect_bub00.wl *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
SetDirectory[$LIECPPPath <> "/verify/bub00/"];
Block[{Print = Identity}, data = Import["PrepareCheckpoint-bub00.wdx", "WDX"]];
Print["Regions head: ", Head[data["Regions"]]];
If[AssociationQ[data["Regions"]],
  Print["First key: ", First[Keys[data["Regions"]]]];
  Print["First key head: ", Head[First[Keys[data["Regions"]]]]];
];
