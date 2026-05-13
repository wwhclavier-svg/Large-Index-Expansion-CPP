(* inspect_box2.wl *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
SetDirectory[$LIECPPPath <> "/verify/Box/"];
data = Import["PrepareCheckpoint-Box.wdx", "WDX"];

Print["data head: ", Head[data]];
Print["Regions head: ", Head[data["Regions"]]];

reg = data["Regions"];
If[AssociationQ[reg],
  k = Keys[reg];
  Print["Number of keys: ", Length[k]];
  Print["First key: ", k[[1]]];
  Print["First key head: ", Head[k[[1]]]];
  Print["First value head: ", Head[reg[k[[1]]]]];
  Print["First value length: ", Length[reg[k[[1]]]]];
];
