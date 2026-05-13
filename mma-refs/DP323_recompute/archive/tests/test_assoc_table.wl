(* Test Association[Table] after loading packages *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
SetDirectory[$LIECPPPath <> "/verify/VerifyUtility/"];
Get["LIEWorkflow.wl"];

(* Simple test *)
Print["Test 1: Association[Table[i -> i^2, {i, 3}]]"];
Print["  Result: ", Association[Table[i -> i^2, {i, 3}]]];

(* Test with Module locals *)
Print["\nTest 2: Module with Table inside Association"];
result = Module[{x, y},
  Association[Table[
    x = i;
    y = i^2;
    x -> y,
    {i, 3}
  ]]
];
Print["  Result: ", result];

(* Test with complex body like regionsBySectors *)
Print["\nTest 3: Complex Table body"];
sectorlist = {{0, 1, 0, 1, 1}};
result2 = Module[{sec, i, expregsub},
  Association[Table[
    sec = sectorlist[[i]];
    expregsub = {sec};
    If[Length[expregsub] > 0,
      sec -> expregsub,
      Nothing
    ],
    {i, Length[sectorlist]}
  ]]
];
Print["  Result: ", result2];
Print["  AssociationQ: ", AssociationQ[result2]];
