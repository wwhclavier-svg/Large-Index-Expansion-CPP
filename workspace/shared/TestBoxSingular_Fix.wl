Needs["LIEFamilyDefine`"];
Get["LIEWorkflow.wl"];
Get["SingularInterface.wl"];
Get["SingularCoordinateRing.wl"];

(* Box family *)
Box = FamilyDatabase`Box;
Alist = Box["Alist"];
vlist = Box["vlist"];

char = 179424673;

(* Sector {1,1,1,1} *)
sector = {1, 1, 1, 1};
Print["=== Box sector ", sector, " ==="];

time = AbsoluteTiming[
  res = regionsBySectors[Box["ibpeqs"], {Alist}, {vlist}, {sector},
    Modulus -> char, "ZeroDim" -> True, Verbose -> True];
][[1]];
Print["Time: ", time];

Do[
  Print["Region ", i, ": VarDeg=", res[[i]]["VarDeg"], 
    ", VarRule length=", Length[res[[i]]["VarRule"]],
    ", MinPoly length=", Length[res[[i]]["MinPoly"]]];
  If[Length[res[[i]]["VarRule"]] > 0,
    Print["  VarRule[[1]]=", res[[i]]["VarRule"][[1]]];
  ];
, {i, Length[res]}];

(* Check VarRule FullForm *)
Do[
  If[!FreeQ[FullForm[res[[i]]["VarRule"]], "Part"],
    Print["ERROR: Region ", i, " has Part expression in VarRule"];
  ];
, {i, Length[res]}];

Print["=== DONE ==="];
