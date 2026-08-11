(* LIE-vs-Kira-Compare — MG1: 端到端运行时间比较 *)
cppRoot = ExpandFileName["~/Large-Index-Expansion-CPP"];
outDir = FileNameJoin[{cppRoot, "results", "comparison"}];
If[FileType[outDir] === None, CreateDirectory[outDir]];

readLIETiming[family_] := Module[
  {f = FileNameJoin[{cppRoot, "results", "lie_e2e", family <> ".json"}], j, ts},
  If[!FileExistsQ[f], Return[<|"Total_s"->"N/A"|>]];
  j = Import[f, "JSON"];
  ts = Lookup[j, "timing_s", <||>];
  <|"Total_s" -> ToString[NumberForm[Lookup[ts, "total", "N/A"], {5,2}]],
    "Breakdown" -> KeyDrop[ts, "total"]|>
];

readKiraTiming[family_] := Module[
  {d = FileNameJoin[{cppRoot, "kira_tests", family}], f, lines, i},
  f = FileNameJoin[{d, "kira.log"}];
  If[!FileExistsQ[f], Return["N/A"]];
  lines = ReadList[f, String];
  For[i = 1, i <= Length[lines], i++,
    If[StringMatchQ[lines[[i]], "*Total time:*"] && i < Length[lines],
      Return[ToString[NumberForm[
        ToExpression[StringTrim[
          StringReplace[lines[[i+1]], {"("->"", ")"->"", "s"->""}]]],
      {4,1}]]]
    ]
  ];
  "N/A"
];

families = {"bub00", "bub10", "bub11", "Tri", "Penta1L", "SR212"};
Print["== LIE vs Kira 端到端时间比较 =="];
Print[""];
Print["family    LIE(s)    Kira(s)   speedup"];
Print["-------------------------------------"];

results = Table[
  lieT = readLIETiming[families[[i]]];
  kiraT = readKiraTiming[families[[i]]];
  lieNum = ToExpression[lieT["Total_s"]];
  kiraNum = ToExpression[kiraT];
  spd = If[NumericQ[lieNum] && NumericQ[kiraNum] && kiraNum > 0,
    ToString[NumberForm[kiraNum / lieNum, {4,2}]] <> "x", "N/A"];
  Print[StringPadRight[families[[i]], 9], lieT["Total_s"], "    ", kiraT, "    ", spd];
  <|"Family"->families[[i]],
    "LIE"->lieT,
    "Kira"-><|"Total_s"->kiraT|>,
    "Speedup"->If[NumericQ[lieNum] && NumericQ[kiraNum] && kiraNum > 0,
      Round[kiraNum / lieNum, 0.01], "N/A"]|>,
  {i, Length[families]}
];

Export[FileNameJoin[{outDir, "all_results.json"}],
  <|"Generated"->DateString[], "Families"->Association[#Family -> # & /@ results]|>,
  "JSON"];
Print[""];
Print["-> ", FileNameJoin[{outDir, "all_results.json"}]];
Print["Done."]
