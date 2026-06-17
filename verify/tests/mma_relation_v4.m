(* mma_relation_v4.m *)
(* 方案A：修复 rank=2 后重新运行，导出关系系数 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];

Get["LIEWorkflow.wl"];
Get["LIEReconstruct.wl"];

modulus = Prime[10000000];
numericRules = {s -> 3, msq -> 0, "d" -> 1/3};

bub00Config = <|
    "Propagators" -> ({-k1^2 + msq, -(k1 - p1)^2 + msq} /. numericRules),
    "LoopMomenta" -> {k1},
    "ExternalMomenta" -> {p1},
    "KinematicRules" -> ({p1^2 -> s} /. numericRules),
    "TopSector" -> {1, 1},
    "Numeric" -> numericRules,
    "Modulus" -> modulus
|>;

data = LIEDefineFamily[
    bub00Config["Propagators"], bub00Config["LoopMomenta"],
    bub00Config["ExternalMomenta"], bub00Config["KinematicRules"],
    bub00Config["TopSector"],
    "Numeric" -> bub00Config["Numeric"],
    Modulus -> modulus, Verbose -> False
];

data = LIESolveRegions[data, Verbose -> False];

data = LIEExpandSeries[data,
    "Order" -> 4, Modulus -> modulus, "Increment" -> 2,
    "LayerByLayer" -> True, Verbose -> False
];

(* 关键修复：传入 "Rank" -> 2 限制层数 *)
data = LIEGetRelations[data,
    "Rank" -> 2,           (* 之前丢失的参数！ *)
    "MaxCoefDeg" -> 2,
    Verbose -> True
];

relations = data["Relations", "Relations"];
Print["\n=== MMA Relation Dimensions ==="];
Print["Dims: ", Dimensions[relations]];

(* 打印每层每度的关系数 *)
Do[
    Do[
        relAtDeg = relations[[lev, deg + 1]];
        nRel = If[ListQ[relAtDeg], Length[relAtDeg], 0];
        Print["(lev=" <> ToString[lev] <> ", deg=" <> ToString[deg] <> ") -> " <> ToString[nRel] <> " relations"];
        If[nRel > 0 && nRel <= 10,
            Do[
                relExpr = relAtDeg[[r]];
                Print["  rel[" <> ToString[r] <> "] = ", Short[InputForm[relExpr], 2]];
            , {r, nRel}];
        ];
    , {deg, 0, 2}];
, {lev, 1, Length[relations]}];

(* 导出系数矩阵 *)
outPath = "/root/Large-Index-Expansion-CPP/verify/bub00/Compare-MMARelation-lev2.m";
stream = OpenWrite[outPath];
WriteString[stream, "(* MMA LIE Relations - rank=2, maxDeg=2 *)\n"];
WriteString[stream, "$MMARelations = " <> ToString[InputForm[relations]] <> ";\n"];
WriteString[stream, "$MMAConfig = <|\"Rank\" -> 2, \"MaxDeg\" -> 2, \"Modulus\" -> " <> ToString[modulus] <> "|>;\n"];
Close[stream];
Print["\nExported to: ", outPath];
Print["\n=== Done ==="];
