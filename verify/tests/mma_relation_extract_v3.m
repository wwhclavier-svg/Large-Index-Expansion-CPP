(* mma_relation_extract_v3.m *)
(* 运行完整 LIE 工作流，提取关系系数向量，导出用于对比 *)

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

data = LIEGetRelations[data, Verbose -> False, "MaxCoefDeg" -> 2];

relations = data["Relations", "Relations"];
Print["relationList type: ", Head[relations]];
Print["relationList dims: ", Dimensions[relations]];

(* The relations output is: relationList[[level_idx, degree_idx+1]] *)
(* Each entry is {dependent_var -> coeff_vector, ...} or a matrix *)

Do[
    level = lev;
    Do[
        deg = d;
        relAtLevel = relations[[level]];
        relAtDeg = relAtLevel[[deg + 1]];
        nRel = If[Head[relAtDeg] === List, Length[relAtDeg], 0];
        Print["\n(lev=" <> ToString[level] <> ", deg=" <> ToString[deg] <> ") -> nRel=" <> ToString[nRel]];
        
        If[nRel > 0,
            Do[
                Print["  rel[" <> ToString[r] <> "]: ", Short[InputForm[relAtDeg[[r]]], 3]];
            , {r, nRel}];
        ];
    , {deg, 0, 2}];
, {lev, 1, Length[relations]}];

(* Export to file *)
outPath = "/root/Large-Index-Expansion-CPP/verify/bub00/Compare-MMARelations_raw.m";
stream = OpenWrite[outPath];
WriteString[stream, "(* MMA LIE Relations - Raw Export *)\n"];
WriteString[stream, "$MMARelations = " <> ToString[InputForm[relations]] <> ";\n"];
Close[stream];

Print["\nExported to: ", outPath];
Print["\n=== Done ==="];
