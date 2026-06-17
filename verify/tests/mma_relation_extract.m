(* mma_relation_extract.m *)
(* 从 MMA 展开数据中提取关系并导出为可比较格式 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];

(* 解析 v1,v2 为符号（非字符串）*)
strToSymRules = {
    "\"v1\"" -> v1, "\"v2\"" -> v2,
    "\"v1^2\"" -> v1^2, "\"v2^2\"" -> v2^2,
    "\"v1^3\"" -> v1^3, "\"v2^3\"" -> v2^3,
    "\"v1^4\"" -> v1^4, "\"v2^4\"" -> v2^4
};
strToSym[str_String] := ToExpression[str];
strToSym[other_] := other;

(* 加载 LIEReconstruct *)
Get["LIEReconstruct.wl"];

(* 读取展开数据 *)
expansionData = Get["/root/Large-Index-Expansion-CPP/verify/bub00/VerifyExpansion-MMAExpansion.m"];

(* 提取数据 *)
reg = expansionData[[1]];
kmax = reg["Kmax"];
ne = reg["NE"];
nb = reg["NB"];
nimax = reg["Nimax"];
topsec = reg["TopSector"];
modulus = reg["Modulus"];
vlist = {v1, v2};

(* 提取 A 矩阵 *)
aregList = {};
Do[
    key = "A" <> ToString[i];
    If[KeyExistsQ[reg["QRingVarRule"], key],
        AppendTo[aregList, {1, reg["QRingVarRule", key]}],
        Break[]
    ],
    {i, 1, 10}
];

(* 提取 h_k 多项式 - 转换字符串为符号 *)
hlistRaw = reg["Solutions", 1, "H"];
hlist = Table[
    If[Head[hlistRaw[[k+1]]] === String,
        ToExpression[hlistRaw[[k+1]]],
        hlistRaw[[k+1]]
    ],
    {k, 0, kmax}
];

Print["=== MMA 展开数据 ==="];
Print["Kmax = ", kmax, ", NE = ", ne, ", NB = ", nb];
Print["Nimax = ", nimax, ", TopSector = ", topsec];
Print["Modulus = ", modulus];
Print["ARegList = ", aregList];
Print["H[0] = ", hlist[[1]]];
Print["H[1] = ", hlist[[2]]];
Print["H[2] = ", hlist[[3]]];
Print["H[3] = ", hlist[[4]]];

(* 对每个 (lev, deg) 运行 LIEReconstruct *)
pairs = {{2, 0}, {2, 1}, {2, 2}};

Do[
    {lev, deg} = pair;
    Print["\n========================"];
    Print["LIEReconstruct[lev=", lev, ", deg=", deg, "]"];
    Print["========================"];
    
    result = LIEReconstruct[lev, deg, kmax, hlist, aregList, topsec, vlist,
        Modulus -> modulus,
        Verbose -> True
    ];
    
    Print["HasSolution: ", result["HasSolution"]];
    
    If[result["HasSolution"],
        (* extract coefficient vectors *)
        coeffs = result["Coefficients"];
        alphas = result["Alphas"];
        betas = result["Betas"];
        
        Print["Number of relations: ", Dimensions[coeffs]];
        Print["Alphas: ", alphas];
        Print["Betas: ", betas];
        
        (* Export unified format *)
        (* Build: each row of coeffs = coefficients for one (alpha, beta) pair *)
        (* coeffs[[a,b,sol]] = c_{alpha_a, beta_b, sol} *)
        (* Flatten to vector: vector[sol][a*nBeta + b] = coeffs[[a,b,sol]] *)
        
        (* Create a simple export *)
        outPath = "/root/Large-Index-Expansion-CPP/verify/bub00/Compare-MMARelation-bub00_lev" <> ToString[lev] <> "_deg" <> ToString[deg] <> ".m";
        
        stream = OpenWrite[outPath];
        WriteString[stream, "(* MMA Relation Export: Unified Format *)\n"];
        WriteString[stream, "(* Family: bub00, Lev=" <> ToString[lev] <> ", Deg=" <> ToString[deg] <> " *)\n"];
        WriteString[stream, "$MMARelationResult = <|\n"];
        WriteString[stream, "  \"Family\" -> \"bub00\",\n"];
        WriteString[stream, "  \"Lev\" -> " <> ToString[lev] <> ",\n"];
        WriteString[stream, "  \"Deg\" -> " <> ToString[deg] <> ",\n"];
        WriteString[stream, "  \"NE\" -> " <> ToString[ne] <> ",\n"];
        WriteString[stream, "  \"Modulus\" -> " <> ToString[modulus] <> ",\n"];
        
        (* Export alphas *)
        WriteString[stream, "  \"Alphas\" -> {\n"];
        Do[
            WriteString[stream, "    {" <> StringRiffle[ToString /@ alphas[[i]], ","] <> "}"];
            If[i < Length[alphas], WriteString[stream, ","]];
            WriteString[stream, "\n"],
            {i, Length[alphas]}
        ];
        WriteString[stream, "  },\n"];
        
        (* Export betas *)
        WriteString[stream, "  \"Betas\" -> {\n"];
        Do[
            WriteString[stream, "    {" <> StringRiffle[ToString /@ betas[[j]], ","] <> "}"];
            If[j < Length[betas], WriteString[stream, ","]];
            WriteString[stream, "\n"],
            {j, Length[betas]}
        ];
        WriteString[stream, "  },\n"];
        
        (* Export coefficient matrix *)
        WriteString[stream, "  \"Coefficients\" -> {\n"];
        Do[
            coeffRow = Table[
                coeffs[[a, b, sol]],
                {a, Length[alphas]}, {b, Length[betas]}
            ] // Flatten;
            WriteString[stream, "    {" <> StringRiffle[ToString /@ coeffRow, ", "] <> "}"];
            If[sol < Dimensions[coeffs][[3]], WriteString[stream, ","]];
            WriteString[stream, "\n"],
            {sol, Dimensions[coeffs][[3]]}
        ];
        WriteString[stream, "  },\n"];
        
        (* Export symbolic relations for reference *)
        WriteString[stream, "  \"Relations\" -> {\n"];
        Do[
            (* Build symbolic expression *)
            terms = Table[
                coeffs[[a, b, sol]] * 
                (Times @@ (vlist^betas[[b]])) * 
                Apply[g, vlist - alphas[[a]]],
                {a, Length[alphas]}, {b, Length[betas]}
            ];
            rel = Total[Flatten[terms]];
            WriteString[stream, "    \"" <> ToString[InputForm[rel], InputForm] <> "\""];
            If[sol < Dimensions[coeffs][[3]], WriteString[stream, ","]];
            WriteString[stream, "\n"],
            {sol, Dimensions[coeffs][[3]]}
        ];
        WriteString[stream, "  },\n"];
        
        WriteString[stream, "  \"HasSolution\" -> True,\n"];
        WriteString[stream, "  \"NumVariables\" -> " <> ToString[Length[alphas] * Length[betas]] <> ",\n"];
        WriteString[stream, "  \"NumSolutions\" -> " <> ToString[Dimensions[coeffs][[3]]] <> "\n"];
        WriteString[stream, "|>;\n"];
        Close[stream];
        
        Print["Exported to: ", outPath];
    ,
        Print["No non-trivial relations found."]
    ],
    {pair, pairs}
];

Print["\n=== Done ==="];
