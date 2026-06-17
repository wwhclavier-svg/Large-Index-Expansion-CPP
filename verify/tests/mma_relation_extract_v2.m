(* mma_relation_extract_v2.m *)
(* 从已有展开数据提取 MMA 关系 *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];

Print["=== 加载包 ==="];
Get["LIEReconstruct.wl"];

(* 读取展开数据 *)
expansionRaw = Get["/root/Large-Index-Expansion-CPP/verify/bub00/VerifyExpansion-MMAExpansion.m"];

(* 访问 Solutions/H 列表 -- 正确路径 *)
reg = expansionRaw[[1]];
sol = reg[["Solutions"]][[1]];
hlistRaw = sol[["H"]];

kmax = reg[["Kmax"]];
ne = reg[["NE"]];
nb = reg[["NB"]];
topsec = reg[["TopSector"]];
modulus = reg[["Modulus"]];

(* 提取 A 矩阵 *)
qRingVarRules = reg[["QRingVarRule"]];
aregList = {1, qRingVarRules[["A"[2]]]};
aregList = {aregList};

(* 转换字符串为符号 *)
strToSymRules = Dispatch[{
    "v1" -> v1, "v2" -> v2,
    "v1^2" -> v1^2, "v2^2" -> v2^2,
    "v1^3" -> v1^3, "v2^3" -> v2^3,
    "v1^4" -> v1^4, "v2^4" -> v2^4,
    "v1^5" -> v1^5, "v2^5" -> v2^5,
    "v1^6" -> v1^6, "v2^6" -> v2^6,
    "v1^7" -> v1^7, "v2^7" -> v2^7,
    "v1^8" -> v1^8, "v2^8" -> v2^8
}];

strToPoly[str_String] := ToExpression[str] /. strToSymRules;
strToPoly[other_] := other;

vlist = {v1, v2};
hlist = Table[strToPoly[hlistRaw[[k+1]]], {k, 0, kmax}];

Print["Kmax = ", kmax, ", NE = ", ne, ", NB = ", nb];
Print["TopSector = ", topsec, ", Modulus = ", modulus];
Print["ARegList = ", aregList];
Print["H[0] = ", hlist[[1]]];
Print["H[1] = ", Short[hlist[[2]]]];
Print["H[2] short = ", Short[hlist[[3]]]];

(* 对 lev=2 运行 LIEReconstruct *)
pairs = {{2, 0}, {2, 1}, {2, 2}};
vlist = {v1, v2};

Do[
    {lev, deg} = pair;
    Print["\n=== LIEReconstruct[lev=", lev, ", deg=", deg, "] ==="];
    
    result = LIEReconstruct[lev, deg, kmax, hlist, aregList, topsec, vlist,
        Modulus -> modulus,
        Verbose -> True
    ];
    
    relList = result;
    Print["result type: ", Head[result]];
    Print["result dimensions: ", Dimensions[result]];
    
    If[Length[relList] > 0 && Length[relList[[1]]] > 0,
        Print["Has non-trivial relations: ", Length[relList[[1, deg+1]]]];
        
        (* Export relations to file *)
        outPath = "/root/Large-Index-Expansion-CPP/verify/bub00/Compare-MMARelation-bub00_lev" <> ToString[lev] <> "_deg" <> ToString[deg] <> ".m";
        
        stream = OpenWrite[outPath];
        WriteString[stream, "(* MMA Relation Export *)\n"];
        WriteString[stream, "(* Family: bub00, Lev=" <> ToString[lev] <> ", Deg=" <> ToString[deg] <> " *)\n"];
        WriteString[stream, "$MMARelationResult = <|\n"];
        WriteString[stream, "  \"HasSolution\" -> True,\n"];
        WriteString[stream, "  \"Relations\" -> \n"];
        
        relations = relList[[1, deg+1]];
        WriteString[stream, "    " <> ToString[InputForm[relations]] <> "\n"];
        WriteString[stream, "|>;\n"];
        Close[stream];
        Print["Exported to: ", outPath];
    ,
        Print["No non-trivial relations found"]
    ],
    {pair, pairs}
];

Print["\n=== Done ==="];
