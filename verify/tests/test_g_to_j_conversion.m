(* ::Package:: *)
(* test_g_to_j_conversion.m *)
(* 精确定位 bub00 KiraVerify FAIL 的根因: g→j 转换和 Kira 规则覆盖 *)
(* 独立运行，不依赖外部工作流 *)

AppendTo[$Path, "/root/Large-Index-Expansion-MMA-Mini"];
Get["LIEWorkflow.wl"];

Print["========================================"];
Print[" g→j 转换与 Kira 规则覆盖诊断"];
Print["========================================"];

modulus = 179424673;
numericRules = {s -> 3, msq -> 0, d -> 1/3};

(* -------------------------------------------------- *)
(* Step 1: 生成 MMA g-form 关系 *)
(* -------------------------------------------------- *)
Print["\n[Step 1] 生成 g-form 关系..."];

bub00Config = <|
    "Propagators" -> {-k1^2 + msq, -(k1 - p1)^2 + msq} /. numericRules,
    "LoopMomenta" -> {k1},
    "ExternalMomenta" -> {p1},
    "KinematicRules" -> {p1^2 -> s} /. numericRules,
    "TopSector" -> {1, 1},
    "Numeric" -> numericRules,
    "Modulus" -> modulus
|>;

data = LIEDefineFamily[bub00Config, "Verbose" -> False];
data = LIESolveRegions[data, "Verbose" -> False];
data = LIEExpandSeries[data, "Order" -> 4, "Increment" -> 2, "Verbose" -> False];
data = LIEGetRelations[data, "Verbose" -> False, "MaxCoefDeg" -> 2];

relations = data["Relations", "Relations"];
Print["  生成关系: ", Length /@ relations // Flatten // Total];

(* -------------------------------------------------- *)
(* Step 2: 加载 Kira 规则 *)
(* -------------------------------------------------- *)
Print["\n[Step 2] 加载 Kira 规则..."];

kiraPath = "/root/kira_tests/bub00_new/bub00";

If[!DirectoryQ[kiraPath],
    Print["  [SKIP] Kira 路径不存在: ", kiraPath];
    Print["  需要先运行 KiraInputGenerate + kira jobs.yaml"];
,
    (* 查找 .rules 文件 *)
    rulesFiles = FileNames["*.rules", {kiraPath}, 2];
    Print["  找到规则文件: ", rulesFiles];
    
    If[Length[rulesFiles] > 0,
        kiraRaw = Import[rulesFiles[[1]], "Text"];
        
        (* 解析 Kira 规则: j[a,b] -> expr *)
        (* Kira 规则格式: j[2,1] -> (d-3)*j[1,1] *)
        kiraRules = StringCases[kiraRaw,
            RegularExpression[".*j\\[(\\d+),\\s*(\\d+)\\].*->.*"] :> 
                Rule @@ ToExpression /@ {"j[" <> $1 <> ", " <> $2 <> "]"}
        ];
        
        (* 实际解析: 从规则文件提取 j[...] 索引 *)
        allKiraIndices = Union @ Flatten @ StringCases[kiraRaw,
            "j[" ~~ DigitCharacter .. ~~ "," ~~ DigitCharacter .. ~~ "]":>
                StringCases[#, DigitCharacter ..]& /@ {"j[" <> $& <> "]"}
        ];
        
        Print["  Kira 规则文件存在，长度: ", StringLength[kiraRaw]];
        Print["  规则文件前500字符: ", StringTake[kiraRaw, 500]];
    ,
        Print["  [WARN] 无 .rules 文件"];
    ];
];

(* -------------------------------------------------- *)
(* Step 3: 分析 g-form 关系的 j 索引覆盖 *)
(* -------------------------------------------------- *)
Print["\n[Step 3] 分析 g-form 关系的 j 索引覆盖..."];

(* 收集所有关系中的 j 索引（假设 g[α,β] → j[{α,β}]）*)
allGIndices = {};
ForEach[levRels, relations,
    ForEach[rel, levRels,
        gCases = Cases[rel, g[a_, b_] :> {a, b}, Infinity];
        allGIndices = Join[allGIndices, gCases];
    ];
];

uniqueGIndices = Union[allGIndices];
Print["  g[...] 索引种类: ", Length[uniqueGIndices]];
Print["  最大 α: ", Max[allGIndices[[All, 1]]]];
Print["  最大 β: ", Max[allGIndices[[All, 2]]]];
Print["  α,β 范围示例: ", Take[uniqueGIndices, UpTo[10]]];

(* -------------------------------------------------- *)
(* Step 4: g → Kira j 格式转换验证 *)
(* -------------------------------------------------- *)
Print["\n[Step 4] g → j 转换验证..."];

(* Kira j 格式: j[{a, b}] — 注意是 List 索引 *)
sampleG = g[-1 + v1, v2];  (* 来自 lev=0,deg=1 关系 *)
kiraEquivalent = j[{-1 + v1, v2}];  (* 转换后的 Kira 格式 *)

Print["  g 格式: ", InputForm[sampleG]];
Print["  Kira j 格式: ", InputForm[kiraEquivalent]];

(* 代入具体 ν 值验证 *)
testNu = {1, 1};
gVal = sampleG /. {v1 -> testNu[[1]], v2 -> testNu[[2]]};
kiraVal = kiraEquivalent /. {v1 -> testNu[[1]], v2 -> testNu[[2]]};
Print["  ν={1,1} 时 g = ", gVal, ", Kira j = ", InputForm[kiraVal]];

(* -------------------------------------------------- *)
(* Step 5: Kira 规则实际约化测试（如果规则可用）*)
(* -------------------------------------------------- *)
Print["\n[Step 5] Kira 规则约化测试..."];

If[Length[rulesFiles] > 0,
    (* 构造简单的测试用例 *)
    testJExpr = j[{1, 1}] + j[{2, 1}];  (* 符号 j *)
    Print["  测试表达式: ", InputForm[testJExpr]];
    
    (* 用原始规则文件做替换 *)
    (* 注意: Kira 规则是真正的 rewrite rules *)
    (* Kira 规则格式: j[2,1] -> (d-3)*j[1,1] *)
    
    (* 从规则文件提取具体规则 *)
    specificRules = StringCases[kiraRaw,
        "j[" ~~ a__ ~~ ", " ~~ b__ ~~ "] -> " ~~ rest__ :>
            With[{aa = ToExpression[a], bb = ToExpression[b]},
                j[{aa, bb}] :> 
                    ToExpression["(" <> rest <> ")"]
            ]
    ];
    Print["  提取到具体规则数: ", Length[specificRules]];
    Print["  规则示例: ", Take[specificRules, UpTo[3]]];
,
    Print["  [SKIP] Kira 规则不可用"]
];

(* -------------------------------------------------- *)
(* Step 6: 诊断报告 *)
(* -------------------------------------------------- *)
Print["\n========================================"];
Print[" 诊断结论"];
Print["========================================"];
Print["  g 索引范围: α∈[", Min[allGIndices[[All, 1]]], ",",
      Max[allGIndices[[All, 1]]], "], β∈[",
      Min[allGIndices[[All, 2]]], ",",
      Max[allGIndices[[All, 2]]], "]"];
Print["  g 索引种类: ", Length[uniqueGIndices];

If[Length[rulesFiles] > 0,
    Print["  Kira 规则: 可用"];
,
    Print["  Kira 规则: 不可用"];
    Print["  [建议] 运行 Kira 生成规则"];
];

Print[""];
Print["  待验证假设:"];
Print["  1. g[α,β] → j[{α,β}] 转换是否正确"];
Print["  2. Kira 规则是否覆盖所有出现的 {α,β}"];
Print["  3. Kira 约化后是否正确代入 d=1/3"];
Print["========================================"];
