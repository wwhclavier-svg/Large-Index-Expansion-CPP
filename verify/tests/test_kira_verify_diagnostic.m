(* ::Package:: *)
(* test_kira_verify_diagnostic.m *)
(* 诊断 bub00 (d=1/3, s=3) MMA-KiraVerify 0/10 FAIL 的根因 *)
(* 定位是关系本身问题、g→j 转换问题、还是 Kira 规则问题 *)

AppendTo[$Path, "/root/Large-Index-Expansion-MMA-Mini"];
Get["LIEWorkflow.wl"];

Print["========================================"];
Print[" KiraVerify FAIL 根因诊断"];
Print["========================================"];

modulus = 179424673;
numericRules = {s -> 3, msq -> 0, d -> 1/3};

(* -------------------------------------------------- *)
(* Step 1: 生成 g-form 关系（MMA LIEReconstruct）*)
(* -------------------------------------------------- *)
Print["\n[Step 1] 生成 g-form 关系..."];

config = <|
    "Propagators" -> {-k1^2 + msq, -(k1 - p1)^2 + msq} /. numericRules,
    "LoopMomenta" -> {k1},
    "ExternalMomenta" -> {p1},
    "KinematicRules" -> {p1^2 -> s} /. numericRules,
    "TopSector" -> {1, 1},
    "Numeric" -> numericRules,
    "Modulus" -> modulus
|>;

data = LIEDefineFamily[config, "Verbose" -> False];
data = LIESolveRegions[data, "Verbose" -> False];
data = LIEExpandSeries[data, "Order" -> 4, "Increment" -> 2, "Verbose" -> False];
data = LIEGetRelations[data, "Verbose" -> False, "MaxCoefDeg" -> 2];

relations = data["Relations", "Relations"];
Print["  关系数: ", Length /@ relations // Flatten // Total];

(* -------------------------------------------------- *)
(* Step 2: 加载 Kira 规则 *)
(* -------------------------------------------------- *)
Print["\n[Step 2] 加载 Kira 规则..."];

kiraRulePath = "/root/kira_tests/bub00_new/bub00";
kiraRulesFile = KiraFindRulesFile[kiraRulePath];
Print["  Kira rules file: ", kiraRulesFile];

If[kiraRulesFile === $Failed || ! FileExistsQ[kiraRulesFile],
    Print["  [WARN] Kira 规则文件不存在，跳过 KiraVerify"];
    Print["  [WARN] 运行: KiraInputGenerate[...] + kira jobs.yaml"];
,
    Print["  Loading Kira rules..."];
    kiraRaw = Import[kiraRulesFile, "Text"];
    
    (* 解析 Kira 规则 *)
    kiraRules = KiraParseRules[kiraRaw];
    Print["  Kira 规则数: ", Length[kiraRules]];
];

(* -------------------------------------------------- *)
(* Step 3: g → j 转换验证（关键诊断点）*)
(* -------------------------------------------------- *)
Print["\n[Step 3] g→j 转换验证..."];

(* 取第一个关系（lev=0, deg=1）的第一个 g-term *)
sampleRel = relations[[1, 1]];  (* lev=0, deg=1 的第一个关系 *)
Print["  示例关系: ", InputForm[sampleRel] // Short];

(* 提取 g[...] 项 *)
gTerms = Cases[sampleRel, g[___], Infinity];
Print["  g[...] 项数: ", Length[gTerms]];
Print["  g[ ] 索引示例: ", gTerms[[1;;Min[3, Length[gTerms]]]]];

(* 验证: g[α,β] 代入 ν 后是否产生 j[α'] 格式 *)
(* Kira 期望 j[{a, b}] 格式 — 检查转换是否匹配 *)
testConversions = Table[
    gTerm = gTerms[[i]];
    alpha = Last[gTerm];
    kiraIdx = alpha;  (* 假设 j[{alpha}] 格式 *)
    {InputForm[gTerm], j[{kiraIdx}]} // InputForm,
    {i, 1, Min[5, Length[gTerms]]}
];
Print["  g→j 转换示例: "];
Scan[Print["    ", #]&, testConversions];

(* -------------------------------------------------- *)
(* Step 4: Kira 约化验证（逐 term）*)
(* -------------------------------------------------- *)
Print["\n[Step 4] Kira 约化验证..."];

If[kiraRulesFile === $Failed || ! FileExistsQ[kiraRulesFile],
    Print["  [SKIP] Kira 规则不可用"];
,
    (* 对 sampleRel 做 Kira 约化 *)
    jExpr = sampleRel /. {g[a_, b_] :> j[{a, b}]};
    Print["  j-form: ", InputForm[jExpr] // Short];
    
    reduced = jExpr /. kiraRules;
    Print["  Kira 约化后: ", InputForm[reduced] // Short];
    
    isZero = (reduced === 0);
    Print["  结果: ", If[isZero, "[PASS] 约化为 0", "[FAIL] residual = "<>ToString[InputForm[reduced]]]];
];

(* -------------------------------------------------- *)
(* Step 5: 零指标规则覆盖检查 *)
(* -------------------------------------------------- *)
Print["\n[Step 5] 零指标规则覆盖检查..."];

allJIndices = Union @ Flatten @ Cases[sampleRel, j[idx_List] :> idx, Infinity];
Print["  j[...] 索引种类: ", Length[allJIndices], " 个"];
Print["  示例索引: ", allJIndices[[1;;Min[5, Length[allJIndices]]]]];

(* 找出哪些索引不被 Kira 规则覆盖 *)
coveredByKira = Union @ Flatten @ Cases[kiraRules, j[idx_] :> idx, Infinity];
uncoveredIndices = Complement[allJIndices, coveredByKira];
Print["  Kira 覆盖: ", Length[coveredByKira], " 个"];
Print["  未覆盖索引: ", Length[uncoveredIndices], " 个"];
If[Length[uncoveredIndices] > 0,
    Print["  未覆盖索引示例: ", uncoveredIndices[[1;;Min[5, Length[uncoveredIndices]]]]]
];

(* 零指标规则: 所有 ≤0 索引应直接为 0 *)
zeroIndices = Select[allJIndices, Max[#] <= 0 &];
positiveIndices = Select[allJIndices, Max[#] > 0 &];
Print["  零指标 j[≤0]: ", Length[zeroIndices], " 个"];
Print["  正指标 j[>0]: ", Length[positiveIndices], " 个"];

(* -------------------------------------------------- *)
(* Step 6: 输出诊断结论 *)
(* -------------------------------------------------- *)
Print["\n========================================"];
Print[" 诊断结论"];
Print["========================================"];

If[Length[uncoveredIndices] > 0,
    Print["  [ROOT CAUSE] Kira 规则不完整"];
    Print["  未覆盖索引数: ", Length[uncoveredIndices]];
    Print["  建议: 增加 Kira 的 RMax/SMax 以覆盖所有索引"];
,
    Print["  [OK] Kira 规则覆盖了所有出现的索引"]
];

If[Length[positiveIndices] > 0,
    Print["  [NOTE] 正指标 j[>0] 无法被零规则归零"];
    Print["  这些项必须被 Kira 规则约化"];
    Print["  正指标数: ", Length[positiveIndices]];
,
    Print["  [OK] 所有 j[>0] 都有 Kira 规则覆盖"]
];

Print["\n  下一步: "];
Print["  1. 检查 Kira RMax/SMax 是否足够大"];
Print["  2. 确认 g→j 转换公式正确（索引偏移量）"];
Print["  3. 在 d=1/13 配置下运行相同测试作为参考"];
Print["========================================"];
