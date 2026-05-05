(* ::Package:: *)
(* Compare-Reconstruct-bub00.wl *)
(* 关系重构验证工作流 *)
(* Step 1: Family 定义 -> Step 2: 区域求解 -> Step 3: 展开 -> Step 4: 关系重建 -> Step 5: Kira验证 *)

SetDirectory[DirectoryName[$InputFileName]];
Print["Working directory: ", Directory[]];

Get["LIEWorkflow.wl"];
Get["LIEReconstruct.wl"];

(* ========================================== *)
(* 配置: 必须与 Compare-FamilyGenerate-bub00.wl 一致 *)
(* ========================================== *)

(* ---- Kira 路径配置 ---- *)
$KiraBaseDir = If[ValueQ[$KiraBaseDir], $KiraBaseDir, "/root/kira_tests/bub00_new/bub00"];

modulus = Prime[10000000];  (* 179424673 *)
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

Print["========================================"];
Print["Compare-Reconstruct-bub00"];
Print["  s = 3, msq = 0, d = 1/3"];
Print["  Modulus = ", modulus];
Print["========================================"];

(* ========================================== *)
(* Step 1: Define Family *)
(* ========================================== *)
Print["\n[Step 1] Defining IBP family..."];
data = LIEDefineFamily[
    bub00Config["Propagators"],
    bub00Config["LoopMomenta"],
    bub00Config["ExternalMomenta"],
    bub00Config["KinematicRules"],
    bub00Config["TopSector"],
    "Numeric" -> bub00Config["Numeric"],
    Modulus -> modulus,
    Verbose -> False
];
Print["  NE = ", data["Family", "NE"]];
Print["  NIBP = ", data["Family", "NIBP"]];

(* ========================================== *)
(* Step 2: Solve Regions *)
(* ========================================== *)
Print["\n[Step 2] Solving regions..."];
data = LIESolveRegions[data, Verbose -> False];
Print["  Regions solved."];

(* ========================================== *)
(* Step 3: Expand Series *)
(* ========================================== *)
Print["\n[Step 3] Expanding series (order=4, LayerByLayer=True)..."];
data = LIEExpandSeries[data,
    "Order" -> 4,
    Modulus -> modulus,
    "Increment" -> 2,
    "LayerByLayer" -> True,
    Verbose -> False
];
Print["  Expansion completed."];

(* 提取数据 *)
hlist = data["Expansion", "HList"];
alist = data["Regions"];
topsec = bub00Config["TopSector"];
ne = data["Family", "NE"];  (* 从 family 数据获取 ne *)
vlist = Table[ToExpression["v" <> ToString[i]], {i, ne}];  (* 生成 v1, v2, ... *)

(* 确保 v1, v2 作为全局符号存在 *)
If[ne >= 1, v1 = v1];
If[ne >= 2, v2 = v2];

Print["  HList shape: ", Dimensions[hlist]];
Print["  Regions: ", Length[alist]];
Print["  VList: ", vlist];
Print["  NE: ", ne];

(* ========================================== *)
(* Step 4: Relation Reconstruction via LIEGetRelations *)
(* ========================================== *)
Print["\n[Step 4] Running LIEGetRelations..."];

rank = 2;       (* max level |alpha| *)
maxDeg = 2;     (* max degree |beta| *)

(* 使用 LIEGetRelations 而不是直接调用 LIEReconstruct *)
data = LIEGetRelations[data, Verbose -> False, "MaxCoefDeg" -> maxDeg];
Print["  Status: ", data["Status"]];

relations = data["Relations", "Relations"];
Print["  Relations computed."];
Print["  Output dimensions: ", Dimensions[relations]];

(* ========================================== *)
(* Step 5: Kira 验证关系 *)
(* ============================================================ *)
(* Step 5: Kira 验证关系（使用 KiraRuleLoader）*)
(* ============================================================ *)
Print["\n[Step 5] Verifying relations using Kira..."];

(* 使用 KiraRuleLoader 加载 Kira 规则 *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    $KiraBaseDir,
    "d" -> 1/3,
    "Modulus" -> modulus,
    "NProp" -> 2
];

(* 从 LIE 关系中提取 g 系数并转换为 Kira j 格式 *)
extractGCoefficientsFromRelation[expr_] := Module[
    {terms, result, termList, coef, gfound, alpha, beta, numG},
    result = {};
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    Do[
        term = terms[[i]];
        termList = List @@ term;
        numG = Count[termList, _?((ToString[Head[#]] === "g")&)];
        Which[
            Head[term] === Times && numG > 0,
            gfound = Cases[termList, _?((ToString[Head[#]] === "g")&)][[1]];
            coef = Coefficient[term, gfound];
            alpha = gfound[[1]] - v1;
            beta = gfound[[2]] - v2;
            AppendTo[result, {coef, alpha, beta}],
            Head[term] === Times,
            Do[
                If[ToString[Head[termList[[j]]]] === "g",
                    gterm = termList[[j]];
                    coef = Times @@ DeleteCases[termList, gterm];
                    alpha = gterm[[1]] - v1;
                    beta = gterm[[2]] - v2;
                    AppendTo[result, {coef, alpha, beta}]
                ]
            , {j, Length[termList]}],
            ToString[Head[term]] === "g",
            alpha = term[[1]] - v1;
            beta = term[[2]] - v2;
            AppendTo[result, {1, alpha, beta}]
        ]
    , {i, Length[terms]}];
    result
];

Print["\n[Step 5] Verifying relations using Kira..."];

Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    $KiraBaseDir,
    "d" -> 1/3,
    "Modulus" -> modulus,
    "NProp" -> 2
];

(* 从 LIE 关系中提取 g 系数并转换为 Kira j 格式 *)
extractGCoefficientsFromRelation[expr_] := Module[
    {terms, result, termList, coef, gfound, alpha, beta, numG},
    result = {};
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    Do[
        term = terms[[i]];
        termList = List @@ term;
        numG = Count[termList, _?((ToString[Head[#]] === "g")&)];
        Which[
            Head[term] === Times && numG > 0,
            gfound = Cases[termList, _?((ToString[Head[#]] === "g")&)][[1]];
            coef = Coefficient[term, gfound];
            alpha = gfound[[1]] - v1;
            beta = gfound[[2]] - v2;
            AppendTo[result, {coef, alpha, beta}],
            Head[term] === Times,
            Do[
                If[ToString[Head[termList[[j]]]] === "g",
                    gterm = termList[[j]];
                    coef = Times @@ DeleteCases[termList, gterm];
                    alpha = gterm[[1]] - v1;
                    beta = gterm[[2]] - v2;
                    AppendTo[result, {coef, alpha, beta}]
                ]
            , {j, Length[termList]}],
            ToString[Head[term]] === "g",
            alpha = term[[1]] - v1;
            beta = term[[2]] - v2;
            AppendTo[result, {1, alpha, beta}]
        ]
    , {i, Length[terms]}];
    result
];

(* 使用 kiraReduce 验证单个 LIE 关系 *)
verifyLIERelation[lieRel_] := Module[
    {gCoeffs, jExpr, kiraReduced, isPass, result},
    gCoeffs = extractGCoefficientsFromRelation[lieRel];
    If[gCoeffs === {}, Return[{True, "No g terms"}]];
    jExpr = Plus @@ (#1 * j[#2, #3] & @@@ gCoeffs);
    kiraReduced = kiraReduce[jExpr];
    isPass = (kiraReduced === 0);
    {isPass, If[!isPass, InputForm[kiraReduced] // Short, "OK"]}
];

(* 对所有 LIE 关系进行验证 *)
Print["  Verifying LIE relations using Kira rules..."];

passCount = 0;
failCount = 0;
totalCount = 0;

Do[
    numRel = Length[relations[[lev+1, deg+1]]];
    If[numRel > 0,
        Print["  Checking (lev=", lev, ", deg=", deg, "): ", numRel, " relations"];
        Do[
            lieRel = relations[[lev+1, deg+1]][[r]];
            totalCount++;
            {isPass, result} = verifyLIERelation[lieRel];
            If[isPass,
                passCount++,
                failCount++;
                Print["    [FAIL] Relation ", r, ": ", result]
            ]
        , {r, numRel}]
    ]
, {lev, 0, rank}, {deg, 0, maxDeg}];

Print[""];
Print["========================================"];
Print["Verification Summary:"];
Print["  Total relations checked: ", totalCount];
Print["  Passed: ", passCount];
Print["  Failed: ", failCount];
Print["========================================"];

(* 退出码 *)
If[failCount > 0, Exit[1]];

(* 打印一些示例关系 *)
Print["\n[Step 6] Sample LIE relations..."];
Do[
    numRel = Length[relations[[lev+1, deg+1]]];
    If[numRel > 0,
        Print["  Relations at (lev=", lev, ", deg=", deg, "): ", numRel, " relations"];
        firstRel = relations[[lev+1, deg+1]][[1]];
        gCoeffs = extractGCoefficientsFromRelation[firstRel];
        Print["    First relation has ", Length[gCoeffs], " g-terms"];
        Print["    First relation: ", InputForm[firstRel]];
    ];
, {lev, 0, rank}, {deg, 0, maxDeg}];

Print["\n========================================"];
Print["Compare-Reconstruct-bub00 complete!"];
Print["========================================"];
Print[""];
Print["工作流数据: /tmp/MMA_bub00_reconstruct_workflow.wdx"];
Print[""];
Print["下一步: 运行 C++ test_relationFF 生成 Compare-CPPRelation-bub00.m"];
Print["  cd ../../build"];
Print["  ./test_relationFF bub00 4 2 2"];
Print["========================================"];