(* ============================================================ *)
(* GenerateBinFiles.wl — 从 FamilyDatabase 批量生成 .bin 文件       *)
(* ============================================================ *)
(* 用法:                                                           *)
(*   在 Mathematica 中:                                             *)
(*     Get["scripts/GenerateBinFiles.wl"]                          *)
(*   或命令行:                                                      *)
(*     math -script scripts/GenerateBinFiles.wl                    *)
(* ============================================================ *)

(* --- 0. 路径解析 ----------------------------------------------------- *)
(* projectDir: C++ 项目根目录 *)
Which[
    (* 模式1: math -script — $InputFileName 是脚本的绝对路径 *)
    StringQ[$InputFileName] && $InputFileName =!= "",
        projectDir = ParentDirectory[DirectoryName[$InputFileName]],
    (* 模式2: Notebook 交互 — 使用当前目录 *)
    True,
        projectDir = Directory[]
];

(* mmaDir: Mathematica 包目录 (固定路径) *)
mmaDir = FileNameJoin[{$HomeDirectory, "文档",
    "Wolfram Mathematica", "[Work] Program_LargeIndexExpansion_2508"}];
If[!DirectoryQ[mmaDir],
    Print["ERROR: MMA package directory not found: ", mmaDir];
    Quit[1]
];

Print["Project directory: ", projectDir];
Print["MMA package directory: ", mmaDir];

(* --- 1. 加载依赖包 ------------------------------------------------- *)
Get[FileNameJoin[{mmaDir, "SingularInterface.wl"}]];
Get[FileNameJoin[{mmaDir, "QuotientAlgebra.wl"}]];
Get[FileNameJoin[{mmaDir, "LargeIndexExpansion.wl"}]];
Get[FileNameJoin[{mmaDir, "BlockRecursion.wl"}]];
Get[FileNameJoin[{mmaDir, "SymbolicBTForm.wl"}]];
Get[FileNameJoin[{mmaDir, "ExportBinary_IBPMatrix.wl"}]];

(* --- 2. 加载 FamilyDatabase --------------------------------------- *)
familyDbPath = FileNameJoin[{projectDir, "verify", "FamilyDatabase", "FamilyDatabase.wl"}];
If[!FileExistsQ[familyDbPath],
    Print["ERROR: FamilyDatabase.wl not found at: ", familyDbPath];
    Quit[1]
];
Get[familyDbPath];

(* --- 3. 目标族列表 -------------------------------------------------- *)
(* 一圈族 + SR *)
targetFamilies = {"bub00", "bub10", "bub11", "Tri", "Box", "SR"};

(* --- 4. 主循环：逐族生成 .bin 文件 -------------------------------- *)
SetDirectory[projectDir];

Do[
    Print["\n========================================"];
    Print["Processing family: ", fam];
    Print["========================================"];

    (* 获取族配置 *)
    If[!KeyExistsQ[$FamilyDatabase, fam],
       Print["WARNING: Family '", fam, "' not found in database, skipping."];
       Continue[]
    ];
    config = $FamilyDatabase[fam];
    Print["  Description: ", config["Description"]];

    (* 提取参数 *)
    pdlist    = config["Propagators"];
    loopmom   = config["LoopMomenta"];
    extmom    = config["ExternalMomenta"];
    spsRep    = config["KinematicRules"];
    topsector = config["TopSector"];
    char      = config["Modulus"];

    (* 分离 Numeric: "d" -> numericD, 其余 -> numeric *)
    numericRules = config["Numeric"];
    numericD  = Cases[numericRules, ("d" -> _)];
    numeric   = DeleteCases[numericRules, ("d" -> _)];

    Print["  Propagators: ", pdlist];
    Print["  TopSector: ", topsector];
    Print["  Numeric: ", numeric];
    Print["  NumericD: ", numericD];
    Print["  Modulus: ", char];

    (* Step A: 定义 IBP 族 *)
    {ne, nl, nE, nibp, Alist, vlist, ibpeqs1, sectorlist} =
        AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topsector];
    Print["  ne=", ne, " nl=", nl, " nE=", nE, " nibp=", nibp];

    (* Step B: 生成展开区域数据 (有限域) *)
    ibpeqsSub = ibpeqs1 /. numeric /. numericD;
    expregdataFF = regionsBySectors[
        ibpeqsSub, sectorlist, Alist, vlist,
        Modulus -> char
    ];

    (* Step C: 导出 .bin 文件 *)
    ibpFile  = projectDir <> "/data/IBPMat_"  <> fam <> ".bin";
    ringFile = projectDir <> "/data/RingData_" <> fam <> ".bin";

    Print["  Exporting: ", ibpFile];
    ExportBinaryIBPMatrix[ibpFile, expregdataFF, char];

    Print["  Exporting: ", ringFile];
    ExportBinaryRingData[ringFile, expregdataFF, Alist, ne, char];

    Print["  Done: ", fam];

, {fam, targetFamilies}];

(* --- 5. 输出摘要 ---------------------------------------------------- *)
Print["\n========================================"];
Print["All families processed. Output directory: ", projectDir];
Print["========================================"];
Print["Generated files in data/:"];
Do[
    fam = targetFamilies[[i]];
    ibpFile  = projectDir <> "/data/IBPMat_"  <> fam <> ".bin";
    ringFile = projectDir <> "/data/RingData_" <> fam <> ".bin";
    If[FileExistsQ[ibpFile],
       Print["  ", ibpFile, " (", Floor[FileByteCount[ibpFile]/1024], " KB)"],
       Print["  ", ibpFile, " — MISSING!"]];
    If[FileExistsQ[ringFile],
       Print["  ", ringFile, " (", Floor[FileByteCount[ringFile]/1024], " KB)"],
       Print["  ", ringFile, " — MISSING!"]];
, {i, Length[targetFamilies]}];
