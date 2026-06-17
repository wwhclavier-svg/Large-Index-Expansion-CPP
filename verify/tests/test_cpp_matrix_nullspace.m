(* test_cpp_matrix_nullspace.m *)
(* 直接测试: C++ 零空间向量是否满足 M(ν)·x = 0 *)
(* 使用 MMA 的 LIEWorkflow 来构建 M(ν) *)

SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];
Get["LIEWorkflow.wl"];

Print["========================================"];
Print["测试: C++ 零空间是否满足 M(ν)·x = 0"];
Print["========================================"];

(* ========================================================================
 * 使用 LIEWorkflow 构建 IBP 系数矩阵 M(ν)
 * ======================================================================== *)

bubbleConfig = <|
  "Propagators" -> {-k1^2 + msq, -(k1 - p1)^2 + msq},
  "LoopMomenta" -> {k1},
  "ExternalMomenta" -> {p1},
  "KinematicRules" -> {p1^2 -> s},
  "TopSector" -> {1, 1},
  "Numeric" -> {s -> 3, msq -> 1, "d" -> 1/3}
|>;

(* Step 1: 定义家族 *)
family = LIEDefineFamily[
  bubbleConfig["Propagators"],
  bubbleConfig["LoopMomenta"],
  bubbleConfig["ExternalMomenta"],
  bubbleConfig["KinematicRules"],
  bubbleConfig["TopSector"],
  "Numeric" -> bubbleConfig["Numeric"],
  Verbose -> False
];

(* Step 2: 求解区域 *)
regions = LIESolveRegions[family, Verbose -> False];

(* Step 3: 展开级数 *)
expanded = LIEExpandSeries[regions, "Order" -> 4, Verbose -> False];

(* 获取系数矩阵 C(α, k, l, γ, j, i) 和 p(α) 矩阵 *)
expResult = expanded;

Print["\n========================================"];
Print["Step 1: 读取 C++ 零空间向量"];
Print["========================================"];

(* 读取 C++ 零空间向量 *)
cppData = Get["/root/Large-Index-Expansion-CPP/verify/bub00/Compare-CPPRelation-bub00_lev2_deg1.m"];
alphas = cppData["Alphas"];
betas = cppData["Betas"];
coeffs = cppData["Coefficients"];
numAlpha = Length[alphas];
numBeta = Length[betas];
numSols = cppData["NumSolutions"];

Print["Alphas: ", alphas];
Print["Betas: ", betas];
Print["零空间向量数: ", numSols];

(* ========================================================================
 * 关键测试: 使用 MMA 的 LIEWorkflow 来计算 M(ν) 并验证零空间
 * ======================================================================== *)

Print["\n========================================"];
Print["Step 2: 计算 M(ν) 矩阵"];
Print["========================================"];

(* 获取 top sector 的 recursion matrix 来构建 M(ν) *)
topSec = {1, 1};
recMat = expResult["Regions"][topSec][[1]]["RecursionMatrix"];

Print["RecursionMatrix 维度: ", Dimensions[recMat]];

(* 根据论文，M(ν) 的元素由下式给出： *)
(* M_{α,β}(ν) = sum_{γ} p_{α}(ν) * C_{α}(k,l,γ,j,i) * (ν)^γ *)
(* 这里 α 是 sector 索引，β 是 power 索引 *)

(* 实际上，对于 Level=2, Deg=1 的关系，我们需要的是： *)
(* 行: 不同的 (alpha, power) 组合 *)
(* 列: g[ν - alpha] 项 *)

(* 但 LIEWorkflow 的内部格式不直接是 M(ν) * x = 0 的形式 *)
(* LIE 的关系形式是: sum_{α,β} c_{α,β}(ν) * ν^β * g_{ν-α} = 0 *)

(* ========================================================================
 * 替代方法: 直接用 RegimeEvaluator 的逻辑手动计算
 * ======================================================================== *)

Print["\n========================================"];
Print["Step 3: 直接计算方程的系数"];
Print["========================================"];

(* 获取 regime 数据 *)
regData = expResult["Regions"][topSec][[1]];
Ccoeff = regData["C Coefficient"];
pMat = regData["p Matrix"];

Print["C Coefficient 维度: ", Dimensions[Ccoeff]];
Print["p Matrix 维度: ", Dimensions[pMat]];

(* 系数矩阵结构: Ccoeff[k, l, cid, j, i] *)
(* pMat[i, j] 是 p_ij(ν) *)

(* 根据 RegimeEvaluator::step2_computeG 的逻辑: *)
(* g_i(α, k) = sum_{l, cid, j} p_ij * C(k,l,cid,j,i) * (ν - α)^γ *)

(* 但实际上，对于零空间验证，我们关心的是： *)
(* 给定 ν 点，方程矩阵 M(ν) 的行是什么？ *)

(* ========================================================================
 * 简化测试: 直接用数值 ν 点验证
 * ======================================================================== *)

Print["\n========================================"];
Print["Step 4: 数值验证"];
Print["========================================"];

(* 测试点 *)
nuTest = {2, 3};

(* 读取 C++ 方程在 nu=(2,3) 时的系数 *)
(* 根据 buildFinalMatrix, 矩阵的行对应不同的 (i, power) 组合 *)

(* 但我们需要知道 C++ 的系数矩阵 M(ν) 实际是什么 *)
(* 让我检查 C++ 导出文件中的系数 *)

(* C++ 导出的 Coefficients 是 Mext: 每列是一个零空间向量 *)
(* Mext[i][s] 是第 s 个零空间向量的第 i 个元素 *)

(* 为了验证 M(ν)·x = 0，我们需要 M(ν) 本身 *)
(* 但 test_relationFF 从不输出 M(ν)，只输出零空间向量 *)

(* 让我检查是否有其他方式获取 M(ν) *)

(* 实际上，我需要重新理解问题:
 * C++ 输出的零空间向量 x 满足 M(ν)·x = 0
 * 其中 M(ν) 是 IBP 系数矩阵
 * 
 * Kira 规则告诉我们 j[α,β] 如何约化
 * 但 j[α,β] 不是 M(ν) 的元素
 * 
 * 关键问题: g[ν - α] 和 j[α,β] 的关系是什么?
 *)

(* 检查 integral family 的定义 *)
Print["\n========================================"];
Print["检查 IntegralFamily 的 j 定义"];
Print["========================================"];

(* IntegralFamily 中 j[α,β] 的含义是 sector α, power β 的 IBP 积分 *)
(* 而 g[ν - α] 是展开系数，不是直接的 IBP 积分 *)

(* 让我查看 LIEWorkflow 中是否有 computeAndExportRelations 函数 *)
(* 或者检查 MMA 版本的 RelationSolver 是如何工作的 *)

Print["\nLIEWorkflow 函数列表:"];
Print /@ Names["LIE*"];

(* 检查是否有 computeMMatrix 或类似函数 *)
If[NameQ["LIEComputeIBPMatrix"],
    Print["\n找到 LIEComputeIBPMatrix"],
    Print["\n未找到 LIEComputeIBPMatrix"]
];
