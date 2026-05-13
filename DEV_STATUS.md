# DEV_STATUS — Large Index Expansion C++ 项目工程状态

> 更新于 2026-05-10，覆盖 2026-05-05 9:00 至今的改动。

---

## 目录

1. [主要管线：C++ Region 计算管线 + Singular 接入](#1-主要管线c-region-计算管线--singular-接入)
2. [验证测试：验证工作流的整理](#2-验证测试验证工作流的整理)
3. [性能建模：Singular 求解 Region 的时间消耗](#3-性能建模singular-求解-region-的时间消耗)
4. [Ansatz 设计与积分约化：AnsatzMode](#4-ansatz-设计与积分约化ansatzmode)
5. [其他工作区的整理](#5-其他工作区的整理)
6. [效率优化：Strategy 设计](#6-效率优化strategy-设计)
7. [约化规则产生：SymbolicRule 管线设计](#7-约化规则产生symbolicrule-管线设计)

---

## 1. 主要管线：C++ Region 计算管线 + Singular 接入

### 1.1 架构概览

全管线四阶段：

```
families/*.json → [B] IBPEqGenerator → [C] RegionSolver → .bin 输出 (IBPMat + RingData)
                  (Singular 子进程)     (Singular 子进程)
                                                           ↓
                                              [D] LayerRecursion → resCache
                                                           ↓
                                              [E] RelationSolver → AllRelation_*.m
                                                           ↓
                                              [F] SymbolicRule (planned) → 完备约化规则
```

Stage A–C 为 header-only C++17，通过 `SingularRunner` 调用 Singular 子进程。管线总代码量约 8500 行（不含 LayerRecursion / RelationSolver / SymbolicRule）。

### 1.2 模块清单

| 模块 | 文件 | 行数 | 功能 |
|------|------|------|------|
| FamilyConfig | `include/FamilyConfig.hpp` | 72 | JSON 解析 → FamilyDef |
| IBPEqGenerator | `include/IBPEqGenerator.hpp` | 1371 | SP2PD → 导数 → IBP 组装 (mon2F) |
| IBPAnalyzer | `include/IBPAnalyzer.hpp` | 368 | A/B 方程构建 + FTable 提取 |
| RegionSolver | `include/RegionSolver.hpp` | 792 | Groebner → 准素分解 → 单项式基底 |
| RecursionBuilder | `include/RecursionBuilder.hpp` | 511 | FractionRule → 递归矩阵 |
| RingBuilder | `include/RingBuilder.hpp` | 193 | A/Ainv 矩阵计算 |
| PolyArith | `include/PolyArith.hpp` | 437 | 有限域多项式运算 |
| BinaryIBPWriter | `include/BinaryIBPWriter.hpp` | 264 | IBPMat_*.bin 写出 |
| BinaryRingWriter | `include/BinaryRingWriter.hpp` | 191 | RingData_*.bin 写出 |
| SingularRunner | `include/SingularRunner.hpp` | 499 | Singular 子进程通用执行器 |
| CLI | `tools/family_generate.cpp` | 421 | 主入口，支持 --diff / --output |

### 1.3 算法流程

1. **SP2PD**：通过 Singular 计算 `C^{-1}`，将标量积映射到 z-基底
2. **导数计算**：`∂D_i/∂l_j · q_k` 用上一步的变换矩阵代入 z-多项式
3. **IBP 组装**：`d·δ_jk·∏z_m + Σ_i n_i·(-deriv[i,j,k]·∏_{m≠i}z_m)` → mon2F → g-算子
4. **A/B 方程**：`buildABEquations` 将 g-算子转为 A/B 多项式
5. **Region 求解**：Singular 计算 Groebner 基 → `minAssPrimes` 准素分解 → 提取单项式基底 nb
6. **递归矩阵**：`RecursionBuilder` 从 FractionRule 构建 M1/N1/K1/F0/F2 等
7. **Ring 矩阵**：`RingBuilder` 计算 A/Ainv
8. **写出**：IBP1 格式的 IBPMat_*.bin + RingData_*.bin

### 1.4 已完成 Family（C++ 管线）

| Family | L | ne | Regions | JSON | 状态 |
|--------|---|----|---------|------|------|
| bub00 | 1 | 2 | 1 | ✅ | ✅ |
| bub10 | 1 | 2 | 3 | ✅ | ✅ |
| bub11 | 1 | 2 | 5 | ✅ | ✅ |
| Tri | 1 | 3 | 2 | ✅ | ✅ |
| Box | 1 | 4 | 15 | ✅ | ✅ |
| SR212 | 2 | 5 | 6 | ✅ | ✅ |
| SR212-3m | 2 | 5 | 8 | ✅ | ✅ |
| SR212-5m | 2 | 5 | 45 | ✅ | ✅ |
| TB123 | 2 | 7 | 37 | ✅ | ⏳ Singular 超时 |
| TB123m | 2 | 7 | — | ❌ | 无 JSON |
| DB313 | 2 | 7 | — | ❌ | 无 JSON |
| NP222 | 2 | 7 | — | ❌ | 无 JSON |
| NP322 | 2 | 8 | — | ❌ | 无 JSON |

### 1.5 关键问题：A/B 方程符号不一致

C++ 和 MMA 的 IBP 组装采用不同中间约定：
- **C++**：乘 ∏z_m 后取负 `-deriv[i,j,k]`
- **MMA**：LargeIndexIBP（n_i→n+v_i）→ Coefficient["n"] 提取

两套约定数学等价，但 A/B 方程部分项符号不同 → Singular 选出等价但不同的单项式基底 → .bin 逐字节不同。功能验证通过（膨胀系数逐字节一致，sol_dim 匹配），但 .bin 文件不与 MMA 参考文件逐字节一致。

### 1.6 待完成

- [ ] Lorentz Invariance (LI) 方程生成
- [ ] A/B 方程符号对齐（达成 .bin 逐字节一致）
- [ ] TB123 超时问题（需优化 Groebner 基计算）
- [ ] DB313 / NP222 / NP322 / NP322m / TB123m 的 JSON config

---

## 2. 验证测试：验证工作流的整理

### 2.1 MMA 参考实现（VerifyUtility）

`commit 33b3e14` — 将 MMA 参考管线从原工作区迁入版本库。经后续清理，当前共 **20 个 .wl 脚本**（约 6500 行），按功能分组：

**核心管线（LIE 系列，7 个）：**

| 文件 | 行数 | 功能 |
|------|------|------|
| `LIEFamilyDefine.wl` | 119 | Family 定义 → 内部数据 |
| `LIEUtility.wl` | 147 | sectorLimitIBP、导数、系数替换 |
| `LIECoreAlgebra.wl` | 295 | expRegSolve2、Groebner/准素分解 |
| `LIERegions.wl` | 231 | regionsBySectors 主入口（含 sector 级计时） |
| `LIEExpand.wl` | 519 | 膨胀系数计算 |
| `LIEReconstruct.wl` | 1022 | 关系重构 |
| `LIEWorkflow.wl` | 570 | 完整管线编排 |

**Singular 接口与导出（3 个）：**

| 文件 | 功能 |
|------|------|
| `SingularInterface.wl` | Singular 子进程调用（新增 `SingularTimeout` 超时选项） |
| `ExportBinary_IBPMatrix.wl` | .bin 导出 |
| `M2Kira.wl` | Kira 格式互转 |

**膨胀系数验证（VerifyExpand 系列，4 个）：**

| 文件 | 功能 |
|------|------|
| `VerifyExpand-Prepare.wl` | C++ .bin → MMA 加载准备 |
| `VerifyExpand-Compare.wl` | C++ vs MMA 膨胀系数逐项对比 |
| `VerifyExpand-MMAExpand.wl` | MMA 独立膨胀计算 |
| `VerifyExpand-SeriesVerify.wl` | 级数自洽性验证 |

**Kira 交叉验证与通用工具（6 个）：**

| 文件 | 功能 |
|------|------|
| `KiraRuleLoader.wl` | Kira IBP 规则导入 |
| `Verify-Blade.wl` | Blade 张量验证 |
| `Verify-Kira.wl` | Kira 约化对比验证 |
| `Verify-MMACompare.wl` | C++ vs MMA 通用对比框架 |
| `Verify-Series.wl` | 级数验证工具 |
| `Debug-DumpFamilyGenerate-v2.wl` | 中间步骤调试 dump（IBP 方程 / FTable / A/B 方程） |
| `Compare-VerifyLog.wl` | 验证日志对比 |

> **已删除**（commit 99ec194）：`Compare-Expand.wl`、`Compare-FamilyGenerate.wl`、`Compare-Results.wl`、`Compare-Reconstruct-bub00.wl`、`VerifyRelation-SeriesVerify.wl`。功能已收敛到上述统一验证框架。

### 2.2 Kira 验证缓存

`VerifyUtility/cache/` 目录包含 bub00、SR212、Tri 的 Kira 约化缓存数据（sector mappings + IBP rule ID），用于 Kira 交叉验证流程：

```
VerifyUtility/cache/
├── bub00_reduce_1/       # bub00: ibps/nids/zidpos 映射 + topology
├── bub00_sectormappings/  # bub00: maximal cut masters
├── SR212_reduce_1/        # SR212: 同上 + Map/SR/SubSym 规则
├── SR212_sectormappings/
├── Tri_reduce_1/          # Tri: 同上 + LI 规则
└── Tri_sectormappings/
```

### 2.3 MMA 参考 .bin 文件与元数据

`verify/<family>/` 目录下 MMA 生成的参考 .bin 和时序/区域信息。当前覆盖 **12 个 family**：

| Family | L | ne | .bin | Timing | RegionInfo | Checkpoint |
|--------|---|---|------|--------|------------|------------|
| bub00 | 1 | 2 | ✅ | ✅ | ✅ | ✅ |
| bub10 | 1 | 2 | ✅ | ✅ | ✅ | ✅ |
| bub11 | 1 | 2 | ✅ | ✅ | ✅ | ✅ |
| Tri | 1 | 3 | ✅ | ✅ | ✅ | ✅ |
| Box | 1 | 4 | ✅ | ✅ | ✅ | ✅ |
| SR212 | 2 | 5 | ✅ | ✅ | ✅ | ✅ |
| SR212-3m | 2 | 5 | ✅ | ✅ | ✅ | ✅ |
| SR212-5m | 2 | 5 | ✅ | ✅ | ✅ | ✅ |
| TB123 | 2 | 7 | ✅ | ✅ | ✅ | ✅ |
| TB123m | 2 | 7 | ✅ | ✅ | ✅ | ✅ |
| DB313 | 2 | 7 | ✅ | ✅ | ✅ | ✅ |
| NP222 | 2 | 7 | ✅ | ✅ | ✅ | ✅ |

> bub00 目录额外包含 `resCache_Expansion_bub00.bin`（膨胀系数缓存）和 `Compare-CPPResult-bub00.m`。

### 2.4 test_relationFF 验证

- **bub00**：膨胀系数逐字节一致，6 个 (lev,deg) 配置 sol_dim 完全匹配
- **Box / SR212 / SR212-3m / SR212-5m / TB123 / NP222**：性能基准完成，relation 重构 sol_dim 验证通过（见 [§3](#3-性能建模singular-求解-region-的时间消耗)）
- 其余 family 待覆盖

### 2.5 验证文档

`verify/docs/` 经大幅整理，保留 3 个核心文档：

| 文件 | 内容 |
|------|------|
| `Test-Expand.md` | 膨胀系数验证结果 |
| `Test-Relation.md` | 关系重构验证结果（大幅精简重写） |
| `CPP-KiraVerify-Debugger.md` | C++ / Kira 交叉调试指南 |

> 已删除：`IBPVerification.md`、`Verify-MMA-KIRA-Guide.md`（commit 99ec194）。

---

## 3. 性能建模：Singular 求解 Region 的时间消耗

### 3.1 测试结果

`docs/Benchmark_Results.md` (712+ 行) 记录了三类基准测试：

**Region 求解性能（C++ FamilyGenerate vs MMA）：**

| Family | ne | C++ 耗时 | MMA 耗时 | 加速比 |
|--------|----|---------|----------|--------|
| bub00 | 2 | 0.11s | 3.4s | 31x |
| Box | 4 | 0.62s | 4.5s | 7x |
| Tri | 3 | 0.25s | 3.7s | 15x |
| SR212 | 5 | 1.2s | — | — |

**test_relationFF 算法优化效果（Incremental RREF + nimax 削减，commit 8391228）：**

| Family | 优化前 | 优化后 | 加速比 |
|--------|--------|--------|--------|
| Box | 71.6s | 7.50s | 9.5x |
| SR212 | 53.6s | 9.85s | 5.4x |
| SR212-3m | 80.2s | 12.75s | 6.3x |
| SR212-5m | OOM | 69.77s | ∞ |
| TB123 | 38min+ | 300.52s | ~7.5x |

**编译优化影响（-O3 vs -O0）：** bub00 加速 2.2x~5.0x，bub11 加速 2.3x~6.0x。

### 3.2 关键发现

- **C++ FamilyGenerate 比 MMA 快 7-34x**（wall clock）：C++ 内联计算 IBP 组装，MMA 走符号计算 + 字符串解析
- **Singular Groebner 基是瓶颈**：TB123 (ne=7) 在 Groebner 基步骤超时
- **Expand-vs-Solve crossover 分析**已记录在 Benchmark_Results 中（含 ne=7 family 的 lev=3 全扇区结果）
- **MMA 端计时细化**（`LIERegions.wl`）：拆分为 `expRegSolve` 时间 + `buildRecursionMatrix` 时间 + 每 sector 总时间

### 3.3 Singular 超时控制

`SingularInterface.wl` 新增 `SingularTimeout` 选项，支持 `timeout --signal=KILL Ns Singular` 机制，防止大 family 的 Groebner 基计算无限挂起。

---

## 4. Ansatz 设计与积分约化：AnsatzMode

### 4.1 实现

`commit 8c98ba2` + `docs/AnsatzModes.md` (214 行)，在 `include/RelationSolver.hpp` (2755 行) 中实现 4 种 ansatz 模式：

| 模式 | 名称 | Kernel | α 生成 |
|------|------|--------|--------|
| 0 | Pyramid (默认) | `g(ν − α)` | level 内所有分量 |
| 1 | DotPyramid | `g(ν − α)` | 点积约束的 level 内分量 |
| 2 | Star | `g(ν + α)` | negated kernel，level 内所有分量 |
| 3 | ExtendedPyramid | `g(ν − α)` | 跨 sector 扩展分量 |

### 4.2 CLI 使用

```bash
./test_relationFF bub 4 1 2 2                  # 默认 Pyramid
./test_relationFF bub 4 1 2 2 --mode 1          # DotPyramid
./test_relationFF bub 4 0 2 2 --mode 3 --sector 110  # ExtendedPyramid
```

### 4.3 设计要点

- 所有模式的 α 列表按 MMA GenerateSeeds 排序（graded by L1 then colexicographic descending）
- RemoveSolvedVariables 跨模式通用
- 导出 .m 文件包含 `"AnsatzMode"` 元数据
- `generateAlphaSeeds()` 与 `alphaLevel()` 解耦，便于扩展新模式

---

## 5. 其他工作区的整理

### 5.1 已清理

| 操作 | 详情 |
|------|------|
| 移除过时验证文档 | `IBPVerification.md`、`Verify-MMA-KIRA-Guide.md`（verify/docs/） |
| 移除冗余对比脚本 | `Compare-Expand.wl`、`Compare-FamilyGenerate.wl`、`Compare-Reconstruct-bub00.wl`、`Compare-Results.wl`（VerifyUtility/） |
| 移除过时关系验证 | `VerifyRelation-SeriesVerify.wl`（VerifyUtility/） |
| 移除旧状态文档 | `CODE_STATUS_260504.md`（被本 DEV_STATUS.md 取代） |

### 5.2 新增文档与工具

| 文件 | 内容 |
|------|------|
| `docs/Plan-FamilyGenerate-CPP-Rewrite.md` | 架构设计计划（已更新至 2026-05-07 状态） |
| `docs/Plan-Ansatz-Redesign.md` | Ansatz 重设计方案 |
| `docs/AnsatzModes.md` | Ansatz 模式规范 |
| `docs/Benchmark_Results.md` | 性能基准测试结果（已更新 ne=7 + crossover） |
| `docs/Strategy-Comparison.md` | Strategy 0 vs 4 对比分析 |
| `docs/SymbolicRuleAlgorithm.md` | 约化规则产生：模 Gröbner 基三锥算法设计与 C++ 实现计划（T005） |
| `docs/RelationSolver_Documentation_Hub.md` | 文档索引 |
| `verify/VerifyUtility/README.md` | 验证工具使用说明 |
| `families/README.md` | Family 数据库说明 |
| `families/FamilyDatabase.wl` | 统一 Family 定义（19 个 family，按 L/E 排序） |
| `verify/README.md` | 验证目录总说明 |
| `DEV_STATUS.md` | 本文件 |

### 5.3 当前目录结构

```
project/
├── families/              # 9 个 JSON family 定义
├── include/               # 11 个管线头文件（header-only）
├── tools/                 # family_generate CLI + test_singular_runner / test_region_solver
├── tests/                 # test_relationFF / test_family_config / test_expandFF 等
├── src/                   # 非模板实现（LayerRecursionCore.cpp）
├── data/                  # 工作二进制数据（IBPMat_*.bin, RingData_*.bin, ExpansionCache_*.bin）
├── relations/             # 导出的关系结果（AllRelations_*.m, RelationMeta_*.m）
├── scripts/               # Shell脚本（compare_relation_*.sh）
├── checkpoint/            # 会话存档和检查点文件
├── mma/reference/         # MMA参考输出（ExpansionMMA_*.m, VerifyAllRelations.wl等）
├── docs/                  # 设计文档 + 性能报告（11 个 .md/.tex）
├── build/                 # CMake 构建产物
├── verify/
│   ├── FamilyDatabase/    # 统一 Family 定义 (MMA) + README
│   ├── VerifyUtility/     # 20 个 MMA 参考管线脚本 + README
│   │   └── cache/         # Kira 验证缓存（bub00/SR212/Tri）
│   ├── docs/              # 验证结果文档（3 个）
│   └── <family>/          # 12 个 family 的 MMA 参考 .bin + 元数据
└── DEV_STATUS.md          # 本文件
```

---

## 6. 效率优化：Strategy 设计

### 6.1 IncrementalNullspaceSolver 重构

`commit 8391228` — 将 RelationSolver 的线性求解器从批量模式重构为增量 RREF 维护：
- 每次 `addRows()` 对新行做增量消元，而非每次从头构建矩阵
- 简化 `nimax_lists`（缩减为 `{0}`，消除冗余解空间维度）
- 效果：Box (9.5x)、TB123 (7.5x+)、SR212-5m 从 OOM 降到 69s

### 6.2 Strategy 0 vs Strategy 4

`docs/Strategy-Comparison.md` 记录了两种方程构建策略的对比：

- **Strategy 0 (Default)**：每个 ν 点一次性评估所有 regime，构建完整矩阵 → RREF
- **Strategy 4 (PerNuRegimeRankCheck)**：每个 ν 点按 |θ| 降序逐 regime 评估，仅保留提升 rank 的行

策略 4 模拟 Laporta 层级消去：|θ| 更大的 regime 优先加入，冗余子扇区行被过滤。

### 6.3 AdaptiveSamplingConfig

可配置的采样参数（`include/RelationSolver.hpp`）：
- `min_nu` / `max_nu`：采样点范围
- `safety_factor`：过采样因子
- `nullity_stable_threshold`：收敛稳定阈值
- `verification_points`：额外验证点

### 6.4 C++ vs MMA 性能对比

C++ FamilyGenerate 管线（含 Singular 子进程）比 MMA 快 **7-34 倍**，主要收益来自：
- IBP 方程组装直接在 C++ 中完成（避免 MMA 符号计算开销）
- 仅 Groebner 基/准素分解通过 Singular 子进程（与 MMA 共享同一后端）

---

## 7. 约化规则产生：SymbolicRule 管线设计

### 7.1 背景

当前 §4 (`RelationSolver`) 使用纯高斯消元重构线性关系，等价于 **TopFlag1** 完备性。但 MMA `setupGeneratingCone` (SymbolicBTForm.wl) 揭示了更深层问题：约化系数的分母可能依赖于 ISP 变量 ($\nu$)，在分母零点上规则退化。

需要实现基于 **模 Gröbner 基**（position-over-term 排序）的完备约化规则生成，区分五级完备性（TopFlag1/2, FullFlag1/2, GlobalFlag），通过三锥结构（B1/B2/B3）逐层级消去目标积分。

### 7.2 设计文档

`docs/SymbolicRuleAlgorithm.md` — 完整算法设计文档，涵盖：

- 输入/输出接口定义（AllRelation + 目标积分 S → ν-符号规则 / 逐积分具体规则）
- 两种完备性检验对比（高斯消元 vs 模 Gröbner 基）
- 三锥结构（B1 高秩锥 → B2 中秩锥 → B3 近角落锥）递归流程
- 核心子程序 `eliminatedRule` + `sectorBlockReducible` 步骤拆解
- C++ 数据结构设计（`ConeRule`, `GeneratingConeResult`, `CompletenessFlags`）
- 7 步实现计划（P0–P7）、Singular 接口要点、验证策略

### 7.3 算法管道

```
RegionSolverAlgorithm  → Groebner → 准素分解 → RingData (.bin)
        ↓
LayerRecursionAlgorithm → IBP 矩阵级数展开 → resCache (.bin)
        ↓
ReconstructAlgorithm   → 高斯消元重构关系 → AllRelation (.m)
        ↓
SymbolicRuleAlgorithm  → 模 Gröbner 三锥规则 → 完备约化规则
```

### 7.4 当前状态 (2026-05-12)

- [x] `include/SymbolicRule.hpp` — 三锥结构主逻辑（~900 行）
  - [x] `buildBlockMatrix` 简化版（高斯消元用）
  - [x] `eliminatedRule` 高斯消元 + 模 Gröbner 基 (ν-多项式生成元)
  - [x] `setupGeneratingCone` B1/B2/B3 锥体迭代
  - [x] `generateSymbolicRules` 从 ConeRule 提取约化规则
  - [x] ν-多项式模 GB：β 向量编码 ν 指数，gen(i) 稀疏格式
  - [x] 无 ν 结构时跳过 Singular（β 全零 → 全对角元系数=1）
- [x] `include/RelationLoader.hpp` — AllRelation_*.m 解析器
- [x] `tests/test_generating_cone.cpp` — 5/5 PASS
- [ ] `tools/generating_cone.cpp` — CLI 入口
- [ ] SR212 全族 B1/B2/B3 与 MMA 对比验证
- [ ] 参数扫描：增大 coefdeg/rank 探索 TopFlag2 通过条件

### 7.5 相关任务

| Task | 内容 | 状态 |
|------|------|------|
| T005 | C++ SymbolicRule 实现 | 骨架完成，5/5 PASS |

---

## 附录 A：已知问题摘要

| 问题 | 影响 | 状态 |
|------|------|------|
| A/B 方程符号不一致 | .bin 不与 MMA 逐字节一致 | 已分析，T004 进行中 |
| LI 方程缺失 | 缺少 Lorentz Invariance 约束 | 待实现 |
| TB123 Singular 超时 | 大 family Groebner 基无法完成 | 待优化 |
| 5 个家族缺 JSON config | DB313/NP222/NP322/NP322m/TB123m | 待补充 |
| 临时调试打印残留 | IBPEqGenerator (ne≤4 导数) / RegionSolver (ne≤4 A/B) | 待清理 |
| SymbolicRule 管线缺失 | 约化规则仅 TopFlag1，无 Gröbner 完备性检验和三锥规则生成 | T005 设计中 |

## 附录 B：commit 记录 (2026-05-05 9:00 — 2026-05-07)

```
99ec194  Add DEV_STATUS.md: comprehensive project status for May 5-7 work
8c98ba2  Add multi-mode ansatz support to RelationSolver
9d5330e  Update Benchmark_Results: ne=7 timing, expand-vs-solve analysis
8391228  Refactor IncrementalNullspaceSolver to incremental RREF
baa0534  Fix .gitignore and add Compare-* utility scripts
8f0df6b  Fix VerifyExpand-Prepare: binary export shadowed by package context
33b3e14  Add VerifyUtility: Mathematica reference implementation
638ea9d  Add FamilyGenerate C++ rewrite architecture plan
b110432  Add FamilyGenerate C++ rewrite Phase 0-2: IBP equation gen
```
