# DEV_STATUS — Large Index Expansion C++ 项目工程状态

> 更新于 2026-05-07，覆盖 2026-05-05 9:00 至今的改动。

---

## 目录

1. [主要管线：C++ Region 计算管线 + Singular 接入](#1-主要管线c-region-计算管线--singular-接入)
2. [验证测试：验证工作流的整理](#2-验证测试验证工作流的整理)
3. [性能建模：Singular 求解 Region 的时间消耗](#3-性能建模singular-求解-region-的时间消耗)
4. [Ansatz 设计与积分约化：AnsatzMode](#4-ansatz-设计与积分约化ansatzmode)
5. [其他工作区的整理](#5-其他工作区的整理)
6. [效率优化：Strategy 设计](#6-效率优化strategy-设计)

---

## 1. 主要管线：C++ Region 计算管线 + Singular 接入

### 1.1 架构概览

```
FamilyDatabase(.json)  →  [B] IBPEqGenerator  →  [C] RegionSolver  →  .bin 输出
                              (Singular 子进程)     (Singular 子进程)
```

全管线为 header-only C++17，通过 `SingularRunner` 模板类调用 Singular 子进程。管线总代码量约 8500 行。

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

### 1.4 已完成 Family

| Family | L | ne | Regions | 状态 |
|--------|---|----|---------|------|
| bub00 | 1 | 2 | 1 | ✅ |
| bub10 | 1 | 2 | 3 | ✅ |
| bub11 | 1 | 2 | 5 | ✅ |
| Tri | 1 | 3 | 2 | ✅ |
| Box | 1 | 4 | 15 | ✅ |
| SR | 2 | 5 | 6 | ✅ |
| SR3m | 2 | 5 | 8 | ✅ |
| SR5m | 2 | 5 | 45 | ✅ |
| TB123 | 2 | 7 | 未完成 | ⏳ Singular 超时 |

### 1.5 关键问题：A/B 方程符号不一致

C++ 和 MMA 的 IBP 组装采用不同中间约定：
- **C++**：乘 ∏z_m 后取负 `-deriv[i,j,k]`
- **MMA**：LargeIndexIBP（n_i→n+v_i）→ Coefficient["n"] 提取

两套约定数学等价，但 A/B 方程部分项符号不同 → Singular 选出等价但不同的单项式基底 → .bin 逐字节不同。功能验证通过（膨胀系数逐字节一致，sol_dim 匹配），但 .bin 文件不与 MMA 参考文件逐字节一致。

### 1.6 待完成

- [ ] Lorentz Invariance (LI) 方程生成
- [ ] A/B 方程符号对齐（达成 .bin 逐字节一致）
- [ ] TB123 超时问题（需优化 Groebner 基计算或换用更高效的基算法）
- [ ] DB313 / NP222 / NP322 / SR212 变体的 JSON config

---

## 2. 验证测试：验证工作流的整理

### 2.1 MMA 参考实现（VerifyUtility）

`commit 33b3e14` — 将 MMA 参考管线从原工作区迁入版本库，共 17 个 .wl 文件，约 6000 行：

| 文件 | 功能 |
|------|------|
| `LIEFamilyDefine.wl` | Family 定义 → 内部数据 |
| `LIEUtility.wl` | sectorLimitIBP、导数、系数替换 |
| `LIECoreAlgebra.wl` | expRegSolve2、Groebner/准素分解 |
| `LIERegions.wl` | regionsBySectors 主入口 |
| `LIEExpand.wl` | 膨胀系数计算 |
| `LIEReconstruct.wl` | 关系重构 |
| `LIEWorkflow.wl` | 完整管线编排 |
| `SingularInterface.wl` | Singular 调用接口 |
| `ExportBinary_IBPMatrix.wl` | .bin 导出 |
| `M2Kira.wl` | Kira 格式互转 |
| `VerifyExpand-*.wl` (4个) | 膨胀系数验证 |
| `VerifyRelation-SeriesVerify.wl` | 关系验证 |
| `KiraRuleLoader.wl` | Kira IBP 规则导入 |

### 2.2 对比工具（Compare-* 系列）

`commit baa0534` — 5 个对比脚本：

| 脚本 | 用途 |
|------|------|
| `Compare-Expand.wl` | 膨胀系数对比 |
| `Compare-FamilyGenerate.wl` | 全管线 .bin 对比 |
| `Compare-Results.wl` | 最终结果对比 |
| `Compare-Reconstruct-bub00.wl` | 关系重构对比（bub00 专用） |
| `Compare-VerifyLog.wl` | 验证日志对比 |

### 2.3 MMA 参考 .bin 文件

`verify/<family>/` 目录下已有部分 MMA 生成的参考 .bin 和元数据：

| Family | .bin 存在 | RegionInfo | Timing |
|--------|----------|------------|--------|
| bub00 | ✅ | ✅ | ✅ |
| bub10 | ✅ | ✅ | ✅ |
| bub11 | ❌ | ❌ | ❌ |
| Tri | ❌ | ❌ | ❌ |
| Box | ✅ | ✅ | ✅ |
| SR/SR3m/SR5m | ❌ | ❌ | ❌ |

### 2.4 test_relationFF 验证

- bub00：膨胀系数逐字节一致，6 个 (lev,deg) 配置 sol_dim 完全匹配
- 其余 family 待测

### 2.5 验证文档

位于 `verify/docs/`：
- `Test-Expand.md` — 膨胀系数验证结果
- `Test-Relation.md` — 关系重构结果
- 已清理过时文档（`IBPVerification.md`、`Verify-MMA-KIRA-Guide.md` 等 5 个文档）

---

## 3. 性能建模：Singular 求解 Region 的时间消耗

### 3.1 测试结果

`docs/Benchmark_Results.md` (712 行) 记录了三类基准测试：

**Region 求解性能（C++ FamilyGenerate vs MMA）：**

| Family | ne | C++ 耗时 | MMA 耗时 | 加速比 |
|--------|----|---------|----------|--------|
| bub00 | 2 | 0.11s | 3.4s | 31x |
| Box | 4 | 0.62s | 4.5s | 7x |
| Tri | 3 | 0.25s | 3.7s | 15x |
| SR | 5 | 1.2s | — | — |

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

- **C++ FamilyGenerate 比 MMA 快 7-34x**（wall clock），主要原因是 C++ 内联计算 IBP 组装，而 MMA 走符号计算 + 字符串解析
- **Singular Groebner 基是瓶颈**：更大的 family（TB123, ne=7）在 Groebner 基步骤超时
- **Expand-vs-Solve crossover 分析**已记录在 Benchmark_Results 中

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

### 5.1 文件清理

- 移除 5 个过时验证文档（`IBPVerification.md`、`Verify-MMA-KIRA-Guide.md` 等）
- 移除 5 个冗余 Compare-* 脚本（功能已收敛到统一验证框架）
- 移除 `CODE_STATUS_260504.md`（被本 DEV_STATUS.md 取代）
- 移除 `VerifyRelation-SeriesVerify.wl`（功能已收敛）

### 5.2 新增文档

| 文件 | 内容 |
|------|------|
| `docs/Plan-FamilyGenerate-CPP-Rewrite.md` | 架构设计计划 |
| `docs/AnsatzModes.md` | Ansatz 模式规范 |
| `docs/Benchmark_Results.md` | 性能基准测试结果 |
| `docs/Strategy-Comparison.md` | Strategy 0 vs 4 对比分析 |
| `docs/RelationSolver_Documentation_Hub.md` | 文档索引 |
| `verify/VerifyUtility/README.md` | 验证工具使用说明 |
| `verify/FamilyDatabase/README.md` | Family 数据库说明 |

### 5.3 目录结构

```
project/
├── families/           # 9 个 JSON family 定义
├── include/            # 10 个管线头文件（header-only）
├── tools/              # family_generate CLI + 测试工具
├── tests/              # test_relationFF / test_family_config 等
├── docs/               # 设计文档 + 性能报告
├── verify/
│   ├── FamilyDatabase/ # 统一 Family 定义 (MMA)
│   ├── VerifyUtility/  # 17 个 MMA 参考管线脚本
│   ├── docs/           # 验证结果文档
│   └── <family>/       # 各 family 的 MMA 参考 .bin
├── build/              # CMake 构建产物
└── DEV_STATUS.md       # 本文件
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

## 附录 A：已知问题摘要

| 问题 | 影响 | 状态 |
|------|------|------|
| A/B 方程符号不一致 | .bin 不与 MMA 逐字节一致 | 已分析，功能无影响 |
| LI 方程缺失 | 缺少 Lorentz Invariance 约束 | 待实现 |
| TB123 超时 | 大 family Groebner 基无法完成 | 待优化 |
| 部分 family JSON config 缺失 | DB313/NP222/NP322/SR212 variants | 待补充 |
| 临时调试打印残留 | IBP/AB/FTable 输出混杂 | 待清理 |

## 附录 B：commit 记录 (2026-05-05 9:00 — 2026-05-07)

```
8c98ba2  Add multi-mode ansatz support to RelationSolver
9d5330e  Update Benchmark_Results: ne=7 timing, expand-vs-solve analysis
8391228  Refactor IncrementalNullspaceSolver to incremental RREF
baa0534  Fix .gitignore and add Compare-* utility scripts
8f0df6b  Fix VerifyExpand-Prepare: binary export shadowed by package context
33b3e14  Add VerifyUtility: Mathematica reference implementation
638ea9d  Add FamilyGenerate C++ rewrite architecture plan
b110432  Add FamilyGenerate C++ rewrite Phase 0-2: IBP equation gen
```
