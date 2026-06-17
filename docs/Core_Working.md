# Core_Working.md — 项目核心模块索引

本文档汇总项目的三类核心资源：理论基础/数学背景、工作流程、完整测试。

---

## 1. 理论基础 / 数学背景

### 核心理论文档

| 文件 | 内容 |
|------|------|
| `LIENotes/Notes LargeIndexExpansion 250810.tex` | LIE 核心理论——费曼积分大指标展开的完整数学框架 |
| `LIENotes/Notes_LIE_260414_Theory.pdf` | 理论笔记 PDF (2026-04-14) |
| `LIENotes/Notes_LIE_260414_CPP.pdf` | C++ 实现设计笔记 PDF (2026-04-14) |
| `LIENotes/Symbolic Block Triangular Form 251209.tex` | IBP 矩阵分块三角结构 (beamer) |
| `LIENotes/ThematicReport260107.tex` | 专题报告：Large Index Expansion for Feynman Integral Reduction |
| `LIENotes/report2509.tex` ~ `report2603.tex` | 月度总结 (2025.09–2026.03) |
| `docs/Reconstruct_Algorithm.md` | 关系重构算法理论框架：数学推导 (§1)、MMA 实现 (§2)、C++ 实现 (§3)、两种方法的符号 vs 数值差异 (§4)、ν 采样验证的数学原理 (§5)。**理论，不含现时测试结果** |
| `docs/LayerRecursion_Algorithm.md` | 层递归算法：多重指标、种子生成、非齐次项组装 |
| `docs/ReconstructReductionRelation_Documentation.md` | Mathematica 包 API 文档 (`ReconstructReductionRelation`) |

### 算法架构文档

| 文件 | 内容 |
|------|------|
| `docs/RelationSolver_ComponentGuide.md` | RelationSolver 组件：类型系统、自适应多点采样 |
| `docs/RelationSolver_QuickReference.md` | RelationSolver 快速参考 |
| `docs/RelationSolver_Documentation_Hub.md` | 文档中心/索引 |
| `CLAUDE.md` | 项目总览、构建、架构、数据流 |
| `AGENTS.md` | 双语详尽指南（中文/English） |

---

## 2. 工作流程

> **验证工作流总览**：[verify/README.md](../verify/README.md) — 展开验证 (3 法) + 关系验证 (4 法) 的全景图、参数配置表、文件分布。

### 完整数据处理管线

```
Integral Family Definition → LIE Expansion → Result Comparison
        (MMA)                   (MMA + C++)     (MMA + C++)

Step 1: Generate .bin       Step 2: Expand    Step 3: Compare
Compare-FamilyGenerate →    test_expandFF /   Compare-Results
                            Compare-Expand
                            
Step 4: Reconstruct         Step 5: Verify
test_relationFF /           Kira ± Verify
Compare-Reconstruct         CPP-KiraVerify
```

### MMA 管线脚本 (`/root/Large-Index-Expansion-MMA-Mini/`)

| 阶段 | 文件 | 说明 |
|------|------|------|
| 族定义 | `Compare-FamilyGenerate-*.wl` | 生成 `IBPMat_*.bin`, `RingData_*.bin` |
| 展开 | `Compare-Expand-*.wl` | MMA LIE 展开计算 |
| 结果对比 | `Compare-Results-*.wl` | C++ vs MMA 结果对比 |
| 关系重构 | `Compare-Reconstruct-*.wl` | 调用 `LIEReconstruct` 重构关系 |
| Kira 验证 | `NuSamplingVerify-MMA.wl` | ν 采样 + Kira 约化验证 |
| 统一验证 | `Unified-NuSampling-Verify.wl` | 统一 ν 采样验证框架 |

### MMA 核心库 (`/root/Large-Index-Expansion-MMA-Mini/`)

| 文件 | 功能 |
|------|------|
| `LIEReconstruct.wl` | 主重构算法 (`LIEReconstruct` 函数) |
| `LIEWorkflow.wl` | 工作流包装 (`LIESolveRegions`, `LIEExpandSeries`, `LIEGetRelations`) |
| `LIEFamilyDefine.wl` | 积分族定义 |
| `LIEExpand.wl` | 级数展开引擎 |
| `LIECoreAlgebra.wl` | 核心代数运算 |
| `LIERegions.wl` | 区域管理 |
| `LIEUtility.wl` | 工具函数 |
| `ExportBinary_IBPMatrix.wl` | 导出 .bin 格式 IBP 矩阵 |
| `KiraRuleLoader.wl` | Kira 约化规则加载器 |
| `/root/M2Kira.wl` | Kira 输入生成接口 |

### C++ 管线 (本仓库)

| 阶段 | 文件 | 说明 |
|------|------|------|
| 矩阵加载 | `include/IBPMatrixLoader_Binary.hpp` | 加载 MMA 生成的 `.bin` 文件 |
| 展开 | `tests/test_expandFF.cpp` | C++ 层递归展开计算 |
| 关系重构 | `tests/test_relationFF.cpp` | C++ 自适应采样关系重构 |
| 核心求解器 | `include/RelationSolver.hpp` | `reconstructReductionRelation<T>()` |
| 线性求解 | `include/LinearSolver.hpp` + `LinearSolver_FF.hpp` | 有限域线性系统求解 |

### 工作流文档

| 文件 | 描述 |
|------|------|
| `verify/docs/Test-Expand.md` | 四步展开测试：Generate → Expand(CPP) → Expand(MMA) → Compare |
| `verify/docs/Test-Relation.md` | 关系重构工作流与现时测试：Kira 验证全流程 (§3–§6)、验证结果总结、MMA vs C++ 差异分析、待解决问题。**结果随测试更新** |
| `verify/docs/Verify-MMA-KIRA-Guide.md` | Kira 交叉验证完整工作流 |
| `verify/docs/IBPVerification.md` | 三种验证方法：CompareVerify, EquationVerify, SeriesVerify |
| `verify/docs/CPP-KiraVerify-Debugger.md` | C++ 到 Kira 验证的调试报告 |

---

## 3. 完整测试

### C++ 测试可执行文件

| 测试 | 用途 | 命令 |
|------|------|------|
| `test_expandFF` | 有限域展开系数计算 | `./build/test_expandFF [family]` |
| `test_relationFF` | 有限域线性关系重构 | `./build/test_relationFF [family] [order] [lev] [deg]` |
| `test_IBPVerification` | IBP 矩阵方程验证 | `./build/test_IBPVerification` |
| `test_expand_family` | 展开族测试（命令行参数） | `./build/test_expand_family [family]` |
| `test_load_bub` | bub 格式矩阵加载测试 | `./build/test_load_bub` |
| `test_ff_verify` | FireFly 库基本功能验证 | `./build/test_ff_verify` |

### MMA 验证脚本

| 文件 | 位置 | 验证内容 |
|------|------|---------|
| `Compare-Results-*.wl` | MMA-Mini | C++ vs MMA 展开系数对比 |
| `Compare-Reconstruct-*.wl` | MMA-Mini | LIE 关系重构 + 内置验证 |
| `NuSamplingVerify-MMA.wl` | MMA-Mini | ν 采样 Kira 验证 (MMA 关系) |
| `Unified-NuSampling-Verify.wl` | MMA-Mini | 统一 ν 采样验证 (MMA + C++) |
| `NuVerify-Relations.wl` | `verify/scripts/` (本仓库) | ν 采样验证 C++ 关系 |

### 独立验证工具

| 文件 | 用途 |
|------|------|
| `tools/test_ff_verify.cpp` | FireFly 有限域运算验证 |
| `tools/test_firefly_simple.cpp` | FireFly 最简功能测试 |
| `tests/archive/` | 历史存档测试（双精度展开、加载器测试等） |

### 验证结果摘要

| 验证类型 | 结果 | 参考文档 |
|---------|------|---------|
| C++ MatrixBuild: M(nu) × coeffs = 0 | **110/110 通过** | `verify/docs/Test-Relation.md` |
| C++ Kira 约化 (单项基底) | 2/10 抵消, 8/10 部分抵消 | `verify/docs/Test-Relation.md` |
| MMA Kira 约化 (d=1/3,s=3) | 0/4 失败 (当前) | `verify/docs/Test-Relation.md` |
| MMA ν={0,0} 平凡验证 | 4/4 通过 (不检验 IBP) | `docs/Reconstruct_Algorithm.md` |
| 历史: d=1/13,s=1 | 60/60 通过 (2026-04-28) | `verify/docs/Test-Relation.md` |
| C++ e2e: test_expandFF + test_relationFF | residual=0 全通过 | `verify/docs/Test-Relation.md` |

---

## 跨模块关联

```
理论基础:  LIENotes/*.tex  →  docs/*.md (算法文档)  →  CLAUDE.md/AGENTS.md
    ↓
工作流:   MMA-Mini/*.wl (MMA)  +  tests/*.cpp (C++)  +  include/*.hpp (C++)
    ↓
测试:     Compare-Results (对比)  +  KiraVerify (约化)  +  MatrixBuild (矩阵)
```

### 关键外部依赖

| 组件 | 路径 | 用途 |
|------|------|------|
| FireFly | `/root/firefly-2.0.3/` | C++ 有限域运算库 |
| Eigen3 | 系统安装 | C++ 线性代数 |
| GMP | 系统安装 | 多精度算术 |
| Mathematica | 系统安装 | `.wl`/`.m` 脚本执行 |
| Kira | `/root/kira/` | IBP 约化规则生成 |
| M2Kira | `/root/M2Kira.wl` | Mathematica → Kira 接口 |
| Kira 测试数据 | `/root/kira_tests/bub00_new/` | bub00 约化规则 |
