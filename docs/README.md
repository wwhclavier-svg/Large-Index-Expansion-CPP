# Large Index Expansion — 文档索引

## 目录结构

```
docs/
├── README.md                        ← 本文档
├── restructure-plan.md              ← 项目重组记录
├── algorithms/                      ← 核心算法理论与实现对比
│   ├── RegionSolverAlgorithm.md     — Region Solver 管线比较（MMA vs C++）
│   ├── ReconstructAlgorithm.md      — 关系重构算法比较（含 MMA API 附录）
│   ├── SymbolicRuleAlgorithm.md     — 符号规则生成与完整性标志
│   └── LayerRecursion_Algorithm.md  — Layer Recursion 展开算法
├── components/                      ← 组件 API 与使用指南
│   ├── RelationSolver_Documentation_Hub.md  ← 导航中心（推荐起点）
│   ├── RelationSolver_ComponentGuide.md     — 完整组件参考
│   ├── RelationSolver_QuickReference.md     — 快速参考
│   ├── AnsatzModes.md                       — Ansatz 模式详解
│   └── Strategy-Comparison.md               — 求解策略对比
├── plans/                           ← 调试与测试计划
│   └── RegionSolver-Debug-Plan.md   — C++ RegionSolver 调试计划
├── reports/                         ← 基准与进度报告
│   ├── Benchmark_Results.md         — 全族性能基准数据
│   └── Progress_Report_April_May_2026.tex
└── verify/ → workspace/shared/verify-docs/  ← 验证工作流
    ├── Test-Expand.md               — C++ vs MMA 展开一致性验证
    ├── Test-Relation.md             — 关系输出验证
    ├── Test-RingData.md             — RingData 二进制验证
    └── CPP-KiraVerify-Debugger.md   — Kira 验证 g 格式 bug 记录
```

## 阅读路径

### 新读者
`Hub.md` → `ComponentGuide.md` → `QuickReference.md` → `ReconstructAlgorithm.md`

### 算法开发
`RegionSolverAlgorithm.md` → `SymbolicRuleAlgorithm.md` → `LayerRecursion_Algorithm.md`

### 调试/测试
`RegionSolver-Debug-Plan.md` → 对应 `docs/verify/` 工作流 → 回查 `algorithms/`

### 性能分析
`Benchmark_Results.md` → `Strategy-Comparison.md` → `RegionSolverAlgorithm.md Appendix B`

### 数据流
```
families/*.json
  → RegionSolver (Stage 1)        — RegionSolverAlgorithm.md
    → IBPMat_*.bin + RingData_*.bin
      → LayerRecursion (Stage 2)  — LayerRecursion_Algorithm.md
        → 展开系数
          → Reconstruct (Stage 3) — ReconstructAlgorithm.md
            → AllRelations_*.m
              → SymbolicRule (Stage 4) — SymbolicRuleAlgorithm.md
                → 完整约化规则
```

## 可执行文件与管线对应

| 可执行文件 | 覆盖管线 | 源码 | 类型 |
|---|---|---|---|
| `build/family_generate` | Stage 1 (RegionSolver) | `tools/family_generate.cpp` | CLI: `family_generate <family.json>` → `.bin` |
| `tools/family_scheduler.py` | Stage 1 (调度包装) | `tools/family_scheduler.py` | Python 调度器: `family_scheduler <family> --workers N` |
| `build/test_relationFF` | Stage 2+3 (展开+重构) | `tests/test_relationFF.cpp` | CLI: `test_relationFF <family> <order> <lev_min> <lev_max> <deg_max>` |
| `build/test_expandFF` | Stage 2 (仅展开) | `tests/test_expandFF.cpp` | CLI: `test_expandFF <family>` |
| `build/generating_cone` | Stage 4 (符号规则) | `tools/generating_cone.cpp` | CLI: `generating_cone <family> <order>` |

> **Stage 2/3 合并**：`test_relationFF` 在一个进程中依次调用 LayerRecursion + Reconstruct，通过 `seriesCoefficient<T>` 在内存中直接传递展开系数，无文件中间态。

## 编译

```bash
mkdir -p build && cd build && cmake .. && make -j$(nproc)
```

CMakeLists.txt 位于项目根目录。`tests/test_relationFF.cpp` 等测试入口编译为独立二进制文件，输出到 `build/`。`tools/` 中的 `.cpp` 文件为独立工具入口，编译为对应的独立二进制文件。`tools/family_scheduler.py` 为 Python 脚本，无需编译。

## 任务模板与 Skill 对应

| Stage | CXXX 任务模板 | 实现 Skill | VXXX 验证模板 | 验证 Skill |
|-------|--------------|-----------|--------------|-----------|
| Stage 1 (RegionSolver) | `region-solve-mma-template` | `region-solver-mma-scheduled` | `verify-ringdata-template` | `verify-ringdata` |
| Stage 1 (C++ 无调度) | `region-solve-cpp-template` | `large-index-expansion` | — | — |
| Stage 1 (C++ 调度) | `region-solve-cpp-scheduled-template` | `region-solver-cpp-scheduled` | — | — |
| Stage 2 (LayerRecursion) | `relation-test-cpp-template` (Stage 2+3) | — | `verify-expand-template` | `verify-expand` |
| Stage 3 (Reconstruct) | ↑ 同上 (合并于 test_relationFF) | — | `verify-relation-template` | `verify-relation` |
| Stage 4 (SymbolicRule) | (T005 开发中) | — | — | — |

> **VXXX 验证流程**: 通过 skill 触发（如 `verify-relation`），使用对应模板创建验证任务，记录 pass/fail。
