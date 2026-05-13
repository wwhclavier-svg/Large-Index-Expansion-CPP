# SymbolicRuleAlgorithm: 基于模 Gröbner 基的完备约化规则生成

> 对应任务 T005 | 参考: `SetupGeneratingCone_Deep_Analysis.md` | MMA 源码: `SymbolicBTForm.wl`
>
> **管线位置**: Stage 4 — 在 ReconstructReductionRelation（Stage 3）输出的关系上验证完备性，生成完整约化规则。

---

## 1. 问题定义

### 1.1 输入

| | 类型 | 说明 |
|---|---|---|
| **输入** | **文件接口** | `AllRelations_<family>_k<order>.m`（来自 Reconstruct Stage 3 的 `.m` 文件）<br>形如 $\sum_i c_i(\nu) \; g(\nu - \alpha_i) = 0$ 的线性关系，$c_i(\nu)$ 是 $\nu$ 的多项式系数 |

### 1.2 输出（双模式）

**模式 A — ν-符号约化规则：**

将目标积分集合 $S$ 划分为子集 $S_1, S_2, \ldots$（$\bigcup S_i \supseteq S$），每个子集附带一组带符号 $\nu$ 的约化规则：

$$j(\nu - \alpha) \;\to\; \sum_k c_k(\nu)\; j(\text{更简单积分}), \quad c_k(\nu) \in \mathbb{Q}(\nu_1,\ldots,\nu_{ne})$$

当 $\nu$ 取特定值时，这些规则覆盖子集 $S_i$ 内所有积分的约化。

**模式 B — 逐积分具体约化规则：**

对每个目标积分 $j(n_1, n_2, \ldots, n_{ne})$，给出具体的数值系数约化规则：

$$j(n_1, n_2, \ldots) \;\to\; \sum_k c_k\; j(\text{主积分})$$

满足**递归代入可达性**：任意目标积分通过反复代入规则，最终可用主积分（master integrals）的线性组合唯一表示。

> **输出接口由 T005 P6 定义**：文件名 `RuleSet_<family>_k<order>.m`，格式待 T005 实现完成后固定。

---

## 2. 算法背景：两种完备性检验

### 2.1 为什么需要模 Gröbner 基

当前 C++ `RelationSolver::reconstructAllRelations()` 使用纯高斯消元（`reduceSolve`），等价于 **TopFlag1** 级别：检验在 generic ISP ($\nu$) 取值下，代数秩是否足够消去目标积分。

但 `setupGeneratingCone` 揭示了更深层的问题：**系数的分母可能依赖于 ISP 变量**。在这些分母的零点上，约化规则退化（分母为零），无法完成约化。

`eliminatedRule` 使用 **模 Gröbner 基**（position-over-term 排序，通过 Singular CAS）检测此问题：将系数视为 $\nu$ 的多项式，检查 Gröbner 基对角元（分母）的解是否依赖于 ISP 变量。

### 2.2 五级完备性标志

| 标志 | 含义 | 检验方式 |
|------|------|---------|
| **TopFlag1** | Generic 运动学点下，高秩积分可约化 | 高斯消元求代数秩 |
| **TopFlag2** | **任意**运动学点下，足够高秩的积分可约化 + 奇异性不依赖 ISP | Gröbner 基检查对角元 ISP 依赖 |
| **FullFlag1** | Generic 运动学点下，**所有**广义 sector 积分可约化 | 高斯消元（全 sector 种子） |
| **FullFlag2** | **任意**运动学点下，所有广义 sector 积分可约化 | Gröbner 基检查对角元 ISP 依赖 |
| **GlobalFlag** | 一切积分、一切秩、一切点都可约化（无奇异性） | 所有对角元 = 1 |

TopFlag1/FullFlag1 用 `reduceSolve`（C++ 已有），TopFlag2/FullFlag2 用 `eliminatedRule`（T005 待实现）。

**SR212 当前状态**：TopFlag1 通过但 TopFlag2 失败——24 个关系代数上足够，但约化系数的分母依赖 $\nu$，在某些 ISP 取值下退化。

---

## 3. 算法结构：三锥生成

`setupGeneratingCone` (SymbolicBTForm.wl:281-365) 构建三个层级锥，从高秩到近角落逐步生成完备规则。

### 3.1 数据流

```
重构关系 rel[g(ν-α_i)] = 0
     │
     ├── B1（高秩锥）: sectorBlockReducible[TopRankOnly→True]
     │     level ← rgen-1 递增至 TopFlag2=True
     │     产出: B1rule (消去所有高秩积分的约化规则)
     │
     ├── B2（中秩锥）: sectorBlockReducible[TopRankOnly→False]
     │     level ← B1level-1 递增至 FullFlag2=True
     │     产出: B2rule (消去中秩+高秩全部积分的约化规则)
     │
     └── B3（近角落锥）: sectorBlockReducible + 超扇区种子
           sectorup ← 超扇区深度 (含 supersector 种子)
           mistab 收敛: 不可约化独立积分列表不再变化
           产出: B3rule (覆盖近角落积分的约化规则)
```

### 3.2 B1: 高秩锥

```
B1Flag = False; level = rgen - 1
WHILE B1Flag == False AND level < levelbound:
    level += 1
    block1data = sectorBlockReducible[TopRankOnly → True](rel, sector, vlist, level)
    B1Flag = block1data["BlockElim"]["TopFlag2"]
    
B1rule = GroupBy[GBdiagonal, firstNonZero]  → 每个对角元位置对应一个约化规则
```

- 只检验高秩积分（`rnk ≥ level` 的 dot-seed 积分）
- 要求 `TopFlag2 == True`：高秩部分在任意 ISP 值下可约化
- 如果 B1 失败 → 返回 `{}`（无法生成规则）

### 3.3 B2: 中秩锥

```
B2Flag = False; level = B1level - 1
WHILE B2Flag == False AND level < levelbound:
    level += 1
    block2data = sectorBlockReducible[TopRankOnly → False](rel, sector, vlist, level,
        NumericReduce → {ν_isp → 0})
    B2Flag = block2data["BlockElim"]["FullFlag2"]
```

- 检验所有秩（高秩 + 中秩）的积分
- `NumericReduce → {ν_isp → 0}`: ISP 变量替换为数值（有限域），加速 Gröbner 基计算
- 要求 `FullFlag2 == True`
- 记录 `SingularLocus`：导致规则退化的 ISP 子流形

### 3.4 B3: 近角落锥

```
level = rgen - 1; sectorup = 1
B3mistab0 = {}
WHILE True:
    level += 1
    block3data = sectorBlockReducible[...,
        SeedLevelMinusRank → sectorup,  (* 超扇区种子 *)
        NumericReduce → {ν_pd → 1, ν_isp → 0}]
    
    B3redtab = block3data["BlockElim"]["ReducedVarTable"]
    B3mistab = { g(seed) | g 不在 B3redtab 中 }  (* 不可约化积分列表 *)
    
    IF Most[B3mistab] == B3mistab0:  BREAK  (* 收敛 *)
    ELSE:  {B3mistab0, block3data0} = {B3mistab, block3data}
```

- `sectorup`：超扇区深度，引入父扇区种子扩大搜索范围
- 收敛条件：`B3mistab`（无法约化的独立积分列表）随 level 不再变化
- `NumericReduce`: pd 变量 → 1, isp 变量 → 0，加速计算

### 3.5 超扇区种子展开 (Seed Generation)

对给定 sector $s = (s_1, \ldots, s_{ne})$ 和 level $L$，生成种子积分 $g(\nu + rk)$ 或 $g(\nu - rk)$：

- **dot-seed**: $g(\nu - rk)$，其中 $rk \ge 0$ 且 $\sum rk_i =$ 指定秩
- **minus-seed**: $g(\nu + rk)$，$rk$ 分量可跨 sector 边界（超扇区展开）
- **generating set**: `circminus[L+1, Range[ne], ne]` 中 $\sum \le L$ 的子集 + 边界项

C++ 端已有 `generateAlphaSeeds()` (RelationSolver.hpp:104-329) 实现了 4 种 ansatz 模式的种子生成，可复用并扩展。

---

## 4. 核心子程序

### 4.1 eliminatedRule (SymbolicBTForm.wl:142-204)

```
Input:  relmat (系数矩阵), relvar (积分变量), sector, level, vlist (ν变量)
Output: Association[GB, Lift, GSDenom, GSVar, BlockStart, BlockEnd,
                    SingularLocus, ReducedVarTable,
                    TopFlag1, TopFlag2, FullFlag1, FullFlag2, GlobalFlag]
```

**步骤**：

1. **排序积分变量** `relvar` 按 sector 层级 → 秩 ℓ → 种子 dot-rank 排序
2. **构建 Block 矩阵**：按 level 分块，行 = 重构关系，列 = 积分变量
3. **计算模 Gröbner 基**：
   ```
   SingularGroebnerBasis[
     relmat[:, 1:nvar] . Array[e, nvar],
     vlist, Array[e, nvar],
     "PositionOverTerm" → True,   (* POT 排序: 积分位置优先 *)
     "LiftMatrix" → True           (* 升维矩阵: GB = Lift . relmat *)
   ]
   ```
4. **对角元分析**：对每个对角元的分子 `gsblockdenom[pos]`，调用 `Solve[denom == 0, ISP变量]`
   - 解不依赖 ISP → 无 ISP 奇异点
   - 解依赖 ISP → 该位置对应关系在 ISP 子流形上退化
5. **五级标志判定**：综合 TopRank/FullRank + ISP 奇异检测

### 4.2 sectorBlockReducible (SymbolicBTForm.wl:212-279)

对给定 sector 和 level，构建 Block 消去系统并调用 `eliminatedRule`。

关键选项：
- `TopRankOnly`: 只展开高秩种子（rnk ≥ level）vs 全部种子
- `SeedLevelDot` / `SeedLevelMinusRank`: 控制 dot-seed 和 minus-seed 的展开层级
- `MinusRankFlat`: 生成 minus-seed 时是否限制分量上限
- `NumericReduce`: 将 ISP 变量替换为数值（有限域加速）

### 4.3 globalBlockReducible (SymbolicBTForm.wl:81-140)

快速诊断接口：纯高斯消元（不调 Gröbner），检验 TopFlag1 / FullFlag1。

C++ 端已有对应实现 `reconstructAllRelations()`，可作为快速预检。

---

## 5. C++ 实现计划

### 5.1 新增模块

| 模块 | 文件 | 行数估计 | 功能 |
|------|------|---------|------|
| SymbolicRule | `include/SymbolicRule.hpp` | ~600 | 三锥结构主逻辑 |
| RelationLoader | `include/RelationLoader.hpp` | ~150 | AllRelation_*.m 解析 |
| GeneratingConeTest | `tests/test_generating_cone.cpp` | ~200 | 集成测试 |
| SymbolicRuleCLI | `tools/generating_cone.cpp` | ~150 | CLI 入口 |

### 5.2 核心数据结构

```cpp
// 输入：从 AllRelation_*.m 解析的重构关系
struct RelationData {
    int ne;                                  // 外部变量数
    std::vector<std::vector<FFInt>> coeffs;  // 系数矩阵 [nrel][nvars]
    std::vector<std::vector<int>> shifts;    // α 向量 [nvars][ne]
    std::vector<std::vector<int>> seeds;     // g(ν-α) 的种子指标 [nvars][ne]
};

// 目标积分集合
struct TargetIntegrals {
    std::vector<int> sector;                 // 目标 sector
    std::vector<std::vector<int>> integrals; // 具体积分指标列表
};

// 完备性标志
struct CompletenessFlags {
    bool TopFlag1, TopFlag2;     // 高秩
    bool FullFlag1, FullFlag2;   // 全秩
    bool GlobalFlag;             // 全局无奇异
};

// 单个锥的输出
struct ConeRule {
    int level;                               // 锥高度
    std::vector<std::vector<FFInt>> denom;    // 对角元（分母）
    std::vector<std::vector<int>> reducedVar; // 可约化积分指标
    std::vector<std::vector<int>> singularISP; // ISP 奇异轨迹
    CompletenessFlags flags;
};

// 顶层输出 (模式 A / 模式 B)
struct GeneratingConeResult {
    std::vector<ConeRule> B1, B2, B3;          // 三锥规则
    std::map<std::vector<int>,                  // 目标积分指标
             std::vector<std::pair<FFInt, std::vector<int>>>> concreteRules; // 模式 B
    bool fullCoverage;                          // S 是否被完全覆盖
};
```

### 5.3 实现步骤

```
P0  接口定义精化: AllRelation 解析 + S 表示 + 规则输出格式
P1  算法分析: 通读 eliminatedRule / sectorBlockReducible / setupGeneratingCone
P2  数据结构: 实现 §5.2 核心 struct + AllRelation 加载器
P3  eliminatedRule: Singular POT Gröbner → 对角元 ISP 分析 → 五级标志
P4  sectorBlockReducible: 种子展开 + Block 矩阵 + eliminatedRule 封装
P5  三锥结构: B1 → B2 → B3 + 覆盖验证
P6  输出生成: 模式 A (ν-符号规则) + 模式 B (逐积分规则) + 递归代入验证
P7  集成验证: SR212 MMA 对比 + 覆盖率 + 性能基准
```

### 5.4 与现有代码的关系

```
现有 C++ 管线                       T005 新增
─────────────                      ──────────
IBPEqGenerator → IBPMat_*.bin
RegionSolver   → RingData_*.bin ──→ A/Ainv 用于 g(ν-α) 移位
LayerRecursion → resCache_*.bin ──→ 膨胀系数 g(ν) 的数值
RelationSolver → AllRelation_*.m ─→ eliminatedRule 的输入关系
                                     setupGeneratingCone → 约化规则
```

### 5.5 Singular 接口要点

`SingularRunner.hpp` 已有 subprocess 调用框架。`eliminatedRule` 需要的新功能：

1. **Position-over-term 排序**: `groebner(M, "POT")` → Singular 原生支持
2. **Lift matrix**: `lift(M)` → 返回 `(GB, Lift)` 满足 `GB = Lift · M`
3. **对角元提取**: 解析 GB 各多项式 → 提取首项系数 → 检查 ISP 变量是否出现
4. **有限域加速**: `ring r = (prime, vlist, dp)` → 与现有 `SingularRunner` 的 `Z/pZ` 一致

---

## 6. 参数策略

### 6.1 影响完备性的参数

| 参数 | 增大后效果 | 代价 |
|------|----------|------|
| **rank**（种子秩上限） | 更多种子 g → 更多代数约束 → 可能消除 ISP 依赖 | 种子数按 $O(\text{rank}^{ne})$ 增长 |
| **coefdeg**（系数多项式最高次数） | $b_0$ 系数自由度增大 → 更多参数消去分母 | 未知数增加 → 线性系统可能欠定 |
| **order**（1/n 展开阶数） | 每 regime/solution 更多方程约束 $b_0$ | 计算量增大 |
| **AnsatzMode** | 不同模式捕获不同代数结构 | 效果依赖具体积分族 |

### 6.2 推荐操作流程

1. **快速诊断**：`rank=2, coefdeg=1, Modulus→p`，看 `globalBlockReducible` 是否通过（≈ 已有 C++ 管线）
2. **TopFlag1 过但 TopFlag2 不过**：首选增大 **coefdeg**（2→3→4），次选增大 **rank**（2→3）
3. **有限域快速探索**：`Modulus → Prime[10^7]` 快速扫参数
4. **精确验证**：确认参数后 `Modulus → 0`

---

## 7. 验证策略

### 7.1 正确性验证

- **SR212**：当前已知 TopFlag1 通过但 TopFlag2 失败。验证增大 coefdeg 后 TopFlag2 是否通过
- **bub00 / Box**：已知 trivial 族，B1/B2/B3 应全部通过且 SingularLocus 为空
- **逐字节对比**：C++ B1/B2/B3 规则表与 MMA `setupGeneratingCone` 输出对比

### 7.2 覆盖率验证

- 对目标积分集合 $S$，验证每条规则链最终可到达主积分
- 检查所有 ISP 采样点（有限域随机采样）上规则不退化

### 7.3 性能基准

- 比较 Groebner 基计算 vs 高斯消元的耗时（per sector, per family）
- 测量 B1/B2/B3 各级 level 的迭代次数和收敛速度

---

## 8. 执行方式

| 入口 | 覆盖范围 | 命令 |
|---|---|---|
| `build/generating_cone` | Stage 4 符号规则生成 | `./build/generating_cone <family> <order>` |
`generating_cone` 读取 `relations/AllRelations_<family>_k<order>.m` 作为输入，通过 Singular CAS 计算模块 Gröbner 基，输出完备性验证结果。

## 9. 参考文献

- `SetupGeneratingCone_Deep_Analysis.md` — 核心问题分析（三锥结构 vs 代数秩检验）
- `SymbolicBTForm.wl:81-365` — MMA 参考实现
- `ReconstructAlgorithm.md` — MMA vs C++ 关系重构对比
- `ReconstructAlgorithm.md Appendix A` — 关系重构 MMA 包 API 参考
- `DEV_STATUS.md §4` — Ansatz 设计与积分约化（当前状态）
- `dev-status.md §6` — RemoveSolvedVariables 效率优化

## 相关文档

- `RegionSolverAlgorithm.md` — Region Solver 算法实现比较（管线上下游）
- `RelationSolver_ComponentGuide.md` — RelationSolver 组件完整指南（RemoveSolvedVariables 等）
