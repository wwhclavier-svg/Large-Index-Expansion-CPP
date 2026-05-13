# Layer Recursion Algorithm Documentation

## 1. 算法概述

层递归（Layer Recursion）是一种用于计算 IBP（Integration By Parts）矩阵级数展开系数的算法。核心思想是**按阶数（order）递增计算，每阶内部按层（layer）递减处理**。

**管线位置**: Stage 2 — 读取 RegionSolver 输出的 `.bin` 文件，计算展开系数，供 ReconstructReductionRelation 进行关系重构。

### 输入/输出接口

| | 类型 | 说明 |
|---|---|---|
| **输入** | **文件接口** | `IBPMat_<family>.bin`, `RingData_<family>.bin`（参见 RegionSolverAlgorithm.md §1.2）<br>加载器: `IBPMatrixLoader_Binary.hpp` → `IBPMatrix`, `RingDataLoader.hpp` → `RingData` |
| **处理** | 内存中 | 每 region 按阶数 order 递增计算，每阶内部按层 layer 递减处理 |
| **输出** | **内存接口** | `seriesCoefficient<T>` 密集数组 `(k,l,seed,j,i)` — 5 维展开系数<br>与下游 Reconstruct 在同一进程中直接传递（`test_relationFF` 内函数调用），无文件序列化 |
| **下游** | 内存调用 | `ReconstructAlgorithm.md` — 接收 `seriesCoefficient` 数组，转换为 `hexpnList` + `aregList` 调用 `reconstructAllRelations()` |

## 2. 数学背景

### 2.1 多重指标与种子（Multi-index and Seed）

- **维度 `ne`**: 外部变量个数
- **种子 `seed = {a1, a2, ..., ane}`**: 非负整数向量，满足 `sum(seed) = l`（层级）
- **层级 `l`**: 种子分量之和，决定多项式的总次数
- **容量 `getCapacity(ne, l)`**: 层级 l 下种子个数，即 `C(l + ne - 1, ne - 1)`

### 2.2 系数存储

系数 `C(k, l, seed, j, i)` 是一个五维数组：
- `k`: 展开阶数（0 到 kmax）
- `l`: 层级别（0 到 incre*k）
- `seed`: 多重指标向量
- `j`: 基索引（0 到 nb-1）
- `i`: 解索引（0 到 nimax），i=0 为特解，i>0 为齐次解

## 3. 算法流程

### 3.1 主循环结构

```
for k = 1 to kmax:
    for C in CTable (branching):
        for l = incre*k - 1 downto 0:
            for seed in seedsAtLevel(l):
                // 组装线性方程组
                // 求解
                // 更新系数
        // 分支处理
```

### 3.2 核心计算：inhomogTerms::buildAll

对于每个 (k, l, seed)，组装线性方程 `M1 * x = -inhomog`，其中：

```
inhomog = NMinus + NZero + NPluMi + NPlus + M1_contrib + MPlus
```

#### NMinus（来自低阶的贡献）
```cpp
// C++: LayerRecursionCore.tpp:341-354
if (l > incre * (k - 1) + 1) return;
for each j where seed[j] > 0:
    NMinus += N1[m][j] * C(k-1, l, seed - unit(j))
```
**物理含义**: 从下层借用到当前层

#### NZero（当前阶的齐次项）
```cpp
// C++: LayerRecursionCore.tpp:359-369
if (l > incre * (k - 1)) return;
opMat = F0 + sum_j seed[j] * K1[j]
NZero = opMat * C(k-1, l, seed)
```
**物理含义**: 同阶次的直接贡献

#### NPluMi（F2 耦合项）
```cpp
// C++: LayerRecursionCore.tpp:377-393
if (l > incre * (k - 1)) return;
for each i, j where seed[j] > 0:
    factor = -(seed[i] + 1)
    NPluMi += F2[m][i][j] * factor * C(k-1, l, seed + unit(i) - unit(j))
```
**物理含义**: F2 算子的梯度耦合

#### NPlus（来自低阶高层的贡献）
```cpp
// C++: LayerRecursionCore.tpp:410-459
lennplus = incre * (k - 1) - l
for l1 = 1 to lennplus:
    // K1 贡献
    for i where seed[i] >= 1:
        factor = seed[i] * BINOM[seed[i]+l1][l1+1]
        NPlus += K1[m][i] * factor * C(k-1, l, seed + l1*unit(i))
    // F2 耦合
    for l2 = 1 to l1+1 (l2 != l1):
        for i, j where seed[j] >= 1:
            factor = sgn(l2) * BINOM[seed[i]+l2][l2] * BINOM[seed[j]+l1-l2][l1-l2+1]
            NPlus += F2[m][i][j] * factor * C(k-1, l, seed + l2*unit(i) + (l1-l2)*unit(j))
```

#### M1_contrib（当前阶的 M1 贡献）
```cpp
// C++: LayerRecursionCore.tpp:396-407
for i = 0 to ncurr-1:
    factor = seed[i] + 1
    M1_contrib += M1[m][i] * factor * C(k, l, seed + unit(i))
```
**物理含义**: M1 算子作用于高一层

#### MPlus（当前阶的高层贡献）
```cpp
// C++: LayerRecursionCore.tpp:462-492
lenmplus = incre * k - l - 1
for l1 = 2 to lenmplus:
    // K1s, K2s 贡献
    for i:
        factor = BINOM[seed[i]+l1][l1]
        mat = K1s + sgn(l1) * K2s
        MPlus += mat * factor * C(k, l, seed + l1*unit(i))
    // F2s 耦合
    for l2 = 1 to l1-1:
        for i, j:
            factor = sgn(l2) * BINOM[seed[i]+l2][l2] * BINOM[seed[j]+l1-l2][l1-l2]
            MPlus += F2s[m][i][j] * factor * C(k, l, seed + l2*unit(i) + (l1-l2)*unit(j))
```

### 3.3 线性系统求解

对于每个 (k, l, seed)，组装方程：
```
M1 * x + sum_{active} C(k,l,seed+unit) = -inhomog
```

其中 `active` 是 `seed[lastNonZero(seed)+1 ... ne-1]` 的范围。

### 3.4 回代与分支

求解后：
1. **回代（Back-substitution）**: 用新解更新低层系数
2. **分支（Branching）**: 如果产生新独立元，则分裂为多个分支

## 4. C++ 与 MMA 实现对比

### 4.1 数据结构

| 方面 | C++ | MMA |
|------|-----|-----|
| 系数存储 | `seriesCoefficient<T>` 密集数组 | `csol` 关联表 |
| 索引方式 | (k,l,cid,j,i) → offset 计算 | `c[k, seed, j]` 符号访问 |
| 默认值 | 数组初始化为 0 | 关联表默认返回 0 |

### 4.2 关键差异点

#### 差异 1: NZero 计算中的 Lookup 行为

**C++** (`LayerRecursionCore.tpp:359-369`):
```cpp
// 直接从 seriesCoefficient 数组获取，已计算的值为实际值，未计算的为 0
T* src = getValuePtrOffSet(C, k - 1, l, seed);  // C(k-1, l, seed)
addProductTo(dest, mat_buf1, src, nindep_p1);
```

**MMA** (`LIEExpand.wl:310-327`):
```mathematica
(* 从关联表 csol 中查找，如果不存在返回默认值（通常为 0） *)
(F0 + sumK1) . LookupOrSelf[csol, c[i - 1, seed]]
```

#### 差异 2: 种子排序

两者都使用 **Degree-reverse-lexicographic 序**（按 Total 降序，Total 相同时按 Reverse 降序）：
```
level=0: {0,0}
level=1: {1,0}, {0,1}
level=2: {2,0}, {1,1}, {0,2}
```

### 4.3 符号处理

`sgn(l)` 函数在两者中等价：
```cpp
// C++: utility::sgn(l) = (l % 2 == 0) ? 1 : -1
// MMA: sgn[l_] := If[EvenQ[l], 1, -1]
```

## 5. 正确性验证

### 5.1 理论不变式

1. **初始条件**: `C(0,0,{0,...,0}, j, 0) = 1`（当 j=0 时），其他为 0
2. **线性系统**: 每个方程组必须相容（有解）
3. **能量守恒**: 同阶多项式系数总和应满足约束

### 5.2 调试方法

1. **打印种子顺序**: 验证 C++ 和 MMA 处理相同顺序的种子
2. **打印非齐次项**: 对比 `inhomogTerms::buildAll` 的中间结果
3. **打印线性系统**: 对比方程组系数矩阵和右端项
4. **打印解向量**: 对比求解结果

## 6. 执行方式

LayerRecursion 作为内嵌模块，不提供独立 CLI。通过以下方式调用：

| 入口 | 覆盖范围 | 命令 |
|---|---|---|
| `build/test_relationFF` | Stage 2+3 串行 | `./build/test_relationFF <family> <order> <min_lev> <max_lev> <max_deg>` |
| `build/test_expandFF` | Stage 2 仅展开 | `./build/test_expandFF <family> [--verbose]` |

两入口均在 process 内调用 `LayerRecursionCore::expandFamily()`，从 `include/LayerRecursion.hpp` 导入 `IBPMatrixLoader_Binary.hpp` + `RingDataLoader.hpp` 加载输入 `.bin` 文件。`test_relationFF` 的展开结果通过 `seriesCoefficient<T>` 数组直接传入 `reconstructAllRelations()`，无文件中间态。

## 7. 已知差异

当 C++ 和 MMA 出现系数不匹配时，排查步骤：

1. **验证 IBP 矩阵数据**: 检查 M1, N1, F0, F2 等算子是否一致
2. **验证种子顺序**: 打印 seedlist 对比
3. **验证非齐次项分解**: 分别计算 NMinus, NZero, NPluMi, NPlus, M1_contrib, MPlus
4. **验证线性系统**: 打印完整的方程组
5. **验证求解器**: 使用相同的测试矩阵验证求解器

---

## 相关文档

- `RegionSolverAlgorithm.md` — Region Solver 管线上游（生成 LayerRecursion 的输入 .bin 文件）
- `ReconstructAlgorithm.md` — 关系重构算法管线下游（使用 LayerRecursion 的展开系数）
- `Benchmark_Results.md` — LayerRecursion 展开性能数据

*文档版本: 1.0*
*最后更新: 2026-04-28*