# Test-RingData.md

C++ vs MMA RingData 二进制逐字节验证工作流。

---

## 概述

`verify_ring_builder` 是一个独立的 C++ 程序，用于将 C++ 管线生成的 `RingData_*.bin` 与 MMA/Singular 参考实现进行逐字节对比。每个 cell 代表一个 sector 的素理想分量，包含 Groebner 基计算 → 单项式基 → A/Ainv 矩阵。

测试入口：`tests/verify_ring_builder.cpp`，通过 `compareRingData()` 函数进行对比。

---

## 测试的积分族

| 族 | ne | Cells | 参考文件 |
|------|----|-------|---------|
| bub00 | 2 | 1 | `RingData_bub00.bin` |
| Tri | 3 | 4 | `RingData_Tri.bin` |
| bub10 | 2 | 2 | `RingData_bub10.bin` |
| bub11 | 2 | 5 | `RingData_bub11.bin` |
| Box | 4 | 15 | `verify/Box/RingData_Box-MMA.bin` |
| SR212 | 5 | 1 | `verify/SR212/RingData_SR212_MMA.bin` |
| TB123 | 7 | 37 | `RingData_TB123.bin` |

---

## 运行

```bash
cd /home/ykm/Large-Index-Expansion-CPP

# 单个族
./build/verify_ring_builder TB123

# 全部默认族（bub00, Tri, bub10, bub11）
./build/verify_ring_builder
```

---

## Cell 匹配逻辑

每个 sector 可分解为多个素理想分量（prime components），因此同一 sector 可能有多个 cell。`compareRingData()` 使用三层匹配：

1. **Sector 分组**：`std::map<sector, std::vector<index>>` 将 MMA cell 按 sector 分组
2. **nb 匹配**：对每个 C++ cell，在同 sector 的未匹配 MMA cell 中优先找相同 `nb` 的
3. **A[0] 矩阵精确匹配**：若同 sector 存在多个相同 nb 的 cell 且 A[0] 矩阵匹配，则确认为同一分量

用 `std::set<int>` 追踪已匹配的 MMA cell，避免重复匹配。

### 注意：不可使用单值 map

早期实现使用 `std::map<sector, int>`，由于 map 对同一 key 只保留最后一个值，导致多分量 sector 中的前几个 cell 全部匹配到错误的 MMA cell。**必须使用 multimap 或 `vector<int>` 值。**

---

## Sort 比较器的确定性

Cell 排序的三级键：

1. **Weight**（sector 元素求和）—— 升序
2. **Sector 二进制值**（逐元素降序比较）—— 降序
3. **nb**（单项式基大小）—— 升序

第三级键（nb）必须存在：同一 sector 的多个素分量具有相同的 weight 和 sector 向量，无 nb 键则 `std::sort` 对等价元素不保证稳定顺序，导致输出文件在不同运行间可能不同。

---

## 参考文件管理

### 何时需要重新生成参考文件

参考文件若在 T004 对齐（`groebner(I)` → `std(I)`，`SingularRunner.hpp`）之前生成，则为过期文件。过期参考与当前 C++ 输出存在矩阵级差异（Groebner 基生成元不同导致 A/Ainv 矩阵不同），不可直接对比。

### 重新生成流程

```bash
# 1. 运行验证（自动生成 RingData_<fam>-CPP.bin）
./build/verify_ring_builder <fam>

# 2. 备份旧参考
cp RingData_<fam>.bin RingData_<fam>.bin.bak

# 3. 用 C++ 输出替换
cp RingData_<fam>-CPP.bin RingData_<fam>.bin

# 4. 重新运行确认逐字节一致
./build/verify_ring_builder <fam>
# 预期：=== All comparisons MATCH ===
```

### 参考文件来源选择

| 来源 | 适用场景 | 注意 |
|------|---------|------|
| MMA 参考 (`*_MMA.bin`) | 使用 `primdecGTZ` 的 Mathematica 管线 | 与 C++ 算法一致，优先使用 |
| Singular 参考 (`*_Singular.bin`) | 使用 `primdecGTZ` 的 Singular 管线 | 素分量数量可能与 MMA/C++ 不同，**不建议直接使用** |

SR212 是典型例子：Singular 参考有 3 个 cell，而 MMA 参考和 C++ 输出均只有 1 个 cell。

---

## 常见问题

### 1. `groebner(I)` vs `std(I)` 差异

T004 对齐将 `SingularRunner::groebnerBasis()` 中的 `groebner(I)` 改为 `std(I)`。两者计算的是同一理想，但生成元集合不同：
- **单分量 sector**：nb=1 的 cell 矩阵完全相同（只有一个基元素，不涉及 NF 计算差异）
- **多分量 sector**：nb 相同但 A/Ainv 矩阵不同（基元素相同但坐标变换矩阵受 Groebner 基生成元影响）

### 2. 同一 sector 的 cell 数量不一致

C++ 和 MMA 对同一 sector 可能产生不同数量的 cell。常见原因：
- MMA 使用 `primdecGTZ` 而 Singular 参考使用 `singularBivarPrimeInfo`
- `groebner(I)` vs `std(I)` 导致素理想分解结果不同

遇到此情况，优先以 MMA 参考为准。

### 3. `hasContent` 过滤

`compareRingData()` 在生成 C++ regions 前会检查 AB 方程是否有实际内容。空方程直接跳过，不调用 `solveRegion`。如果 MMA 参考中某 sector 的 cell 在 C++ 侧全部被跳过，说明该 sector 的 AB 方程为空（通常是高权重 sector 退化的结果）。

## 相关文档

- `../../docs/algorithms/RegionSolverAlgorithm.md` — Region Solver 算法（RingData 来源）
- `../../docs/plans/RegionSolver-Debug-Plan.md` — C++ RegionSolver 调试计划（含 verify_ring_builder）
