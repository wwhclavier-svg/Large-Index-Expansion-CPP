# RegionSolver Debug Plan

## 目的
全面测试 C++ RegionSolver 模块与 MMA 参考实现的一致性，覆盖从单区域（bub00）到多区域（Box）再到多维生成元（Tri）的各种场景。

---

## 当前测试状态

### ✅ 已通过测试

| 测试 | 族 | nb | VarIndep | 状态 |
|------|-----|----|----------|------|
| `test_region_solver_full` | bub00 | 1 | 空 | ✅ 1 region, VarRule/FR 正确 |
| `test_recursion_builder` | bub00 | 1 | 空 | ✅ M1+K2s=K1s, N1=F2D+K1 通过 |
| `test_recursion_builder_tri` | Tri | 1 | 空 | ✅ 2 regions, 矩阵恒等式通过 |
| `test_ring_builder` | bub00+Tri | 1 | 空 | ✅ A·A⁻¹=I 全部通过 |
| `test_ibp_analyzer` | bub00 | - | - | ✅ FTable 提取正确 |

### ⚠️ 待验证

| 测试 | 族 | 预期 | 状态 |
|------|-----|------|------|
| `test_ring_builder` (verify) | bub00/Tri | RingData 与 MMA .bin 对比 | 🔲 |
| Box 区域求解 | Box | 多 region / 多维 VarIndep | 🔲 |
| classifyVariablesPrimeA | 通用 | MMA 变量分类一致 | 🔲 |

---

## Debug 执行计划

### 阶段 1：单元验证（不依赖 Singular）

#### 1.1 PolyArith 单元测试
**文件**: `tests/test_poly_arith.cpp`
**命令**: `./build/test_poly_arith`
**验证点**:
- `polyAdd`: 合并同类项、模运算
- `polyMul`: 指数相加、系数相乘
- `polySubstitute`: 变量替换
- `canonicalize`: 排序+合并
- `parseSingularPolynomial`: 解析 "A1^2*B2+3" 格式
- `polyToSingularString`: 格式化回字符串

#### 1.2 IBPAnalyzer 独立测试
**文件**: `tests/test_ibp_analyzer.cpp`
**命令**: `./build/test_ibp_analyzer`
**验证点**:
- `buildABEquations`: A/B 多项式字符串格式
- `extractFTable`: F0/F1/F2D/f2 系数完整性
- gShift → A/B 转换正确性（v=all-1 约定）

#### 1.3 classifyVariablesPrimeA 独立测试
**目的**: 验证新变量分类函数的正确性
```cpp
// 手动构造测试用例
std::vector<std::string> Avars = {"A1", "A2", "A3"};
std::vector<std::string> Bvars = {"B1", "B2", "B3"};
std::vector<std::string> compGb = {"A1^2 - 2", "A2*A3 + 1", "B1*A1 - 1"};
// 预期: primeA = {"A1^2 - 2", "A2*A3 + 1"}
// 预期: vargen = {"A1"}, varpar = {"A2", "A3"}
// 预期: ximinpoly = {"A1^2 - 2"}
```
**关键验证**: 输出与 MMA `classifyVariablesPrimeA` 对比

---

### 阶段 2：RegionSolver 集成测试

#### 2.1 bub00 单区域测试
**文件**: `tests/test_region_solver_full.cpp`
**命令**: `./build/test_region_solver_full`
**预期输出**:
```
Found 1 regions
  nb=1, VarIndep=[], VarDep=[A1,A2]
  VarRule: A1=A1_const, A2=A2_const
  FractionRule: B[j,i] = const
  MonomialBasisMatrix: [1] (恒等矩阵)
```
**与 MMA 对比项**:
- VarRule 值（常数）
- FractionRule 值
- nb = 1

#### 2.2 Tri 多区域测试
**文件**: `tests/test_recursion_builder_tri.cpp`
**命令**: `./build/test_recursion_builder_tri`
**预期**: 2 regions，nb=1 各区域
**与 MMA 对比项**:
- region 数量
- VarRule 值（应不同）
- FractionRule 值
- M1/K1/K2 矩阵值

#### 2.3 Box 4-propagator 测试
**文件**: 需新建 `tests/test_region_solver_box.cpp`
**命令**: TBD
**背景**: Box 有 4 个传播子（ne=4），可能产生多 region 且 VarIndep 非空
**验证点**:
- region 数量（预期 MMA 有 3 regions）
- 每个 region 的 nb, VarIndep, VarDep
- VarRule 和 FractionRule
- MonomialBasisMatrix 维度

---

### 阶段 3：RecursionBuilder 验证

#### 3.1 递归矩阵恒等式
**所有族测试都验证**:
```
M1[m][i] + K2s[m][i] = K1s[m][i]  (mod p)
N1[m][i] = F2D[m][i] * I + K1[m][i]
```
**通过条件**: 所有 m,i, matrix element 相等

#### 3.2 RingBuilder A·A⁻¹=I 验证
**文件**: `tests/test_ring_builder.cpp`
**通过条件**: 对每个 A_i，矩阵乘积等于单位阵

#### 3.3 C++ vs MMA 二进制对比
**文件**: `tests/verify_ring_builder.cpp`
**命令**: `./build/verify_ring_builder`
**对比内容**:
- `RingData_bub00.bin` 中的 A_list / Ainv_list
- `RingData_Tri.bin` 中的 A_list / Ainv_list
**通过条件**: 字节级完全一致

---

### 阶段 4：Singular 集成问题排查

#### 4.1 Singular 进程调用
**检查点**:
- `SingularRunner::groebnerBasis` 返回非空
- `SingularRunner::minimalAssPrimes` 返回 0 维分量
- `reducePolynomialsSingular` 返回约化多项式

**常见错误**:
1. **Singular 未安装**: `which Singular` 应返回路径
2. **超时**: 增大 ` SingularRunner.hpp` 中的超时限制
3. **变量名冲突**: A/B 变量名不能是 Singular 关键字

#### 4.2 MonomialBasisMatrix 计算
**问题**: 当 VarIndep 非空时，`computeMonomialBasisMatrix` 可能返回全零
**原因**: `reducePolynomialsSingular` 使用 vargen ring，但 GB 引用了 VarDep 变量
**检查**: 调用前 VarDep 已被 VarRule 替换为 VarIndep 表达式

#### 4.3 VarRule 求解
**问题**: `solveVarRule` 对某些族返回空
**检查**: GB 中是否有形如 `A_i + P(vargen)` 的多项式
**备选**: 使用 MMA 计算 VarRule 作为参考

---

### 阶段 5：MMA 交叉验证

#### 5.1 generate AllRelations 对比
**文件**: `workspace/shared/VerifyUtility/VerifyExpand-Compare.wl`
**目的**: 对比 C++ 层递归展开与 MMA 展开结果
**命令**: MMA 中运行 `VerifyExpand[family, order]`

#### 5.2 RingData 二进制对比
**目的**: 验证 RingBuilder 输出的 A_i 和 Ainv_i 矩阵与 MMA 一致
**工具**: `tests/verify_ring_builder.cpp`

#### 5.3 EquationVerify 自洽性检查
**目的**: 将 C++ 系数代回 IBP 方程验证 M1*C + Total = 0
**注意**: 自洽性检查不能替代独立验证（因为生成器和验证器可能共用错误代码）

---

## 测试矩阵

| 族 | ne | nl | 预期 Regions | nb 范围 | VarIndep | 关键测试 |
|-----|----|----|-------------|---------|----------|---------|
| bub00 | 2 | 1 | 1 | 1 | 空 | 基线 |
| Tri | 3 | 1 | 2 | 1 | 空 | 多 region |
| Box | 4 | 1 | 3 | 1-2+ | 部分非空 | 多维 VarIndep |
| bub10 | 2 | 2 | ? | ? | ? | 高 nl |
| SR212 | 3 | 2 | ? | ? | ? | 高 nl, ne |

---

## 执行顺序

```
1. test_poly_arith          (5s)
   ↓
2. test_ibp_analyzer        (10s)
   ↓
3. test_region_solver_full  (30s)
   ↓
4. test_recursion_builder    (30s)
   ↓
5. test_recursion_builder_tri (60s)
   ↓
6. test_ring_builder         (60s)
   ↓
7. verify_ring_builder       (120s) — 需要 RingData_*.bin 文件
   ↓
8. test_region_solver_box    (NEW, 300s) — 需要 Box .bin 文件
```

---

## 已知问题

### 问题 1：VarDep 为空时 VarRule 解为空
- **现象**: bub00 的 VarDep = [A1, A2]，但 VarRule 有解
- **解释**: A_i 在 GB 中有线性形式 `A_i + const = 0`

### 问题 2：MonomialBasisMatrix 全零
- **现象**: VarIndep 非空时，矩阵全零
- **原因**: `reducePolynomialsSingular` 使用的 ring 与 GB 不匹配
- **解决**: 确保 `ximinpoly` 只含 VarIndep 变量

### 问题 3：ParallelSolver 无独立测试
- **现状**: `ParallelSolver.hpp` 只有类定义，没有单元测试
- **建议**: 添加 test_parallel_solver.cpp 验证并行高斯消元

### 问题 4：C++ vs MMA RingData 不一致（已修复部分）
- **根因1**: `parseSingularPolynomial` 在 `vargen=[]` 且 `varNames=[]` 时无法正确解析常数
  - **修复**: `PolyArith.hpp` 特殊处理 `nVars==0` 情况，确保非数字常数返回 `{}`
- **根因2**: `solveVarRule` 对空 `vargen` 生成的 VarRule 字符串格式错误
  - **修复**: `RegionSolver.hpp` 在 `vargen` 为空时直接输出常数而非调用 `polyToSingularString`
- **根因3**: `MonomialBasisMatrix` 在 `nb==0` 且 `VarIndep` 空时未正确初始化
  - **修复**: `solveRegion` 在 `nb==0` 且 `VarIndep` 空时显式设置 `MonomialBasisMatrix = [[[1]]]`
- **剩余差异**: C++ A_i 值仍与 MMA 不同（`A1=119616450` vs `MMA=59808223`），这是 VarRule 提取逻辑与 MMA 的系统性差异，需要进一步调查

### 问题 5：Tri region 数量不匹配
- **C++**: 4 regions，`solveAllSectors` 对所有 subsector 调用 `solveRegion`
- **MMA**: 3 regions，`regionsBySectors` 使用 `RegionWise` 策略，部分 subsector 可能被合并
- **需要**: 对比 C++ 和 MMA 的 subsector 列表，确认哪些被合并

---

## 成功标准

| 级别 | 标准 |
|------|------|
| P0 | 所有现有测试通过，无 regression | ✅ (test_ring_builder PASS) |
| P1 | `verify_ring_builder` bub00/Tri 与 MMA .bin 字节一致 | ⚠️ (A_i 值系统性差异) |
| P2 | Box 区域求解结果与 MMA 对比通过 | 🔲 |
| P3 | ParallelSolver 有独立测试覆盖 | 🔲 |
| P4 | Tri region 数量与 MMA 对齐 | 🔲 |

---

## 已修复问题

| 日期 | 问题 | 文件 | 修复内容 |
|------|------|------|---------|
| 2026-05-10 | `parseSingularPolynomial` 空变量列表无法处理非数字常数 | `PolyArith.hpp` | 特殊处理 `nVars==0` 情况 |
| 2026-05-10 | `solveVarRule` 空 `vargen` 生成错误 VarRule 字符串 | `RegionSolver.hpp` | 空 `vargen` 时直接输出常数 |
| 2026-05-10 | `MonomialBasisMatrix` 在 `nb>0` 且 `VarIndep` 空时未正确初始化 | `RegionSolver.hpp` | 显式设置 `MonomialBasisMatrix = [[[1]]]` |
| 2026-05-10 | `solveRegion` 中 `ximinpoly` 空的边界情况 | `RegionSolver.hpp` | 注释确认 `nb=0` 时跳过 MBM 计算 |

---

## 附录：关键文件位置

```
include/
  RegionSolver.hpp          — 主模块
  SingularRunner.hpp        — Singular 进程管理
  PolyArith.hpp             — 多项式工具
  IBPAnalyzer.hpp           — g→A/B 转换 + FTable
  RecursionBuilder.hpp      — 递归矩阵构建
  RingBuilder.hpp           — 环矩阵构建

tests/
  test_region_solver_full.cpp   — bub00 集成
  test_recursion_builder.cpp    — bub00 递归矩阵
  test_recursion_builder_tri.cpp — Tri 多 region
  test_ring_builder.cpp         — 环矩阵验证
  verify_ring_builder.cpp       — vs MMA 二进制

verify/{Family}/
  Compare-RegionInfo-{Family}.m  — MMA RegionData 对比
```
## 相关文档

- `../algorithms/RegionSolverAlgorithm.md` — Region Solver 算法实现比较
- `../verify/Test-RingData.md` — RingData 二进制验证工作流
- `../verify/Test-Expand.md` — C++ vs MMA 展开一致性验证
