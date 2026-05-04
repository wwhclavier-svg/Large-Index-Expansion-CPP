# CODE STATUS — 2026-05-04

主干程序（`Large-Index-Expansion-CPP`）相对最初版本（初始提交 `5facc13`, 2026-03-19）的改动总结。

时间跨度：2026-03-19 ~ 2026-05-04，共 12 个 commit，涉及 36+ 个文件，+5454/-1458 行。

---

## 1. 核心算法修复 (commit `204d650`)

### RelationSolver.hpp — 自适应方程构建器重构

- 每次迭代都调用 `solver.update()` 跟踪零空间维度稳定性
- **关键改进**：验证失败时将失败的验证方程重新加入求解系统，使零空间维度单调递减，确保最终收敛到正确解
- `Mext` 格式从 `[cols] x [nullity]` 改为 `[cols] x [1 + nullity]`，第一列留作特解位置
- 新增 `exportRelationCoefficientToMMA()` 函数，将关系系数导出为 Mathematica `.m` 格式
- `LinearSystemResult` 结构新增 `pivot_cols` 字段

### LayerRecursionCore.tpp — 二项式系数溢出修复

- 将 `BINOM[a][b] * BINOM[c][d]` 拆成两个 `long long` 分别取值再乘，避免中间结果溢出

### LinearSolver_FF.hpp — 有限域求解器

- 归一化和消元操作改用 `= a * b` 显式写法替代 `*=`
- 新增 `pivot_cols` 追踪

### LinearSolver_Eigen.hpp — 解验证修复

- `if constexpr` 分支处理浮点类型的 `isApprox()` 和有限域类型的 `==` 精确比较

### LinearSolver.hpp

- 结果中转发 `pivot_cols`

---

## 2. 缓冲区溢出修复 (commit `0f57166`)

### LayerRecursion.hpp

- `nimax` 不再硬编码为 `4`，改为动态计算
- 遍历所有 k/l 组合，取 `BINOM[l+ne-1][ne-1] * nb` 的最大值

---

## 3. 新增功能

### SeriesCoefficientIO.hpp

- 新增 `exportAllResultsToMMA()` — 将级数展开结果导出为 Mathematica `.m` 格式

### include/firefly/FFInt.hpp（新文件，78 行）

- 最小化 FFInt 实现，用于编译测试

### UnifiedStorage.hpp

- 新增 `size()` 方法

### IBPMatrixLoader.hpp / IBPMatrixLoader_Binary.hpp

- `IBPMatrixE` 结构新增 `incre` 字段存储
- 清理死代码和无用注释

---

## 4. 测试变更

| 操作 | 文件 |
|------|------|
| 重写 | `test_relationFF.cpp` (341→504行) |
| 新增 | `test_IBPVerification.cpp`, `test_expand_family.cpp`, `test_load_bub.cpp`, `test_ff_verify.cpp`, `test_firefly_simple.cpp`, `IBPVerification.hpp` |
| 删除 | `test_RelationNew.cpp`, `test_recons.cpp`, `test_expand.cpp`, `IBPMatrixLoader_Test.cpp`, `RingDataLoader_test.cpp` |

---

## 5. 文档新增

- `CLAUDE.md` (282行)
- `ReconstructAlgorithm.md` (605行)
- `Test-Expand.md` (232行)
- `Test-Relation.md` (472行)
- `docs/LayerRecursion_Algorithm.md` (205行)
- `AGENTS.md` 大幅更新

---

## 6. 其他

- `CMakeLists.txt` 更新（新增 `test_expand_family`, `test_load_bub` 等构建目标）
- 符号修正 (commit `ebf89df`)：将 j[α] 替换为 g[ν-α] 避免歧义
- `.gitignore` 更新

---

## 7. 未提交 bug 修复 (commit `057d511`, 2026-05-04)

编译和测试验证通过 (`test_expandFF bub00`, `test_relationFF bub00`)。

### seed 索引偏移修复（5 个文件）

- `include/LayerRecursion.hpp`, `include/LayerRecursionCore.tpp`, `src/LayerRecursionCore.cpp`, `src/layerRecursion.cpp`:
  `max(ncurr, 0)` → `max(ncurr - 1, 0)` — 修正方程组变量枚举的 off-by-one 错误

### f2 卷积存储维度修复 (RelationSolver.hpp)

- `f2_len` / `total_k` / `rowsPerNu()`: `k_max_ + 1` → `deg_ + k_max_ + 1`
- `f2_store_` 分配大小同步更新
- 卷积索引从 `target_pow = l + k` 改为 `idx = deg_ - l + k`，正确处理 g 函数的负指数项
- 修复前：超出 `k_max_` 范围的项被静默丢弃，矩阵缺行导致结果错误但不崩溃

### sgn() 有限域溢出修复 (Utilities.hpp)

- `static_cast<T>(-1)` → `-static_cast<T>(1)` — 避免无符号有限域类型的溢出

### 测试配套更新

- `test_expandFF.cpp`: 新增 `exportMetaToMMA()`；order 改为命令行参数；路径改到 `verify/<fam>/`
- `test_relationFF.cpp`: `equations_per_sample` 公式同步更新；j[α]→g[ν-α] 符号修正

### 文件重组

- 旧脚本移动到 `archive/`、`docs/`、`mma/`、`tools/`、`verify/docs/`
