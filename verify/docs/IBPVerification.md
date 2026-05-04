# IBP 验证测试套件

## 三种单元测试

### 1. CompareVerify - C++ vs MMA 展开结果比较

**目的**：验证 C++ 和 MMA 两种实现产生完全相同的展开系数。

**文件**：
- `Compare-Results-[famname].wl` - MMA 脚本
- `Compare-CPPResult-[famname].m` - C++ 输出

**运行**：
```bash
# Step 1: 生成 .bin 文件
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-FamilyGenerate-[famname].wl

# Step 2: C++ 展开
cd /root/Large-Index-Expansion-CPP/build
./test_expandFF [famname]

# Step 3: MMA 展开
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-Expand-[famname].wl

# Step 4: 比较
cd /root/Large-Index-Expansion-CPP
wolframscript -file Compare-Results-[famname].wl
```

**验证内容**：
- 对每个阶 k=0..4，比较 h^(k) 多项式
- 结果：`[MATCH] k=0,1,2,3,4` 表示完全一致

---

### 2. EquationVerify - IBP 方程验证

**目的**：验证展开系数满足递归关系 `M1 * C + Total = 0`

**文件**：
- `tests/IBPVerification.hpp` - 验证器类
- `tests/test_IBPVerification.cpp` - 主程序

**运行**：
```bash
cd build
make test_IBPVerification
./test_IBPVerification [famname]
```

**验证内容**：
- 对每个 (k, l, seed) 组合
- 构建 `inhomogTerms`（包含 NMinus, NZero, NPluMi, NPlus, M1, MPlus）
- 检查 Total[m][j][i] 是否为零

**算法**：
```
For each (k, l, seed):
    buildAll(ibpmat, C, k, l, seed, nindep, ncurr)
    if Total != 0: FAIL
```

---

### 3. SeriesVerify — IBP 方程级数代入验证

**目的**：将展开级数代回 IBP 方程，验证在给定截断阶以下各阶 1/n 系数归零。

**原理**：
IBP 方程（Large-Index 形式）为 Σ_α c_α(n, ν) · g[ν+α] = 0，其中 g[ν+α] = G(ν+α)。

1. **系数展开**：c_α(n, ν) 直接展开为 n 的多项式（不移动 ν，因为算符因子 (n+ν_a) 已包含大参数 n）：c_α(n, ν) = Σ_j c_α^(j)(ν) · n^j
2. **积分展开**：g(ν+α) = A^(+(ν+α)) · Σ_k h^(k)(ν+α) / n^k，其中 A^(+(ν+α)) = ∏ A_i^(ν_i+α_i) 为主导渐近因子（匹配 eigen-substitution g[ν+α] → A^α）
3. **代入归零**：因子 A^(+ν) ≠ 0 提出后，对 m = 0,1,...,K，1/n^m 的系数必须为零：

```
Σ_{α, j}  A^(+α) · c_α^(j)(ν) · h^(j+m)(ν+α) ≡ 0  (mod p)
```

其中 A^(+α) = ∏ A_i^(α_i) 为有限域常数，α_i 为 g[ν+α] 中的整数偏移。

**文件**：
- `VerifyExpand-SeriesVerify.wl` — MMA 脚本（已通用化）

**运行**：
```bash
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-SeriesVerify.wl <famname> [order]
```

**验证内容**：对每条 IBP 方程、每个 1/n^m 阶 (m=0..K)，检查代入展开后的 residual 是否 ≡ 0 (mod p)。

**与 EquationVerify (1b) 的区别**：
- EquationVerify 验证**生成**展开系数的层递归内部一致性（M1·C + Total = 0）
- SeriesVerify 验证**展开产物**是否真正满足原始 IBP 方程的平移形式（外部一致性）
- 两者互补，覆盖"递归算法正确"和"最终结果满足方程"两个层面

---

## 测试覆盖

| 验证类型 | 验证内容 | 通过标准 |
|---------|---------|---------|
| CompareVerify | C++ vs MMA | 所有 k=0..4 MATCH |
| EquationVerify | 递归方程 | Total == 0 |
| SeriesVerify | IBP 方程代入 | 所有 m=0..K residual ≡ 0 |

---

## 运行完整测试流程

```bash
# 1. CompareVerify (展开系数比较)
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-Prepare.wl bub00
cd /root/Large-Index-Expansion-CPP && ./build/test_expandFF bub00
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-MMAExpand.wl bub00
wolframscript -file VerifyExpand-Compare.wl bub00

# 2. EquationVerify (IBP递归方程验证)
cd /root/Large-Index-Expansion-CPP/build
./test_IBPVerification bub00

# 3. SeriesVerify (IBP方程级数代入验证)
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-SeriesVerify.wl bub00

# 4. RelationVerify (关系重建比较)
cd /root/Large-Index-Expansion-CPP/build
./test_relationFF bub00 4 2 2
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-ExportRelation-bub00.wl
cd /root/Large-Index-Expansion-CPP
cp build/Compare-CPPRelation-bub00.m .
cp Compare-MMARelation-bub00.m .
wolframscript -file Compare-Relation-Results.wl
```

---

## 当前测试状态

### 1-loop 积分族

| 积分族 | NE | NB (Top) | Ring Type | CompareVerify | EquationVerify | SeriesVerify |
|--------|----|---------|-----------|-------------|----------------|--------------|
| bub00  | 2 | 1 | F_p (NB=1) | k=0..4 MATCH | PASS | PASS (m=0..3) |
| bub10  | 2 | 2 | F_p[A]/(MinPoly) | k=0..4 MATCH | PASS | PASS (m=0..3) |
| bub11  | 2 | 1 | F_p (NB=1) | k=0..4 MATCH | PASS | PASS (m=0..3) |
| Tri    | 3 | 1 | F_p (NB=1) | k=0..4 MATCH | PASS | PASS (m=0..3) |
| Box    | 4 | 2 | F_p[A]/(MinPoly) | k=0..4 MATCH, 12矩阵 | PASS | PASS (m=0..3) |

### SeriesVerify 实现细节

- **NB=1 路径**：A_i 为 F_p 中数值常量，直接进行模运算。h^(k) 为标量多项式。
- **NB>1 路径**：A_i 属于域扩张 F_p[A_free]/(MinPoly)，表示为单项式基 {1, A, ..., A^{NB-1}} 中的 NB-向量。环运算通过 companion matrix 方法实现：
  - M1 = A 的 companion matrix（次对角线为1，最后一列为 MinPoly 系数的负值）
  - basisMatrices[[k+1]] = M1^k 表示乘以 A^k 的矩阵
  - ringMultiply(v1, v2) = mat(v1) · v2（其中 mat(v) = Σ v_k · M1^k）
  - ringInverse(v) = LinearSolve[mat(v), e1, Modulus->p]
  - ringPower(v, n) = 重复平方/求逆
  - 该方法适用于任意 NB ≥ 1，不限于 NB=2

### 已知问题

- **bub10/bub11**：`bub10` 和 `bub11` 有相同的 propagators 但不同的 Numeric 设置。C++ 端 `bub10` 配置可能不准确（使用 msq=0 而非 msq=1），需在 C++ 端进行数值测试以确认一致性。
- **Box checkpoint**：`PrepareCheckpoint-Box.wdx` 由于 `AbsoluteTiming` 中的 `CompoundExpression` 漏洞而存储 Null，导致每次运行 MMAExpand 或 SeriesVerify 都需要重新生成。

---

## 重要教训：FFInt 符号转换 Bug

### 问题描述

C++ 展开结果对所有积分族都产生错误系数。具体症状：
- **EquationVerify (自洽性检查)** 通过 — `M1*C + Total == 0`
- **SeriesVerify (独立验证)** 失败 — 将展开代入原始 IBP 方程 residual ≠ 0
- **CompareVerify (MMA对比)** 失败 — C++ 系数与 MMA 不匹配

例如 bub00 的 h^(1) 线性项：C++ 输出 `105689379`，正确值应为 `14952056`（≡ MMA 的 `-164472617 mod 179424673`）。

### 为什么 EquationVerify 通过但结果错误？

EquationVerify 检查的是 `M1*C + Total = 0`，这是**层递归算法自身的方程**（生成系数时使用的同一组方程）。当 `sgn()` 返回错误值时：
1. 系数生成时，`buildAll()` 用错误的 sgn 值构建 Total
2. 然后用 `M1*C + Total = 0` 求解出系数
3. EquationVerify 复用同一套 `buildAll()` 逻辑，所以 Total 仍为零 → **自洽但不正确**

**核心教训**：自洽性检查不能替代独立验证。必须用独立的数学关系（原始 IBP 方程代入）或参考实现（MMA）进行交叉验证。

### 根因定位过程

**Step 1 — 排除数据问题**：加 F0/N1 debug print → 确认 .bin 输入数据与 MMA 一致 ✓

**Step 2 — 定位错误方程**：在 `assembleLinearSystem` 后加 debug print → 发现 `Total[0] = 115516503`（应为 `119616449`），RHS 偏差 `4099946`

**Step 3 — 逐分量追踪**：在 `buildAll()` 中分别打印 NMinus, NZero, NPluMi, M1, MPlus 分量 → 定位到 `MPlus[2]` 贡献了 `175324727`（应为 `0`）

**Step 4 — 逐项追踪 MPlus**：在 `add_MPlus_contribution` 中打印每个 K1s/K2s 和 F2s 项的 mat, factor, src → 发现 F2s 耦合项的 `factor=4099945`（应为 `179424672 ≡ -1`）

**Step 5 — 追溯 factor 来源**：`factor = l2_sign * static_cast<T>(bin1 * bin2)`，其中 `l2_sign = sgn(1) = 4099945`（应为 `-1 ≡ 179424672`）

**Step 6 — 定位 sgn() 缺陷**：
```cpp
// Utilities.hpp:13 (原实现)
inline T sgn(int l) {
    return static_cast<T>((l % 2 == 0) ? 1 : -1);
}
```
`FFInt` 只有 `FFInt(uint64_t)` 构造函数。`static_cast<FFInt>(-1)` 中 `int(-1)` 被隐式转换为 `uint64_t(2^64-1)`，再 `mod p = 4099945`。[Python验证：`(2^64-1) % 179424673 = 4099945`]

### 修复

```cpp
// Utilities.hpp:13 (修复后)
inline T sgn(int l) {
    if (l % 2 == 0) return static_cast<T>(1);
    return -static_cast<T>(1);  // 用 operator-() 而非从 int(-1) 构造
}
```

### 经验总结

1. **`static_cast<T>(负整数)` 对无符号构造函数的类型是陷阱** — 当 T 的构造函数接受 `uint64_t` 时，负 `int` 会静默转换为巨大的无符号值（C++ 隐式转换规则），再 mod-p 产生看似"合法"但完全错误的值（无编译警告）
2. **自洽性验证不等于正确性** — 当 bug 存在于方程构建代码中时，生成器和验证器使用相同的错误逻辑，自洽性检查必然通过。必须使用**独立于算法实现**的验证方法
3. **调试方法**：从输出差异反向追溯，逐层加打印定位到具体的数值偏差（Total → MPlus → F2s factor → sgn），每一步用 Python 独立计算期望值进行对比
4. **Python 辅助验证**：对于有限域计算，Python 的 `pow(x, -1, p)` 可以快速验证模逆、解方程等，是独立于 C++ 实现的黄金参考

---

## 4. RelationVerify - 关系重建比较

**目的**：验证 C++ 和 MMA 两种实现的关系重建结果一致。

**文件**：
- `tests/test_relationFF.cpp` - C++ 关系重建（含统一格式导出）
- `Compare-ExportRelation-[famname].wl` - MMA 关系重建导出
- `Compare-Relation-Results.wl` - 比较脚本

**统一格式** (`$RelationResult`)：
```mathematica
$RelationResult = <|
  "Family" -> "bub00",
  "Lev" -> 2,
  "Deg" -> 2,
  "NE" -> 2,
  "Modulus" -> 179424673,
  "Alphas" -> {{0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {2,0}},
  "Betas" -> {{0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {2,0}},
  "Coefficients" -> {{...}, {...}, ...},  (* 36 x 31 矩阵 *)
  "Relations" -> {"Particular" -> "...", "Basis1" -> "..."},
  "HasSolution" -> True,
  "NumVariables" -> 36,
  "NumSolutions" -> 30
|>
```

**运行**：
```bash
# Step 1: 生成家族数据
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-FamilyGenerate-bub00.wl

# Step 2: C++ 关系重建
cd /root/Large-Index-Expansion-CPP/build
./test_relationFF bub00 4 2 2

# Step 3: MMA 关系重建
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-Expand-bub00.wl
wolframscript -file Compare-ExportRelation-bub00.wl

# Step 4: 比较
cd /root/Large-Index-Expansion-CPP
cp build/Compare-CPPRelation-bub00.m .
cp Compare-MMARelation-bub00.m .
wolframscript -file Compare-Relation-Results.wl
```

**验证内容**：
- 元数据一致：Family, Lev, Deg, NE, Modulus
- 索引一致：Alphas, Betas 列表相同
- 系数一致：Coefficients 矩阵相同（模 modulus）

---

## 重要教训：KinematicRules 中的 dot product 规则

### 问题描述

运行 Compare-FamilyGenerate-Tri.wl 时，SINGULAR 报错：
```
? `p1` is not defined
? error occurred in or before STDIN line 5: `...p1*p2...`
```
导致 `#Non-triv Sectors = 0`，无法生成有效的 binary 文件。

### 根本原因

在三角积分族（Tri）和 Box 积分族的定义中，propagators 展开后包含 `p1*p2` 类型的交叉项。例如：
```
-(k1 + p2)^2 + msq = -k1^2 - 2*k1*p2 - p2^2 + msq
```

MMA 的 `LIEDefineFamily` 会将这些多项式传递给 SINGULAR 进行计算，但 `p1`、`p2` 等符号变量没有在 KinematicRules 中被消除。SINGULAR 在处理时遇到了未定义的 `p1` 变量。

### 解决方案

**对于 Tri (3点积分)**：
```mathematica
KinematicRules -> ({p1^2 -> s1, p2^2 -> s2, (p1 + p2)^2 -> s3, p1*p2 -> (s3 - s1 - s2)/2} /. numericRules)
```

**对于 Box (4点积分)**：
```mathematica
KinematicRules -> ({
    p1^2 -> s, p2^2 -> t, p3^2 -> u,
    (p1 + p2)^2 -> 0, (p2 + p3)^2 -> 0, (p1 + p2 + p3)^2 -> s + t + u,
    p1*p2 -> (0 - s - t)/2,
    p2*p3 -> (0 - t - u)/2,
    p1*p3 -> (s + t + u - s - t - u)/2 - (0 - s - t)/2 - (0 - t - u)/2
} /. numericRules)
```

### 公式说明

对于任意两个动量的点乘，在质壳条件下有：
```
p_i · p_j = (p_i + p_j)^2 - p_i^2 - p_j^2 / 2
```

因此：
- `p1*p2 = ((p1+p2)^2 - p1^2 - p2^2) / 2`
- 对于 Tri：`s3 - s1 - s2` 除以 2
- 对于 Box：需要根据具体数值计算

### 经验总结

1. **检查 propagators 中是否包含交叉项**（如 `p1*p2`、`p2*p3`）
2. **如果包含，必须在 KinematicRules 中添加对应的 dot product 规则**
3. **规则需要在 numericRules 应用后仍然有效**，确保 SINGULAR 能正确处理

---

## 文件清单

| 文件 | 作用 |
|------|------|
| `Compare-Results-[famname].wl` | CompareVerify |
| `tests/test_IBPVerification.cpp` | EquationVerify |
| `VerifyExpand-SeriesVerify.wl` | SeriesVerify |
| `tests/test_relationFF.cpp` | RelationVerify (C++) |
| `Compare-ExportRelation-[famname].wl` | RelationVerify (MMA) |
| `Compare-Relation-Results.wl` | RelationVerify (比较) |