# Test-Expand.md

C++ 和 MMA 展开计算对比测试工作流文档。工作流已通用化：通过 famname 参数切换积分族，无需为每个族复制脚本。

---

## 三条并行验证方法

展开系数的正确性可通过三种独立方法交叉确认。**推荐 1a + 1b 组合**达到最佳覆盖。

| 方法 | 验证内容 | 独立性 | 状态 |
|------|---------|:---:|:---:|
| **1a. CompareVerify** | C++ vs MMA 展开系数逐阶 diff | 需 MMA 参考 | k=0..4 MATCH |
| **1b. EquationVerify** | IBP 递归方程 M1*C+Total=0 | 纯 C++ 自洽 | PASS |
| **1c. SeriesVerify** | IBP 方程级数代入 (ν→ν+θn) | 自洽 (需除 `A^(-ν)`) | 待测试 |

> 参数配置与工作流见 [verify/README.md](../README.md)。

---

## 概述

本项目使用统一的工作流来验证 C++ 程序和 Mathematica (MMA) 对积分族展开计算的一致性。

所有 MMA 脚本已通用化：通过命令行参数 `famname` 切换积分族，配置从 [FamilyDatabase.wl](../FamilyDatabase/FamilyDatabase.wl) 集中读取。**添加新族只需在数据库中注册，无需复制脚本。**

### 关键文件

| 文件 | 说明 |
|------|------|
| `VerifyExpand-Prepare.wl` | MMA: 定义家族 + 求解区域 → .bin + checkpoint (参数: famname) |
| `VerifyExpand-MMAExpand.wl` | MMA: 加载 checkpoint，执行展开计算 (参数: famname [order]) |
| `VerifyExpand-Compare.wl` | MMA: 比较 C++ vs MMA 结果 + 生成 VerifyLog.md (参数: famname) |
| `families/FamilyDatabase.wl` | 所有积分族配置的集中数据库 |
| `test_expandFF` | C++: 执行展开计算 (参数: famname) |
| `test_IBPVerification` | C++: EquationVerify 验证 IBP 方程 |
| `VerifyExpand-SeriesVerify.wl` | MMA: IBP 方程级数代入验证 (参数: famname [order]) |

---

## 测试不同积分族

C++ 和 MMA 两侧均通过命令行参数切换积分族：

```bash
# C++ 侧（无需重新编译）
cd build
./test_expandFF bub00     # 测试 bub00
./test_expandFF bub10     # 测试 bub10
./test_expandFF Tri       # 测试 Tri
./test_expandFF Box       # 测试 Box

# MMA 侧（统一脚本 + famname 参数）
wolframscript -file VerifyExpand-Prepare.wl bub00
wolframscript -file VerifyExpand-MMAExpand.wl bub00
wolframscript -file VerifyExpand-Compare.wl bub00
```

可用族名可通过以下命令查看：
```bash
wolframscript -file Compare-FamilyGenerate.wl    # 不带参数列出所有族
```

---

## 完整测试工作流

假设要测试积分族 `[famname]`：

### Step 1: Prepare — 生成积分族二进制文件 + checkpoint

```bash
cd workspace/shared/VerifyUtility
wolframscript -file VerifyExpand-Prepare.wl [famname]
```

生成 `verify/<fam>/IBPMat_[famname].bin`、`RingData_[famname].bin` 和 `PrepareCheckpoint-[famname].wdx`。

### Step 2: 运行 C++ 测试

```bash
cd /path/to/Large-Index-Expansion-CPP
./build/test_expandFF [famname]
```

生成 `Compare-CPPResult-[famname].m`。

### Step 3: MMA 展开

```bash
cd workspace/shared/VerifyUtility
wolframscript -file VerifyExpand-MMAExpand.wl [famname]
```

从 checkpoint 加载已求解的区域，执行展开。可选第二个参数指定阶数（默认 4）：
```bash
wolframscript -file VerifyExpand-MMAExpand.wl [famname] 6   # order=6
```

### Step 4: 比较 + 生成验证日志

```bash
cd workspace/shared/VerifyUtility
wolframscript -file VerifyExpand-Compare.wl [famname]
```

输出比较结果并生成 `verify/<fam>/VerifyLog-[famname].md`。

---

## 积分族配置

所有积分族的配置集中定义在 `families/FamilyDatabase.wl` 中，不再在脚本中硬编码。

每个族包含以下字段：`Propagators`、`LoopMomenta`、`ExternalMomenta`、`KinematicRules`、`TopSector`、`Numeric`（含 `"d"` 维度参数）、`Modulus`。

已注册的积分族列表见 [FamilyDatabase/README.md](../FamilyDatabase/README.md)。

---

## 添加新的积分族

只需一步：在 `families/FamilyDatabase.wl` 的 `$FamilyDatabase` 中添加新条目。

例如添加一个名为 `MyFamily` 的新族：

```mathematica
"MyFamily" -> <|
    "Description" -> "Description of this family",
    "Propagators" -> {...},      (* 传播子列表 *)
    "LoopMomenta" -> {...},       (* 圈动量 *)
    "ExternalMomenta" -> {...},   (* 独立外动量 *)
    "KinematicRules" -> {...},    (* 标量积替换规则 *)
    "TopSector" -> {...},         (* 1=主传播子, 0=辅助 *)
    "Numeric" -> {...},           (* 数值替换，含 "d" 维度参数 *)
    "Modulus" -> Prime[10000000]
|>
```

添加后，上述 4 步工作流自动适用于新族名，无需创建任何新脚本。

---

## 已知问题排查

### 1. "No matrices loaded"

检查 .bin 文件是否存在于工作目录：

```bash
ls -la IBPMat_[famname].bin
```

### 2. 比较结果不一致

- 确认 C++ 和 MMA 使用**完全相同的配置**（s, msq, d, modulus）
- 确认 `.bin` 文件是用正确的配置生成的（Step 1 和 Step 3 使用相同 famname）
- 用以下命令验证数据库中的配置：
  ```bash
  wolframscript -code 'Get["families/FamilyDatabase.wl"]; PrintFamilyInfo["bub00"]'
  ```

### 3. MMA 展开失败

检查 MMA 环境是否正常，特别是 SINGULAR 接口。

### 4. Unknown family

检查 famname 是否在 FamilyDatabase 中注册：

```bash
wolframscript -code 'Get["families/FamilyDatabase.wl"]; PrintAllFamilies[]'
```

### 5. C++ 与 MMA 不一致，但 EquationVerify 通过

这是典型的**自洽性陷阱**：EquationVerify (1b) 复用系数生成时的方程构建代码，若 bug 在共享代码中（如 `sgn()` 返回错误符号），生成器与验证器 agree on the wrong answer。此时 CompareVerify (1a) 和 SeriesVerify (1c) 会失败。

排查步骤：
1. 先用 `CompareVerify` 确认差异的具体系数和阶数
2. 在 `buildAll()` 中加 debug print，逐分量（NMinus, NZero, NPluMi, NPlus, M1, MPlus）对比期望值
3. 对有限域计算用 Python `pow(x, -1, p)` 独立验证模逆和符号

详见下方 §FFInt 符号转换 Bug 的完整调试过程。

### 6. KinematicRules 中的 dot product 规则

**问题**：SINGULAR 报错 `p1 is not defined`，导致 `#Non-triv Sectors = 0`。

**原因**：propagators 展开后包含 `p1*p2` 交叉项，但 KinematicRules 未定义这些 dot product。

**解决方案**：在 KinematicRules 中添加 dot product 规则：
```mathematica
(* 对于 Tri（3点） *)
KinematicRules -> {p1^2 -> s1, p2^2 -> s2, (p1+p2)^2 -> s3,
                   p1*p2 -> (s3 - s1 - s2)/2} /. numericRules

(* 对于 Box（4点） *)
KinematicRules -> {
    p1^2 -> s, p2^2 -> t, p3^2 -> u,
    (p1+p2)^2 -> 0, (p2+p3)^2 -> 0, (p1+p2+p3)^2 -> s+t+u,
    p1*p2 -> (0 - s - t)/2, p2*p3 -> (0 - t - u)/2} /. numericRules
```

通用公式：`p_i · p_j = ((p_i+p_j)^2 - p_i^2 - p_j^2) / 2`

### 7. SeriesVerify NB>1 实现细节

- **NB=1 路径**：A_i 为 F_p 中数值常量，h^(k) 为标量多项式
- **NB>1 路径**：A_i 属于域扩张 F_p[A_free]/(MinPoly)，通过 companion matrix 实现环运算
  - M1 = A 的 companion matrix；basisMatrices[[k+1]] = M1^k
  - ringMultiply(v1, v2) = mat(v1) · v2
  - ringInverse(v) = LinearSolve[mat(v), e1, Modulus->p]
  - 该方法适用于任意 NB ≥ 1

---

## 文件命名约定

```
VerifyExpand-Prepare.wl               # Step 1: 定义家族 + 求解区域 → .bin + checkpoint
VerifyExpand-MMAExpand.wl             # Step 2: MMA 展开计算
VerifyExpand-Compare.wl               # Step 3: 比较 + 生成验证日志
Compare-CPPResult-[famname].m         # C++ 展开结果
Compare-CPPMeta-[famname].m           # C++ 元数据 (增量、参数)
VerifyExpansion-MMAExpansion.m        # MMA 展开结果 (新命名)
Compare-RegionInfo-[famname].m        # Region 特征方程摘要
Compare-MMATiming-[famname].m         # MMA 族定义 + Region 求解耗时
Compare-MMAExpandTiming-[famname].m   # MMA 展开耗时
PrepareCheckpoint-[famname].wdx       # Prepare → MMAExpand 检查点
VerifyLog-[famname].md                # 统一验证日志 (最终输出)
IBPMat_[famname].bin                  # IBP 矩阵二进制数据
RingData_[famname].bin                # 环数据二进制数据
```

> 旧脚本 `Compare-FamilyGenerate.wl`、`Compare-Expand.wl`、`Compare-Results.wl`、`Compare-VerifyLog.wl` 仍可用。

---

## 示例：测试 bub00

```bash
# Step 1: Prepare
cd workspace/shared/VerifyUtility
wolframscript -file VerifyExpand-Prepare.wl bub00

# Step 2: C++ 展开
cd /home/ykm/Large-Index-Expansion-CPP
./build/test_expandFF bub00

# Step 3: MMA 展开
cd workspace/shared/VerifyUtility
wolframscript -file VerifyExpand-MMAExpand.wl bub00

# Step 4: 比较 + 生成日志
wolframscript -file VerifyExpand-Compare.wl bub00
```

期望输出：
```
[MATCH] k=0
[MATCH] k=1
[MATCH] k=2
[MATCH] k=3
[MATCH] k=4
[PASS] All orders match! C++ and MMA are consistent.
```

---

## 工程教训：FFInt 符号转换 Bug

### 问题描述

C++ 展开结果对所有积分族都产生错误系数。具体症状：
- **EquationVerify (自洽性检查)** 通过 — `M1*C + Total == 0`
- **CompareVerify (MMA对比)** 失败 — C++ 系数与 MMA 不匹配

例如 bub00 的 h^(1) 线性项：C++ 输出 `105689379`，正确值应为 `14952056`。

### 为什么 EquationVerify 通过但结果错误？

EquationVerify 检查 `M1*C + Total = 0`，这是**层递归算法自身的方程**。当 `sgn()` 返回错误值时：
1. 系数生成时，`buildAll()` 用错误的 sgn 值构建 Total
2. 然后用 `M1*C + Total = 0` 求解出系数
3. EquationVerify 复用同一套 `buildAll()` 逻辑 → **自洽但不正确**

**核心教训**：自洽性检查不能替代独立验证。

### 根因定位

```
Step 1 — 排除数据问题：加 debug print → .bin 与 MMA 一致 ✓
Step 2 — 定位错误方程：assembleLinearSystem 后 debug → Total[0] 偏差 4099946
Step 3 — 逐分量追踪：buildAll() 中分别打印 NMinus/NZero/NPluMi/M1/MPlus
           → MPlus[2] 贡献了 175324727（应为 0）
Step 4 — 逐项追踪 MPlus：F2s 耦合项的 factor=4099945（应为 -1 ≡ 179424672）
Step 5 — 追溯 factor：factor = sgn(1) = 4099945
Step 6 — 定位 sgn() 缺陷
```

### 根因

```cpp
// Utilities.hpp:13 (原实现)
inline T sgn(int l) {
    return static_cast<T>((l % 2 == 0) ? 1 : -1);
}
```

`FFInt` 只有 `FFInt(uint64_t)` 构造函数。`static_cast<FFInt>(-1)` 中 `int(-1)` 被隐式转换为 `uint64_t(2^64-1)`，再 `mod p = 4099945`。

### 修复

```cpp
inline T sgn(int l) {
    if (l % 2 == 0) return static_cast<T>(1);
    return -static_cast<T>(1);  // 用 operator-() 而非从 int(-1) 构造
}
```

### 经验

1. `static_cast<T>(负整数)` 对无符号构造函数的类型是陷阱 — C++ 隐式转换不会产生编译警告
2. 自洽性验证 ≠ 正确性 — 必须使用独立于算法实现的验证方法（MMA 对比、原始 IBP 方程代入）
3. 调试方法：从输出差异反向追溯，逐层缩小范围（Total → MPlus → F2s factor → sgn），每步用 Python 独立计算期望值

---

## 验证记录

| 日期 | 积分族 | 结果 | 验证的 Commits |
|------|--------|------|----------------|
| 2026-05-01 | bub00 | [PASS] k=0..4 全部匹配 | ebf89df (通用脚本化) |
| 2026-04-30 | bub00 | [PASS] k=0..4 全部匹配 | 0f57166, 204d650, e6e6121 |

**验证命令：**
```bash
./build/test_expandFF bub00
wolframscript -file Compare-Results.wl bub00
```

## 相关文档

- `../../docs/algorithms/RegionSolverAlgorithm.md` — Region Solver 算法实现比较（输出用于展开验证）
- `../../docs/algorithms/ReconstructAlgorithm.md` — MMA vs C++ 关系重构算法对比
- `../../docs/plans/RegionSolver-Debug-Plan.md` — C++ RegionSolver 调试计划
