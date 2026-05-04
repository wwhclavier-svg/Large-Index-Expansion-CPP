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

> 参数配置与工作流见 [verify/README.md](../README.md)。三种方法的原理详解见 [IBPVerification.md](IBPVerification.md)。

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
| `verify/FamilyDatabase/FamilyDatabase.wl` | 所有积分族配置的集中数据库 |
| `test_expandFF` | C++: 执行展开计算 (参数: famname) |
| `test_IBPVerification` | C++: EquationVerify 验证 IBP 方程 |
| `VerifyExpand-SeriesVerify.wl` | MMA: IBP 方程级数代入验证 (参数: famname [order]) |

> 旧脚本 `Compare-FamilyGenerate.wl`、`Compare-Expand.wl`、`Compare-Results.wl`、`Compare-VerifyLog.wl` 仍可用。

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
cd /path/to/Large-Index-Expansion-MMA-Mini
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
cd /path/to/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-MMAExpand.wl [famname]
```

从 checkpoint 加载已求解的区域，执行展开。可选第二个参数指定阶数（默认 4）：
```bash
wolframscript -file VerifyExpand-MMAExpand.wl [famname] 6   # order=6
```

### Step 4: 比较 + 生成验证日志

```bash
cd /path/to/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-Compare.wl [famname]
```

输出比较结果并生成 `verify/<fam>/VerifyLog-[famname].md`。

---

## 积分族配置

所有积分族的配置集中定义在 `verify/FamilyDatabase/FamilyDatabase.wl` 中，不再在脚本中硬编码。

每个族包含以下字段：`Propagators`、`LoopMomenta`、`ExternalMomenta`、`KinematicRules`、`TopSector`、`Numeric`（含 `"d"` 维度参数）、`Modulus`。

已注册的积分族列表见 [FamilyDatabase/README.md](../FamilyDatabase/README.md)。

---

## 添加新的积分族

只需一步：在 `verify/FamilyDatabase/FamilyDatabase.wl` 的 `$FamilyDatabase` 中添加新条目。

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
  wolframscript -code 'Get["verify/FamilyDatabase/FamilyDatabase.wl"]; PrintFamilyInfo["bub00"]'
  ```

### 3. MMA 展开失败

检查 MMA 环境是否正常，特别是 SINGULAR 接口。

### 4. Unknown family

检查 famname 是否在 FamilyDatabase 中注册：

```bash
wolframscript -code 'Get["verify/FamilyDatabase/FamilyDatabase.wl"]; PrintAllFamilies[]'
```

### 5. C++ 与 MMA 不一致，但 EquationVerify 通过

这是典型的**自洽性陷阱**：EquationVerify (1b) 复用系数生成时的方程构建代码，若 bug 在共享代码中（如 `sgn()` 返回错误符号），生成器与验证器 agree on the wrong answer。此时 CompareVerify (1a) 和 SeriesVerify (1c) 会失败。

排查步骤：
1. 先用 `CompareVerify` 确认差异的具体系数和阶数
2. 在 `buildAll()` 中加 debug print，逐分量（NMinus, NZero, NPluMi, NPlus, M1, MPlus）对比期望值
3. 对有限域计算用 Python `pow(x, -1, p)` 独立验证模逆和符号

详见 [IBPVerification.md §FFInt 符号转换 Bug](IBPVerification.md#重要教训ffint-符号转换-bug) 的完整调试过程。

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
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-Prepare.wl bub00

# Step 2: C++ 展开
cd /root/Large-Index-Expansion-CPP
./build/test_expandFF bub00

# Step 3: MMA 展开
cd /root/Large-Index-Expansion-MMA-Mini
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
