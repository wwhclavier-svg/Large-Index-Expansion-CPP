# Test-Expand.md

C++ 和 MMA 展开计算对比测试工作流文档。

---

## 概述

本项目使用统一的工作流来验证 C++ 程序和 Mathematica (MMA) 对积分族展开计算的一致性。

### 关键文件

| 文件 | 说明 |
|------|------|
| `Compare-FamilyGenerate-[famname].wl` | MMA: 生成积分族的 .bin 文件 |
| `Compare-Expand-[famname].wl` | MMA: 执行展开计算 |
| `Compare-Results-[famname].wl` | MMA: 比较 C++ 和 MMA 结果 |
| `test_expandFF.cpp` | C++: 执行展开计算 |
| `test_IBPVerification` | C++: EquationVerify 验证 IBP 方程 |
| `Compare-IBPVerify-[famname].wl` | MMA: SeriesVerify 级数自洽性验证 |
| `IBPVerification.md` | 三种验证方法详细文档 |

---

## 测试不同积分族

### 方式一：命令行参数（无需重新编译）

```bash
cd build
make test_expandFF
./test_expandFF [famname]
```

例如：
```bash
./test_expandFF bub00   # 测试 bub00
./test_expandFF bub     # 测试 bub
./test_expandFF Tri     # 测试 Tri
```

### 方式二：修改源码重新编译

修改 `tests/test_expandFF.cpp` 中的 `famname` 变量后重新编译。

---

## 完整测试工作流

假设要测试积分族 `[famname]`：

### Step 1: 生成积分族二进制文件

在 MMA 目录下运行：

```bash
cd /path/to/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-FamilyGenerate-[famname].wl
```

这会生成：
- `IBPMat_[famname].bin`
- `RingData_[famname].bin`

并自动复制到 `Large-Index-Expansion-CPP/` 和 `Large-Index-Expansion-CPP/build/` 目录。

### Step 2: 运行 C++ 测试

```bash
cd /path/to/Large-Index-Expansion-CPP/build
make test_expandFF
./test_expandFF [famname]
```

这会生成 `Compare-CPPResult-[famname].m` 文件。

### Step 3: 运行 MMA 展开

```bash
cd /path/to/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-Expand-[famname].wl
```

这会生成 `Compare-MMAResult-[famname].m` 文件。

### Step 4: 比较结果

```bash
cd /path/to/Large-Index-Expansion-CPP
# 如果 MMA 结果文件在 MMA 目录，先复制过来
cp ../Large-Index-Expansion-MMA-Mini/Compare-MMAResult-[famname].m .

wolframscript -file Compare-Results-[famname].wl
```

如果所有 k=0,1,2,3,4 都显示 `[MATCH]`，则测试通过。

---

## 积分族配置说明

每个积分族的配置（s, msq, d 值）定义在生成脚本中：

| famname | Propagators | s | msq | d |
|---------|-------------||------|-----|
| bub00 | {-k1^2+msq, -(k1-p1)^2+msq} | 3 | 0 | 1/3 |
| bub | {-k1^2+msq, -(k1-p1)^2+msq} | 3 | 1 | 1/3 |
| Tri | {-k1^2+msq, -(k1-p1)^2+msq, -(k1-p1-p2)^2+msq} | 3 | 0 | 1/3 |
| Box | {-k1^2+msq, -(k1-p1)^2+msq, -(k1-p1-p2)^2+msq, -(k1-p1-p2-p3)^2+msq} | s=3,t=2 | 0 | 1/3 |

---

## 添加新的积分族

### 1. 创建生成脚本

在 MMA 目录创建 `Compare-FamilyGenerate-[famname].wl`：

```mathematica
(* Compare-FamilyGenerate-[famname].wl *)
modulus = Prime[10000000];  (* 179424673 *)
numericRules = {s -> 3, msq -> 0, "d" -> 1/3};  (* 根据需要修改 *)

config = <|
    "Propagators" -> ({-k1^2 + msq, -(k1 - p1)^2 + msq} /. numericRules),
    "LoopMomenta" -> {k1},
    "ExternalMomenta" -> {p1},
    "KinematicRules" -> ({p1^2 -> s} /. numericRules),
    "TopSector" -> {1, 1},
    "Numeric" -> numericRules,
    "Modulus" -> modulus
|>;
```

### 2. 创建展开脚本

复制并修改 `Compare-Expand-[famname].wl`，使用相同的配置。

### 3. 创建比较脚本

复制并修改 `Compare-Results-[famname].wl`。

### 4. 运行测试

按照上面的工作流执行。

---

## 已知问题排查

### 1. "No matrices loaded"

检查 .bin 文件是否存在于当前工作目录：

```bash
ls -la IBPMat_[famname].bin
```

### 2. 比较结果不一致

- 确认 C++ 和 MMA 使用**完全相同的配置**（s, msq, d, modulus）
- 确认 MMA 展开脚本和生成脚本的配置一致
- 确认 .bin 文件是用正确的配置生成的

### 3. MMA 展开失败

检查 MMA 环境是否正常，特别是 SINGULAR 接口。

---

## 文件命名约定

```
Compare-FamilyGenerate-[famname].wl  # 生成积分族 .bin 文件
Compare-Expand-[famname].wl         # MMA 展开计算
Compare-Results-[famname].wl         # 比较结果
Compare-CPPResult-[famname].m        # C++ 展开结果
Compare-MMAResult-[famname].m        # MMA 展开结果
```

---

## 示例：测试 bub00

```bash
# Step 1: 生成 .bin 文件
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-FamilyGenerate-bub00.wl

# Step 2: 运行 C++ 测试
cd /root/Large-Index-Expansion-CPP/build
./test_expandFF bub00

# Step 3: 运行 MMA 展开
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file Compare-Expand-bub00.wl

# Step 4: 比较结果
cd /root/Large-Index-Expansion-CPP
cp ../Large-Index-Expansion-MMA-Mini/Compare-MMAResult-bub00.m .
wolframscript -file Compare-Results-bub00.wl
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
| 2026-04-30 | bub00 | [PASS] k=0..4 全部匹配 | 0f57166, 204d650, e6e6121 |

**验证命令：**
```bash
./build/test_expandFF bub00
wolframscript -file Compare-Results-bub00.wl
```

**验证内容：**
- `nimax` 动态计算修复 (0f57166)
- LayerRecursionCore BINOM 溢出修复 (204d650)
- RelationSolver 收敛改进 (204d650)
- LinearSolver_FF FFInt 验证修复 (204d650)
