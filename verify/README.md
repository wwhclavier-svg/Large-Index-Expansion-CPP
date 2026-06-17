# Verification Workflow

C++ 程序 与 Mathematica (MMA) 程序交叉验证展开系数和线性关系重构的正确性。

## 前置检查：环境烟雾测试 (test_ff_verify)

在运行任何验证之前，先确认 FireFly 有限域库在当前环境下正常工作：

```bash
./build/test_ff_verify
```

**期望输出**：`✓ FireFly is working correctly!`

### 何时需要运行

| 情境 | 说明 |
|------|------|
| **首次使用** | 确认 FireFly 编译链接正确 |
| **环境变更后** | 系统更新、GMP/FireFly 重编译、编译器升级后 |
| **切换素数模数后** | 验证新模数下基本运算（+、-、*、/、模逆）正确 |
| **调试有限域问题时** | 其他 FF 测试（test_expandFF、test_relationFF）失败时，先跑此测试排除 FireFly 本身的问题 |
| **CI/CD 流水线** | 作为构建后的冒烟测试，快速判定构建环境是否健康 |

### 测试内容

对模数 p=2147483647 (2^31-1) 执行加减乘除及模逆运算，验证 `(b/a) * a == b`。该模数独立于实际计算所用模数，专门测试 FireFly 基础功能。

## 验证工作流全景

```
                          ┌─ 1a. CompareVerify ─── C++ vs MMA 展开系数逐阶 diff
                          │
  展开验证 (Test-Expand) ─┼─ 1b. EquationVerify ─── IBP 递归方程自洽性 M1*C+Total==0
                          │
                          └─ 1c. SeriesVerify ──── IBP 方程级数代入 (ν→ν+θn)
                          
                          │
                          ├─ 2a. MatrixBuild ───── M(ν)×coeffs=0 (110/110 通过)
                          │
  关系验证 (Test-Relation)┼─ 2b. CPP-KiraVerify ─── C++ 关系→Kira 约化→0
                          │
                          ├─ 2c. CPP-SeriesVerify ─ 关系代入展开系数验证
                          │
                          └─ 2d. NuVerify ──────── ν 采样 Kira 验证
```

## 工作流参数（所有测试的入口）

修改以下变量以切换测试例。**同一轮验证中确保所有脚本使用一致的参数。**

> **积分族定义数据库**：[verify/FamilyDatabase/](FamilyDatabase/) — 所有族的传播子、圈动量、外动量、运动学替换规则、数值参数集中在此。脚本可通过 `Get["verify/FamilyDatabase/FamilyDatabase.wl"]` 加载，避免参数散落。

| 参数 | 含义 | 当前值 | 设置位置 |
|------|------|--------|---------|
| `famname` | 积分族名称 | `bub00` | C++ CLI arg / MMA: `$FamilyDatabase` 键名 |
| `modulus` | 有限域特征 p 或模数 | `179424673` (Prime[10^7]) | `$FamilyDatabase[fam]["Modulus"]` |
| `dVal` | 时空维度参数 | `1/3` | `$FamilyDatabase[fam]["Numeric"]` 中的 `"d"` |
| `s`, `t`, `u` 等 | 运动学不变量 | 见各族定义 | `$FamilyDatabase[fam]["Numeric"]` |
| `order` | 展开阶数上限 | `4` | C++ CLI / MMA: `"Order"->4` |
| `lev_max` | max \|α\|待重构关系中积分指标的范围 | `2` | C++ CLI / MMA: `"MaxCoefDeg"` |
| `deg_max` | max \|β\|待重构关系中系数作为指标多项式的阶数 | `2` | C++ CLI / MMA: `"MaxCoefDeg"` |

> **设计缺陷（已修复）**：MMA 脚本已通用化并重构为 `VerifyExpand-Prepare.wl`、`VerifyExpand-MMAExpand.wl`、`VerifyExpand-Compare.wl` 三步管线，共享工具集中在 Prepare 的 Utils section。旧 per-family 脚本保留向后兼容。

## 展开验证 (Test-Expand)

三种平行方法，可独立运行，也可组合使用。

### 1a. CompareVerify — C++ vs MMA 交叉对比

两个独立实现（C++ 层递归 + MMA LIEExpand）跑相同输入，逐阶 diff 展开系数。

```bash
# Step 1: Prepare — 生成 .bin 文件 + checkpoint
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-Prepare.wl bub00

# Step 2: C++ 展开 → Compare-CPPResult-bub00.m
cd /root/Large-Index-Expansion-CPP
./build/test_expandFF bub00

# Step 3: MMA 展开 → VerifyExpansion-MMAExpansion.m
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-MMAExpand.wl bub00

# Step 4: 对比 + 生成 VerifyLog → verify/bub00/VerifyLog-bub00.md
wolframscript -file VerifyExpand-Compare.wl bub00
```

期望: `[MATCH] k=0..4` → `[PASS]`

> 旧脚本 `Compare-FamilyGenerate.wl`、`Compare-Expand.wl`、`Compare-Results.wl` 仍可用，推荐使用新的 `VerifyExpand-*.wl` 三步管线。

### 1b. EquationVerify — IBP 递归方程自洽性 (C++)

不依赖 MMA。验证每个 `(k, l, seed)` 的展开系数满足 `M1 * C + Total == 0`。

```bash
./build/test_IBPVerification bub00
```

### 1c. SeriesVerify — IBP 方程级数代入 (MMA)

将展开系数代入平移后的 IBP 方程，验证在给定截断阶以下各阶 1/n 系数 ≡ 0。

```bash
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file VerifyExpand-SeriesVerify.wl bub00
```

### 展开验证状态

| 积分族 | CompareVerify | EquationVerify | SeriesVerify |
|--------|:---:|:---:|:---:|
| bub00 | k=0..4 MATCH | PASS | — (待测试新脚本) |
| bub | k=0..4 MATCH | PASS | — |
| Tri | k=0..4 MATCH | PASS | — |
| Box | k=0..4 MATCH | PASS | — |

## 关系验证 (Test-Relation)

### 2a. MatrixBuild — 矩阵零空间自洽 (C++)

在多个 ν 点构造系数矩阵 M(ν)，验证 `M(ν) × coeffs = 0`。**110/110 通过**，证明 C++ 零空间计算正确。

```bash
./build/test_relationFF bub00 4 2 2
```

### 2b. CPP-KiraVerify — Kira 约化验证

将 C++ 关系代入 Kira IBP 约化规则验证归零。2/10 部分抵消 — C++ 输出的是零空间**基底向量**（单项形式），不是 MMA 的多项式完整关系。

### 2c. CPP-SeriesVerify — 展开代入

C++ 关系代入展开系数验证等式。对基底关系单独不适用 — 需要线性组合成完整关系。

### 2d. NuVerify — ν 采样 Kira 验证

在任意 ν 点验证 C++ 关系是否被 Kira 约化为 0。

```bash
wolframscript -file verify/scripts/NuVerify-Relations.wl
```

加载 `verify/results/bub00/Compare-CPPRelation-*.m`。

### 关系验证状态

| 方法 | 结果 | 说明 |
|------|------|------|
| MatrixBuild | **110/110 PASS** | 零空间计算正确 |
| Kira 约化 (单项基底) | 2/10 部分抵消 | 基底 vs 完整关系的差异 |
| MMA Kira (d=1/3,s=3) | 0/4 FAIL | 当前调试中 |
| 历史 (d=1/13,s=1) | **60/60 PASS** | 2026-04-28 |

## 注意事项 / 经验教训

以下来自通用化改造和 bub00 验证过程中踩过的坑。

### 1. FamilyDatabase 与 per-family 脚本的配置一致性

**FamilyDatabase 创建晚于 per-family 脚本，两者存在配置漂移（Configuration Drift）。**

- 脚本习惯将 `KinematicRules` 预替换为数值（如 `{p1^2 -> s} /. {s -> 3}` → `{p1^2 -> 3}`），再传给 `LIEDefineFamily`
- 数据库保持符号形式（`{p1^2 -> s}` + `Numeric` 中 `s -> 3`），依赖 `LIEDefineFamily` 内部 `ibpeqsfull/.numeric` 做替换
- **两种方式经全流程验证结果等价**，但不保证所有族都如此。新族加入数据库后，必须跑一轮 diff 对比
- 验证方法：写脚本分别用两种配置跑 `LIEDefineFamily → LIESolveRegions → LIEExpandSeries`，比较 HList

### 2. wolframscript 参数传递格式

`$ScriptCommandLine` 的格式是 `{scriptPath, arg1, arg2, ...}`，**不包含** `wolframscript` 或 `-file` 前缀。

```mathematica
(* 正确 *)
famname = $ScriptCommandLine[[2]];  (* 第一个参数 *)

(* 错误 — 不存在 wolframscript, -file *)
famname = $ScriptCommandLine[[4]];
```

验证方法：`wolframscript -file test.wl arg1 2>&1` 打印 `$ScriptCommandLine` 确认。

### 3. 不要对 wolframscript 输出用 aggressive grep

Piping wolframscript 输出到 `grep` 可能匹配不到任何行，让失败看起来像成功（空输出无法区分 "没匹配" 和 "脚本崩溃"）。

```bash
# 错误 — 失败不可见
wolframscript -file script.wl 2>&1 | grep -E '^ExpectedPattern'

# 正确 — 先捕获再筛选
wolframscript -file script.wl > /tmp/out.wl 2>&1
grep 'ExpectedPattern' /tmp/out.wl
```

或者至少先检查退出码：`wolframscript ... || echo "FAILED with $?"`。

### 4. SINGULAR 版本兼容性警告

SINGULAR 4.3.2 对空变量列表的 `ring r = p, (), (lp(0),C)` 语法报 error，但仅影响平凡 sector（如 `{1,0}`、`{0,1}`），顶 sector 不受影响。这些报错可以忽略。

### 5. FFInt 类型安全：禁止 `static_cast<FFInt>(负整数)`

`FFInt` 仅有 `uint64_t` 构造函数。`static_cast<FFInt>(-1)` 会导致 `int(-1)` 隐式转换为 `uint64_t(2^64-1)` 再 mod p，产生垃圾值（如 p=179424673 时得 4099945），**无编译警告**。正确写法是 `-FFInt(1)`。

此 bug 曾导致所有积分族展开系数错误，而 EquationVerify (1b) 因复用同一错误逻辑仍通过。详见 [IBPVerification.md](docs/IBPVerification.md#重要教训ffint-符号转换-bug)。

### 6. 自洽性检查不可替代独立验证

EquationVerify (1b) 验证的是 `M1*C + Total == 0`，即**层递归自身的方程**。当 bug 在方程构建代码中时（如上述 sgn bug），生成器和验证器使用相同的错误逻辑 → 自洽但结果错误。必须搭配 CompareVerify (1a) 或 SeriesVerify (1c) 等独立方法。详见 [IBPVerification.md](docs/IBPVerification.md#为什么-equationverify-通过但结果错误)。

## 目录结构

```
verify/
├── README.md                    # 本文件 — 验证工作流总览
├── FamilyDatabase/
│   ├── README.md                # 积分族数据库文档
│   └── FamilyDatabase.wl        # 所有族的集中定义 (MMA 可加载)
├── docs/
│   ├── Test-Expand.md           # 展开验证详情
│   ├── Test-Relation.md         # 关系验证详情
│   ├── IBPVerification.md       # 三种验证方法原理
│   ├── Verify-MMA-KIRA-Guide.md # Kira 交叉验证指南
│   └── CPP-KiraVerify-Debugger.md
├── scripts/
│   └── NuVerify-Relations.wl    # ν 采样验证 C++ 关系
└── results/                     # 已验证结果快照 (.m)
    └── bub00/
        ├── Compare-CPPRelation-bub00_lev*_deg*.m
        └── Relations_bub00_lev*_deg*.m
```

## 中间文件分布

| 文件模式 | 生成者 | 位置 |
|---------|--------|------|
| `IBPMat_*.bin`, `RingData_*.bin` | MMA `VerifyExpand-Prepare.wl <fam>` | `verify/<fam>/` |
| `Compare-CPPResult-*.m` | C++ `test_expandFF <fam>` | CWD (通常项目根) |
| `VerifyExpansion-MMAExpansion.m` | MMA `VerifyExpand-MMAExpand.wl <fam>` | `verify/<fam>/` |
| `Compare-CPPRelation-*_lev*_deg*.m` | C++ `test_relationFF` | CWD |
| `Relations_*_lev*_deg*.m` | C++ `test_relationFF` | CWD |
| `ExpansionCache_*.bin` | C++ `test_relationFF` | CWD |
| `resCache_Expansion_*.bin` | C++ `test_expandFF` | CWD |
| `PrepareCheckpoint-*.wdx` | MMA `VerifyExpand-Prepare.wl <fam>` | `verify/<fam>/` |

> C++ 测试用相对路径写 CWD。从 `build/` 运行时输出落 `build/`。`verify/results/` 存已验证快照供 MMA 脚本读取。

## 添加新积分族

1. 在 `verify/FamilyDatabase/FamilyDatabase.wl` 的 `$FamilyDatabase` 中添加条目
2. `wolframscript -file VerifyExpand-Prepare.wl <family>` 生成 `.bin` + checkpoint
3. `./build/test_expandFF <family>` + `./build/test_relationFF <family>` 跑 C++ 测试
4. `wolframscript -file VerifyExpand-MMAExpand.wl <family>` MMA 展开
5. `wolframscript -file VerifyExpand-Compare.wl <family>` 对比 + 日志
6. 更新本文件中的验证状态表

> **重要**：新族加入 FamilyDatabase 后，必须先跑一轮全流程比对（与旧 per-family 脚本生成的结果 diff），确认配置等价。特别是 KinematicRules —— 脚本侧习惯预替换为数值，数据库侧保持符号，虽经 `LIEDefineFamily` 内部 `ibpeqsfull/.numeric` 补偿后等价，但需逐族验证。
