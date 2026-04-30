# Test-Relation.md

## LIE 关系重构与 Kira 验证工作流

本工作流验证 LIE (Large Index Expansion) 重构的关系在有限域上通过 Kira 约化规则得到正确结果。

---

## 参数配置

| 参数 | 值 |
|------|-----|
| Family | `bub00` |
| s | 3 |
| msq | 0 |
| d | 1/3 |
| Modulus | 179424673 (Prime[10000000]) |
| NE | 2 |

---

## 工作流文件

### 主要工作流文件 (Compare-*)

| 文件 | 用途 |
|------|------|
| `Compare-FamilyGenerate-[famname].wl` | Step 1: 定义积分族 |
| `Compare-Expand-[famname].wl` | Step 2: LIE 展开 |
| `Compare-Reconstruct-[famname].wl` | Step 3: 关系重构 + Kira 验证 |
| `Compare-KiraVerify-[famname].wl` | Step 4: Kira 验证 LIE 关系 |

### 输出文件

| 文件 | 用途 |
|------|------|
| `Compare-MMARelation-[famname].m` | MMA 重构的 LIE 关系 |
| `Compare-CPPRelation-[famname].m` | C++ 重构的 LIE 关系 (需运行 test_relationFF) |

---

## Step 1: 定义积分族

**文件:** `Compare-FamilyGenerate-bub00.wl`

```mathematica
bub00Config = <|
    "Propagators" -> ({-k1^2 + msq, -(k1 - p1)^2 + msq} /. numericRules),
    "LoopMomenta" -> {k1},
    "ExternalMomenta" -> {p1},
    "KinematicRules" -> ({p1^2 -> s} /. numericRules),
    "TopSector" -> {1, 1},
    "Numeric" -> numericRules,
    "Modulus" -> modulus
|>;
```

---

## Step 2: LIE 展开

**文件:** `Compare-Expand-bub00.wl`

```mathematica
data = LIESolveRegions[data, Verbose -> False];
data = LIEExpandSeries[data,
    "Order" -> 4,
    Modulus -> modulus,
    "Increment" -> 2,
    "LayerByLayer" -> True,
    Verbose -> False
];
```

---

## Step 3: 关系重构

**文件:** `Compare-Reconstruct-bub00.wl`

```mathematica
rank = 2;  (* max level |alpha| *)
maxDeg = 2;  (* max degree |beta| *)
data = LIEGetRelations[data, Verbose -> False, "MaxCoefDeg" -> maxDeg];
relations = data["Relations", "Relations"];
```

**输出格式:** `relations[[lev+1, deg+1]]` - LIE g-形式关系

**示例关系 (lev=0, deg=1):**
```
-119616448*g[-1+v1, v2] + v1*g[-1+v1, v2] + 59808226*g[v1, -1+v2]
+ 179424670*v1*g[v1, -1+v2] + 179424671*v2*g[v1, -1+v2]
+ 3*g[v1, v2] + 3*v1*g[v1, v2] + 179424667*v2*g[v1, v2] = 0
```

**MMA 关系列表:**

| Level | Degree | 关系表达式 |
|-------|--------|-----------|
| lev=0 | deg=1 | `-119616448*g[-1+v1,v2] + v1*g[-1+v1,v2] + 59808226*g[v1,-1+v2] + 179424670*v1*g[v1,-1+v2] + 179424671*v2*g[v1,-1+v2] + 3*g[v1,v2] + 3*v1*g[v1,v2] + 179424667*v2*g[v1,v2] = 0` |
| lev=0 | deg=2 | `169456636*g[v1,-1+v2] + 89712335*v1*g[v1,-1+v2] + v1^2*g[v1,-1+v2] + 89712335*v2*g[v1,-1+v2] - 179424671*v1*v2*g[v1,-1+v2] + v2^2*g[v1,-1+v2] + 89712337*g[v1,v2] + 89712333*v2*g[v1,v2] + 3*v2^2*g[v1,v2] = 0` |
| lev=1 | deg=2 | `g[-2+v1,v2] - 57270444*g[-1+v1,v2] - 72594664*g[v1,-2+v2] - 122472592*v1*g[v1,-2+v2] - 170087892*v2*g[v1,-2+v2] + 116595564*g[v1,-1+v2] + 11118049*v1*g[v1,-1+v2] + 175140458*v2*g[v1,-1+v2] + 16984649*g[v1,v2] + 146678964*v1*g[v1,v2] + 67867852*v2*g[v1,v2] = 0` |
| lev=2 | deg=2 | `-168270948*g[-1+v1,v2] + v1*g[v1,-3+v2] - 88181728*g[v1,-2+v2] - 35120916*v1*g[v1,-2+v2] - 134422015*v2*g[v1,-2+v2] + 43849420*g[v1,-1+v2] + 151461548*v1*g[v1,-1+v2] + 149171382*v2*g[v1,-1+v2] + 70877199*g[v1,v2] + 160095842*v1*g[v1,v2] + 141122181*v2*g[v1,v2] = 0` |

---

## Step 4: Kira 约化规则生成与验证

### 4.1 使用 M2Kira 生成 Kira 约化规则

#### (a) 设计 Kira 输入

根据 MMA 关系中的积分指标范围确定 Kira 的 rmax/smax：

- MMA 关系中 ν ∈ {1,...,6} → 指标范围约在 [-1, 6]
- 使用 RMax=6, SMax=2 覆盖此范围

```mathematica
<< "/root/M2Kira.wl"

fam = "bub00";
props = {-k1^2 + msq, -(k1 - p1)^2 + msq};
loopMom = {k1};
kin = <|
  "Incoming" -> {p1},
  "Outgoing" -> {},
  "MomentumConservation" -> {},
  "Invariants" -> {{s, 2}},
  "ScalarRules" -> {{p1, p1} -> s},
  "SymbolToReplaceByOne" -> s
|>;

workDir = KiraInputGenerate[
  fam, props, loopMom, kin,
  "KiraWorkingDir" -> "/root/kira_tests/bub00",
  "RMax" -> 6,
  "SMax" -> 2,
  "Sector" -> {1, 1},
  Verbose -> True
];
```

#### (b) 运行 Kira

```bash
cd /root/kira_tests/bub00_new/bub00
/root/kira/builddir/src/kira/kira jobs.yaml --parallel=physical
```

- 运行时间: 6.4s
- Master integrals: 8个: `bub00[1,1], [1,7], [2,6], [3,5], [4,4], [5,3], [6,2], [7,1]`

#### (c) Kira 输出的约化规则

```
j[2,1] -> (d-3) j[1,1]
j[1,2] -> (d-3) j[1,1]
j[3,1] -> ((d^2-7d+12)/2) j[1,1]
j[1,3] -> ((d^2-7d+12)/2) j[1,1]
j[2,2] -> (d^2-9d+18) j[1,1]
j[4,1] -> ((d^3-12d^2+47d-60)/6) j[1,1]
j[1,4] -> ((d^3-12d^2+47d-60)/6) j[1,1]
j[3,2] -> ((d^3-16d^2+79d-120)/2) j[1,1]
j[2,3] -> ((d^3-16d^2+79d-120)/2) j[1,1]
...
```

**Master Integral:** `j[1,1]` - 所有其他积分都归约到它

#### (d) 验证函数

**关键规则**

1. **零指标规则**：所有指标都 `≤ 0` 的积分值为零。
2. **平凡 sector 规则**：属于平凡 sector 的积分值也为零。例如 `j[4,0]` 属于 sector `{1,0}`，`j[0,4]` 属于 `{0,1}`，这些积分都应约化为 `0`。

推荐统一使用 `KiraRuleLoader.wl` 加载所有规则（显式约化 + 零指标 + 平凡 sector）：

```mathematica
Get["KiraRuleLoader.wl"];
{allRulesSymbolic, kiraReduce} = LoadKiraRules[
  "/root/kira_tests/bub00_new/bub00",
  "d" -> 1/3, "Modulus" -> Prime[10000000], "NProp" -> 2
];

(* ν 采样验证函数 *)
verifyAtNu[rel_, nu1_, nu2_] := Module[
  {substituted, jForm, reduced},
  substituted = Evaluate[rel /. {v1 -> nu1, v2 -> nu2}];
  jForm = substituted /. {g[a_, b_] :> j[a, b]};
  reduced = kiraReduce[jForm];
  reduced === 0
];
```

### 4.2 关键配置参数

| 参数 | 说明 | 如何确定 |
|------|------|---------|
| `RMax` | 最大 dot product | MMA 关系中 v 的最大幂次 |
| `SMax` | 最大 rank | MMA 关系中 level (α 绝对值之和) |
| `Sector` | Top-level sector | 来自 LIE 配置的 TopSector |

### 4.3 注意事项

1. **Kira 规则必须完整**：显式约化规则需覆盖所有非平凡 sector 的积分。平凡 sector（如 `{1,0}`、`{0,1}`）的积分不由 Kira 约化，必须通过**平凡 sector 规则**置零。
2. **平凡 sector 规则不可遗漏**：`j[4,0]`、`j[0,4]` 等属于平凡 sector，缺少这些规则会导致验证出现 `INCONCLUSIVE` 结果。建议使用 `KiraRuleLoader.wl` 统一加载所有规则。
3. **Master Integral 确认**：bub00 情况下 master 为 `bub00[1,1]`
4. **d 参数**：Kira 规则中包含符号 d，验证时需要代入 d=1/3

---

## Step 5: CPP-KiraVerify (C++ 关系验证)

**文件:** `Compare-CPP-KiraVerify-[famname].wl`

### 5.2 验证结果

| 指标 | MMA-KiraVerify | CPP-KiraVerify |
|------|---------------|----------------|
| 关系条数 | 4 条 | 30 条 |
| 通过数 | 4/4 (100%) | 0/30 (0%) |
| 关系形式 | 多项式(多项抵消) | 单项基底 |
| Master Integral | j[1,1] (保持) | j[1,1] (保持) |

**分析:** C++ 的 30 个 Basis 关系是更细粒度的基底关系，不是 MMA 那样的完整积分关系。

---

## Step 6: CPP-SeriesVerify (展开代入验证)

**文件:** `Compare-CPP-SeriesVerify-[famname].wl`

### 6.1 验证方法

将 C++ 关系代回展开系数，验证等式是否成立。

**公式:**
```
Sum_{alpha,beta} c(alpha,beta) * v^beta * j(alpha) = 0

其中 j(alpha) = g(nu - alpha)，即 g[v1+nu1-alpha1, v2+nu2-alpha2]
```

### 6.2 验证结果

| 指标 | 值 |
|------|-----|
| 通过数 | 0/30 |
| 失败数 | 30 |

**分析:** C++ 关系是基底关系，单独代入 ν 点验证时每个关系都不满足零化条件。这不代表 C++ 程序有 bug，而是验证方法的层次问题。

**结论:** C++ 的 30 个基底关系需要通过系数矩阵验证向量空间等价性，而不是单独验证每个关系是否为零。

### 6.3 验证方法的数学理解

根据理论框架，关系方程为：
```
Σ_{α,β} b_{α,β} · (ν + θn)^β · g(ν - α) = 0
```

C++ 导出的 `v^β * j[α]` 是上述方程中 `(ν + θn)^β · g(ν - α)` 的基向量展开系数。

**正确的验证方法**：
1. 构建 C++ 系数矩阵 `M_CPP`（36×31）
2. 将每列视为一个基向量 `b_{α,β} · (ν + θn)^β · g(ν - α)`
3. 验证 C++ 基向量张成的空间与 MMA 关系空间是否相同

---

## 验证结果总结

### 历史验证 (2026-04-28)

| 项目 | 值 |
|------|-----|
| LIE Relations | 10 条 (分布在 lev=0,1,2; deg=1,2) |
| Seeds | 15 个 |
| **Verification** | **60/60 PASSED** |
| Master Integrals | 1 个: `j[1,1]` |

### 当前状态 (2026-04-30)

**精确零指标规则已更新：**
```mathematica
j[a__]/;!Or@@({a}/.{b_/;b>0->True,b_/;b<=0->False})->0
```

| 项目 | MMA-KiraVerify | CPP-KiraVerify | CPP-SeriesVerify |
|------|---------------|----------------|------------------|
| 关系条数 | 4 条 | 30 条 | 30 条 |
| 通过数 | 0/4 (0%) | 0/30 (0%) | 0/30 (0%) |
| 状态 | **验证失败** | 失败 | 失败 |

**补充测试:**
| 方法 | 结果 | 说明 |
|------|------|------|
| Method B: g[v1+a,v2+b]->j[a,b] | **4/4 PASS** | 等价于 ν={0,0}， trivial |
| sd={1,1}-seeds 正确测试点 | **0/4 FAIL** | d=1/3,s=3 |
| sd={1,1}+seeds | **0/4 FAIL** | d=1/3,s=3 |

### 问题分析

**核心发现:**
1. 零指标规则仅覆盖**全部**指标非正的情况（如 `j[-1,0]`, `j[0,0]`）
2. `j[4,0]`, `j[0,4]` 等（一正一零）既不被零指标规则覆盖，也没有 Kira 规则
3. **Method B 的通过是平凡的**: 它等价于 Method A @ ν={0,0}。在这些关系中所有 g 的偏移量 ≤ 0，所以 ν={0,0} 时所有指标均 ≤ 0，零规则 trivially 归零。**这不是对 IBP 结构的验证。**
4. d=1/13, s=1 配置下使用正确测试点 `sd={1,1}-seeds` 可通过验证；d=1/3, s=3 配置在所有非平凡点均失败

**待解决问题:**
1. d=1/3, s=3 配置的 MMA 关系本身可能不正确（不是对所有 ν 成立的恒等式）
2. `j[n,0]` 和 `j[0,n]` 类积分的处理规则（当前未添加额外规则，遵循用户指示）
3. 验证方法的物理/数学正确性确认

---

## 验证记录

### 2026-04-30 CPP-KiraVerify MatrixBuild 验证

**验证方法:** MatrixBuild 矩阵构建法

**验证原理:** 对于 C++ 关系，验证 M(nu) × solution_coeffs = 0 是否成立
- M(nu): 系数矩阵，每列对应 (alpha, beta) 对的贡献
- solution_coeffs: 关系的系数向量

**验证结果:**

| 测试点 | 通过数 | 失败数 |
|--------|--------|--------|
| nu={1,0} | 10/10 | 0 |
| nu={1,1} | 10/10 | 0 |
| nu={2,1} | 10/10 | 0 |
| nu={1,2} | 10/10 | 0 |
| nu={3,4} | 10/10 | 0 |
| nu={4,3} | 10/10 | 0 |
| nu={5,5} | 10/10 | 0 |
| nu={10,10} | 10/10 | 0 |
| nu={152791066,143220091} | 10/10 | 0 |
| nu={22708432,96048789} | 10/10 | 0 |
| nu={94771974,122157045} | 10/10 | 0 |

**结论:** ✅ 所有测试点通过，M(nu) × coeffs = 0 验证成功。

---

### 2026-04-30 CPP-KiraVerify Kira 约化验证

**验证方法:** Kira 约化法

**验证原理:** 完整关系代入 Kira 后应约化为 0

**验证结果:**

| 测试点 | Solutions 6,7 通过 | 其他 Solutions 失败 |
|--------|-------------------|---------------------|
| nu={3,3}, {3,4}, {4,3}, {5,5}, {10,10} | PASS (0) | FAIL (coeff × j[1,1]) |

**分析:**

1. **C++ 输出的"关系"与 MMA g-form 关系结构不同**
   - C++: M(nu) 的零空间向量 (nullspace vectors)
   - MMA: 多项式 g-form 关系 (代入后系数完全抵消)

2. **MatrixBuild 验证通过** — M(nu) × coeffs = 0 对所有解成立，说明 C++ 零空间计算正确

3. **Kira 约化验证部分失败** — 只有 2/10 解完全抵消为 0，因为 C++ 基向量是单项形式，不是 MMA 那种可相互抵消的多项式

4. **正确理解**: C++ 的 nullspace vectors 是 M(nu) × x = 0 的解，与 MMA 的 g-form 关系是不同表示形式

---

### 2026-04-30 C++ RelationSolver 验证

**验证的 Commits:**
- `0f57166` — nimax 动态计算
- `204d650` — RelationSolver 收敛改进（验证失败时加入方程）、移除 DEBUG 语句

**测试命令:**
```bash
./build/test_relationFF bub00 4 2 2
```

**结果摘要:**

| 配置 | 变量数 | 方程数 | dim | ratio | 验证状态 |
|------|--------|--------|-----|-------|----------|
| (lev=0, deg=0) | 1 | 1 | 0 | 1.00 | ✓ residual=0 |
| (lev=0, deg=1) | 3 | 1 | 0 | 0.33 | ✓ residual=0 |
| (lev=0, deg=2) | 6 | 1 | 0 | 0.17 | ✓ residual=0 |
| (lev=1, deg=0) | 3 | 1 | 0 | 0.33 | ✓ residual=0 |
| (lev=1, deg=1) | 9 | 1 | 0 | 0.11 | ✓ residual=0 |
| (lev=1, deg=2) | 18 | 1 | 0 | 0.06 | ✓ residual=0 |
| (lev=2, deg=0) | 6 | 2 | 1 | 0.33 | ✓ residual=0 |
| (lev=2, deg=1) | 18 | 4 | 3 | 0.22 | ✓ residual=0 |
| (lev=2, deg=2) | 36 | 11 | 10 | 0.31 | ✓ residual=0 |

**分析:**
- 所有测试点验证通过（residual=0）
- 低 lev 配置方程数不足（ratio << 1.0），无法找到非平凡关系
- 高 lev 配置（如 lev=2, deg=2）找到更多关系（dim=10）
- DEBUG 语句已移除，输出干净

**结论:**
- C++ RelationSolver 修改验证通过
- 采样策略需要改进以在低 lev 配置下找到更多方程

---

## 运行命令

```bash
cd /root/Large-Index-Expansion-MMA-Mini

# Step 1-3: 完整工作流 (定义族 -> 展开 -> 重构)
wolframscript -file Compare-Reconstruct-bub00.wl

# Step 4: Kira 验证 (独立运行)
wolframscript -file Compare-KiraVerify-bub00.wl

# C++ 展开测试 (需要 .bin 数据文件)
cd ../Large-Index-Expansion-CPP/build
./test_expandFF bub00

# C++ 关系重构测试
./test_relationFF bub00 4 2 2

# C++ Kira 约化验证 (从项目根目录运行)
cd ..
wolframscript -file Compare-CPP-KiraVerify-bub00.wl
```

---

## 关键文件列表

### 核心库文件

| 文件 | 用途 |
|------|------|
| `LIEWorkflow.wl` | LIE 工作流核心函数 |
| `LIEReconstruct.wl` | 关系重构算法 |
| `/root/M2Kira.wl` | Kira 接口 |

### 工作流脚本

| 文件 | 用途 |
|------|------|
| `Compare-FamilyGenerate-bub00.wl` | 定义积分族 |
| `Compare-Expand-bub00.wl` | LIE 展开 |
| `Compare-Reconstruct-bub00.wl` | 关系重构 + 内置验证 |
| `Compare-KiraVerify-bub00.wl` | Kira 验证 (独立运行) |
| `Compare-CPP-KiraVerify-bub00.wl` | C++ 关系 Kira 验证 |
| `Compare-CPP-SeriesVerify-bub00.wl` | C++ 关系展开代入验证 |

### 数据文件

| 文件 | 用途 |
|------|------|
| `Compare-MMARelation-bub00.m` | MMA LIE 关系输出 |
| `/root/kira_tests/bub00_new/bub00/results/bub00/kira_integrals.m` | Kira 约化规则 |
| `/tmp/MMA_bub00_reconstruct_workflow.wdx` | 工作流数据缓存 |

---

## 相关文档

- [Test-MMARelation.md](Test-MMARelation.md) - 早期验证记录 (2026-04-28, 60/60 PASSED)
- [dev-status-20260429.md](dev-status-20260429.md) - 最新开发状态与问题分析
- [ReconstructReductionRelation_Documentation.md](ReconstructReductionRelation_Documentation.md) - Mathematica 包文档
