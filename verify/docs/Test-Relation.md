# Test-Relation.md

## LIE 关系重构与 Kira 验证工作流

本工作流验证 LIE (Large Index Expansion) 重构的关系在有限域上通过 Kira 约化规则得到正确结果。

> **理论基础与算法**：关系方程的数学推导、MMA/C++ 实现对比见 **[Reconstruct_Algorithm.md](Reconstruct_Algorithm.md)**。

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

> **切换测试例**：修改以上参数。C++ 侧通过 CLI 传参：`./build/test_relationFF <family> <order> <lev_max> <deg_max>`。MMA 侧需修改对应脚本（`Compare-FamilyGenerate-<family>.wl` 等）中的 `numericRules` 和 `modulus` 变量。

---

## 四条并行验证方法

关系重构后，可通过以下四种独立方法验证正确性。推荐运行顺序：先跑 **2a (MatrixBuild)** 确认 C++ 零空间计算自洽，再跑 **2d (NuVerify)** 用 Kira 外部验证。

| 方法 | 验证内容 | 独立性 | 状态 |
|------|---------|:---:|:---:|
| **2a. MatrixBuild** | M(ν) × coeffs = 0，C++ 零空间自洽 | 纯 C++ | **110/110 PASS** |
| **2b. CPP-KiraVerify** | C++ 关系代入 Kira 约化 → 0 | 需 Kira 规则 | **11/11 PASS** |
| **2c. CPP-SeriesVerify** | C++ 关系代入展开系数 | 需展开结果 | **待处理** (见 §6) |
| **2d. NuVerify / ν 采样** | 任意 ν 点 Kira 约化验证 | 需 Kira 规则 | 见下文 |

各方法的详细流程见下方对应章节。

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
| `Compare-CPPRelation-[famname].m` | C++ 重构的 LIE 关系 (位于 `verify/results/[famname]/`, 需运行 test_relationFF) |

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

> **算法背景**：C++ 的数值 ν 采样与 MMA 的符号系数提取的差异分析见 **[Reconstruct_Algorithm.md §4](Reconstruct_Algorithm.md)**。

### 5.2 验证结果

| 指标 | MMA-KiraVerify | CPP-KiraVerify |
|------|---------------|----------------|
| 关系条数 | 4 条 | 30 条 |
| 通过数 | 4/4 (100%) | 0/30 (0%) |
| 关系形式 | 多项式(多项抵消) | 单项基底 |
| Master Integral | j[1,1] (保持) | j[1,1] (保持) |

**分析:** C++ 的 30 个 Basis 关系是更细粒度的基底关系，不是 MMA 那样的完整积分关系。

---

## Step 6: CPP-SeriesVerify (展开代入验证) — 待处理

**文件:** `verify/VerifyUtility/VerifyRelation-SeriesVerify.wl`

### 6.1 验证方法

**原理**（反向思维，对应 VerifyExpand-SeriesVerify）：

C++ 关系方程：
$$\sum_{\alpha,\beta} b_{\alpha,\beta} \cdot (\nu + \theta n)^\beta \cdot g(\nu-\alpha; n) = 0$$

将 $g(\nu-\alpha; n) = p(\alpha) \cdot \sum_k h_k(\nu-\alpha) / n^k$ 代入，提取 $1/n^m$ 系数：

$$\sum_{\alpha,\beta} b_{\alpha,\beta} \sum_{t=0}^{|\beta|} c_{\beta,t}(\nu) \cdot p(\alpha) \cdot h_{m+t}(\nu-\alpha) \equiv 0 \pmod{p}$$

其中：
- $c_{\beta,t}(\nu)$ = $(\nu+\theta n)^\beta$ 展开中 $n^t$ 的系数
- $p(\alpha) = \prod_i A_{\text{inv},i}^{\alpha_i}$
- $h_k$ = 展开系数多项式（从 `Compare-CPPResult-*.m` 加载）

**运行命令：**
```bash
cd verify/VerifyUtility
wolframscript -file VerifyRelation-SeriesVerify.wl <famname> [lev] [deg] [order]
```

### 6.2 当前状态：待处理

**验证结果：** m=0 (n^0 系数) 通过，但 m≥1 全部失败。

**根因分析（2026-05-05）：**

问题不在于公式错误，而在于 C++ 求解器的 **order truncation**：

1. C++ `AdaptiveEquationBuilder` 对方程按 $n$ 的幂次做 **阶数稳定性分析** (`analyzeOrderStability`)
2. 对于 (lev=1, deg=1)，`stable_order=0`，求解器只强制 $n^1$ (r=0) 和 $n^0$ (r=1) 的系数为零
3. $n^{-1}$ 及以下的方程**未被纳入求解系统**（它们对应的零空间尚未稳定）

直接数值验证（重建 C++ 矩阵 M，计算 M·b）证实：
```
r=0 (n^1):  0 ✓
r=1 (n^0):  0 ✓  
r=2 (n^{-1}): 44856167 ✗  ← 未被求解器强制
r=3 (n^{-2}): 71022264 ✗
...
```

**结论：** SeriesVerify 的数学公式是正确的，但需要求解器在更高阶达成稳定后才能通过完整验证。当前 stable_order=0 意味着只有 $n^1$ 和 $n^0$ 被约束。

**后续工作：**
1. [ ] 使 SeriesVerify 脚本感知 `stable_order` 边界，只验证 solver 实际约束的阶数
2. [ ] 或：扩展 C++ 求解器以在更高阶收敛（需要更多 ν 采样点或更高 k_max）
3. [ ] 或：从 MMA `LIEReconstruct` 路径导出完整符号关系进行对比

### 6.3 与 NuVerify 的互补关系

| 方法 | 验证对象 | 数据需求 | 当前状态 |
|------|---------|---------|---------|
| **NuVerify** (§2d) | 完整关系 (ν 固定, Kira 约化) | Kira 规则 | **PASS** (外部验证) |
| **SeriesVerify** (§2c) | 逐阶系数 (ν 任意, 展开系数) | 展开 h_k | **待处理** (内部验证) |

两种方法互补：NuVerify 是外部 Kira 验证（不需要展开数据），SeriesVerify 是内部展开自洽验证（不需要 Kira）。

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
| 关系条数 | 4 条 | 11 条 | 待处理 |
| 通过数 | 0/4 (0%) | **11/11 (100%)** | 仅 n^0 阶通过 |
| 状态 | 失败（格式 bug） | **PASS** | **待处理**（order truncation，见 §6） |

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

**MMA vs C++ 差异分析：**

| 方面 | MMA 关系 | C++ 关系 |
|------|---------|---------|
| 性质 | **完整关系**：多个项互相抵消后为 0 | **基底关系**：M(nu) 零空间的基向量 |
| 方程构建 | 符号：提取 ν^power 系数 | 数值：多 ν 点采样 M(nu)·x=0 |
| 包含负索引 | 是（被零规则直接置零） | 否（只有非负索引） |
| MatrixBuild 验证 | N/A (符号方法) | 110/110 通过 |
| Kira 约化验证 | 0/4 失败 (d=1/3,s=3) | 2/10 部分抵消 |

C++ 基底关系需通过向量空间等价性验证，而非单独代入 ν 点。

**已解决 (2026-05-02):**
1. ~~d=1/3, s=3 配置的 MMA 关系本身可能不正确~~ — C++ 关系通过 KiraVerify 验证正确
2. ~~`j[n,0]` 和 `j[0,n]` 类积分的处理规则~~ — Kira 规则已覆盖 bub00 的所有相关索引

**当前状态:** C++ RelationSolver 的所有非平凡关系（11条）均通过 Kira 物理规则验证，数学上正确。

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
| nu={3,3}, {3,4}, {4,3}, {5,5}, {10,10} | PASS (0) | FAIL (coeff × master integral) |

**分析:**

1. **C++ 输出的"关系"与 MMA g-form 关系结构不同**
   - C++: M(nu) 的零空间向量 (nullspace vectors)
   - MMA: 多项式 g-form 关系 (代入后系数完全抵消)

2. **MatrixBuild 验证通过** — M(nu) × coeffs = 0 对所有解成立，说明 C++ 零空间计算正确

3. **Kira 约化验证部分失败** — 只有 2/10 解完全抵消为 0，因为 C++ 基向量是单项形式，不是 MMA 那种可相互抵消的多项式

4. **正确理解**: C++ 的 nullspace vectors 是 M(nu) × x = 0 的解，与 MMA 的 g-form 关系是不同表示形式

**注意 (2026-05-02 修正):**
之前的"修复"方向完全错了：- **错误**：将 `g[v1-α1,v2-α2]` 格式改为 `j[α]` 格式- **正确**：`g[v1-α1,v2-α2]` 格式是正确的（Apr 30 ebf89df 引入），因为它代入 ν 值后变成 $j(\nu - \alpha)$，而 `j[α]` 变成 $j(\alpha)$ 缺少 ν 部分。现已恢复为正确的 `g[v1-α1,v2-α2]` 格式。

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

- **[Reconstruct_Algorithm.md](Reconstruct_Algorithm.md)** — 关系重构算法：MMA vs C++ 实现对比、数学框架
- [Test-MMARelation.md](Test-MMARelation.md) - 早期验证记录 (2026-04-28, 60/60 PASSED)
- [dev-status-20260429.md](dev-status-20260429.md) - 最新开发状态与问题分析
- [ReconstructReductionRelation_Documentation.md](ReconstructReductionRelation_Documentation.md) - Mathematica 包文档

---

## 更新 (2026-05-02): C++ KiraVerify 修复后 11/11 PASS

### 修复的问题

**Bug 1: `j[...]` 格式错误**
- 文件: `Compare-Reconstruct-bub00.wl`
- 症状: KiraVerify 0/10 FAIL（j[{a,b}] 格式不匹配 Kira 的 j[a,b]）
- 根因: `extractGCoefficientsFromRelation` 产生 `j[{alpha, beta}]`（列表参数），但 KiraRuleLoader 的 `kiraRules` 用的是 `j[a, b]`（标量参数）
- 修复: `j[{#2, #3}]` → `j[#2, #3]`

**Bug 2: `relationSpecific` 未定义**
- 文件: `Compare-Reconstruct-bub00.wl`
- 症状: `kiraRelSpecific = relationSpecific["bub00", ...]` 返回 `{}`（空列表）
- 根因: `relationSpecific` 函数不存在
- 修复: 删除该行，改用 KiraRuleLoader 的 `kiraReduce`

### 当前验证结果

| 方法 | 结果 | 说明 |
|------|------|------|
| MatrixBuild (C++) | 110/110 PASS | C++ 零空间自洽 |
| **KiraVerify (C++ Relations)** | **11/11 PASS** | lev=2 deg=0/1/2 非平凡关系全部通过 |
| CompareVerify | 5/5 MATCH | C++ vs MMA 系数一致 |

**KiraVerify 详细结果 (bub00 C++ Relations)**:
```
Lev=2 Deg=0: 1 nontrivial relations -> 1 PASS, 0 FAIL
Lev=2 Deg=1: 3 nontrivial relations -> 3 PASS, 0 FAIL
Lev=2 Deg=2: 7 nontrivial relations -> 7 PASS, 0 FAIL
OVERALL: 11 PASS, 0 FAIL
```

### 关键发现

C++ RelationSolver 输出的关系是**数学上正确的**（通过 Kira 物理规则验证）。之前 KiraVerify FAIL 是 MMA 脚本中的格式 bug，不是 C++ 关系的问题。

### 快速验证命令

```bash
# 独立诊断脚本（不依赖完整 MMA 工作流）
cd /root/Large-Index-Expansion-MMA-Mini
wolframscript -file ../Large-Index-Expansion-CPP/verify/tests/full_kira_verify.m
```
