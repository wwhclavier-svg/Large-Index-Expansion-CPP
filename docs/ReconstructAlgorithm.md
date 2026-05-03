# LIE 关系重构算法：MMA vs C++ 实现对比

## 1. 基本数学框架

### 1.1 核心展开公式

IBP 系数在  $\underline{\tilde{\nu}}=\underline{\nu} + \underline{\theta} n$  处关于 $n$ 的展开：
   
$$j(\underline{\theta} n + \underline{\nu}) = h(\underline{\nu}; n) = \sum_{k=0}^{\infty} \sum_{\alpha_i \leq 0, |\underline{\alpha}|\leq incre*k} c_{k,\underline{\alpha}} \underline{\nu}^{\underline{\alpha}} n^{-k}$$
   

其中：
- $j(\tilde{\underline{\nu}}) = j(\underline{\theta} n + \underline{\nu})=h(\underline{\nu};n)$ 是费曼积分在变量替换 $\underline{\nu} \to \underline{\nu} + \underline{\theta} n$ 后对 $n \to \infty$ 的级数展开.
- $j(\nu)$ 表示原本定义的费曼积分，$g(\nu;n)$ 表示 $n$ -平移之后的费曼积分，$h(\nu;n)$表示 $g$ 对 $n$ 的渐进展开.
- $\underline{\nu} = (\nu_1, \nu_2, ..., \nu_e)$  是费曼积分的指标向量
-  $\underline{\theta} = (\theta_1, \theta_2, ..., \theta_e)$  是 sector 相关的 $\theta$ 向量
-  $n$  是渐进展开参数
-  $c_{k,\underline{\alpha}}$  是展开系数

### 1.2 关系方程的一般形式

在渐进展开框架下，假设不同指标 $\nu$ 处费曼积分满足关系方程为：

$$
\sum_{\underline{\alpha},\underline{\beta}} b_{\underline{\alpha},\underline{\beta}} \underline{\tilde{\nu}}^{\underline{\beta}}  j(\underline{\tilde{\nu}} - \underline{\alpha}) = 
\sum_{\underline{\alpha},\underline{\beta}} b_{\underline{\alpha},\underline{\beta}} (\underline{\theta} n + \underline{\nu})^{\underline{\beta}}  g(\underline{\nu} - \underline{\alpha};n) = 0 
$$


其中：
- $b_{\underline{\alpha},\underline{\beta}}$ 是待求的 relation 系数
- $(\underline{\theta} n + \underline{\nu})^{\underline{\beta}} = \prod_i (\nu_i + \theta_i n)^{\beta_i}$ 是带平移的多项式

### 1.3 关键变换

将 $g(\underline{\nu} - \underline{\alpha};n)$ 用级数展开代入：

$$ 
g(\underline{\nu} - \underline{\alpha};n) = \sum_{k,\underline{\gamma}} c_{k,\underline{\gamma}}  (\underline{\nu}-\underline{\alpha})^{\underline{\gamma}} / n^k
$$


因此：
$$
\sum_{\underline{\alpha},\underline{\beta},\underline{\gamma},k} b_{\underline{\alpha},\underline{\beta}} \, c_{k,\underline{\gamma}}  (\underline{\theta} n + \underline{\nu})^{\underline{\beta}}  (\underline{\nu}-\underline{\alpha})^{\underline{\gamma}} / n^k = 0
$$

这是一个关于 $\nu$ 的多项式方程。我们可以从中提取中单项式方程，或直接取 $\nu$ 的不同数值采样得到不同的标量方程，最后解出 $b_{\alpha,\beta}$.

---

## 2. MMA 实现

### 2.1 核心流程

```
LIEReconstruct[rankLevel, degree, order, hexpnList, aregList, topsec, vlist]
```

### 2.2 关键步骤

#### Step 1: GenerateAnsatz - 创建 ansatz 多项式

```mathematica
GenerateAnsatz[mode_, rank_, maxDeg_, limitSector_, nVars_, vlist_] := Module[
  {seeding, fullAnsatz, levelList, getPoly, monomials},

  (* 生成种子 *)
  seeding[n_, lev_] := Join @@ Table[
    SortBy[Flatten[Permutations /@ (IntegerPartitions[l + n, {n}] - 1), 1], -Reverse[#] &],
    {l, 0, lev}
  ];
  seedList = seeding[nVars, maxDeg];
  monomials[d_] := (Times @@ Thread@Power[vlist, #]) & /@ seeding[nVars, d];

  (* 构造多项式: Sum b[rk, power] * v^power *)
  getPoly[rk_, deg_] := Dot[("b"[rk, #] & /@ seeding[nVars, deg]), monomials[deg]];

  (* Pyramid mode: g[v - rk] -> b[rk,deg] * v^deg *)
  base = Association@Table[
    "g" @@ (vlist - rk) -> getPoly[rk, maxDeg],
    {rk, rkG}
  ];
  ...
];
```

**结果**：ansatz 是 $g(\underline{\nu} - \underline{\alpha}) = \sum b_{\underline{\alpha}, \underline{\beta}}  \underline{\nu}^{\underline{\beta}}$ 形式的符号多项式

#### Step 2: ComputeBasisPower - 计算 p(α) = A^(-α)

```mathematica
ComputeBasisPower[alphaList_, areg_, modulus_] := Module[
  {powerTable, getPower},

  (* 递归计算: B^0 = I *)
  powerTable[ConstantArray[0, ne]] = IdentityMatrix[nb];

  (* idx[i] > 0: B[i] * B^(idx - e_i) *)
  (* idx[i] < 0: A[i] * B^(idx + e_i) *)
  getPower[idx_] := Module[{i, res},
    i = FirstNonZero[idx];
    If[idx[[i]] > 0,
      res = varBMatrix[i] . getPower[idx - unit],
      res = varAMatrix[i] . getPower[idx + unit]
    ];
    powerTable[idx] = res
  ];
  ...
];
```

**结果**：对每个 α，计算 $p(\underline{\alpha}) = \underline{A}^{-\underline{\alpha}}$ 矩阵

#### Step 3: UpdateExpansionTable - 计算 g(v - α)

```mathematica
UpdateExpansionTable[expnTable_, basePower_, hexpnList_, indexNew_, vlist_, maxOrder_, modulus_] := Module[
  {expnNew},

  expnNew = Table[
    alpha -> Dot[
      basePower[[j]][alpha],
      hSeries[[i]] /. Dispatch@Thread[vlist -> (vlist - alpha)]
    ],
    {alpha, indexNew}
  ];
  ...
];
```

**关键**：`Thread[vlist -> (vlist - alpha)]` 替换后，$ h^{(i)}(\underline{\nu} - \underline{\alpha})$ 仍是 $\underline{\nu}$ 的**符号多项式**

**结果**：`expntable[reg, sol, k][α]` 是 $\underline{\nu}$ 的符号多项式

#### Step 4: SetupEquationsAll - 构建方程组

```mathematica
SetupEquationsAll[coeftable_, expntable_, alphaList_, ...] := Module[
  {relansatz},

  (* 构建关系 ansatz: 系数 × 展开式 的卷积 *)
  relansatz = Table[
    Sum[
      coeftable[[j, kk + 1, l]] * expntable[[j, m, k - kk + 1]][alphaList[[l]]],
      {l, Length[alphaList]}, {kk, 0, Min[k, deg]}
    ],
    {m, Length[expntable[[j]]]}
  , {j, nreg}];

  (* 提取方程: 获取 v^power 的系数 *)
  eqs = CoefficientEquations[Flatten[relansatz /. numeric], vlist];
  ...
];
```

**关键**：
- `coeftable[[j, kk+1, l]]` 是 `(n+v)^β` 展开后 `n^k` 的系数（数值）
- `expntable[[j, m, k-kk+1]][alphaList[[l]]]` 是 `g(v - α)`（符号多项式）
- `relansatz` 合成后仍是 `v` 的符号多项式

#### Step 5: CoefficientEquations - 提取系数得到方程

```mathematica
CoefficientEquations[expr_, vars_] := Module[
  {eqs, varsSym, exprConverted, convRule},

  (* 转换为符号列表 *)
  varsSym = Table[Unique["v"], {Length[vars]}];
  convRule = Thread[vars -> varsSym];
  exprConverted = expr /. convRule;

  (* 提取系数数组 *)
  eqs = CoefficientArrays[exprConverted, varsSym];

  (* 转换为方程列表 *)
  eqs = Join @@ (If[Head[#] === SparseArray, Most@ArrayRules[#][[;;,2]], {#}] & /@ eqs);
  eqs // DeleteCases[#, 0] &
];
```

**结果**：`v¹, v², ...` 每个幂次对应一个方程

#### Step 6: Solve - 解线性系统

```mathematica
ReduceSolve[eqs_, vars_, p_] := Module[
  {mat, dep, indep, rule},

  (* 系数矩阵 *)
  mat = CoefficientArrays[eqs, varRev];
  mat = Join[mat[[2]], List /@ mat[[1]], 2];

  (* 行简化 *)
  mat = If[p =!= 0, RowReduce[mat, Modulus -> p], RowReduce[mat]];

  (* 提取依赖/独立变量 *)
  dep = FirstNonZero /@ mat;
  indep = Complement[Range[Length[vars]], dep];

  (* 回代得到解 *)
  rule = Thread[varRev[[dep]] -> (-mat[[;;, indep]] . varRev[[indep]] - mat[[;;, -1]])];
  rule
];
```

### 2.3 MMA 算法总结

| 步骤 | 输入 | 输出 | 性质 |
|------|------|------|------|
| GenerateAnsatz | rank, deg | `g(v-α) = Σ b_{α,β}·v^β` | 符号 |
| ComputeBasisPower | α list | `p(α) = A^(-α)` | 数值矩阵 |
| UpdateExpansionTable | `h^(i)(v)` | `g(v-α) = p(α) * h^(i)(v-α)` | 符号 |
| SetupEquationsAll | coeftable, expntable | `relansatz = Σ coeff·g(v-α)` | 符号 |
| CoefficientEquations | 符号表达式 | 方程组 `M · b = 0` | 数值 |
| ReduceSolve | 方程组 | 系数 `b_{α,β}` | 数值 |

---

## 3. C++ 实现

### 3.1 核心流程

```
RelationSolver::reconstructReductionRelation<FFInt>(...)
```

### 3.2 关键步骤

#### Step 1: RegimeData::prepare - 预计算 p(α)

```cpp
template<typename T>
void RegimeData<T>::prepare(int lev, const std::vector<std::vector<int>>& alphas) {
    int ne = theta.size();
    p_store = std::make_unique<IndexStorage<T>>(ne, nb * nb);

    // 递归计算所有 p(alpha)
    for (const auto& alpha : alphas) {
        computePRecursive(alpha, alphas, A_ops, A_inv_ops);
    }
    is_prepared = true;
}
```

**结果**：与 MMA 相同，`p(α) = A^(-α)` 矩阵

#### Step 2: step2_computeG - 计算 g(nu - alpha)

```cpp
template<typename T>
void RegimeEvaluator<T>::step2_computeG(const std::vector<T>& nu, int nimax_idx) {
    for (const auto& alpha : alphas_) {
        T* p_ptr = reg_->getP(alpha);
        std::vector<T> g_val(nb_ * (k_max_ + 1), T(0));

        for (int k = 0; k <= k_max_; ++k) {
            std::vector<T> h_jk(nb_, T(0));

            for (int l = 0; l <= k; ++l) {
                long long states = BINOM[l + ne_ - 1][ne_ - 1];
                for (long long cid = 0; cid < states; ++cid) {
                    std::vector<int> gamma = readIndex(static_cast<int>(cid), l, ne_);

                    // weight = (nu - alpha)^gamma
                    T weight = T(1);
                    for (int m = 0; m < ne_; ++m) {
                        T diff = nu[m] - static_cast<T>(alpha[m]);
                        weight *= detail::power(diff, gamma[m]);
                    }

                    // 累加 h_jk
                    for (int j = 0; j < nb_; ++j) {
                        h_jk[j] += (*C_)(k, l, static_cast<int>(cid), j, nimax_idx) * weight;
                    }
                }
            }

            // 矩阵乘: g_i = sum_j p_ij * h_j
            for (int i = 0; i < nb_; ++i) {
                for (int j = 0; j < nb_; ++j) {
                    g_val[i * (k_max_ + 1) + k] += p_ptr[i * nb_ + j] * h_jk[j];
                }
            }
        }
        g_store_->insert(alpha, g_val);
    }
}
```

**关键**：`nu[m] - alpha[m]` 是**数值**，不是符号

**结果**：`g_store_[alpha]` 是 `nb × (k_max+1)` 的**数值数组**

#### Step 3: step3_computeF1 - 计算 (theta + nu)^beta

```cpp
template<typename T>
void RegimeEvaluator<T>::step3_computeF1(const std::vector<T>& nu) {
    for (const auto& alpha : alphas_) {
        for (const auto& beta : betas_) {
            std::vector<T> poly(k_max_ + 1, T(0));
            poly[0] = T(1);

            for (int i = 0; i < ne_; ++i) {
                int b = beta[i];
                if (b == 0) continue;

                std::vector<T> term_poly(b + 1, T(0));
                T th = static_cast<T>(theta_[i]);
                T n_val = nu[i];  // <-- 采样点数值

                // 二项式展开: (theta + nu)^b = sum C(b,m) * theta^(b-m) * nu^m
                for (int m = 0; m <= b; ++m) {
                    T coef = static_cast<T>(BINOM[b][m]) *
                             detail::power(th, b - m) *
                             detail::power(n_val, m);
                    term_poly[m] = coef;
                }
                // 多项式乘法: poly = poly * term_poly
                ...
            }
            f1_store_->insert(alpha, poly, beta);
        }
    }
}
```

**关键**：`n_val = nu[i]` 是**采样点数值**

**结果**：`f1_store_[alpha,beta]` 是 `(k_max+1)` 的**数值数组**，表示多项式 `(theta+nu)^beta` 在各阶 k 的系数

#### Step 4: step4_computeF2 - 计算 f2 = f1 * g 卷积

```cpp
template<typename T>
void RegimeEvaluator<T>::step4_computeF2(const std::vector<int>& alpha,
                                         const std::vector<int>& beta) {
    T* g_ptr = g_store_->retrieve(alpha);
    T* f1_ptr = f1_store_->retrieve(alpha, beta);

    int f2_len = k_max_ + 1;
    std::vector<T> f2_val(nb_ * f2_len, T(0));

    for (int i = 0; i < nb_; ++i) {
        for (int l = 0; l <= k_max_; ++l) {
            T f1_val = f1_ptr[l];
            if (f1_val == T(0)) continue;

            for (int k = 0; k <= k_max_ - l; ++k) {
                T g_val = g_ptr[i * f2_len + k];
                if (g_val == T(0)) continue;

                int target_pow = l + k;
                f2_val[i * f2_len + target_pow] += f1_val * g_val;
            }
        }
    }
    f2_store_->insert(alpha, f2_val, beta);
}
```

**结果**：`f2_store_[alpha,beta]` 是 `nb × (k_max+1)` 的**数值数组**

#### Step 5: buildFinalMatrix - 组装最终矩阵

```cpp
template<typename T>
std::vector<std::vector<T>> RegimeEvaluator<T>::buildFinalMatrix() {
    // 遍历所有 (alpha, beta) 组合
    for (size_t a = 0; a < alphas_.size(); ++a) {
        for (size_t b = 0; b < betas_.size(); ++b) {
            T* f2_ptr = f2_store_->retrieve(alphas_[a], betas_[b]);
            if (!f2_ptr) continue;

            // 将 f2_val 按 k 阶展开到矩阵行
            for (int k = 0; k <= k_max_; ++k) {
                for (int i = 0; i < nb_; ++i) {
                    size_t row = ...;
                    size_t col = a * betas_.size() + b;
                    matrix[row][col] = f2_ptr[i * (k_max_+1) + k];
                }
            }
        }
    }
    return matrix;
}
```

**结果**：得到 `M(nu) · x = 0` 的系数矩阵 `M`

#### Step 6: AdaptiveEquationBuilder - 自适应采样求解

```cpp
BuildResult AdaptiveEquationBuilder<T>::build(
    const std::vector<RegimeData<T>>& regimes,
    const std::vector<std::vector<int>>& nimax_lists,
    int ne)
{
    // 主循环：自适应采样
    while (result.nu_count < config_.max_nu) {
        // 获取下一个采样点
        std::vector<T> nu = sampler.next();

        // 在该点评估所有 regime 和 nimax
        auto rows = assembler.evaluateAtNu(nu);

        // 添加到求解器
        solver.addRows(rows);

        // 检查收敛
        if (sampler.shouldCheck()) {
            int current_nullity = solver.getNullity();
            sampler.update(current_nullity);

            if (sampler.hasConverged()) {
                result.converged = true;
                result.nullspace = solver.getNullspace();
                return result;
            }
        }
    }
}
```

**关键**：需要在**多个采样点** `nu` 上评估方程，然后求解

### 3.3 C++ 算法总结

| 步骤 | 输入 | 输出 | 性质 |
|------|------|------|------|
| prepare | α list | `p(α) = A^(-α)` | 数值矩阵 |
| step2_computeG | nu (采样点) | `g(nu-α)` | 数值 |
| step3_computeF1 | nu (采样点) | `(theta+nu)^beta` 系数 | 数值 |
| step4_computeF2 | g, f1 | `f2 = f1 * g` | 数值 |
| buildFinalMatrix | f2 | `M(nu) · x = 0` | 数值 |
| build() | 多个 nu | `M · x = 0` 方程组 | 数值 |

---

## 4. MMA vs C++ 核心差异

### 4.1 符号 vs 数值

| 方面 | MMA | C++ |
|------|-----|-----|
| `g(v-α)` | 符号多项式 `P(v)` | 数值 `g(nu-α)` |
| `(theta+v)^beta` | 符号多项式 | 数值 `(theta+nu)^beta` |
| 方程来源 | `CoefficientEquations[P(v), v]` | `M(nu) · x = 0` 在多个 nu 点 |

### 4.2 方程构建方式

**MMA**:
```
P(v) = Σ b_{α,β} · v^β · g(v-α)  [符号多项式]
eqs = CoefficientEquations[P(v), v]  →  每个 v^power 一个方程
```

**C++**:
```
M(nu) = Σ b_{α,β} · (theta+nu)^β · g(nu-α)  [标量]
在多个 nu 点采样 → M1·x=0, M2·x=0, ...
```

### 4.3 关键公式对比

**MMA 的系数提取**:
```
coeftable[[j, kk+1, l]]  -- (n+v)^beta 展开后 n^{kk} 的系数
expntable[[j, m, k-kk+1]][alphaList[[l]]]  -- g(v-alpha) 符号多项式
relansatz = Sum_{l,kk} coeftable * expntable  -- v 的符号多项式
eqs = CoefficientEquations[relansatz, v]  -- 提取 v^power 系数
```

**C++ 的采样求解**:
```
f1[k] = coefficient of (theta+nu)^beta at order k  -- 数值
g[k] = g(nu-alpha) at order k  -- 数值
f2[k] = Sum_{l} f1[l] * g[k-l]  -- 数值卷积
M(nu)[row] = f2  -- 标量
在 nu1, nu2, ... 上采样得到 M1, M2, ...
求解 M · x = 0
```

---

## 5. 验证方法

### 5.1 ν 采样验证框架

根据理论框架（1.3节），关系方程为：
$$
\sum_{\underline{\alpha},\underline{\beta},\underline{\gamma},k} b_{\underline{\alpha},\underline{\beta}} \, c_{k,\underline{\gamma}}  (\underline{\theta} n + \underline{\nu})^{\underline{\beta}}  \underline{\nu}^{\underline{\gamma}} / n^k = 0
$$

**ν 采样验证流程**：

对于关系 `rel = Σ b_{α,β} · v^β · g(v-α)`：

1. **代入 ν 值**：将符号变量 v1,v2 替换为具体的 ν 点值 (ν1, ν2)
   ```
   substituted = rel /. {v1 -> ν1, v2 -> ν2}
   ```
   例如：`g[-1+v1,v2]` 变为 `g[ν1-1, ν2]`

2. **Kira 约化**：对所有出现的 g[...] 应用 Kira 约化规则
   - 负索引 `g[负数, *]` 或 `g[*, 负数]` → 0
   - 正索引按 Kira 规则约化到主积分 `j[1,1]`

3. **取模验证**：代入 d=1/3 并对素数取模
   ```
   result = PolynomialMod[reduced /. d -> 1/3, modulus]
   result = 0 表示关系在该 ν 点成立
   ```

### 5.2 MMA 关系 ν 采样验证

**MMA 关系格式**：`rel1 = -119616448*g[-1+v1,v2] + v1*g[-1+v1,v2] + ... + 3*g[v1,v2] + ...`

**验证步骤**：
```mathematica
(* 1. 对每个 ν 点代入 *)
For ν = (ν1, ν2):
  substituted = rel1 /. {v1 -> ν1, v2 -> ν2}
  (* 例: g[-1+v1,v2] -> g[ν1-1, ν2] *)

(* 2. Kira 约化（负索引置零，正索引约化到 j[1,1]） *)
reduced = substituted //. {
  g[a_?Negative, _] :> 0,
  g[_, b_?Negative] :> 0,
  g[a_, b_] /; a >= 0 && b >= 0 :> KiraReduce[g[a, b]]
};

(* 3. 代入 d=1/3 并取模 *)
result = PolynomialMod[reduced /. d -> 1/3, modulus];
```

**注意**：MMA 关系中的 `v1*g[-1+v1,v2]` 项：
- 代入 ν 后变成 `ν1 * g[ν1-1, ν2]`
- Kira 约化：`g[ν1-1, ν2]` 中若 `ν1-1 < 0` 则为 0
- 所以负索引项对任意 ν 都为 0

**MMA 验证结果**（4/4 通过）：因为 MMA 关系包含负索引项，这些项对任意 ν 值都为 0。

### 5.3 C++ 关系验证

**C++ 关系格式**：`rel = Σ coeff · v^β · j[α]`，其中 `j[α] = g(v-α)`

**验证步骤**：
```mathematica
(* 1. 解析关系，提取 coeff, β, α *)
terms = parseRelation[rel];  (* {coeff, {α1,α2}, {β1,β2}] *)

(* 2. 对每个 ν 点计算 *)
For each ν = (ν1, ν2):
  value = Sum_{terms} coeff · ν1^β1 · ν2^β2 · g(ν1-α1, ν2-α2)

  (* 3. Kira 约化 *)
  value = KiraReduce[value]

  (* 4. 验证 *)
  If value ≠ 0: 关系在 ν 点不成立
```

**C++ 验证结果**（0/30 失败）：C++ 导出的 `j[α]` 主要是正索引（如 `j[0,0]`, `j[1,1]`），这些不是主积分，会被 Kira 规则约化到 `j[1,1]`。单独验证每个基底关系时，它们不满足零化条件。

### 5.4 差异原因分析

| 方面 | MMA 关系 | C++ 关系 |
|------|---------|---------|
| 性质 | **完整关系**：多个项互相抵消后为 0 | **基底关系**：描述个别项的行为 |
| 包含负索引 | 是（被 Kira 直接置零） | 否（只有非负索引） |
| 单独验证 | 通过（负索引项抵消） | 失败（正索引项不归零） |

**注意**：C++ 的 30 个基底关系需要一起理解，它们张成的向量空间才对应 MMA 的完整关系。

### 5.5 正确的 C++ 关系验证方法

验证 C++ 关系应该检查：
1. **向量空间等价性**：C++ 基底张成的空间与 MMA 关系空间是否相同
2. **代入验证**：将 C++ 关系代入原始方程（包含 (ν+θn)^β 因子），验证是否满足

当前验证失败是因为单独验证每个基底关系，而不是验证完整的方程 `Σ b_{α,β} · (ν+θn)^β · g(ν-α) = 0`。

---

## 6. 待解决问题

1. **采样点选取**：C++ 的自适应采样是否收敛到正确解？
2. **方程完整性**：在多个 nu 点采样是否等价于 MMA 的系数提取？
3. **数值精度**：有限域算术是否存在奇异性？

---

## 7. 参考文件

### MMA
- `LIEReconstruct.wl` - 主重构算法
- `LIEWorkflow.wl` - 工作流包装
- `Compare-Reconstruct-bub00.wl` - 测试脚本

### C++
- `RelationSolver.hpp` - 核心求解器
- `test_relationFF.cpp` - 测试入口
- `Compare-CPPRelation-bub00.m` - 输出格式
