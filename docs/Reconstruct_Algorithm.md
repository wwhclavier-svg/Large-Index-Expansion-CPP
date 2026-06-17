# LIE 关系重构算法：MMA vs C++ 实现对比

## 1. 基本数学框架

### 1.1 核心展开公式

在 LIE (Large Index Expansion) 中我们考虑费曼积分 $j(\underline{\tilde{\nu}})$ 在 $\underline{\tilde{\nu}}=\underline{\nu} + \underline{\theta} n$  处关于 $n$ 的展开：
   
$$j(\underline{\theta} n + \underline{\nu}) 
= g(\underline{\nu}; n) = \underline{A}^{\underline{\nu}}\sum_{k=0}^{\infty} \left( \sum_{\alpha_i \leq 0, |\underline{\alpha}|\leq \mathrm{incre} \times k} c_{k,\underline{\alpha}} \underline{\nu}^{\underline{\alpha}} \right) n^{-k}$$
   

其中：
- $g(\underline{\nu};n) = j(\underline{\theta} n + \underline{\nu})$ 是费曼积分 $j(\underline{\nu})$ 在变量替换 $\underline{\nu} \to \underline{\nu} + \underline{\theta} n$ 后对 $n \to \infty$ 的级数展开.
- $j(\nu)$ 表示原本定义的费曼积分，$g(\nu;n)$ 表示 $n$ -平移之后的费曼积分，$h_k(\nu;n)$ 表示 $g$ 对 $n^{-k}$ 渐进展开系数的多项式部分.
- $\underline{\nu} = (\nu_1, \nu_2, ..., \nu_e)$  是费曼积分的指标向量
-  $\underline{\theta} = (\theta_1, \theta_2, ..., \theta_e)$  是 sector 相关的 $\theta$ 向量
- 向量的向量幂次 (如 $\underline{\nu}^{\underline{\alpha}}$) 表示相应分量的分量幂次乘积 (如$\nu_1^{\alpha_1} \ldots \nu_n^{\alpha_n}$).
-  $n$  是渐进展开参数
-  $c_{k,\underline{\alpha}}$ 是 $h_{k}$ 数值系数

### 1.2 关系方程的一般形式

在渐进展开框架下，假设不同指标 $\nu$ 处费曼积分满足关系方程为：

$$
\sum_{\underline{\alpha},\underline{\beta}} b_{\underline{\alpha},\underline{\beta}} \underline{\nu}^{\underline{\beta}}  j(\underline{\nu} - \underline{\alpha}) = 0 $$
 
平移并展开得到的式子是级数展开需要满足的：
$$
\sum_{\underline{\alpha},\underline{\beta}} b_{\underline{\alpha},\underline{\beta}} (\underline{\theta} n + \underline{\nu})^{\underline{\beta}}  g(\underline{\nu} - \underline{\alpha};n) = 0 
$$


其中：
- $b_{\underline{\alpha},\underline{\beta}}$ 是待求的 relation 系数
- $(\underline{\theta} n + \underline{\nu})^{\underline{\beta}} = \prod_i (\nu_i + \theta_i n)^{\beta_i}$ 是带平移的多项式

### 1.3 关键变换

将 $g(\underline{\nu} - \underline{\alpha};n)$ 用级数展开代入：

$$ 
g(\underline{\nu} - \underline{\alpha};n) = \underline{A}^{\underline{\nu}-\underline{\alpha}}\sum_{k,\underline{\gamma}} c_{k,\underline{\gamma}}  (\underline{\nu}-\underline{\alpha})^{\underline{\gamma}} / n^k
$$


因此：
$$
\sum_{\underline{\alpha},\underline{\beta},\underline{\gamma},k} b_{\underline{\alpha},\underline{\beta}} \,\underline{A}^{-\underline{\alpha}} c_{k,\underline{\gamma}}  (\underline{\theta} n + \underline{\nu})^{\underline{\beta}}  (\underline{\nu}-\underline{\alpha})^{\underline{\gamma}} / n^k = 0
$$

这是一个关于 $\underline{\nu}$, $n$ 和 $1/n$ 的多项式方程。我们可以从中提取中单项式方程，或直接取 $\underline{\nu}$ 的不同数值采样得到不同的标量方程，最后解出 $b_{\underline{\alpha},\underline{\beta}}$.

---

## 2. MMA 实现

### 2.1 核心流程

MMA 实现位于 `/root/Large-Index-Expansion-MMA-Mini/LIEReconstruct.wl`，主函数：

```
LIEReconstruct[rankLevel, degree, order, hexpnList, aregList, topsec, vlist, opts]
```

参数与 `ReconstructReductionRelation` 相同（见 `ReconstructReductionRelation_Documentation.md`）。

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

**结果**：ansatz 是 $ \sum b_{\underline{\alpha}, \underline{\beta}}  \underline{\nu}^{\underline{\beta}} g(\underline{\nu} - \underline{\alpha}) $ 形式的符号多项式

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

**结果**：对每个 $\alpha$，计算 $p(\underline{\alpha}) = \underline{A}^{-\underline{\alpha}}$ 矩阵

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

**关键**：`Thread[vlist -> (vlist - alpha)]` 替换后，$ h_i(\underline{\nu} - \underline{\alpha})$ 仍是 $\underline{\nu}$ 的**符号多项式**

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
- `coeftable[[j, kk+1, l]]` 是 $(n\underline{\theta}+\underline{\nu})^{\underline{\beta}}$ 展开后 $n^k$ 的系数（数值）
- `expntable[[j, m, k-kk+1]][alphaList[[l]]]` 是 $g(\underline{\nu} - \underline{\alpha})$（符号多项式）
- `relansatz` 合成后仍是 $\underline{\nu}$ 的符号多项式

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
| GenerateAnsatz | rank, deg | `Σ b_{α,β}·v^β g(v-α)` | 符号 |
| ComputeBasisPower | α list | `p(α) = A^(-α)` | 数值矩阵 |
| UpdateExpansionTable | `h_i(v)` | `g(v-α) = p(α) * h_i(v-α)` | 符号 |
| SetupEquationsAll | coeftable, expntable | `relansatz = Σ coeff·g(v-α)` | 符号 |
| CoefficientEquations | 符号表达式 | 方程组 `M · b = 0` | 数值 |
| ReduceSolve | 方程组 | 系数 `b_{α,β}` | 数值 |

---

## 3. C++ 实现

### 3.1 核心流程

```
RelationSolver::reconstructReductionRelation<FFInt>(CTable, sector, A_list, Ainv_list, ne, lev, deg, config)
  → AdaptiveEquationBuilder<T>::build(regimes, nimaxLists, ne)        // 主引擎：自适应采样 + 收敛检测
      → GlobalEquationAssembler<T>::evaluateAtNu(nu)                  // 跨 regime 组装
          → RegimeEvaluator<T>::evaluate(nu, nimax_idx)               // 单点评估：computeG → computeF1 → computeF2
```

`AdaptiveEquationBuilder<T>` 是顶层引擎（见 `RelationSolver_ComponentGuide.md`），负责随机采样 ν 点、收敛检测、求解。`RegimeEvaluator<T>` 是内部 worker，在单个 ν 点执行 §3.2 的 Step 2–5。

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

#### Step 3: step3_computeF1 - 计算 (theta*n + nu)^beta

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
                T nu_i = nu[i];  // 采样点 ν 的第 i 分量

                // 二项式展开: (theta + nu_i)^b = sum C(b,m) * theta^(b-m) * nu_i^m
                for (int m = 0; m <= b; ++m) {
                    T coef = static_cast<T>(BINOM[b][m]) *
                             detail::power(th, b - m) *
                             detail::power(nu_i, m);
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

**结果**：`f1_store_[alpha,beta]` 是 `(k_max+1)` 的**数值数组**，表示多项式 `(theta*n+nu)^beta` 在 `n` 的各阶 `k` 的系数

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

        // 每次迭代更新零空间追踪
        int current_nullity = solver.getNullity();
        sampler.update(current_nullity);

        // 检查收敛
        if (sampler.shouldCheck() && sampler.hasConverged()) {
            // 额外验证：在独立采样点上验证零空间基
            auto verify_points = sampler.getVerificationPoints(config_.verification_points);
            auto current_info = solver.getNullspace();
            bool verification_passed = true;

            for (const auto& vnu : verify_points) {
                auto vrows = assembler.evaluateAtNu(vnu);
                T residual = solver.verifyBasis(vrows, current_info, T(config_.tolerance));
                if (residual != T(0)) {
                    verification_passed = false;
                    break;
                }
            }

            if (verification_passed) {
                result.converged = true;
                result.nullspace = current_info;
                return result;
            }
            // 验证失败：将不满足的方程加入系统，继续迭代
            // 这确保零空间维度单调递减，最终收敛到正确解
        }
    }
}
```

**关键**：
- 需要在**多个采样点** `nu` 上评估方程，然后求解
- **验证反馈机制**：收敛后，在独立采样点上验证零空间基。若残差非零，将失败方程重新加入系统并继续迭代，确保零空间维度单调递减至正确解

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
| `(theta+v)^beta` | 符号多项式 | 数值 `(theta*n+nu)^beta` |
| 方程来源 | `CoefficientEquations[P(v), v]` | `M(nu) · x = 0` 在多个 nu 点 |

### 4.2 方程构建方式

**MMA**:
```
P(v) = Σ b_{α,β} · v^β · g(v-α)  [符号多项式]
eqs = CoefficientEquations[P(v), v]  →  每个 v^power 一个方程
```

**C++**:
```
M(nu) = Σ b_{α,β} · (theta*n+nu)^β · g(nu-α)  [标量]
在多个 nu 点采样 → M1·x=0, M2·x=0, ... (x=b_{α,β})
```

### 4.3 关键公式对比

**MMA 的系数提取**:
```
coeftable[[j, kk+1, l]]  -- (n*theta+v)^beta 展开后 n^{kk} 的系数
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

### 4.4 步骤对应关系

| 数学操作 | MMA (§2) | C++ (§3) |
|---------|----------|----------|
| 生成 α,β 多指标 | `GenerateAnsatz` 内 `seeding[]` | `generateAllIndices()` (Combinatorics.hpp) |
| 计算 `p(α) = A^(-α)` | `ComputeBasisPower` (递归 memoization) | `RegimeData::prepare` → `computePRecursive` |
| 计算 `g(ν-α)` | `UpdateExpansionTable` (符号多项式) | `RegimeEvaluator::step2_computeG` (数值) |
| 展开 `(θn+ν)^β` | `SetupCoefficientTable` (符号) | `RegimeEvaluator::step3_computeF1` (数值) |
| 卷积 `f1 * g` | `SetupEquationsAll` 内循环 | `RegimeEvaluator::step4_computeF2` |
| 构建方程 | `CoefficientEquations` (提取多项式系数) | `RegimeEvaluator::buildFinalMatrix` (组装数值矩阵) |
| 求解 | `ReduceSolve` (RowReduce, 模数 p) | `AdaptiveEquationBuilder::build` (自适应 ν 采样 + 求解) |

---

## 5. 验证方法

### 5.1 ν 采样验证的数学原理

由 §1.3，关系方程本身是 $\underline{\nu}$ 和 $n$ 的多项式恒等式：

$$
\sum_{\underline{\alpha},\underline{\beta},\underline{\gamma},k} b_{\underline{\alpha},\underline{\beta}} \,\underline{A}^{-\underline{\alpha}} c_{k,\underline{\gamma}}  (\underline{\theta} n + \underline{\nu})^{\underline{\beta}}  (\underline{\nu}-\underline{\alpha})^{\underline{\gamma}} / n^k = 0
$$

**两种等价处理方式**：

| 路径 | 方法 | 本质 |
|------|------|------|
| **符号系数提取** | 将等式视为 ν 的多项式，提取每个单项式的系数 → 每个 ν^power 一个方程 | MMA: `CoefficientEquations[P(v), v]` (§2 Step 5) |
| **数值 ν 采样** | 取不同 ν 值代入，得到标量方程 $M(\nu^{(1)})x=0$, $M(\nu^{(2)})x=0$, ... → 自适应采样求解 | C++: `M(nu)·x = 0` (§3 Step 6) |

**ν 采样验证流程**：

给定关系 `rel = Σ b_{α,β} · v^β · g(v-α)`，检验其在 ν 点是否成立：

1. **代入 ν 值**：`substituted = rel /. {v_i → ν_i}` — 将符号多项式变为标量，`g(v-α)` 变为 `g(ν-α)` = `j(ν-α)` 即一个具体的费曼积分
2. **Kira 约化**：对 `j(ν-α)` 应用约化规则 —
   - 若所有指标 ≤ 0（零 sector）：`j[...] → 0`
   - 若属于 trivial sector：`j[...] → 0`
   - 其余：按 Kira 规则约化到主积分 (master integral)
3. **取模验证**：代入 d 并对模数取模，`result = PolynomialMod[reduced /. d → dVal, modulus]`。`result = 0` 表示该关系在该 ν 点成立

> **实际执行与测试结果**：完整的 Kira 约化规则生成、Mathematica 验证代码、各配置的测试结论和差异分析见 **[Test-Relation.md](../verify/docs/Test-Relation.md)** §3–§6。

---

## 6. 参考文件

### MMA 核心库 (`/root/Large-Index-Expansion-MMA-Mini/`)

| 文件 | 用途 |
|------|------|
| `LIEReconstruct.wl` | 主重构算法包 (`LIEReconstruct` 函数) |
| `LIEWorkflow.wl` | 工作流包装 (LIESolveRegions, LIEExpandSeries, LIEGetRelations) |
| `LIEFamilyDefine.wl` | 积分族定义 |
| `LIEExpand.wl` | 级数展开引擎 |
| `LIECoreAlgebra.wl` | 核心代数运算 |
| `LIERegions.wl` | 区域管理 |
| `LIEUtility.wl` | 工具函数 |
| `ExportBinary_IBPMatrix.wl` | 导出 .bin 格式 IBP 矩阵 |
| `KiraRuleLoader.wl` | Kira 约化规则加载器 |
| `/root/M2Kira.wl` | Kira 输入生成接口 |

### MMA 工作流脚本 (`/root/Large-Index-Expansion-MMA-Mini/`)

| 文件 | 用途 |
|------|------|
| `Compare-FamilyGenerate-bub00.wl` | 定义 bub00 积分族并导出 .bin |
| `Compare-Expand-bub00.wl` | bub00 LIE 展开 |
| `Compare-Reconstruct-bub00.wl` | bub00 关系重构 (调用 LIEReconstruct) |
| `NuSamplingVerify-MMA.wl` | ν 采样验证 MMA 关系 |
| `Unified-NuSampling-Verify.wl` | 统一 ν 采样验证框架 |

### C++ (本仓库)

| 文件 | 用途 |
|------|------|
| `include/RelationSolver.hpp` | 核心求解器 (`RegimeData<T>`, `RegimeEvaluator<T>`, `AdaptiveEquationBuilder<T>`, `RelationCoefficient<T>`) |
| `include/IncrementalRelationSolver.hpp` | 高级多级求解器 |
| `tests/test_relationFF.cpp` | 关系重构测试入口 |
| `tools/test_ff_verify.cpp` | FireFly 有限域验证 |

### C++ 文档 (本仓库)

| 文件 | 用途 |
|------|------|
| `docs/RelationSolver_ComponentGuide.md` | 组件 API、算法步骤、性能考量 |
| `docs/RelationSolver_QuickReference.md` | 快速参考、常见模式、参数调优 |
| `docs/RelationSolver_Documentation_Hub.md` | 文档中心、阅读路线、交叉引用索引 |

### 本仓库 MMA 脚本 (`mma/`)

| 文件 | 用途 |
|------|------|
| `mma/ReconstructReductionRelation.wl` | MMA 包副本 (函数名 `ReconstructReductionRelation`，与 `LIEReconstruct` 一致) |

### 生成输出

| 文件 | 用途 |
|------|------|
| `build/Compare-CPPRelation-bub00*.m` | C++ 导出的关系结果 (Mathematica 可读) |
| `build/ExpansionMMA_*.m` | C++ 导出的展开系数 |
| `IBPMat_*.bin` + `RingData_*.bin` | 测试数据文件 (由 Compare-FamilyGenerate 生成) |
