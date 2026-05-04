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

### 1.4 θ 参数与 limitSector 语义

**θ 向量**（在代码中存储为 `RingCell::limitSector`，运行时赋值给 `RegimeData::theta`）定义了渐进展开的 **sector 方向**。

当取 $n \to \infty$ 极限时，做变量替换 $\underline{\nu} \to \underline{\nu} + \underline{\theta} n$。不同的 $\underline{\theta}$ 值对应 $\underline{\nu}$ 行为的不同区域（sector），同一费曼积分在不同 sector 下有不同的渐进展开形式。

**关键点**：

- `limitSector` / `theta` **不是** "边界标记" 或 "索引范围限制" — 它是展开 sector 的**参数**
- 每个 regime（特征方程组的极小伴随素理想对应的解分支）从它所属的 sector 继承 $\theta$ 值
- 在 C++ 代码中，$\theta$ 通过 `RingDataLoader` 从二进制文件 `RingData_*.bin` 加载，存储在 `RingCell::limitSector`
- 在 `reconstructReductionRelation()` 中，`reg.theta = sector[r]` 将 sector 的 $\theta$ 赋给该 sector 内的所有 regime
- `RegimeEvaluator::step3_computeF1` 使用 $\theta$ 计算 $(\underline{\theta} + \underline{\nu})^{\underline{\beta}}$ 的多项式展开系数

**数据流**：`RingData_*.bin` → `RingCell::limitSector` → `sector[r]` → `RegimeData::theta` → `step3_computeF1`

### 1.5 统一 1/n 符号基准与 k_max 截断

**问题**：$(\nu+\theta n)^\beta$ 产生 $n$ 的正幂次（$n^0,...,n^{|\beta|}$），而 $g = \sum_k h_k/n^k$ 产生 $n$ 的负幂次（$1/n^0,...,1/n^{k_{max}}$）。直接在正/负混合幂次下做多项式乘法容易混淆哪些阶的贡献是完整的。

**方法**：提出 $n^{|\beta|}$ 因子，将所有量统一到 $1/n$ 基准下：

$$(\nu+\theta n)^\beta = n^{|\beta|} \cdot (\nu/n + \theta)^\beta = n^{|\beta|} \cdot \sum_{i=0}^{|\beta|} \tilde{p}_i / n^i$$

$$g = \sum_{k=0}^{k_{max}} q_k / n^k$$

两者现在都是 $1/n$ 的（截断）幂级数。乘积为：

$$n^{|\beta|} \cdot \sum_{r=0}^{|\beta|+k_{max}} \left(\sum_{i=0}^{\min(|\beta|,r)} \tilde{p}_i \cdot q_{r-i}\right) / n^r$$

**截断依据**：对 $r > k_{max}$，$c_r$ 必然缺少 $\tilde{p}_{|\beta|} \cdot q_r$ 的贡献（$q_r$ 因 $g$ 截断到 $k_{max}$ 而不可得）。因此 $c_r$ 对 $r > k_{max}$ 是**不精确的**，必须丢弃。仅保留 $r \in [0, k_{max}]$ 的系数。

乘以 $n^{|\beta|}$ 还原后，使用的 n 指数范围为 $[|\beta|-k_{max},\; |\beta|]$。

**关键点**：
- f1 **不应**在正 n 幂次域截断到 $k_{max}$——这会错误丢弃 $n^{k_{max}+1},...,n^{|\beta|}$ 的贡献，而这些贡献与 $g$ 的低阶项（$q_0, q_1,...$）相乘后产生的是可靠的方程
- f1 应**完整计算**（$i=0..|\beta|$），然后在 $1/n$ 域卷积，仅在乘积阶段截断到 $k_{max}$
- 旧代码的错误：在 `step3_computeF1` 中将 f1 截断到 `k_max`（丢弃正 n 高幂次），同时保留了 $n^{-1}, n^{-2},...$ 的方程行——前者丢失了可靠的高阶方程，后者产生的负幂次方程也不可靠

**C++ 实现**：
- `step3_computeF1`: `poly[i]` = coeff of $1/n^i$，存储大小 `deg_+1`，不做截断
- `step4_computeF2`: $c_r = \sum_{i=0}^{\min(deg, r)} f1[i] \cdot g[r-i]$（$r=0..k_{max}$），存储大小 `k_max+1`
- `buildFinalMatrix`: 行数 `nb·(k_max+1)`，行 $r$ 对应 n 指数 $deg - r$

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

#### Step 3: step3_computeF1 - 计算 (ν/n + θ)^β 的 1/n 展开

```cpp
template<typename T>
void RegimeEvaluator<T>::step3_computeF1(const std::vector<T>& nu) {
    // poly[i] = coefficient of 1/n^i in (ν/n + θ)^β = n^{-|β|} · (ν + θn)^β
    std::vector<T> work_buffer(deg_ + 1, T(0));

    for (const auto& alpha : alphas_) {
        for (const auto& beta : betas_) {
            std::vector<T> poly(deg_ + 1, T(0));  // 完整 deg+1，不截断
            poly[0] = T(1);

            for (int i = 0; i < ne_; ++i) {
                int b = beta[i];
                if (b == 0) continue;

                // term_poly[m] = coeff of 1/n^m in (ν_i/n + θ_i)^b
                std::vector<T> term_poly(b + 1, T(0));
                T th = static_cast<T>(theta_[i]);
                T nu_i = nu[i];  // <-- 采样点数值

                // 二项式展开: (ν_i/n + θ_i)^b = Σ_m C(b,m)·θ_i^{b-m}·ν_i^m·(1/n)^m
                for (int m = 0; m <= b; ++m) {
                    T coef = static_cast<T>(BINOM[b][m]) *
                             detail::power(th, b - m) *
                             detail::power(nu_i, m);
                    term_poly[m] = coef;
                }
                // 多项式乘法 (1/n 域，不做 k_max 截断，仅限 deg_)
                for (int p1 = 0; p1 <= deg_; ++p1) {
                    if (poly[p1] == T(0)) continue;
                    for (int p2 = 0; p2 <= b; ++p2) {
                        int target_idx = p1 + p2;
                        if (target_idx <= deg_) {
                            work_buffer[target_idx] += poly[p1] * term_poly[p2];
                        }
                    }
                }
                poly = work_buffer;
            }
            f1_store_->insert(alpha, poly, beta);
        }
    }
}
```

**关键**：poly[i] 存储 $1/n^i$ 的系数（非 $n^i$），存储大小 `deg_+1`（完整，不截断到 `k_max`）。系数计算公式 $C(b,m) \cdot \theta^{b-m} \cdot \nu^m$ 在两个域下通用。

**结果**：`f1_store_[alpha,beta]` 是 `(deg_+1)` 的**数值数组**，含义为 $(\nu/n+\theta)^\beta$ 在 $1/n$ 下的展开系数

#### Step 4: step4_computeF2 - 1/n 域卷积 f1 ⊛ g

```cpp
template<typename T>
void RegimeEvaluator<T>::step4_computeF2(const std::vector<int>& alpha,
                                         const std::vector<int>& beta) {
    T* g_ptr = g_store_->retrieve(alpha);
    T* f1_ptr = f1_store_->retrieve(alpha, beta);

    // f2[r] = coefficient of 1/n^r, r ∈ [0, k_max_]
    int f2_len = k_max_ + 1;
    std::vector<T> f2_val(nb_ * f2_len, T(0));

    for (int i = 0; i < nb_; ++i) {
        for (int r = 0; r <= k_max_; ++r) {
            T sum = T(0);
            int i_max = (deg_ < r) ? deg_ : r;
            for (int s = 0; s <= i_max; ++s) {
                // f1[s] = coeff of 1/n^s,  g[r-s] = coeff of 1/n^{r-s}
                sum += f1_ptr[s] * g_ptr[i * (k_max_ + 1) + (r - s)];
            }
            f2_val[i * f2_len + r] = sum;  // 对应 n 指数 deg_ - r
        }
    }
    f2_store_->insert(alpha, f2_val, beta);
}
```

**关键**：在 $1/n$ 域做标准卷积 $c_r = \sum_{s=0}^{\min(deg, r)} f1[s] \cdot g[r-s]$，截断到 $r \le k_{max}$。乘以 $n^{deg}$ 还原后对应 n 指数 $deg - r$。

**结果**：`f2_store_[alpha,beta]` 是 `nb × (k_max+1)` 的**数值数组**，`f2[r]` = $n^{deg-r}$ 的系数

#### Step 5: buildFinalMatrix - 组装最终矩阵

```cpp
template<typename T>
std::vector<std::vector<T>> RegimeEvaluator<T>::buildFinalMatrix() {
    int total_k = k_max_ + 1;  // 仅使用可靠阶数 r = 0..k_max
    int rows = nb_ * total_k;
    int cols = alphas_.size() * betas_.size();

    std::vector<std::vector<T>> mat(rows, std::vector<T>(cols, T(0)));

    for (size_t a = 0; a < alphas_.size(); ++a) {
        for (size_t b = 0; b < betas_.size(); ++b) {
            T* f2_ptr = f2_store_->retrieve(alphas_[a], betas_[b]);
            if (!f2_ptr) continue;

            int col_idx = a * betas_.size() + b;
            // f2[kt] = coefficient of n^{deg_ - kt}
            for (int i = 0; i < nb_; ++i) {
                for (int kt = 0; kt < total_k; ++kt) {
                    int row_idx = i * total_k + kt;
                    mat[row_idx][col_idx] = f2_ptr[i * total_k + kt];
                }
            }
        }
    }
    return mat;
}
```

**关键**：矩阵行数为 `nb·(k_max+1)`（曾为 `nb·(deg+k_max+1)`）。每行对应一个可靠的 n 指数 $deg, deg-1, ..., deg-k_{max}$。

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

| 步骤 | 输入 | 输出 | 存储大小 | 性质 |
|------|------|------|----------|------|
| prepare | α list | `p(α) = A^(-α)` | nb×nb | 数值矩阵 |
| step2_computeG | nu (采样点) | `g[1/n^k]` | nb·(k_max+1) | 1/n 形式 |
| step3_computeF1 | nu (采样点) | `f1[1/n^i]` | deg+1 | 1/n 形式，完整 |
| step4_computeF2 | g, f1 | `f2[1/n^r]` | nb·(k_max+1) | 1/n 卷积，截断到 k_max |
| buildFinalMatrix | f2 | `M(nu) · x = 0` | nb·(k_max+1) 行 | 每行对应可靠 n 指数 |
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

## 6. RemoveSolvedVariables 策略：减少多配置方程冗余

### 6.1 问题背景

当需要对多个 `(lev, deg)` 配置求解关系时，高配置的变量空间中包含大量已被低配置求解的变量。如果不加处理，每个配置独立求解会产生大量冗余方程，且解空间的维度膨胀。

例如 bub00 在 `(lev=2, deg=2)` 时有 36 个变量 `(6 alpha × 6 beta)`，但其中约一半已被 `(lev=1, deg=1)` 和 `(lev=1, deg=2)` 的关系覆盖。

### 6.2 MMA 实现 (`ReconstructReductionRelation.wl`)

#### 核心函数：RemoveSolvedVariables（第 419 行）

```mathematica
RemoveSolvedVariables[cfList0_, alphaList0_, gvList0_, levelBindep_, bindep_, deg_] := Module[
  {cfList, ...},
  
  (* 度数截断 *)
  cfList = cfList0 /. {"b"[alpha_, beta_] /; Total[beta] > deg -> 0};
  
  (* 同层低度消除 *)
  cfList = cfList /. {"b"[alpha_, beta_] /; Or @@ (
    alpha == #[[1]] && And @@ Thread[beta >= #[[2]]] & /@ levelBindep) -> 0};
  
  (* 跨层消除 *)
  cfList = cfList /. {"b"[alpha_, beta_] /; Or @@ (
    And @@ Thread[alpha >= #[[1]]] && And @@ Thread[beta >= #[[2]]] & /@ Join @@ bindep) -> 0};
  
  ...
];
```

两条分量级支配（componentwise dominance）规则：

| 规则 | 来源 | 条件 | 含义 |
|------|------|------|------|
| 同层低度 | `levelBindep` | `α == α_s` 且 `β >= β_s` | 同层已求解的独立变量，其高次单项式冗余 |
| 跨层 | `bindep` | `α >= α_s` 且 `β >= β_s` | 低层的独立变量支配高层的更大平移 |

其中 `levelBindep` 是当前层已求解的独立变量对 `{α_s, β_s}`，`bindep` 是所有更低层的累计。

#### 主循环结构（第 889-991 行）

```
for lev in 0..rankLevel:
    levelBindep = {}
    for deg in 0..maxCoefDeg:
        {cfList, bVars} = RemoveSolvedVariables(..., levelBindep, bindep, deg)
        bSol = SolveDegreeEquations(...)
        bVarIndep = extract independent variables from bSol
        levelBindep += bVarIndep
    bindep += levelBindep
```

关键点：
- 外层循环按 seed level 递增，内层循环按 coefficient degree 递增
- 每次求解前先过滤被支配变量，求解后记录独立变量供后续过滤
- `bindep` 跨层累积，`levelBindep` 仅同层累积

### 6.3 C++ 实现 (`RelationSolver.hpp`)

#### 辅助函数

```cpp
// 分量级支配判断
inline bool componentwiseDominates(const vector<int>& a, const vector<int>& b) {
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i] < b[i]) return false;
    return true;
}

// 变量过滤：返回活跃列索引和掩码
inline pair<vector<size_t>, vector<bool>> filterVariablePairs(
    const vector<vector<int>>& alphas,
    const vector<vector<int>>& betas,
    const vector<AlphaBeta>& levelBindep,
    const vector<vector<AlphaBeta>>& bindep);
```

#### 列压缩与展开

由于变量过滤后列数减少，采用 **列压缩 → 求解 → 展开** 的策略：

```
完整矩阵 (full_cols 列)
    │ stripInactiveColumns()  ← 按 active_mask 剥离
    ▼
压缩矩阵 (active_count 列)
    │ IncrementalNullspaceSolver
    ▼
压缩零空间基 (active_count 列)
    │ expandNullspaceToFull()  ← 被消除列填 0
    ▼
完整零空间基 (full_cols 列)
```

这避免了修改 `RegimeEvaluator` 和 `GlobalEquationAssembler` 的内部实现，仅在 `AdaptiveEquationBuilder::build()` 中插入列压缩步骤。

#### 高层 API

```cpp
template<typename T>
std::vector<LevDegResult<T>> reconstructAllRelations(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne, int lev_max, int deg_max,
    const AdaptiveSamplingConfig& config = {});
```

`LevDegResult<T>` 包含：
- `linear_result` — 展开回完整变量空间的 `LinearSystemResult`
- `coeff` — `RelationCoefficient`（可直接用于 MMA 导出）
- `independent_pairs` — 该配置的独立变量对（供验证/调试）
- `active_vars` / `total_vars` — 过滤前后的变量数

#### 独立变量提取

求解后从压缩列空间的零空间信息中提取独立变量：

```
solver.free_cols (压缩空间)
    → active_indices[free_col] → 完整空间列索引
    → (col / nBeta, col % nBeta) → (alpha_idx, beta_idx)
    → (alphas[a], betas[b]) → AlphaBeta 对
```

### 6.4 效果数据（bub00, ne=2, k_max=4）

| (lev, deg) | 总变量 | 活跃变量 | 消除 | 解空间维数 | 独立变量 |
|-----------|--------|---------|------|----------|---------|
| (0, 0) | 1 | 1 | 0 | 0 | 0 |
| (0, 1) | 3 | 3 | 0 | 0 | 0 |
| (0, 2) | 6 | 6 | 0 | 0 | 0 |
| (1, 0) | 3 | 3 | 0 | 0 | 0 |
| (1, 1) | 9 | 9 | 0 | 2 | 2 |
| (1, 2) | 18 | 13 | 5 | 1 | 1 |
| (2, 0) | 6 | 6 | 0 | 0 | 0 |
| (2, 1) | 18 | 12 | 6 | 0 | 0 |
| (2, 2) | 36 | **19** | **17** | 4 | 4 |

最高配置 `(lev=2, deg=2)` 变量减少 **47%**（36→19）。

### 6.5 与 MMA 的对应关系

| 概念 | MMA | C++ |
|------|-----|-----|
| 变量空间 | `b[α, β]` 符号 | `(alpha_idx, beta_idx)` 列索引 |
| 同层消除 | `levelBindep` 列表 | `filterVariablePairs(..., levelBindep, {})` |
| 跨层消除 | `bindep` 列表 | `filterVariablePairs(..., {}, bindep)` |
| 变量过滤 | 模式匹配 `/. b[...]->0` | `active_mask[col]=false` |
| 独立变量 | `bVarIndep = Cases[recIbp, b[...]]` | `solver.free_cols → active_indices → AlphaBeta` |
| 累积机制 | `levelBindep = Join[levelBindep, bVarIndep]` | `levelBindep.push_back(independent_pairs)` |

---

## 7. 待解决问题

1. **采样点选取**：C++ 的自适应采样是否收敛到正确解？
2. **方程完整性**：在多个 nu 点采样是否等价于 MMA 的系数提取？
3. **数值精度**：有限域算术是否存在奇异性？
4. **独立变量提取精度**：当前从压缩列空间映射独立变量，需验证与 MMA 的 `bVarIndep` 提取逻辑是否完全一致
5. ✅ **阶数稳定性检测**：已实现。C++ 两阶段零冗余设计 — Phase 1 ν-采样 + Phase 2 阶数分析，完整复现 MMA `SolveDegreeEquations` plateau 逻辑（详见 [Section 7](#7-阶数稳定性检测)）

---
## 7. 阶数稳定性检测 (Order-Stability Convergence)

### 7.1 问题背景

C++ 求解器使用预计算的展开系数（固定 `k_max`）一次性构建方程。然而，不同 `(lev, deg)` 配置需要的最小展开阶数不同——低配置可能只需 1-2 阶，高配置可能需要 4+ 阶。如果不检测阶数稳定性，可能：
- 对低配置使用过多展开阶数（浪费计算资源）
- 对高配置使用过少展开阶数（解空间未真正收敛）

### 7.2 MMA 实现 (`ReconstructReductionRelation.wl`)

MMA 通过 `SolveDegreeEquations` (L541-678) 实现逐阶递增求解：

```wl
While[k <= maxorder && !converged,
    (* 构建当前阶数 k 的方程 *)
    eqs = SetupEquationsAll[..., k, deg];
    
    (* 求解 *)
    bSolNew = DispatchLinearSolve[eqs, bVarReg, ...];
    
    (* 稳定性检测 *)
    If[Length[bSolAcc] == Length[bVars],   (* Trigger 2: 所有变量已确定 *)
        state["Stable"] = True;
        state["StableLevel"] = k;
    ];
    
    (* Plateau 确认 *)
    If[state["Stable"] && k >= state["StableLevel"] + plateausize,
        converged = True; Break[];
    ];
    k++;
];
```

三个稳定性触发器与 Plateau 确认详见 `ReconstructReductionRelation.wl` L577-671。

### 7.3 C++ 实现 (`RelationSolver.hpp`)

C++ 采用**两阶段零冗余**设计：

**Phase 1 (ν-采样)**: 标准自适应采样循环。每次 `evaluateAtNu()` 返回的行通过 `splitRowsByOrder()` 按展开阶数 r 拆分为 `rows_by_order[r]` 并累积存储。行布局已知：`row_idx = i·(k_max+1) + r`。

**Phase 2 (阶数分析)**: ν-收敛后，`analyzeOrderStability()` 使用 `IncrementalNullspaceSolver` 逐阶添加预存行 (r = 0, 1, ..., k_max)，追踪每阶 nullity。

```
analyzeOrderStability(rows_by_order, k_max, plateau_size, num_vars):
    for cur_order = 0..k_max:
        solver.addRows(rows_by_order[cur_order])
        nullity[cur_order] = solver.getNullity()
        
        if nullity[cur_order] == 0:
            return {cur_order, 0}   // 决定性终止（对应 MMA Trigger 2）
        
        if nullity stable for plateau_size+1 orders:
            return {first_stable_order, final_nullity}  // plateau 确认
    
    return {-2, nullity[k_max]}  // 可用阶数内未稳定
```

**零冗余保证**: `step2_computeG`（最昂贵的操作，访问 5D 张量）、`step3_computeF1`、`step4_computeF2` 在采样时每个 ν 点仅运行一次。Phase 2 仅对预存行做增量高斯消元。

### 7.4 输出格式

```
=== Stability bounds (lev x deg) ===
           0       1       2
  0        0       1       2
  1        1       3      -2
  2        2       4      -2

=== Relation counts (lev x deg) ===
           0       1       2
  0        0       0       0
  1        0       2       1
  2        0       0       4
```

- `stable_order` 非负 = 在该展开阶数稳定；`-2` = 可用阶数不足
- 只有稳定解应被视为可靠

### 7.5 效果数据 (bub00, ne=2, k_max=4)

| (lev, deg) | 活跃变量 | 解维数 | stable_order | 说明 |
|-----------|---------|--------|-------------|------|
| (0, 0) | 1 | 0 | 0 | 零阶确定 |
| (0, 1) | 3 | 0 | 1 | 一阶确定 |
| (0, 2) | 6 | 0 | 2 | 二阶确定 |
| (1, 0) | 3 | 0 | 1 | 一阶确定 |
| (1, 1) | 9 | 2 | 3 | 三阶稳定 |
| (1, 2) | 13 | 1 | **-2** | 阶数不足 |
| (2, 0) | 6 | 0 | 2 | 二阶确定 |
| (2, 1) | 12 | 0 | 4 | 四阶达到 nullity=0 |
| (2, 2) | 19 | 4 | **-2** | 阶数不足 |

对角线模式（更高 `lev + deg` → 需要更高展开阶数）与 MMA 经验一致。

---
## 8. 参考文件

### MMA
- `mma/ReconstructReductionRelation.wl` - 主重构算法（含 RemoveSolvedVariables）
- `LIEWorkflow.wl` - 工作流包装
- `Compare-Reconstruct-bub00.wl` - 测试脚本

### C++
- `include/RelationSolver.hpp` - 核心求解器（含 `reconstructAllRelations`、`filterVariablePairs`、列压缩）
- `tests/test_relationFF.cpp` - 测试入口（使用 `reconstructAllRelations`）
- `AllRelations_bub00_k5.m` - 统一输出文件（所有 lev/deg 的结果）
