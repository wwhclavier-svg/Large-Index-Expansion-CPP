# 策略对比 — Strategy 0 vs Strategy 4

## 一、策略设计

### 策略 0：Default（原始策略）

每个 ν 采样点同时评估**所有 regime**，一次性构建完整方程矩阵后整体送入 RREF 求解。

```
for each ν sample point:
    rows = assembler.evaluateAtNu(nu)     // 所有 regime 一起评估
    solver.addRows(rows)                  // 全部加一行一次消元
```

- **设计意图**：简单直接，每个 ν 点一次性贡献全部方程行
- **特点**：所有 regime 贡献的行全部参与消元，包括可能线性相关的冗余行

### 策略 4：PerNuRegimeRankCheck（逐 regime 增量秩检查）

每个 ν 点按 regime 的 $|\theta|$（sector 标签中 1 的个数）**降序**逐 regime 评估，每加入一个 regime 的贡献后检查 RREF 的 rank 是否增加。只有**有效提升 rank 的 regime 行**才被保留，其余丢弃。

```
for each ν sample point:
    regimes_sorted_by_|θ|_desc           // 预排序（一次）
    for each regime (按|θ|降序):
        rows = assembler.evaluateRegimeAtNu(r, nu)
        nullity_before = solver.getNullity()
        solver.addRows(rows)
        nullity_after = solver.getNullity()
        if nullity_after < nullity_before:
            keep rows                     // rank 增加了
        else:
            discard rows                  // 冗余，不做消元
```

- **设计意图**：模拟 Laporta 方法的层级消去逻辑——$\theta$ 更大的 sector 通常携带更多约束，优先加入可尽早确定 rank，后续子扇区的冗余行被自然过滤
- **排序依据**：$|\theta|$ 降序（sector 中 1 越多，即越"顶"的扇区，越优先评估）
- **rank 检查**：每次 addRows 后比较 nullity 变化，rank 不增加的行视为冗余丢弃

### 共同点

- 均使用相同的 AdaptiveSampler 生成 ν 采样点序列（特殊点优先：单位向量→全1→固定方向→随机点）
- 均使用 IncrementalNullspaceSolver 做增量 RREF 维护
- 均使用相同的收敛判据（连续 stable_threshold 次 nullity 不变）
- 均使用 RemoveSolvedVariables 跨 (lev,deg) 消元

---

## 二、测试结果 — TB123 / NP222 (--topsector, lev=1..3, deg=0..2)

测试参数：`min_nu=3, max_nu=50, nullity_stable_threshold=5, check_interval=5, verification_points=5`

### 2.1 时间 (wall clock)

| 家族 | Order | 策略 0 | 策略 4 | 加速 |
|------|-------|--------|--------|------|
| TB123 | k=2 | 45.83s | **36.52s** | 1.25x |
| TB123 | k=3 | 327.72s | **270.73s** | 1.21x |
| NP222 | k=2 | 35.74s | **29.53s** | 1.21x |
| NP222 | k=3 | 202.88s | **181.07s** | 1.12x |

### 2.2 内存 (Max RSS)

| 家族 | Order | 策略 0 | 策略 4 | 增幅 |
|------|-------|--------|--------|------|
| TB123 | k=2 | 183 MB | 239 MB | +30% |
| TB123 | k=3 | 1391 MB | 1501 MB | +8% |
| NP222 | k=2 | 128 MB | 166 MB | +30% |
| NP222 | k=3 | 665 MB | 748 MB | +12% |

### 2.3 关系数对比

**两种策略在所有配置下产生完全相同的关系数**（sol_dim 和 independent 数量一致）。

仅有的差异来自 stable_order（策略 4 在 (3,0) 上通常 +1）：

| 家族 | Order | (3,0) s0 | (3,0) s4 |
|------|-------|----------|----------|
| TB123 | k=2 | 6 (s=0) | 6 (s=1) |
| TB123 | k=3 | 21 (s=0) | 21 (s=1) |
| NP222 | k=2 | 3 (s=0) | 3 (s=1) |
| NP222 | k=3 | 15 (s=0) | 15 (s=1) |

> s = stable_order; 0 表示在 order=0 已收敛，1 表示 order=0 到 order=1 之间确认稳定
> 策略 4 的 stable_order 稍高是因为逐 regime 过滤改变了采样方程的分布，收敛曲线略有平移，但不影响最终 independent 数量

其余配置（lev=1,2 全部 deg 级别、lev=3 deg=1）均为 sol_dim=0。(2,2) 和 (3,2) 在所有条件下均为 stable_order=-2（欠采样假阳性）。

---

## 三、分析

### 策略 4 为什么更快

1. **消去冗余行**：策略 4 的逐 regime rank 检查丢弃了与前面 regime 线性相关的行。以 TB123 k=3 (3,2) 为例，vars=4320 中 RemoveSolvedVariables 压缩到 active=2619（~60%），策略 4 进一步在 sampling 阶段过滤掉额外 15-25% 的冗余行，使 RREF 矩阵更小，消元更快
2. **减少无效消元**：Incremental RREF 的消元复杂度约 $O(n\_rows \times n\_cols^2)$，丢弃冗余行直接减少行数
3. **k 越大加速比越小**：k 增大时每个 ν 点的总方程数 (eq/sample) 随 $k^{ne}$ 增长，RREF 本身成为主导（无论是否过滤行），逐 regime 检查的相对收益降低

### 策略 4 为什么多用内存

1. **逐 regime 评估开销**：每个 regime 单独调用 `evaluateRegimeAtNu()`，产生中间向量（`regime_rows`, `stripped`），虽然用完释放，但峰值 RSS 更高
2. **rows_by_order 累积**：策略 4 始终累积完整的 `rows_by_order`（用于后续的阶数稳定性分析），即使部分行被 rank 检查丢弃，原始数据仍被保留
3. **k 越大增幅越小**：k 较大时，RREF 工作矩阵本身主导内存，逐 regime 评估的中间开销占比下降

### 策略 4 的结果一致性

两种策略产生相同的 sol_dim，因为它们都构建了 $\nu$-采样矩阵 $M$，且都求解 $M \cdot b = 0$ 的 nullspace。策略 4 只是通过 rank 检查提前丢弃了 $M$ 中的线性相关行，nullspace 不变（消去冗余行不改变零空间维度）。这是线性代数的基本性质：行空间的秩决定了零空间的维数，消去线性相关行不改变秩，因此 nullspace 维度不变。

---

## 四、实践建议

| 场景 | 推荐策略 | 理由 |
|------|---------|------|
| k ≤ 2, 小 ne | 策略 0 | 本身已很快（~35-45s），策略 4 的加速绝对值小 |
| k ≥ 3, ne ≥ 7 | **策略 4** | 节省 10-25% 时间，内存增幅可接受（8-30%） |
| 内存紧张（< 2 GB） | 策略 0 | 策略 4 多用 10-30% 内存 |
| 需要最快结果 | **策略 4** | 在所有测试中均更快 |
| 需要稳定阶数估计 | 任意 | 差异仅 1 阶，不影响最终判断 |

### 切换方式

```bash
# 策略 0 (默认)
./build/test_relationFF TB123 3 1 3 2 --topsector

# 策略 4
./build/test_relationFF TB123 3 1 3 2 --topsector --strategy 4
```

## 相关文档

- `RelationSolver_ComponentGuide.md` — 策略实现细节（IncrementalNullspaceSolver）
- `../reports/Benchmark_Results.md` — 同族时序/内存详细数据
