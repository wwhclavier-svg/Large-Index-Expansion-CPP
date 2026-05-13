# CPP-KiraVerify Debugger Report

> **更正 (2026-05-02)**：之前的分析将正确格式 `g[v1-α,...]` 误判为 bug 并修改为 `j[α]`，这是错误的。`g[v1-α,...]` 格式是正确的（见 Apr 30 commit ebf89df）。现已恢复为正确的 `g[v1-α,...]` 格式。

## 1. 理论框架

### 1.1 LIE 关系方程的正确形式

根据 LIE 理论（LIE Notes, algorithms/ReconstructAlgorithm.md §1.2），IBP 关系方程的形式为：

$$\sum_{\alpha,\beta} b_{\alpha,\beta}(\nu) \cdot \nu^\beta \cdot j(\nu - \alpha) = 0 \quad (1)$$

其中：
- $j(\nu)$ 是原始费曼积分（未平移）
- $b_{\alpha,\beta}(\nu) = \sum_\beta c_\beta \nu^\beta$ 是 $\nu$ 的多项式系数
- 方程 (1) 是一个关于 $\nu$ 和 $n$ 的多项式恒等式

**关键**：在渐进展开框架下，关系方程 (1) 对所有 ν 成立。

### 1.2 g 和 j 的关系

$$g(\nu - \alpha; n) = j(\theta n + \nu - \alpha)$$

- g 是做了 n→∞ 平移后的积分
- j 是原始积分

### 1.3 为什么用 g[ν-α] 格式

MMA 验证脚本（NuVerify-Relations.wl）的工作流程：
1. 代入 ν 值：`g[v1-α1, v2-α2]` → `g[ν1-α1, ν2-α2]`
2. g→j 转换：`g[ν1-α1, ν2-α2]` → `j[ν1-α1, ν2-α2]` = $j(\nu - \alpha)$

因此 `g[v1-α,...]` 格式是**正确**的，Apr 30 的 commit (ebf89df) 已经正确说明了这一点。

**错误分析（已撤销）**：之前误认为 `g[v1-α,...]` 是错误的，并修改为 `j[α]`。这是完全错误的——`j[α]` 表示的是 $j(\alpha)$，缺少 ν 部分，无法用于验证。

## 2. Bug 修复记录

### 2.1 问题描述

**错误修复（已撤销）**：2026-05-02 错误地将 `g[v1-α,...]` 格式改为 `j[α]` 格式。

- **文件**：`tests/test_relationFF.cpp`
- **影响**：两处输出函数 `exportRelationToMMA_Polynomial` 和 `exportRelationToMMA_Unified`
- **状态**：已恢复为正确的 `g[v1-α,...]` 格式（2026-05-02）

### 2.2 恢复的修复内容

**修复后（正确）**：
```cpp
out << "*g[";
for (int i = 0; i < ne; ++i) {
    if (i > 0) out << ",";
    int alpha_i = alphas[a][i];
    if (alpha_i == 0) out << "v" << (i + 1);
    else if (alpha_i == 1) out << "v" << (i + 1) << "-1";
    else out << "v" << (i + 1) << "-" << alpha_i;
}
out << "]";
// 输出: g[v1,v2-1] 等
```

## 3. 验证结果

修复后输出的格式为 `g[v1,v2-1]`（v-offset 形式），与 MMA 验证脚本的期望格式一致：

```
78498295*g[v1,v2] + 89712338*g[v1,v2-1] + 1*g[v1,v2-2]
```

## 4. 关键教训

1. **不要轻易修改已通过验证的代码**：Apr 30 的 commit (ebf89df) 已经正确说明了 `g[ν-α]` 格式的用途
2. **区分 j[α] 和 g[ν-α]**：
   - `j[α]` = $j(\alpha)$（原始积分在 α 处）
   - `g[ν-α]` = $g(\nu - \alpha; n)$（平移后的积分）
3. **验证脚本期望 g 格式**：NuVerify-Relations.wl 中的 `gToJ` 函数期望 `g[...]` 格式

## 5. 相关文件

- `/home/ykm/Large-Index-Expansion-CPP/tests/test_relationFF.cpp` — 已恢复正确的 g-format 输出
- `/home/ykm/Large-Index-Expansion-CPP/verify/bub00/Relations_*.m` — 9 个正确格式的关系文件
- `/home/ykm/Large-Index-Expansion-CPP/verify/bub00/Compare-CPPRelation-*.m` — 9 个正确格式的比较文件
- `/home/ykm/Large-Index-Expansion-CPP/verify/scripts/NuVerify-Relations.wl` — MMA 验证脚本（期望 g 格式）

## 相关文档

- `docs/verify/Test-Relation.md` — 关系输出验证方法（含 KiraVerify）
- `ReconstructAlgorithm.md §1.2` — g 格式数学定义
