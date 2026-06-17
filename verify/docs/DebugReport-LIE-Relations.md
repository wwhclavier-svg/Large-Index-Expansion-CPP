# C++ LIE 零空间问题诊断报告

## 问题确认

**测试结果**：
- Kira 规则本身正确：j[2,1] - (d-3)*j[1,1] = 0 ✓
- C++ 零空间向量不满足 Kira 规则：所有测试点都 FAIL ✗

**文件：**
- `/root/Large-Index-Expansion-CPP/verify/tests/simple_kernel_test.m`
- `/root/Large-Index-Expansion-CPP/verify/tests/debug_lie_relations.m`

## 关键发现

1. **MatrixBuild 测试只验证 M·x = 0，不验证 Kira 规则**
   - MatrixBuild 使用 M(ν)·x = 0 验证零空间
   - 但这个验证不能保证零空间向量的"物理正确性"

2. **KiraVerify 失败（0/9 PASS）**
   - 将 C++ 零空间向量代入 Kira 规则后，结果不为 0

3. **Kira 规则本身是自洽的**
   - j[2,1] - (d-3)*j[1,1] = 0 正确
   - 说明问题在 C++ 输出的关系，不在 Kira 规则

## 需要进一步检查的地方

### 1. step2_computeG 中的 h_jk 计算
```cpp
// 第 600-613 行
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
```

**可能的 bug**：
- `gamma` 的计算是否正确？readIndex 函数的行为是否符合预期？
- `(nu - alpha)^gamma` 的计算顺序和方式是否正确？

### 2. step3_computeF1 中的多项式展开
```cpp
// 第 647-680 行
for (int i = 0; i < ne_; ++i) {
    int b = beta[i];
    if (b == 0) continue;
    
    // 二项式展开: (theta + nu)^b = sum_{m=0}^{b} C(b,m) * theta^{b-m} * nu^m
    for (int m = 0; m <= b; ++m) {
        T coef = static_cast<T>(BINOM[b][m]) *
                 detail::power(th, b - m) *
                 detail::power(nu_i, m);
        term_poly[m] = coef;
    }
    // ...
}
```

**检查**：二项式展开是否正确？

### 3. buildFinalMatrix 中的矩阵组装
```cpp
// 第 770-780 行
// f2[kt] = 系数 of n^{b_sum - kt} → n-指数 e = b_sum - kt
// 行映射：n-指数 e → deg_ - e = deg_ - b_sum + kt
for (int i = 0; i < nb_; ++i) {
    for (int kt = 0; kt < total_k; ++kt) {
        int row_idx = i * rows_per_basis + (deg_ - b_sum + kt);
        if (row_idx >= 0 && row_idx < rows && col_idx < cols) {
            mat[row_idx][col_idx] = f2_ptr[i * total_k + kt];
        }
    }
}
```

**检查**：行索引计算是否正确？

## 建议的调试步骤

1. **在 test_relationFF 中添加 DEBUG 输出**：
   - 打印每个采样点 ν 的 g 函数值
   - 打印 Mext 矩阵的前几行

2. **创建单元测试**：
   - 隔离测试 step2_computeG
   - 隔离测试 step3_computeF1
   - 隔离测试 step4_computeF2

3. **交叉验证**：
   - 用 MMA 的 LIEReconstruct 函数计算相同配置的关系
   - 对比 C++ 和 MMA 的输出

## 已知信息

- **Alphas**: {{0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {2,0}}
- **Betas**: {{0,0}, {0,1}, {1,0}}
- **零空间向量数**: 3
- **系数矩阵维度**: 18 x 4 (18 variables, 4 columns = 1 particular + 3 nullspace)
