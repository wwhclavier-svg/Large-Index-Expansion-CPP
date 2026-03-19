#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include <vector>
#include <cstddef>

// 包含组合数表定义和初始化函数（已使用 inline 避免多重定义）
#include "binomial.hpp"

// 安全的组合数访问（检查边界）
inline long long comb(int n, int k) {
    if (k < 0 || k > n || n >= MAX_VAL) return 0;
    return BINOM[n][k];
}

// 计算给定变量个数和总和下的多重指标状态数（非负整数解个数）
inline long long getCapacity(int n_vars, int sum_val) {
    if (n_vars == 0) return 1;
    if (sum_val < 0) return 0;
    int n = sum_val + n_vars - 1;
    int k = n_vars - 1;
    return comb(n, k);
}

/**
 * 根据多重指标向量 v（长度为 ne）和总层级 level，计算其在降序枚举中的索引。
 * 算法：对前 ne-1 个分量，累加所有可能比当前值大的情况。
 */
inline int getIndex(const std::vector<int>& v, int current_level) {
    int n = static_cast<int>(v.size());
    int index = 0;
    int remaining = current_level;
    for (int i = 0; i < n - 1; ++i) {
        if (remaining > v[i]) {
            int slots = n - i;
            index += comb(remaining - v[i] + slots - 2, slots - 1);
        }
        remaining -= v[i];
    }
    return index;
}

/**
 * 从索引还原多重指标向量（长度为 ne，总和为 level）。
 * 与 getIndex 互逆。
 */
inline std::vector<int> readIndex(int idx, int level, int ne) {
    std::vector<int> seed(ne, 0);
    int remaining = level;
    for (int i = 0; i < ne - 1; ++i) {
        for (int x = remaining; x >= 0; --x) {
            long long cnt = comb((remaining - x) + (ne - i - 2), ne - i - 2);
            if (idx < cnt) {
                seed[i] = x;
                remaining -= x;
                break;
            } else {
                idx -= cnt;
            }
        }
    }
    seed[ne - 1] = remaining;
    return seed;
}

/**
 * 生成所有 level 从 0 到 maxlevel 的多重指标向量。
 * 返回三维向量：第一维索引为 level，第二维为该 level 的所有种子向量。
 */
inline std::vector<std::vector<std::vector<int>>> seedGenerator(int maxlevel, int ne) {
    std::vector<std::vector<std::vector<int>>> result(maxlevel + 1);
    for (int level = 0; level <= maxlevel; ++level) {
        long long cap = getCapacity(ne, level);
        result[level].reserve(cap);
        for (long long idx = 0; idx < cap; ++idx) {
            result[level].push_back(readIndex(idx, level, ne));
        }
    }
    return result;
}

// 查找向量中最后一个非零元素的索引（从0开始），若无则返回 -1
inline int lastNonZero(const std::vector<int>& list) {
    for (int i = static_cast<int>(list.size()) - 1; i >= 0; --i) {
        if (list[i] != 0) return i;
    }
    return -1;
}

// 检查整数 a 是否不在 list 中
inline bool notin(const std::vector<int>& list, int a) {
    for (int x : list) if (x == a) return false;
    return true;
}

/**
 * 以下两个函数用于在给定种子和增量时快速获取新种子在新层级中的索引。
 * 注意：它们会临时修改种子向量，但调用后会恢复原值，因此 seed 参数不能为 const。
 */
inline int getIndexOffSet(int old_level, std::vector<int>& seed, int k, int i) {
    seed[i] += k;
    int new_level = old_level + k;
    int idx = getIndex(seed, new_level);
    seed[i] -= k;
    return idx;
}

inline int getIndexOffSet(int old_level, std::vector<int>& seed, int k, int i, int l, int j) {
    seed[i] += k;
    seed[j] += l;
    int new_level = old_level + k + l;
    int idx = getIndex(seed, new_level);
    seed[i] -= k;
    seed[j] -= l;
    return idx;
}

#endif // COMBINATORICS_HPP