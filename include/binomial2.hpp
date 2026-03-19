#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include <cstddef>
#include <array>

namespace combinatorics {
    constexpr int MAX_VAL = 30;

    // 返回组合数表的引用，并在第一次调用时初始化
    inline const std::array<std::array<long long, MAX_VAL>, MAX_VAL>& getBinomial() {
        static std::array<std::array<long long, MAX_VAL>, MAX_VAL> table{};
        static bool initialized = false;
        if (!initialized) {
            for (int i = 0; i < MAX_VAL; ++i) {
                table[i][0] = 1;
                table[i][i] = 1;
                for (int j = 1; j < i; ++j) {
                    table[i][j] = table[i-1][j-1] + table[i-1][j];
                }
            }
            initialized = true;
        }
        return table;
    }

    // 便捷访问函数
    inline long long C(int n, int k) {
        if (n < 0 || k < 0 || n >= MAX_VAL || k > n) return 0;
        return getBinomial()[n][k];
    }

    // 可选：提前初始化函数，如果需要可以调用
    inline void initBinomial() {
        getBinomial(); // 确保初始化
    }
}
#endif