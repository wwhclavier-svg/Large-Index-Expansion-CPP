#ifndef LINEAR_SOLVER_FF_HPP
#define LINEAR_SOLVER_FF_HPP

#include "firefly/FFInt.hpp"
#include <vector>
#include <algorithm>
#include <iostream>

namespace firefly {

/**
 * 线性方程组解结构体（有限域版本）
 * Mext: 第一列是特解 v，后续列是基础解系（零空间基）
 * S   : 自由变量在原始变量中的列索引
 */
template <typename Field>
struct LinearSystemResult {
    bool hasSolution;                  // 是否有解（非齐次时有效）
    std::vector<std::vector<Field>> Mext; // 拼接矩阵 (v | M)
    std::vector<int> S;                 // 自由变量索引
};

/**
 * 有限域线性求解器
 * @param A_in 系数矩阵 (rows × cols)
 * @param b_in 右端项向量 (长度为 rows)，若为空则解齐次方程组
 * @return 解结构体
 */
template <typename Field>
LinearSystemResult<Field> solveLinearSystem_Firefly(
    const std::vector<std::vector<Field>>& A_in,
    const std::vector<Field>& b_in = {})
{
    using namespace std;

    int rows = A_in.size();
    if (rows == 0) return {false, {}, {}};
    int cols = A_in[0].size();

    // 构建增广矩阵 [A | b]
    vector<vector<Field>> aug(rows, vector<Field>(cols + 1));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j)
            aug[i][j] = A_in[i][j];
        if (b_in.empty())
            aug[i][cols] = Field(0);
        else
            aug[i][cols] = b_in[i];
    }

    int row = 0;
    int col = 0;
    vector<int> pivot_cols;  // 记录主元所在的列

    // 高斯-若尔当消元（化为简化行阶梯形 RREF，仅行交换，不交换列）
    while (row < rows && col < cols) {
        // 在当前列 col 中，从 row 行开始找第一个非零主元
        int pivot_row = -1;
        for (int r = row; r < rows; ++r) {
            if (aug[r][col] != Field(0)) {
                pivot_row = r;
                break;
            }
        }
        if (pivot_row == -1) {
            // 该列全零，跳过
            ++col;
            continue;
        }

        // 交换当前行与主元行
        if (pivot_row != row)
            swap(aug[row], aug[pivot_row]);

        // 归一化当前行
        Field inv = Field(1) / aug[row][col];
        for (int j = col; j <= cols; ++j)
            aug[row][j] *= inv;

        // 消去其他所有行的当前列
        for (int r = 0; r < rows; ++r) {
            if (r == row) continue;
            Field factor = aug[r][col];
            if (factor != Field(0)) {
                for (int j = col; j <= cols; ++j)
                    aug[r][j] -= factor * aug[row][j];
            }
        }

        // 记录主元列，移动到下一行下一列
        pivot_cols.push_back(col);
        ++row;
        ++col;
    }

    // 检查矛盾（非齐次时）
    bool consistent = true;
    if (!b_in.empty()) {
        for (int i = row; i < rows; ++i) {
            if (aug[i][cols] != Field(0)) {
                consistent = false;
                break;
            }
        }
    }
    if (!consistent)
        return {false, {}, {}};

    // 确定自由变量（所有未成为主元的列）
    vector<int> free_cols;
    int p_idx = 0;
    for (int j = 0; j < cols; ++j) {
        if (p_idx < (int)pivot_cols.size() && pivot_cols[p_idx] == j)
            ++p_idx;
        else
            free_cols.push_back(j);
    }

    // 构造特解 v (大小 cols)
    vector<Field> v(cols, Field(0));
    for (size_t i = 0; i < pivot_cols.size(); ++i) {
        int pc = pivot_cols[i];
        v[pc] = aug[i][cols];   // 对应行的最后一列
    }

    // 构造零空间基（每个自由列对应一个基向量）
    vector<vector<Field>> kernel;  // 每个基向量为列向量
    for (int fj : free_cols) {
        vector<Field> vec(cols, Field(0));
        vec[fj] = Field(1);
        // 根据 RREF 方程：主元变量 = - (自由变量系数)
        for (size_t i = 0; i < pivot_cols.size(); ++i) {
            int pc = pivot_cols[i];
            vec[pc] = -aug[i][fj];   // aug[i][fj] 是主元行中自由列 fj 的系数
        }
        kernel.push_back(vec);
    }

    // 组装结果 Mext (cols × (1 + kernel.size()))
    vector<vector<Field>> Mext(cols, vector<Field>(1 + kernel.size()));
    for (int i = 0; i < cols; ++i) {
        Mext[i][0] = v[i];
        for (size_t k = 0; k < kernel.size(); ++k)
            Mext[i][k + 1] = kernel[k][i];
    }

    return {true, Mext, free_cols};
}

/**
 * 打印解结构体（用于调试）
 */
template <typename Field>
void printResult_FF(const LinearSystemResult<Field>& res) {
    using namespace std;
    if (!res.hasSolution) {
        cout << "No solution (inconsistent)" << endl;
        return;
    }
    cout << "Solved (finite field):" << endl;
    for (size_t i = 0; i < res.Mext.size(); ++i) {
        for (size_t j = 0; j < res.Mext[i].size(); ++j) {
            cout << res.Mext[i][j] << "\t";
        }
        cout << endl;
    }
    cout << "Free variables (original columns): ";
    for (int idx : res.S) cout << idx << " ";
    cout << endl;
}

} // namespace firefly

#endif // LINEAR_SOLVER_FF_HPP