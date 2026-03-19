#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP 

#include <vector>
#include <type_traits>

// 包含 Eigen 求解器，它定义了全局的 LinearSystemResult 和 solveLinearSystem_Eigen
#include "LinearSolver_Eigen.hpp"

// 包含 Firefly 求解器，它在 firefly 命名空间中定义了 LinearSystemResult 和 solveLinearSystem_Firefly
#include "LinearSolver_FF.hpp"

// 全局的 LinearSystemResult 已在 LinearSolver_Eigen.hpp 中定义，此处不再重复

/**
 * 统一线性求解器调度函数
 * 根据数据类型 T 自动选择求解后端：
 *   - 若 T 为 firefly::FFInt，调用 Firefly 有限域求解器
 *   - 否则调用 Eigen 浮点求解器
 * @tparam T 矩阵元素类型
 * @param A 系数矩阵 (rows × cols)
 * @param b 右端项向量 (可选，默认为空表示齐次方程组)
 * @return 统一的 LinearSystemResult<T> 结构
 */
template<typename T>
auto solveLinearSystem(const std::vector<std::vector<T>>& A, const std::vector<T>& b = {}) {
    if constexpr (std::is_same_v<T, firefly::FFInt>) {
        // 调用 Firefly 求解器
        auto ff_res = firefly::solveLinearSystem_Firefly(A, b);
        // 转换为全局的 LinearSystemResult
        LinearSystemResult<T> res;
        res.hasSolution = ff_res.hasSolution;
        res.Mext = std::move(ff_res.Mext);
        res.S = std::move(ff_res.S);
        return res;
    } else {
        // 调用 Eigen 求解器
        return solveLinearSystem_Eigen(A, b);
    }
}

#endif // LINEAR_SOLVER_HPP