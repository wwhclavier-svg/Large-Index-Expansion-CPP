#ifndef LAYER_RECURSION_WRAPPER_HPP
#define LAYER_RECURSION_WRAPPER_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>

#include "IBPMatrixLoader_Binary.hpp"
#include "SeriesCoefficient.hpp"
#include "LayerRecursionCore.hpp"
#include "Combinatorics.hpp"
#include "LinearSolver.hpp"
#include "IncrementDetector.hpp"


/**
 * 验证展开结果：对特解部分 (index=0) 逐 seed 检查 M1*C + TotalInhomog == 0。
 * @return true 若所有方程残差为零，否则 false
 */
template<typename T>
bool validateExpansion(const std::vector<seriesCoefficient<T>>& branches,
                       const IBPMatrixE<T>& mat, int order, int incre,
                       int ne, int nb)
{
    if (branches.empty()) return false;
    const auto& C = branches[0];
    int nibp = mat.M1.size();

    // 计算 nimax
    int nimax = 4;
    for (int k = 0; k <= order; ++k) {
        int lmax = incre * k;
        for (int l = 1; l <= lmax; ++l)
            nimax = std::max(nimax, (int)BINOM[l + ne - 1][ne - 1] * nb);
    }

    // 只检查 k=1（低阶相容是必要条件，满足则高阶自动相容）
    for (int k = 1; k <= 1; ++k) {
        for (int l = incre * k - 1; l >= 0; --l) {
            long long nseeds = BINOM[l + ne - 1][ne - 1];
            for (long long cid = 0; cid < nseeds; ++cid) {
                std::vector<int> seed = readIndex(cid, l, ne);
                int ncurr = lastNonZero(seed) + 1;  // MMA 1-indexed

                // 重建 inhomog total，nindep=0 只需特解
                LayerRecursionCore::inhomogTerms<T> terms(k, incre, nibp, nb, nimax, ne);
                // terms.buildAll 需要非 const 引用，但实际只读 k,l,seed
                seriesCoefficient<T>& C_ref = const_cast<seriesCoefficient<T>&>(C);
                terms.buildAll(mat, C_ref, k, l, seed, 0, ncurr, BINOM);

                // 对每个 IBP 方程检查残差
                int n_active = ne - std::max(ncurr - 1, 0);
                for (int m = 0; m < nibp; ++m) {
                    for (int j = 0; j < nb; ++j) {
                        T residual = terms.get_Total_row(m, j)[0]; // 常数部分
                        // M1 * C_active (index=0 特解)
                        for (int j1 = 0; j1 < n_active; ++j1) {
                            int target = j1 + std::max(ncurr - 1, 0);
                            T factor = static_cast<T>(seed[target] + 1);
                            const T* src = getValuePtrOffSet(C, k, l, seed, 1, target);
                            if (src) {
                                for (int j2 = 0; j2 < nb; ++j2) {
                                    residual += mat.M1[m][target][j * nb + j2] * factor * src[j2 * (nimax + 1)];
                                }
                            }
                        }
                        if (residual != T(0)) return false;
                    }
                }
            }
        }
    }
    return true;
}


/**
 * 对单个 IBP 矩阵执行层递归，计算所有解。
 * @tparam T 数据类型（double 或 firefly::FFInt）
 * @param ibpmat  IBP 矩阵数据
 * @param ne      外部变量数（多重指标维度）
 * @param nb      矩阵块大小
 * @param nibp    积分个数（IBP 方程个数）
 * @param order   最大展开阶数
 * @param incre   每阶递增的层数（通常为 1）
 * @return 所有解（每个解是一个 seriesCoefficient<T> 容器）
 */
template<typename T>
auto layerRecursion(const struct IBPMatrixE<T> &ibpmat, int ne, int nb, int nibp, int order, int incre=2)
{
    // auto-detect optimal increment if default
    if (incre == 2) {
        int detected = detectIncrement(ibpmat);
        // only override if we get a clear signal
        if (detected == 1 || detected == 3) incre = detected;
    }

    // initialize coefficient: C[order][level][index][symb][nsol]
    // Dynamically compute nimax to avoid buffer overruns in large problems.
    int nimax = 4;
    for (int k = 0; k <= order; ++k) {
        int lmax = incre * k;
        for (int l = 1; l <= lmax; ++l) {
            nimax = std::max(nimax, static_cast<int>(BINOM[l + ne - 1][ne - 1]) * nb);
        }
    }

    // 1. 初始化系数容器
    seriesCoefficient<T> C0(order,incre,ne,nb,nimax,BINOM);
    std::vector<seriesCoefficient<T>> CTable(1, std::move(C0));
    std::vector<seriesCoefficient<T>> CNew(0);

    // 2. 预生成所有种子
    std::vector<std::vector<std::vector<int>>> seedlist = seedGenerator(incre*order-1,ne);

    // 3. 初始化非齐次项构建器
    LayerRecursionCore::inhomogTerms<T> terms(order, incre, nibp, nb, nimax, ne);
    
    // 4. 初始化独立变量集合和主元列表
    std::vector<std::array<int,4>> indepSet;
    std::vector<std::array<int,2>> pivots = {{0,0}}; 

    // 5. 设置初始条件：C(0,0,0,0)[0] = 1
    CTable[0](0,0,0,0)[0] = T(1);

    auto total_start = std::chrono::high_resolution_clock::now();

    // 6. 主循环：逐阶推进
    for(int k = 1; k <= order; ++k)
    {
        std::cout << "   ======   order: " << k << "   ======"<< std::endl; 
        auto start_k = std::chrono::high_resolution_clock::now();
        
        CNew.clear(); 
        std::vector<std::array<int, 2>> new_pivots;

        for(auto &C : CTable) {
            indepSet.clear();
            for(int l = incre*k - 1; l >= 0; --l) {
                //std::cout << "   =======   level: " << l << "   ======" << std::endl;
                //std::cout << "#seeds: " << seedlist[l].size() << std::endl;
                for(auto seed : seedlist[l]) {
                    int nindep = static_cast<int>(indepSet.size());
                    int ncurr = lastNonZero(seed) + 1;  // MMA 1-indexed
                    int idxcurr = getIndexOffSet(l,seed,1,max(ncurr - 1,0));
                    // 构建非齐次项
                    terms.buildAll(ibpmat, C, k, l, seed, nindep, ncurr, BINOM);
                    // 组装方程组
                    auto [eqnmat, eqnvec] = LayerRecursionCore::assembleLinearSystem(ibpmat, terms, seed, ncurr, ne, nb, nindep, nimax + 1);
                    // 生成变量映射表
                    std::vector<std::array<int, 4>> eqnvar;
                    LayerRecursionCore::equationVariable(eqnvar, k, l, seed, nb, ne, indepSet);
                    // 消除冗余
                    LayerRecursionCore::removeRedundancy(eqnmat, eqnvec, eqnvar, pivots);
                    // 求解
                    auto res = solveLinearSystem(eqnmat, eqnvec); 
                    // 更新系数
                    if (res.hasSolution) {
                        LayerRecursionCore::updateSeriesCoefficient(C, res, eqnvar, indepSet, k, l, incre, idxcurr, ne, nb);
                    }
                } // Seed-loop Over!.
            } // Level-loop Over!
            // 根据 indepSet 迁移分支
            if(!indepSet.empty()) {
                // 情况 1：产生了新分支 (n 个新解)
                LayerRecursionCore::migrateAllTables(C, CNew, indepSet, k, order, incre, ne, nb, nimax, BINOM);
                for(const auto& item : indepSet) {
                    new_pivots.push_back({item[1], item[2]});
                }
            } 
            else {
                // 情况 2：没有产生新分支 (唯一解)
                CNew.push_back(std::move(C)); 
            }
        }
        // 更新主元列表
        if(!new_pivots.empty()) {
            pivots.insert(pivots.end(), new_pivots.begin(), new_pivots.end());
            std::sort(pivots.begin(), pivots.end());
            pivots.erase(unique(pivots.begin(), pivots.end()), pivots.end());
        }
        CTable = std::move(CNew); 

        auto end_k = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_k - start_k;
        std::cout << "order "<< k << "  ......  " << std::fixed << std::setprecision(4) << diff.count() << " s. " << std::endl;
    }
        
    auto total_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_diff = total_end - total_start;
    std::cout << "Total time: " << total_diff.count() << " s." << std::endl;
    std::cout << "#sol = " << CTable.size() << std::endl;

    return CTable;
}


/**
 * 批量处理多个 IBP 矩阵的层递归。
 * @tparam T 数据类型
 * @param allMatrices IBP 矩阵列表
 * @param order 最大阶数
 * @param incre 增量
 * @return 每个矩阵对应的所有解列表
 */
template<typename T>
std::vector<std::vector<seriesCoefficient<T>>> batchProcessRecursion(const std::vector<IBPMatrixE<T>>& allMatrices, int order, int incre_hint=2) 
{
    std::vector<std::vector<seriesCoefficient<T>>> allResults;
    allResults.reserve(allMatrices.size());

    std::cout << "开始批量计算，总计: " << allMatrices.size() << " 组数据..." << std::endl;
    for (size_t i = 0; i < allMatrices.size(); ++i) {
        const auto& mat = allMatrices[i];
        std::cout << "   ......   Reg. " << i+1 << "/"<< allMatrices.size() << "   nb = " << mat.nb << "    ......  " << std::endl;

        // 对每个 matrix 单独检测 incre
        int base = incre_hint;
        bool incompat = false;
        if (base == 2) {
            int detected = detectIncrement(mat);
            base = detected;
            incompat = (detected == 3);
        }

        std::vector<seriesCoefficient<T>> result;
        int trial = (incompat && base < 3) ? 3 : base;
        int max_attempts = incompat ? 10 : 1;
        bool succeeded = false;

        for (int attempt = 0; attempt < max_attempts; attempt++, trial++) {
            result = layerRecursion<T>(mat, mat.ne, mat.nb, mat.nibp, order, trial);

            if (!incompat) {
                succeeded = true;
                break;
            }

            if (validateExpansion<T>(result, mat, order, trial, mat.ne, mat.nb)) {
                succeeded = true;
                std::cout << "  incre=" << trial << " validated OK" << std::endl;
                break;
            }

            std::cout << "!!! Incompatible at current incre. Restarting with incre = " << trial + 1 << std::endl;
        }

        if (!succeeded) {
            std::cerr << "  ERROR: failed to find valid incre for matrix " << i+1 << std::endl;
        }
        allResults.push_back(std::move(result));
    }
    return allResults;
}

#endif