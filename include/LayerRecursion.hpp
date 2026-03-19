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
    // initialize coefficient: C[order][level][index][symb][nsol]
    int nimax = 4;

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
                    int ncurr = lastNonZero(seed);
                    int idxcurr = getIndexOffSet(l,seed,1,max(ncurr,0)); 
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
std::vector<std::vector<seriesCoefficient<T>>> batchProcessRecursion(const std::vector<IBPMatrixE<T>>& allMatrices, int order, int incre=2) 
{
    std::vector<std::vector<seriesCoefficient<T>>> allResults;
    allResults.reserve(allMatrices.size());

    std::cout << "开始批量计算，总计: " << allMatrices.size() << " 组数据..." << std::endl;
    for (size_t i = 0; i < allMatrices.size(); ++i) {
        const auto& mat = allMatrices[i];
        std::cout << "   ......   Reg. " << i+1 << "/"<< allMatrices.size() << "   nb = " << mat.nb << "    ......  " << std::endl;   
        auto result = layerRecursion<T>(mat, mat.ne, mat.nb, mat.nibp, order, incre);
        allResults.push_back(std::move(result));
    }
    return allResults;
}

#endif