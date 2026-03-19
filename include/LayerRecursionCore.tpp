#ifndef LAYER_RECURSION_CORE_TPP
#define LAYER_RECURSION_CORE_TPP

#include "LayerRecursionCore.hpp"
#include "Utilities.hpp"
#include <set>

namespace LayerRecursionCore {

template <typename T>
pair<std::vector<std::vector<T>>, std::vector<T>> assembleLinearSystem(const IBPMatrixE<T> &ibpmat, inhomogTerms<T> &terms, const std::vector<int> &seed, int ncurr, int ne, int nb, int nindep, int nimax_p1) 
{
    int nibp = ibpmat.M1.size();
    int num_eqns = nibp * nb;
    int n_active_vars = ne - max(ncurr, 0);
    int num_vars = n_active_vars * nb + nindep;
    // 1. 初始化矩阵和向量
    std::vector<std::vector<T>> eqnmat(num_eqns, std::vector<T>(num_vars, 0));
    std::vector<T> eqnvec(num_eqns, 0); 
    for (int m = 0; m < nibp; ++m) {
        for (int i = 0; i < nb; ++i) {
            int row_idx = m * nb + i;
            // --- A. 填充未知系数的矩阵部分 (M1 贡献) ---
            for (int j1 = 0; j1 < n_active_vars; ++j1) {
                int target_seed_idx = j1 + std::max(ncurr, 0);
                T factor = static_cast<T>(seed[target_seed_idx] + 1);
                for (int j2 = 0; j2 < nb; ++j2) {
                    // M1[m][target_seed_idx] 是一个 nb*nb 的矩阵
                    // 我们取其第 i 行第 j2 列
                    eqnmat[row_idx][j1 * nb + j2] = ibpmat.M1[m][target_seed_idx][i * nb + j2] * factor;
                }
            }
            // --- B. 填充非齐次项部分 (来自 inhomogTerms) ---
            T* total_row = terms.get_Total_row(m, i);
            // 1. 常数项偏移 (Index 0) -> 移到方程右侧需取负号
            eqnvec[row_idx] = -total_row[0];
            // 2. 独立变量的系数 (Index 1 to nindep) -> 填充到矩阵后部
            for (int j_indep = 0; j_indep < nindep; ++j_indep) eqnmat[row_idx][n_active_vars * nb + j_indep] = total_row[j_indep + 1];
        }
    }
    return {eqnmat, eqnvec};
}

/**
 * @brief 消除冗余：通过追加恒等方程 (1*x = 0) 强制指定变量为零
 * * @tparam T 数据类型 (double, complex等)
 * @param eqnmat 线性方程组左侧矩阵 (将被修改，追加新行)
 * @param eqnvec 线性方程组右侧向量 (将被修改，追加新元素)
 * @param eqnvar 变量映射表，用于查找列索引
 * @param pivot  需要取零的约束列表 {{l, cid}, ...}
 */
template<typename T>
void removeRedundancy(std::vector<std::vector<T>>& eqnmat, 
                      std::vector<T>& eqnvec, 
                      const std::vector<std::array<int, 4>>& eqnvar, 
                      const std::vector<std::array<int, 2>>& pivot) 
{
    if (pivot.empty() || eqnmat.empty()) return;

    // 1. 优化查找：将 pivot 转换为 set，查找复杂度降为 O(log M) 或 O(1)
    // 使用 set 存储 {l, cid} 对，方便快速判定
    std::set<std::pair<int, int>> pivot_set;
    for (const auto& p : pivot) {
        pivot_set.insert({p[0], p[1]});
    }

    int num_cols = eqnmat[0].size();
    
    // 2. 一次性收集所有目标列，避免多次 resize
    std::vector<int> target_cols;
    target_cols.reserve(eqnvar.size()); // 预估最大容量

    for (size_t col = 0; col < eqnvar.size(); ++col) {
        // 构建当前变量的 key
        std::pair<int, int> current_key = {eqnvar[col][1], eqnvar[col][2]};
        
        // 快速查找
        if (pivot_set.count(current_key)) {
            target_cols.push_back(col);
        }
    }

    if (target_cols.empty()) return;

    // 3. 内存优化：一次性分配所有需要的行，避免循环 push_back 导致的反复拷贝
    int old_rows = eqnmat.size();
    int add_rows = target_cols.size();
    
    // 扩展 eqnmat
    eqnmat.resize(old_rows + add_rows);
    // 扩展 eqnvec
    eqnvec.resize(old_rows + add_rows, static_cast<T>(0));

    // 4. 填充新行 (使用移动语义或直接构建)
    // 注意：vector<vector> 的 resize 可能会用默认构造填充，我们需要确保新行大小正确
    for (int i = 0; i < add_rows; ++i) {
        // 直接初始化该行，避免先构造再赋值
        eqnmat[old_rows + i].assign(num_cols, static_cast<T>(0));
        
        // 设置约束位为 1
        int target_col = target_cols[i];
        eqnmat[old_rows + i][target_col] = static_cast<T>(1);
    }
}

template <typename T>
void performBackSubstitution(
    seriesCoefficient<T>& C, int k, int l_start, int l_end, int idxcurr,
    int nindep, int nS, const std::vector<std::vector<T>>& Mext, 
    int row_offset, int ne, int nb, long long (&BINOM)[MAX_VAL][MAX_VAL]
) {
    // --- 1. 单位矩阵检测 (Identity Check) ---
    bool is_identity = (nindep == nS);
    if (is_identity) {
        for (int i = 0; i < nindep; ++i) {
            const std::vector<T>& row = Mext[row_offset + i];
            // 检查常数项是否为 0
            if constexpr (std::is_floating_point_v<T>) {
                if (std::abs(row[0]) > 1e-15) { is_identity = false; break; }
                } else {
                if (row[0] != T(0)) { is_identity = false; break; }
            }
            // 检查对应对角线位置 (i+1) 是否为 1，其他位置是否为 0
            for (int j = 0; j < nS; ++j) {
                T expected = (i == j) ? T(1) : T(0);
                if constexpr (std::is_floating_point_v<T>) {
                    if (std::abs((row[j + 1] - expected)) > 1e-15) { is_identity = false; break; }
                    } else {
                    if ((row[j + 1] - expected) != T(0)) { is_identity = false; break; }
                }
            }
            if (!is_identity) break;
        }
    }

    // --- 2. 如果是单位阵变换，且没有缩减，直接跳过全量回代 ---
    if (is_identity && nS == nindep) {
        return; 
    }

    // --- 3. 执行回代逻辑 ---
    for (int l1 = l_end; l1 >= l_start; --l1) {
        int max_cid = (l1 == l_start) ? idxcurr : (int)BINOM[l1 + ne - 1][ne - 1];
        
        for (int cid = 0; cid < max_cid; ++cid) {
            for (int j = 0; j < nb; ++j) {
                T* sol0_ptr = C(k, l1, cid, j);
                if (!sol0_ptr) continue;

                // 特殊情况优化：如果是单位阵但发生了缩减 (nS < nindep)
                // 此时系数数值不变，只需要抹除末尾
                if (is_identity) {
                    if (nS < nindep) {
                        std::fill(sol0_ptr + nS + 1, sol0_ptr + nindep + 1, T(0));
                    }
                    continue; 
                }

                // 普通情况：执行全量投影计算
                static std::vector<T> sol1;
                sol1.assign(nS + 1, T(0));
                sol1[0] = sol0_ptr[0];

                for (int i1 = 0; i1 < nindep; ++i1) {
                    T old_val = sol0_ptr[i1 + 1];
                    if constexpr (std::is_floating_point_v<T>) {
                        if (std::abs(old_val) < 1e-18) { continue; }
                        } else {
                        if (old_val == T(0)) { continue; }
                    }
                    

                    const std::vector<T>& trans_row = Mext[row_offset + i1];
                    sol1[0] += old_val * trans_row[0];
                    for (int i2 = 0; i2 < nS; ++i2) {
                        sol1[i2 + 1] += old_val * trans_row[i2 + 1];
                    }
                }

                std::copy(sol1.begin(), sol1.end(), sol0_ptr);
                if (nS < nindep) {
                    std::fill(sol0_ptr + nS + 1, sol0_ptr + nindep + 1, T(0));
                }
            }
        }
    }
}

template <typename T>
void updateSeriesCoefficient(
    seriesCoefficient<T>& C,
    const LinearSystemResult<T>& res,
    const std::vector<std::array<int, 4>>& current_eqnvar, // 瞬时映射表
    std::vector<std::array<int, 4>>& indepSet,              // 独立集
    int k, int l, int incre, int idxcurr, int ne, int nb
) {
    if (res.Mext.empty()) return;

    int old_nindep = indepSet.size();
    int new_nindep = res.S.size();
    // 活跃变量占据的列数 = (总列数 - 旧独立元列数)
    int n_active_cols = current_eqnvar.size() - old_nindep;

    // --- 1. 执行 Back Substitution (回代更新) ---
    // 逻辑：更新那些已经存在于 C 中、但基于旧 indepSet 表达的系数
    // 范围：从 l+1 到 incre*k
    if (old_nindep > 0) {
        performBackSubstitution(
            C, k, l + 1, incre * k, idxcurr, 
            old_nindep, new_nindep, res.Mext, n_active_cols, ne, nb, BINOM);
    }

    // --- 2. 写入当前活跃变量的定义 (Injection) ---
    // 这些变量是刚刚通过线性方程组解出来的，其定义直接来自 Mext 的前 n_active_cols 行
    for (int i = 0; i < n_active_cols; ++i) {
        const auto& phys = current_eqnvar[i]; // {k, l+1, cid, j}
        T* ptr = C(phys[0], phys[1], phys[2], phys[3]);
        
        // 写入 [1个特解 + new_nindep 个基础解系数]
        std::copy(res.Mext[i].begin(), res.Mext[i].begin() + new_nindep + 1, ptr);
        
        // 内存深度清理
        if (new_nindep < old_nindep) {
            std::fill(ptr + new_nindep + 1, ptr + old_nindep + 1, T(0));
        }
    }

    // --- 3. 更新 indepSet (基底替换) ---
    indepSet.clear();
    for (int i = 0; i < new_nindep; ++i) {
        int col_idx = res.S[i];
        // 将 eqnvar 中的标识符覆盖到 indepSet 的前部
        indepSet.push_back(current_eqnvar[col_idx]);
    }
}

template<typename T>
void migrateAllTables(
    seriesCoefficient<T>& C,               // 当前 Order 的主容器
    std::vector<seriesCoefficient<T>>& CNew,    // 固化存储容器序列
    const std::vector<std::array<int, 4>>& indepSet, // 当前的独立元集合
    int k, int order, int incre, int ne, int nb, int nimax,
    long long (&BINOM)[MAX_VAL][MAX_VAL]
) {
    int nindep = static_cast<int>(indepSet.size());
    size_t old_size = CNew.size();
    
    // 1. 扩容 CNew。每个新独立元都将获得一个完整的、
    // 以其为基准的 seriesCoefficient 容器，外加一个常数项容器 (i=0)
    CNew.resize(old_size + nindep + 1);

    int nimax_p1 = nimax + 1;

    // 2. 遍历每一个新的分量 (0 到 nindep)
    for (int i = 0; i < nindep + 1; ++i) {
        // 初始化新容器
        CNew[old_size + i] = seriesCoefficient<T>(order, incre, ne, nb, nimax, BINOM);
        seriesCoefficient<T>& target_C = CNew[old_size + i];

        // --- 逻辑 A: 迁移历史阶级 (prev_k < k) ---
        // 只有常数分量 (i=0) 需要承载历史 Order 的完整级数信息
        if (i == 0) {
            for (int prev_k = 0; prev_k < k; ++prev_k) {
                for (int l = 0; l <= incre * prev_k; ++l) {
                    int states = (int)BINOM[l + ne - 1][ne - 1];
                    for (int cid = 0; cid < states; ++cid) {
                        for (int j = 0; j < nb; ++j) {
                            T* src = C(prev_k, l, cid, j);
                            T* dest = target_C(prev_k, l, cid, j);
                            // 历史数据是完整迁移 (包含之前所有的独立元系数)
                            std::copy(src, src + nimax_p1, dest);
                        }
                    }
                }
            }
        }

        // --- 逻辑 B: 处理当前阶级 (order k) ---
        // 将旧容器中第 i 个独立元的系数，提取出来作为新容器的常数项 (index 0)
        for (int l = 0; l <= incre * k; ++l) {
            int states = (int)BINOM[l + ne - 1][ne - 1];
            for (int cid = 0; cid < states; ++cid) {
                for (int j = 0; j < nb; ++j) {
                    T* src = C(k, l, cid, j);
                    T* dest = target_C(k, l, cid, j);
                    
                    // 核心转换：C(k)[i] -> CNew_i(k)[0]
                    dest[0] = src[i];
                    
                    // 此时 dest[1...nimax] 保持为 0，因为在 target_C 
                    // 构造时已初始化。这代表当前 Order k 的结果已被
                    // 固化为常数，不再依赖本层的其他变量。
                }
            }
        }
    }
}

/*template<typename T>
void removeRedundancy(std::vector<std::vector<T>>& eqnmat,
                      std::vector<T>& eqnvec,
                      const std::vector<std::array<int, 4>>& eqnvar,
                      const std::vector<std::array<int, 2>>& pivot)
{
    if (pivot.empty() || eqnmat.empty()) return;

    // 将 pivot 转换为 set 快速查找
    std::set<std::pair<int, int>> pivot_set;
    for (const auto& p : pivot) {
        pivot_set.insert({p[0], p[1]});
    }

    int num_cols = eqnmat[0].size();
    std::vector<int> target_cols;
    target_cols.reserve(eqnvar.size());

    for (size_t col = 0; col < eqnvar.size(); ++col) {
        std::pair<int, int> key = {eqnvar[col][1], eqnvar[col][2]};
        if (pivot_set.count(key)) {
            target_cols.push_back(col);
        }
    }

    if (target_cols.empty()) return;

    int old_rows = eqnmat.size();
    int add_rows = target_cols.size();

    eqnmat.resize(old_rows + add_rows);
    eqnvec.resize(old_rows + add_rows, T(0));

    for (int i = 0; i < add_rows; ++i) {
        eqnmat[old_rows + i].assign(num_cols, T(0));
        eqnmat[old_rows + i][target_cols[i]] = T(1);
    }
}*/

// --- inhomogTerms Implementation (Inline for Template) ---

template <typename T>
void inhomogTerms<T>::add_NMinus_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, std::vector<int>& seed, int k, int l, int nindep_p1) {
    // 边界条件：只有 l 足够小时，k-1 阶才有贡献
    if (l > incre * (k - 1) + 1) return;
    for (int m = 0; m < nibp; ++m) {
        T* dest = get_NMinus(m);
        for (int j = 0; j < ne; ++j) {
            if (seed[j] > 0) {
                // 直接利用 addProductTo 累加： N1[m][j] * C(k-1, l, seed-unit(j))
                T* src = getValuePtrOffSet(C, k - 1, l, seed, -1, j);
                addProductTo(dest, ibpmat.N1[m][j], src, nindep_p1);
            }
        }
    }
}
// ... [包含其他 inhomogTerms 成员函数的实现，照搬原代码] ...
// (此处省略中间实现代码以节省空间，实际文件需包含原文件中 inhomogTerms 类的所有函数体)

template <typename T>
void inhomogTerms<T>::add_NZero_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, std::vector<int>& seed, int k, int l, int nindep_p1) {
    if (l > incre * (k - 1)) return;
    for (int m = 0; m < nibp; ++m) {
        // 利用 mat_buf1 构造临时算子矩阵: F0 + sum(seed[j] * K1[j])
        std::copy(ibpmat.F0[m].begin(), ibpmat.F0[m].end(), mat_buf1.begin());
        for (int j = 0; j < ne; ++j) {
            if (seed[j] >= 1) {
                T s_val = static_cast<T>(seed[j]);
                const std::vector<T>& K1_mj = ibpmat.K1[m][j];
                for (int idx = 0; idx < nb * nb; ++idx) mat_buf1[idx] += K1_mj[idx] * s_val;
            }
        }
        T* src = getValuePtrOffSet(C, k - 1, l, seed);
        addProductTo(get_NZero(m), mat_buf1, src, nindep_p1);
    }
}

template <typename T>
void inhomogTerms<T>::add_NPluMi_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, std::vector<int>& seed, int k, int l, int nindep_p1) {
    if (l > incre * (k - 1)) return;
    for (int m = 0; m < nibp; ++m) {
        T* dest = get_NPluMi(m);
        for (int i = 0; i < ne; ++i) {
            for (int j = 0; j < ne; ++j) {
                if (seed[j] > 0) {
                    // 因子：-(seed[i] + 1)  注意：如果 i==j, seed[i] 还是原来的值
                    T factor = -static_cast<T>(seed[i] + 1);
                    // 获取偏移后的源指针：seed + unit(i) - unit(j)
                    T* src = getValuePtrOffSet(C, k - 1, l, seed, 1, i, -1, j);
                    addScaledProductTo(dest, ibpmat.F2[m][i][j], factor, src, nindep_p1);
                }
            }
        }
    }
}

template <typename T>
void inhomogTerms<T>::add_M1_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, std::vector<int>& seed, int k, int l, int nindep_p1, int ncurr) {
    // M1 贡献通常针对当前 Order k
    for (int m = 0; m < nibp; ++m) {
        T* dest = get_M1(m);
        for (int i = 0; i < ncurr; ++i) {
            T factor = static_cast<T>(seed[i] + 1);
            // 获取当前阶 k 的高层级系数：seed + unit(i)
            T* src = getValuePtrOffSet(C, k, l, seed, 1, i);
            addScaledProductTo(dest, ibpmat.M1[m][i], factor, src, nindep_p1);
        }
    }
}

template <typename T>
void inhomogTerms<T>::add_NPlus_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, std::vector<int>& seed, int k, int l, int nindep_p1, long long (&BINOM)[MAX_VAL][MAX_VAL]) {
    int lennplus = incre * (k - 1) - l;
    if (lennplus < 1) return;
    for (int l1 = 1; l1 <= lennplus; ++l1) {
        T l1_sign = utility::sgn<T>(l1);
        for (int m = 0; m < nibp; ++m) {
            T* dest_ptr = get_NPlus(l1, m);
            // --- 第一部分：处理涉及 K1 和 F2 的单指标 i 贡献 ---
            for (int i = 0; i < ne; ++i) {
                // 重置缓冲区 mat_buf1 用于计算 matNP1
                std::fill(mat_buf1.begin(), mat_buf1.end(), T(0));
                // 1.1 K1 贡献: seed[i] * BINOM[seed[i]+l1][l1+1]
                if (seed[i] >= 1) {
                    T k_factor = static_cast<T>(BINOM[seed[i] + l1][l1 + 1]);
                    const std::vector<T>& K1_mi = ibpmat.K1[m][i];
                    for (int idx = 0; idx < nb * nb; ++idx) mat_buf1[idx] += K1_mi[idx] * k_factor;
                }
                // 1.2 F2 累加贡献: sum_j(F2[m][i][j] * seed[j] * BINOM[seed[i]+l1][l1])
                for (int j = 0; j < ne; ++j) {
                    if (seed[j] >= 1) {
                        T f_factor = l1_sign * static_cast<T>(seed[j] * BINOM[seed[i] + l1][l1]);
                        const std::vector<T>& F2_mij = ibpmat.F2[m][i][j];
                        for (int idx = 0; idx < nb * nb; ++idx) mat_buf1[idx] += F2_mij[idx] * f_factor;
                    }
                }
                // 1.3 执行点积并累加到目标块
                // getValuePtrOffSet 内部会处理 seed[i] += l1 的逻辑
                T* src_ptr = getValuePtrOffSet(C, k - 1, l, seed, l1, i);
                addProductTo(dest_ptr, mat_buf1, src_ptr, nindep_p1);
            }
            // --- 第二部分：处理涉及 F2 的双指标 i, j 耦合项 (l2 循环) ---
            for (int l2 = 1; l2 <= l1 + 1; ++l2) {
                if (l2 == l1) continue; // 排除 l2 == l1 的特殊情况
                T l2_sign = utility::sgn<T>(l2);
                for (int i = 0; i < ne; ++i) {
                    for (int j = 0; j < ne; ++j) {
                        if (seed[j] >= 1) {
                            // 计算耦合因子: sgn(l2) * BINOM[seed[i]+l2][l2] * BINOM[seed[j]+l1-l2][l1-l2+1]
                            T coupling_factor = l2_sign * static_cast<T>(BINOM[seed[i] + l2][l2]) * static_cast<T>(BINOM[seed[j] + l1 - l2][l1 - l2 + 1]);
                            // 直接利用 ibpmat.F2[m][i][j] 和增量 seed 指针进行点积
                            // getValuePtrOffSet 内部处理 seed[i]+=l2, seed[j]+=l1-l2
                            T* src_ptr_coupling = getValuePtrOffSet(C, k - 1, l, seed, l2, i, l1 - l2, j);
                            addScaledProductTo(dest_ptr, ibpmat.F2[m][i][j], coupling_factor, src_ptr_coupling, nindep_p1);
                        }
                    }
                }
            }
        }
    }
}

template <typename T>
void inhomogTerms<T>::add_MPlus_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, std::vector<int>& seed, int k, int l, int nindep_p1, long long (&BINOM)[MAX_VAL][MAX_VAL]) {
    int lenmplus = incre * k - l - 1;
    if (lenmplus < 1) return;
    for (int l1 = 2; l1 <= incre * k - l; ++l1) {
        T sign_l1 = utility::sgn<T>(l1);
        for (int m = 0; m < nibp; ++m) {
            T* dest = get_MPlus(l1, m);
            // 1. K1s, K2s 项
            for (int i = 0; i < ne; ++i) {
                T factor = static_cast<T>(BINOM[seed[i] + l1][l1]);
                // 构造临时矩阵：K1s + sign(l1)*K2s
                for (int idx = 0; idx < nb * nb; ++idx) {
                    mat_buf1[idx] = (ibpmat.K1s[m][i][idx] + ibpmat.K2s[m][i][idx] * sign_l1) * factor;
                }
                T* src = getValuePtrOffSet(C, k, l, seed, l1, i);
                addProductTo(dest, mat_buf1, src, nindep_p1);
            }
            // 2. F2s 耦合项
            for (int l2 = 1; l2 <= l1 - 1; ++l2) {
                T l2_sign = utility::sgn<T>(l2);
                for (int i = 0; i < ne; ++i) {
                    for (int j = 0; j < ne; ++j) {
                        T factor = l2_sign * static_cast<T>(BINOM[seed[i] + l2][l2] * BINOM[seed[j] + l1 - l2][l1 - l2]);
                        T* src = getValuePtrOffSet(C, k, l, seed, l2, i, l1 - l2, j);
                        addScaledProductTo(dest, ibpmat.F2s[m][i][j], factor, src, nindep_p1);
                    }
                }
            }
        }
    }
}

/**
 * 将所有非齐次项分量汇总到 Total 缓冲区
 * @param k 当前 Order
 * @param l 当前 Level
 */
template <typename T>
void inhomogTerms<T>::finalizeTotal(int k, int l) {
    // 1. 获取 Total 缓冲区的起始指针
    T* total_ptr = &data[off_Total];
    // 2. 计算每个方程块的逻辑总长度 (nibp * nb * nimax_p1)
    size_t total_elements = (size_t)nibp * nb * nimax_p1;
    auto aggregate_block = [&](size_t source_offset) {
        T* src_ptr = &data[source_offset];
        for (size_t i = 0; i < total_elements; ++i) total_ptr[i] += src_ptr[i];
    };
    // --- A. 汇总基础项 ---
    aggregate_block(off_NMinus);
    aggregate_block(off_NZero);
    aggregate_block(off_NPluMi);
    aggregate_block(off_M1);
    // --- B. 汇总动态层数项 NPlus ---
    int nplus_limit = incre * (k - 1) - l;
    for (int l1 = 1; l1 <= nplus_limit; ++l1) {
        size_t offset = off_NPlusBase + (size_t)(l1 - 1) * block_len;
        aggregate_block(offset);
    }
    // --- C. 汇总动态层数项 MPlus ---
    int mplus_limit = incre * k - l;
    for (int l1 = 2; l1 <= mplus_limit; ++l1) {
        size_t offset = off_MPlusBase + (size_t)(l1 - 2) * block_len;
        aggregate_block(offset);
    }
}

template <typename T>
inhomogTerms<T>::inhomogTerms(int kmax, int incre, int nibp, int nb, int nimax, int ne) 
    : nibp(nibp), nb(nb), nimax_p1(nimax + 1), incre(incre), ne(ne) {
    // 1. 确定各个算子的最大层数
    nplus_max = incre * (kmax - 1);
    mplus_max = incre * kmax - 1;
    block_len = (size_t)nibp * nb * nimax_p1;
    // 2. 计算各个算子的起始位置 (手动平铺逻辑)
    size_t current_ptr = 0;
    off_NMinus = current_ptr; current_ptr += block_len;
    off_NZero  = current_ptr; current_ptr += block_len;
    off_NPluMi = current_ptr; current_ptr += block_len;
    off_M1     = current_ptr; current_ptr += block_len;
    off_Total  = current_ptr; current_ptr += block_len;
    // 对应 NPlus[l1-1] 和 MPlus[l1-2]
    off_NPlusBase = current_ptr; 
    current_ptr += (nplus_max > 0 ? nplus_max : 0) * block_len;
    off_MPlusBase = current_ptr;
    current_ptr += (mplus_max > 0 ? mplus_max : 0) * block_len;
    // 3. 一次性物理内存分配
    data.assign(current_ptr, T(0));
    mat_buf1.assign(nb * nb, T(0));
    mat_buf2.assign(nb * nb, T(0));
    sum_buf.assign(nb * nb, T(0));
}

template <typename T>
void inhomogTerms<T>::buildAll(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, int k, int l, std::vector<int> &seed, int nindep, int ncurr, long long (&BINOM)[MAX_VAL][MAX_VAL]) {
    this->reset();
    int nindep_p1 = nindep + 1;
    add_NMinus_contribution(ibpmat, C, seed, k, l, nindep_p1);
    add_NZero_contribution(ibpmat, C, seed, k, l, nindep_p1);
    add_NPluMi_contribution(ibpmat, C, seed, k, l, nindep_p1);
    add_NPlus_contribution(ibpmat, C, seed, k, l, nindep_p1, BINOM);
    add_M1_contribution(ibpmat, C, seed, k, l, nindep_p1, ncurr);
    add_MPlus_contribution(ibpmat, C, seed, k, l, nindep_p1, BINOM);
    finalizeTotal(k, l);
}

}

#endif