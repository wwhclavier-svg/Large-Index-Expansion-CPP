#ifndef LAYER_RECURSION_CORE_HPP
#define LAYER_RECURSION_CORE_HPP

#include <vector>
#include <array>
#include <cstddef>

#include "Combinatorics.hpp"          // 组合数、种子生成、索引转换
#include "SeriesCoefficient.hpp"      // seriesCoefficient 类
#include "IBPMatrixLoader_Binary.hpp" // IBPMatrixE 结构
#include "LinearSolver.hpp"           // 统一求解器接口"
#include "Utilities.hpp"

// 核心辅助函数声明（模板实现放在 .tpp 文件中）
namespace LayerRecursionCore {

template<typename T> class inhomogTerms;

/**
 * 组装线性方程组
 * @param ibpmat  IBP 矩阵数据
 * @param terms   已构建的非齐次项
 * @param seed    当前种子向量
 * @param ncurr   最后一个非零种子索引
 * @param ne      维度
 * @param nb      块大小
 * @param nindep  当前独立变量个数
 * @param nimax_p1 通解个数 + 1
 * @return 方程组 (eqnmat, eqnvec)
 */
template<typename T>
std::pair<std::vector<std::vector<T>>, std::vector<T>>
assembleLinearSystem(const IBPMatrixE<T>& ibpmat,
                     inhomogTerms<T>& terms,
                     const std::vector<int>& seed,
                     int ncurr, int ne, int nb,
                     int nindep, int nimax_p1);

/**
 * 回代更新：将新解表达转换到旧独立基下
 * @param C          系数容器（将被修改）
 * @param k          当前阶数
 * @param l_start    起始层（l+1）
 * @param l_end      结束层（incre*k）
 * @param idxcurr    当前种子的索引
 * @param nindep     旧独立变量个数
 * @param nS         新独立变量个数
 * @param Mext       变换矩阵（来自求解结果）
 * @param row_offset 活跃变量列数
 * @param ne         维度
 * @param nb         块大小
 * @param BINOM      组合数表
 */
template<typename T>
void performBackSubstitution(seriesCoefficient<T>& C,
                              int k, int l_start, int l_end,
                              int idxcurr, int nindep, int nS,
                              const std::vector<std::vector<T>>& Mext,
                              int row_offset, int ne, int nb,
                              long long (&BINOM)[MAX_VAL][MAX_VAL]);

/**
 * 更新当前解的系数：将刚解出的活跃变量写入，并更新独立变量集合
 * @param C                系数容器（将被修改）
 * @param res              求解结果
 * @param current_eqnvar   当前方程变量映射表
 * @param indepSet         独立变量集合（将被更新）
 * @param k                当前阶数
 * @param l                当前层
 * @param incre            增量
 * @param idxcurr          当前种子索引
 * @param ne               维度
 * @param nb               块大小
 */
template<typename T>
void updateSeriesCoefficient(seriesCoefficient<T>& C,
                             const LinearSystemResult<T>& res,
                             const std::vector<std::array<int, 4>>& current_eqnvar,
                             std::vector<std::array<int, 4>>& indepSet,
                             int k, int l, int incre,
                             int idxcurr, int ne, int nb);

/**
 * 迁移所有解表：根据新独立变量集生成多个分支容器
 * @param C          当前容器（将被移动）
 * @param CNew       目标容器序列（将被扩充）
 * @param indepSet   新独立变量集合
 * @param k          当前阶数
 * @param order      最大阶数
 * @param incre      增量
 * @param ne         维度
 * @param nb         块大小
 * @param nimax      最大通解个数
 * @param BINOM      组合数表
 */
template<typename T>
void migrateAllTables(seriesCoefficient<T>& C,
                      std::vector<seriesCoefficient<T>>& CNew,
                      const std::vector<std::array<int, 4>>& indepSet,
                      int k, int order, int incre,
                      int ne, int nb, int nimax,
                      long long (&BINOM)[MAX_VAL][MAX_VAL]);

// 生成变量映射表（非模板，普通函数）
void equationVariable(std::vector<std::array<int,4>>& eqnvar,
                      int order, int level,
                      std::vector<int>& seed,
                      int nb, int ne,
                      const std::vector<std::array<int,4>>& indepSet);

// 消除冗余（模板函数）
template<typename T>
void removeRedundancy(std::vector<std::vector<T>>& eqnmat,
                      std::vector<T>& eqnvec,
                      const std::vector<std::array<int, 4>>& eqnvar,
                      const std::vector<std::array<int, 2>>& pivot);

template <typename T> 
class inhomogTerms {
private:
    std::vector<T> data;
    int nibp, nb, ne, nimax_p1, incre;
    int nplus_max, mplus_max;
    size_t off_NMinus, off_NZero, off_NPluMi, off_M1, off_Total;
    size_t off_NPlusBase, off_MPlusBase;
    size_t block_len;
    std::vector<T> mat_buf1; 
    std::vector<T> mat_buf2; 
    std::vector<T> sum_buf;

    inline void addProductTo(T* dest_ptr, const std::vector<T>& mat, T* src_ptr, int nindep_p1) {
        for (int r = 0; r < nb; ++r) {
            T* row_dest = dest_ptr + (size_t)r * nimax_p1;
            for (int c = 0; c < nb; ++c) {
                T scalar = mat[r * nb + c];
                if (scalar == 0) continue; 
                T* col_src = src_ptr + (size_t)c * nimax_p1;
                for (int i = 0; i < nindep_p1; ++i) row_dest[i] += scalar * col_src[i];
            }
        }
    }

    inline void addScaledProductTo(T* dest, const vector<T>& mat, T factor, T* src, int nindep_p1) {
        if (factor == 0) return;
        for (int r = 0; r < nb; ++r) {
            T* row_d = dest + (size_t)r * nimax_p1;
            for (int c = 0; c < nb; ++c) {
                T scalar = mat[r * nb + c] * factor; 
                if (scalar == 0) continue;
                T* col_s = src + (size_t)c * nimax_p1;
                for (int i = 0; i < nindep_p1; ++i) row_d[i] += scalar * col_s[i];
            }
        }
    }

    // 内部函数声明
    void add_NMinus_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, vector<int>& seed, int k, int l, int nindep_p1);
    void add_NZero_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, vector<int>& seed, int k, int l, int nindep_p1);
    void add_NPluMi_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, vector<int>& seed, int k, int l, int nindep_p1);
    void add_M1_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, vector<int>& seed, int k, int l, int nindep_p1, int ncurr);
    void add_NPlus_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, vector<int>& seed, int k, int l, int nindep_p1, long long (&BINOM)[MAX_VAL][MAX_VAL]);
    void add_MPlus_contribution(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, vector<int>& seed, int k, int l, int nindep_p1, long long (&BINOM)[MAX_VAL][MAX_VAL]);
    void finalizeTotal(int k, int l);

public:
    inhomogTerms(int kmax, int incre, int nibp, int nb, int nimax, int ne);
    void reset() { std::fill(data.begin(), data.end(), T(0)); }

    inline T* get_NMinus(int m) { return &data[off_NMinus + (size_t)m * nb * nimax_p1]; }
    inline T* get_NZero(int m)  { return &data[off_NZero  + (size_t)m * nb * nimax_p1]; }
    inline T* get_NPluMi(int m) { return &data[off_NPluMi + (size_t)m * nb * nimax_p1]; }
    inline T* get_M1(int m)     { return &data[off_M1     + (size_t)m * nb * nimax_p1]; }
    inline T* get_Total(int m)  { return &data[off_Total  + (size_t)m * nb * nimax_p1]; }
    inline T* get_NPlus(int l1, int m) { return &data[off_NPlusBase + (size_t)(l1 - 1) * block_len + (size_t)m * nb * nimax_p1]; }
    inline T* get_MPlus(int l1, int m) { return &data[off_MPlusBase + (size_t)(l1 - 2) * block_len + (size_t)m * nb * nimax_p1]; }
    inline T* get_Total_row(int m, int j) { return &data[off_Total + ((size_t)m * nb + j) * nimax_p1]; }
    
    void buildAll(const IBPMatrixE<T>& ibpmat, seriesCoefficient<T>& C, int k, int l, vector<int> &seed, int nindep, int ncurr, long long (&BINOM)[MAX_VAL][MAX_VAL]);
};

} // namespace LayerRecursionCore



// 包含模板实现（若使用分离编译模型）
#include "LayerRecursionCore.tpp"

#endif // LAYER_RECURSION_CORE_HPP