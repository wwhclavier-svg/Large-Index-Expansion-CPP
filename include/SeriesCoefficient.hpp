#ifndef SERIES_COEFFICIENT_HPP
#define SERIES_COEFFICIENT_HPP

#include <vector>
#include "Combinatorics.hpp"

template<typename T>
class seriesCoefficient {
private:
    std::vector<T> data;
    // offsets[k][l] 存储逻辑层 (k, l) 在一维 data 中的起始下标
    std::vector<std::vector<size_t>> offsets;

    int kmax;
    int incre;
    int ne;    // 向量指标维数
    int nb;    // 矩阵 block 维度 (nb*nb，针对 nb=1 优化)
    int nimax; // 通解个数

public:
    seriesCoefficient()  : kmax(0), incre(0), ne(0), nb(0), nimax(0) {}
    /**
     * @param kmax  最大 Order
     * @param incre 每级 Order 对应的 Level 增量 (lmax = incre * k)
     * @param ne    对应组合数公式中的参数 (BINOM[l+ne-1][ne-1])
     * @param nb    对应维度 nb
     * @param nimax 对应通解个数
     * @param BINOM 外部传入的组合数预计算表 (二维数组)
     */
    seriesCoefficient(int kmax, int incre, int ne, int nb, int nimax, long long (&BINOM)[MAX_VAL][MAX_VAL]) 
        : kmax(kmax), incre(incre), ne(ne), nb(nb), nimax(nimax) {

        offsets.resize(kmax + 1);
        size_t total_elements = 0;

        // --- 1. 预计算偏移量表 (一次性完成) ---
        for (int k = 0; k <= kmax; ++k) {
            int lmax = incre * k;
            offsets[k].resize(lmax + 1);
            
            for (int l = 0; l <= lmax; ++l) {
                // 记录当前 (k, l) 组合在 data 数组中的绝对起始位置
                offsets[k][l] = total_elements;
                
                // 计算当前层级的 Seed (cid) 数量
                long long states = BINOM[l + ne - 1][ne - 1];
                
                // 累加总大小：当前层 seeds 数 * nb * (nimax + 1)
                // 使用 nimax + 1 是为了兼容索引从 0 到 nimax 的访问
                total_elements += (size_t)states * nb * (nimax + 1);
            }
        }

        // --- 2. 物理内存分配 ---
        // 使用 assign 确保内存被分配并初始化为 0
        data.assign(total_elements, T(0));

        // 打印内存预警（可选）
        // std::cout << "[Storage] Allocated " << (total_elements * sizeof(T)) / (1024.0 * 1024.0) << " MB." << std::endl;
    }

    /**
     * 5维索引访问：C(k, l, cid, j, i)
     * 时间复杂度：O(1)，仅包含 1 次向量查表和 3 次乘法
     */
    inline T& operator()(int k, int l, int cid, int j, int i) {
        // 核心公式：层起始偏移 + cid偏移 + j偏移 + i偏移
        return data[offsets[k][l] + (size_t)cid * nb * (nimax + 1) + j * (nimax + 1) + i];
    }

    inline const T& operator()(int k, int l, int cid, int j, int i) const {
        return data[offsets[k][l] + (size_t)cid * nb * (nimax + 1) + j * (nimax + 1) + i];
    }
        

    /**
     * 4维层级指针访问：C(k, l, cid, j)
     * 返回指向连续 i 维数据的指针，方便进行批量计算或 memcpy
     */
    inline T* operator()(int k, int l, int cid, int j) {
        return &data[offsets[k][l] + (size_t)cid * nb * (nimax + 1) + j * (nimax + 1)];
    }

    /**
     * 3维层级指针访问：C(k, l, cid)
     * 返回指向该 seed 下所有 j 和 i 数据的起始指针
     */
    inline T* operator()(int k, int l, int cid) {
        return &data[offsets[k][l] + (size_t)cid * nb * (nimax + 1)];
    }

    // 获取底层一维数组的引用，方便整体操作
    std::vector<T>& raw_data() { return data; }
    const std::vector<T>& raw_data() const { return data; }

    // 重置所有系数为 0
    void clear() {
        std::fill(data.begin(), data.end(), T(0));
    }

    // 获取总元素个数
    size_t total_size() const { return data.size(); }
    //int basis_size() { return nb; } 
    int basis_size() const { return nb; }

    int getKmax() const { return kmax; }
    int getIncre() const { return incre; }
    int getNe() const { return ne; }
    int getNb() const { return nb; }
    int getNimax() const { return nimax; }
};

// ==========================================
// 便捷指针访问函数（用于 seriesCoefficient）
// ==========================================

/**
 * 获取偏移后的系数指针（单指标偏移）
 * @tparam T 数据类型
 * @param C          seriesCoefficient 对象
 * @param order      阶数 k
 * @param old_level  原层级 l
 * @param seed       当前种子向量（将被修改后恢复）
 * @param k_inc      种子增量（可为负）
 * @param di         要增加的维度
 * @return 指向新位置 (order, new_level, new_cid) 处整个数据块的指针
 */
template <typename T>
inline T* getValuePtrOffSet(seriesCoefficient<T>& C, int order, int old_level,
                            std::vector<int>& seed, int k_inc, int di) {
    seed[di] += k_inc;
    int new_level = old_level + k_inc;
    int cid = getIndex(seed, new_level);
    T* result_ptr = C(order, new_level, cid);
    seed[di] -= k_inc;
    return result_ptr;
}

/**
 * 获取偏移后的系数指针（双指标偏移）
 * @tparam T 数据类型
 * @param C          seriesCoefficient 对象
 * @param order      阶数 k
 * @param old_level  原层级 l
 * @param seed       当前种子向量（将被修改后恢复）
 * @param k_inc      第一个维度的增量
 * @param di         第一个维度索引
 * @param l_inc      第二个维度的增量
 * @param dj         第二个维度索引
 * @return 指向新位置 (order, new_level, new_cid) 处整个数据块的指针
 */
template <typename T>
inline T* getValuePtrOffSet(seriesCoefficient<T>& C, int order, int old_level,
                            std::vector<int>& seed, int k_inc, int di,
                            int l_inc, int dj) {
    seed[di] += k_inc;
    seed[dj] += l_inc;
    int new_level = old_level + k_inc + l_inc;
    int cid = getIndex(seed, new_level);
    T* result_ptr = C(order, new_level, cid);
    seed[di] -= k_inc;
    seed[dj] -= l_inc;
    return result_ptr;
}

/**
 * 获取当前种子对应的系数指针（无偏移）
 * @tparam T 数据类型
 * @param C          seriesCoefficient 对象
 * @param order      阶数 k
 * @param old_level  层级 l
 * @param seed       种子向量（不会被修改）
 * @return 指向当前位置 (order, old_level, cid) 处整个数据块的指针
 */
template <typename T>
inline T* getValuePtrOffSet(seriesCoefficient<T>& C, int order, int old_level,
                            const std::vector<int>& seed) {
    int cid = getIndex(seed, old_level);
    return C(order, old_level, cid);
}

// 单指标偏移（只读版本）
template <typename T>
inline const T* getValuePtrOffSet(const seriesCoefficient<T>& C, int order, int old_level,
                                  std::vector<int>& seed, int k_inc, int di) {
    seed[di] += k_inc;
    int new_level = old_level + k_inc;
    int cid = getIndex(seed, new_level);
    const T* result_ptr = C(order, new_level, cid);
    seed[di] -= k_inc;
    return result_ptr;
}

// 双指标偏移（只读版本）
template <typename T>
inline const T* getValuePtrOffSet(const seriesCoefficient<T>& C, int order, int old_level,
                                  std::vector<int>& seed, int k_inc, int di,
                                  int l_inc, int dj) {
    seed[di] += k_inc;
    seed[dj] += l_inc;
    int new_level = old_level + k_inc + l_inc;
    int cid = getIndex(seed, new_level);
    const T* result_ptr = C(order, new_level, cid);
    seed[di] -= k_inc;
    seed[dj] -= l_inc;
    return result_ptr;
}


// 无偏移版本（只读）
template <typename T>
inline const T* getValuePtrOffSet(const seriesCoefficient<T>& C, int order, int old_level,
                                  const std::vector<int>& seed) {
    int cid = getIndex(seed, old_level);
    return C(order, old_level, cid);
}

#endif