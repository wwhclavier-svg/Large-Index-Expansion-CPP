#ifndef UNIFIED_STORAGE_HPP
#define UNIFIED_STORAGE_HPP

#include <vector>
#include <map>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <utility>

#include "Combinatorics.hpp"

// ==========================================
// 2. 地址结构定义
// ==========================================

struct TensorAddress {
    uint32_t secNo;
    int dot, rank;
    long long dotNo, rankNo;
    int pd_dim, isp_dim;
};

inline TensorAddress rawToAddress(const std::vector<int>& indices) {
    TensorAddress addr;
    addr.secNo = 0;
    std::vector<int> pdSet, ispSet;
    int n = indices.size();

    for (int i = 0; i < n; ++i) {
        if (indices[i] > 0) {
            addr.secNo |= (1 << i);
            pdSet.push_back(indices[i]);
        } else {
            ispSet.push_back(-indices[i]);
        }
    }

    addr.dot = std::accumulate(pdSet.begin(), pdSet.end(), 0);
    addr.rank = std::accumulate(ispSet.begin(), ispSet.end(), 0);
    addr.pd_dim = pdSet.size();
    addr.isp_dim = ispSet.size();
    addr.dotNo = getIndex(pdSet, addr.dot);
    
    // 保持与原逻辑一致，ispSet 逆序计算 rankNo
    std::vector<int> revIspSet = ispSet;
    std::reverse(revIspSet.begin(), revIspSet.end());
    addr.rankNo = getIndex(revIspSet, addr.rank); 

    return addr;
}

// ==========================================
// 3. 核心存储结构：GenericStorage
// ==========================================

struct BlockKey {
    uint32_t secNo;
    int dot, rank;
    // 用于 DualIndexStorage 的额外维度 (Level), IndexStorage 时设为 -1
    int level; 

    bool operator<(const BlockKey& o) const {
        return std::tie(secNo, dot, rank, level) < std::tie(o.secNo, o.dot, o.rank, o.level);
    }
};

/**
 * @brief 通用数据块
 * 统一管理 IndexStorage 和 DualIndexStorage 的内存
 */
template <typename T>
class GenericBlock {
public:
    long long dot_cap, rank_cap, exp_cap;
    int element_size; // 每个格点储存的数据大小 (Matrix: N*N, std::vector: K)
    std::vector<T> data;

    // 构造函数
    GenericBlock(int pd_dim, int isp_dim, int dot, int rank, 
                 int exp_dim, int level, int elem_size) 
        : element_size(elem_size) {
        
        dot_cap = getCapacity(pd_dim, dot);
        rank_cap = getCapacity(isp_dim, rank);
        
        // 如果 exp_dim <= 0，说明是 IndexStorage模式，exp_cap 视为 1
        exp_cap = (exp_dim > 0) ? getCapacity(exp_dim, level) : 1;

        long long total_items = dot_cap * rank_cap * exp_cap;
        data.resize(total_items * element_size, T(0));
    }

    // 通用 Setter (通过指针写入，避免拷贝)
    void set(long long dNo, long long rNo, long long eNo, const std::vector<T>& val) {
        // Index = ( (d * R_cap + r) * E_cap + e ) * Size
        long long idx = ((dNo * rank_cap + rNo) * exp_cap + eNo) * element_size;
        for(int i=0; i<element_size; ++i) data[idx + i] = val[i];
    }
    
    // 获取数据指针
    T* getPtr(long long dNo, long long rNo, long long eNo = 0) {
        long long idx = ((dNo * rank_cap + rNo) * exp_cap + eNo) * element_size;
        return &data[idx];
    }
};

/**
 * @brief 模板化存储管理类
 * Mode 0: IndexStorage (Key: alpha)
 * Mode 1: DualIndexStorage (Key: alpha, beta)
 */
template <typename T, int Mode>
class UnifiedStorage {
private:
    std::map<BlockKey, GenericBlock<T>> table;
    int total_vars; // 对应 ne
    int default_elem_size; // 默认元素大小

public:
    UnifiedStorage(int ne, int elem_size) : total_vars(ne), default_elem_size(elem_size) {}

    // 获取或创建块 (Index 模式)
    GenericBlock<T>& getBlock(const TensorAddress& addr, int level = -1) {
        int l_key = (Mode == 0) ? -1 : level;
        BlockKey key = {addr.secNo, addr.dot, addr.rank, l_key};
        
        if (table.find(key) == table.end()) {
            // Index模式下 exp_dim=0, Dual模式下 exp_dim=total_vars
            int e_dim = (Mode == 0) ? 0 : total_vars;
            table.emplace(std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(addr.pd_dim, addr.isp_dim, addr.dot, addr.rank, 
                                      e_dim, l_key, default_elem_size));
        }
        return table.at(key);
    }

    // 写入 (重载：支持 Index 和 Dual Index)
    void insert(const std::vector<int>& alpha, const std::vector<T>& val, 
                const std::vector<int>& beta = {}) {
        TensorAddress addr = rawToAddress(alpha);
        
        long long expNo = 0;
        int level = -1;

        if constexpr (Mode == 1) {
            level = std::accumulate(beta.begin(), beta.end(), 0);
            expNo = getIndex(beta, level);
        }

        GenericBlock<T>& block = getBlock(addr, level);
        block.set(addr.dotNo, addr.rankNo, expNo, val);
    }

    // 读取指针
    T* retrieve(const std::vector<int>& alpha, const std::vector<int>& beta = {}) {
        TensorAddress addr = rawToAddress(alpha);
        
        long long expNo = 0;
        int level = -1;
        if constexpr (Mode == 1) {
            level = std::accumulate(beta.begin(), beta.end(), 0);
            expNo = getIndex(beta, level);
        }

        BlockKey key = {addr.secNo, addr.dot, addr.rank, level};
        if (table.find(key) == table.end()) return nullptr;
        
        return table.at(key).getPtr(addr.dotNo, addr.rankNo, expNo);
    }

    // 检查是否存在 (用于递归基)
    bool exists(const std::vector<int>& alpha) {
        TensorAddress addr = rawToAddress(alpha);
        BlockKey key = {addr.secNo, addr.dot, addr.rank, -1};
        // 注意：这里只检查 Block 是否存在，不检查具体 slot 是否填入
        // 实际应用需配合 computed flag 使用
        return table.find(key) != table.end();
    }

    // 将存储的所有 b_{alpha, beta} 对应的系数导出为一个巨大的矩阵
    // 行索引：(i, k_total)，列索引：(alpha_idx, beta_idx)
    std::vector<std::vector<T>> toFinalMatrix(
        const std::vector<std::vector<int>>& all_alphas,
        const std::vector<std::vector<int>>& all_betas,
        int nb, int total_k_max) {
        
        int row_count = nb * total_k_max;
        int col_count = all_alphas.size() * all_betas.size();
        std::vector<std::vector<T>> matrix(row_count, std::vector<T>(col_count, T(0)));

        for (size_t a_idx = 0; a_idx < all_alphas.size(); ++a_idx) {
            for (size_t b_idx = 0; b_idx < all_betas.size(); ++b_idx) {
                T* data_ptr = this->retrieve(all_alphas[a_idx], all_betas[b_idx]);
                if (!data_ptr) continue;

                int col = a_idx * all_betas.size() + b_idx;
                // data_ptr 存储的是 f2_{i, k_total}
                for (int i = 0; i < nb; ++i) {
                    for (int kt = 0; kt < total_k_max; ++kt) {
                        int row = i * total_k_max + kt;
                        matrix[row][col] = data_ptr[i * total_k_max + kt];
                    }
                }
            }
        }
        return matrix;
    }

    void clear() {
        // 遍历所有 table 中的 block，将其 data 向量重置为 0
        // 或者直接清空 table 映射
        table.clear(); 
    }
};

// 类型别名定义
template<typename T> using IndexStorage = UnifiedStorage<T,0>;      // 储存 Matrix (p, g)
template<typename T> using DualIndexStorage = UnifiedStorage<T,1>;  // 储存 std::vector/Matrix (f1, f2)

#endif