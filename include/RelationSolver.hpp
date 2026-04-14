#ifndef RELATION_SOLVER_HPP
#define RELATION_SOLVER_HPP

#include <vector>
#include <map>
#include <random>
#include <memory>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <cstdint>

#include "Combinatorics.hpp"
#include "SeriesCoefficient.hpp"
#include "UnifiedStorage.hpp"
#include "LinearSolver.hpp"
#include "Utilities.hpp"
#include "RingDataLoader.hpp"

// 幂运算辅助函数（支持泛型）
namespace detail {
    template<typename T>
    inline T power(T base, int exp) {
        T result = T(1);
        for (int i = 0; i < exp; ++i) result *= base;
        return result;
    }
}

namespace RelationSolver {

// ==========================================
// a.随机数生成（已适配 double 和 firefly::FFInt）
// ==========================================
template<typename T>
std::vector<std::vector<T>> generateRandomNu(
    int ne, int num_nu,
    T min_val, T max_val,
    uint64_t modulus = 0)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    if constexpr (std::is_floating_point_v<T>) {
        std::uniform_real_distribution<T> dis(min_val, max_val);
        std::vector<std::vector<T>> result(num_nu, std::vector<T>(ne));
        for (int i = 0; i < num_nu; ++i)
            for (int j = 0; j < ne; ++j)
                result[i][j] = dis(gen);
        return result;
    }
    else if constexpr (std::is_same_v<T, firefly::FFInt>) {
        uint64_t low = 0;
        uint64_t high = (modulus > 0) ? (modulus - 1) : static_cast<uint64_t>(max_val.n);
        std::uniform_int_distribution<uint64_t> dis(low, high);
        std::vector<std::vector<T>> result(num_nu, std::vector<T>(ne));
        for (int i = 0; i < num_nu; ++i)
            for (int j = 0; j < ne; ++j)
                result[i][j] = firefly::FFInt(dis(gen));
        return result;
    }
    else {
        static_assert(sizeof(T) == 0, "Unsupported type for generateRandomNu");
    }
}

// ==========================================
// b.多重指标生成器（用于 alpha/beta 枚举）
// ==========================================
inline void generateAllIndices(int dim, int max_sum, std::vector<int>& current,
                               std::vector<std::vector<int>>& result, bool allow_neg) {
    if (current.size() == static_cast<size_t>(dim)) {
        result.push_back(current);
        return;
    }
    int start = allow_neg ? -max_sum : 0;
    for (int i = start; i <= max_sum; ++i) {
        current.push_back(i);
        int current_sum = 0;
        for (int x : current) current_sum += std::abs(x);
        if (current_sum <= max_sum)
            generateAllIndices(dim, max_sum, current, result, allow_neg);
        current.pop_back();
    }
}

// ==========================================
// 1. Regime
// 系统区域结构（扩展版，包含预计算能力）
// ==========================================
template<typename T>
struct RegimeData {
    // ===== 原始数据 =====
    const seriesCoefficient<T>* C;                     // 指向当前分支的系数容器
    std::vector<int> theta;                       // theta 向量
    std::vector<std::vector<T>> A_ops;            // A_i 矩阵组
    std::vector<std::vector<T>> A_inv_ops;        // A_i 逆矩阵组
    int nb;                                       // 矩阵块大小
    
    // ===== 预计算数据 =====
    std::unique_ptr<IndexStorage<T>> p_store;     // 预计算的 p(alpha) = A^{-alpha}
    bool is_prepared = false;                     // 预计算完成标记
    
    // ===== 预计算方法 =====
    // 计算所有 p(alpha)，与 nu 和 nimax 无关，只需执行一次
    void prepare(int lev, const std::vector<std::vector<int>>& alphas);
    
    // 获取 p(alpha) 指针 (matrix)
    T* getP(const std::vector<int>& alpha) const {
        return p_store ? p_store->retrieve(alpha) : nullptr;
    }
    
    // ===== 移动语义支持 =====
    // 默认构造函数（初始化为未准备状态）
    RegimeData() = default;
    
    // 显式删除复制构造和赋值（unique_ptr 不支持复制）
    RegimeData(const RegimeData&) = delete;
    RegimeData& operator=(const RegimeData&) = delete;
    
    // 默认支持移动构造和移动赋值
    RegimeData(RegimeData&&) noexcept = default;
    RegimeData& operator=(RegimeData&&) noexcept = default;
    
private:
    // 递归计算 p(alpha) 的辅助函数
    void computePRecursive(const std::vector<int>& alpha,
                          const std::vector<std::vector<int>>& all_alphas,
                          const std::vector<std::vector<T>>& A_ops,
                          const std::vector<std::vector<T>>& A_inv_ops);
};

// ==========================================
// 2. Evaluator
// 单点计算器：负责单个 regime、单个 nimax、单个 nu 的计算
// step1 (p_alpha) 已在 RegimeData 中预计算，此处只处理 step2-step5
// ==========================================
template<typename T>
class RegimeEvaluator {
public:
    // 默认构造
    RegimeEvaluator() = default;
    
    // 初始化：保存引用，创建临时存储器
    // 此方法必须在 evaluate() 前调用，否则 evaluate() 会抛出异常
    // 参数验证：alphas/betas 不能为空，reg 必须有效，k_max 必须非负
    void init(const RegimeData<T>& reg, 
              int k_max,
              int lev,
              int deg,
              const std::vector<std::vector<int>>& alphas,
              const std::vector<std::vector<int>>& betas);
    
    /**
     * 计算单个采样点在该 Regime 的方程贡献。
     * 
     * 此方法执行以下步骤：
     *   Step2: g(α) = p(α) * h(α, ν, nimax_idx)
     *   Step3: f1(α,β) = (θ + ν)^β 的多项式展开系数
     *   Step4: f2(α,β) = f1 ✕ g 的卷积
     *   Step5: 组装成最终矩阵
     * 
     * @param nu       是采样点向量，维度必须等于 ne
     * @param nimax_idx nimax 索引（对应查找财数 C(上操) 的是哪个展开阶数）
     * 
     * @return 方程矩阵 M，维度 [nb·(k_max+1)] × [|αs|·|βs|]
     *         索引 i*total_k + k 对应（第 i 个基量、第 k 阶）
     *         列索引 a_idx*|betas| + b_idx 对应（特定 (α, β) 組合）
     * 
     * @throw std::runtime_error 如果未初始化（isInitialized() 返回 false）
     * @throw std::invalid_argument 如果 nu 为空、维度不匹配或 alphas/betas 未初始化
     * 
     * @see init() 方法 （必须先调用以初始化评估器）
     * @see clearTemporaryStorage() 方法 （可清理结果）
     */
    std::vector<std::vector<T>> evaluate(const std::vector<T>& nu, int nimax_idx);
    
    // 清理临时存储（每次 evaluate 后或新 nu 前调用）
    void clearTemporaryStorage();
    
    // 获取该 regime 在单个 nu 下的贡献行数
    size_t rowsPerNu() const { return nb_ * (k_max_ + 1); }
    
    // 获取该 regime 支持的最大 nimax
    int maxNimax() const { return C_ ? C_->getNimax() : 0; }
    
    // 检查是否已初始化
    bool isInitialized() const { return C_ != nullptr && nb_ > 0; }

private:
    // step2: g(alpha) = p(alpha) * h(alpha, nu, nimax_idx)
    void step2_computeG(const std::vector<T>& nu, int nimax_idx);
    
    // step3: f1(beta) = (theta + nu)^beta 展开系数
    void step3_computeF1(const std::vector<T>& nu);
    
    // step4: f2 = f1 * g (多项式卷积)
    void step4_computeF2(const std::vector<int>& alpha, const std::vector<int>& beta);
    
    // step5: 组装最终矩阵
    std::vector<std::vector<T>> buildFinalMatrix();
    
    // 引用外部数据
    const RegimeData<T>* reg_;                    // 指向 regime（包含 p_store）
    const seriesCoefficient<T>* C_;               // 系数容器
    std::vector<int> theta_;                      // theta 向量
    
    // 维度参数
    int nb_ = 0;                                  // 矩阵块大小
    int ne_ = 0;                                  // 多重指标维度
    int k_max_ = 0;                               // 最大展开阶数
    int lev_ = 0;                                 // max |alpha|
    int deg_ = 0;                                 // max |beta|
    
    // alpha/beta 列表（全局共享）
    std::vector<std::vector<int>> alphas_;
    std::vector<std::vector<int>> betas_;
    
    // 临时存储器（仅 g, f1, f2，每次 evaluate 前需清理）
    std::unique_ptr<IndexStorage<T>> g_store_;           // 存储 g(alpha)
    std::unique_ptr<DualIndexStorage<T>> f1_store_;      // 存储 f1(alpha, beta)
    std::unique_ptr<DualIndexStorage<T>> f2_store_;      // 存储 f2(alpha, beta)
};

// ==========================================
// 3. Sampler
// 自适应采样器：智能生成采样点序列
// 策略：特殊点（单位向量、全1）优先，随后随机点
// ==========================================
struct AdaptiveSamplingConfig {
    int min_nu = 3;                    // 最小采样数
    int max_nu = 50;                   // 最大采样数
    int check_interval = 1;            // 秩检查间隔（每几个点检查一次）
    int nullity_stable_threshold = 3;  // 零空间维度连续稳定次数
    int verification_points = 3;       // 验证用额外点数
    double tolerance = 1e-10;          // 数值容差
    
    // 采样点生成配置（使用 double，生成时转换为目标类型）
    bool use_special_points = true;    // 是否使用特殊点（单位向量等）
    double random_min = 3.0;           // 随机点最小值（避免0）
    double random_max = 100.0;         // 随机点最大值
    
    // 用于生成 alphas/betas 的参数（必须提供）
    int lev_hint = 2;                  // max |alpha|
    int deg_hint = 2;                  // max |beta|
};

template<typename T>
class AdaptiveSampler {
public:
    struct State {
        int nu_count = 0;              // 已生成点数
        int current_nullity = -1;      // 当前零空间维度
        int stable_count = 0;          // 连续稳定次数
        bool converged = false;        // 是否收敛
    };
    
    // 构造：指定维度 ne 和配置
    AdaptiveSampler(int ne, const AdaptiveSamplingConfig& config);
    
    // 生成下一个采样点
    std::vector<T> next();
    
    // 查看下一个点（不推进）
    std::vector<T> peek() const;
    
    // 更新状态（在求解后调用）
    void update(int nullity);
    
    // 检查是否应该进行求解验证
    bool shouldCheck() const;
    
    // 检查是否已收敛
    bool hasConverged() const { return state_.converged; }
    
    // 获取当前状态
    State getState() const { return state_; }
    
    // 获取用于验证的额外点（不推进主序列）
    std::vector<std::vector<T>> getVerificationPoints(int count);
    
    // 重置状态
    void reset();

private:
    int ne_;                                    // 维度
    AdaptiveSamplingConfig config_;             // 配置
    State state_;                               // 当前状态
    std::vector<std::vector<T>> sequence_;      // 预生成的采样序列
    size_t position_ = 0;                       // 当前位置
    std::mt19937 rng_;                          // 随机数，防止每次重新创建
    
    // 生成采样序列
    void generateSequence();
    
    // 生成单个随机点
    std::vector<T> generateRandomPoint();
};

// ==========================================
// 4. Solver
// 增量零空间求解器：用于齐次方程 Ax = 0
// 支持增量添加方程行，维护零空间基
// ==========================================
template<typename T>
class IncrementalNullspaceSolver {
public:
    struct NullspaceInfo {
        std::vector<std::vector<T>> basis;  // 零空间基向量（行向量）
        int rank = 0;                        // 矩阵秩
        int nullity = 0;                     // 零空间维度 = num_vars - rank
        bool is_valid = false;               // 是否有效
    };
    
    // 构造：指定变量数（列数）
    explicit IncrementalNullspaceSolver(size_t num_vars);
    
    // 添加单方程行
    void addRow(const std::vector<T>& row);
    
    // 批量添加方程行
    void addRows(const std::vector<std::vector<T>>& rows);
    
    // 获取当前零空间信息（按需计算）
    NullspaceInfo getNullspace();
    
    // 快速获取零空间维度（不计算基向量）
    int getNullity();
    
    // 验证给定基向量在新方程行上的残差
    // 返回最大残差绝对值
    T verifyBasis(const std::vector<std::vector<T>>& test_rows, 
                 const NullspaceInfo& info,
                 T tolerance = T(1e-10));
    
    // 获取当前累积的矩阵（调试用）
    const std::vector<std::vector<T>>& getMatrix() const { return matrix_; }
    
    // 获取变量数
    size_t getNumVars() const { return num_vars_; }
    
    // 获取当前行数
    size_t getNumRows() const { return matrix_.size(); }
    
    // 清空所有数据
    void clear();

private:
    size_t num_vars_;                           // 变量数（列数）
    std::vector<std::vector<T>> matrix_;        // 累积的方程矩阵
    
    // 缓存
    NullspaceInfo cached_info_;
    bool cache_valid_ = false;
    
    // 内部计算零空间（使用高斯消元）
    NullspaceInfo computeNullspace();
    
    // 执行高斯消元，返回秩
    int gaussianElimination(std::vector<std::vector<T>>& A);
};

// ==========================================
// 5. Assembler
// 全局方程组装器：管理多个 regime 和所有 nimax
// 对每个 nu 点，组装所有 regime 和 nimax 的贡献
// ==========================================
template<typename T>
class GlobalEquationAssembler {
public:
    struct RegimeConfig {
        const RegimeData<T>* reg;           // 指向 regime 数据
        std::vector<int> nimax_indices;     // 该 regime 的所有 nimax 索引
    };
    
    // 默认构造
    GlobalEquationAssembler() = default;
    
    // 初始化：传入全局参数
    void init(int k_max, int lev, int deg,
              const std::vector<std::vector<int>>& alphas,
              const std::vector<std::vector<int>>& betas);
    
    // 添加 regime（及其 nimax 列表）
    void addRegime(const RegimeData<T>& reg, const std::vector<int>& nimax_indices);
    
    // 核心：在指定 nu 点计算所有贡献（所有 regime，所有 nimax）
    // 返回: 全局方程行（已拼接好，每行长度为 alphas.size() * betas.size()）
    std::vector<std::vector<T>> evaluateAtNu(const std::vector<T>& nu);
    
    // 获取单个 nu 贡献的总行数（所有 regime，所有 nimax）
    size_t totalRowsPerNu() const;
    
    // 获取变量数（列数）= alphas.size() * betas.size()
    size_t numVariables() const { return alphas_.size() * betas_.size(); }
    
    // 获取 regime 数量
    size_t numRegimes() const { return regimes_.size(); }

private:
    std::vector<RegimeConfig> regimes_;             // regime 配置列表
    std::vector<RegimeEvaluator<T>> evaluators_;    // 每个 regime 一个 evaluator
    
    // 全局参数
    int k_max_ = 0;
    int lev_ = 0;
    int deg_ = 0;
    std::vector<std::vector<int>> alphas_;
    std::vector<std::vector<int>> betas_;
};

// ==========================================
// 6. Builder
// 自适应方程构建器：整合采样、组装、求解的主控制器
// ==========================================
template<typename T>
class AdaptiveEquationBuilder {
public:
    struct BuildResult {
        // 求解结果
        std::vector<std::vector<T>> equations;          // 最终方程矩阵
        typename IncrementalNullspaceSolver<T>::NullspaceInfo nullspace;
        
        // 采样信息
        std::vector<std::vector<T>> nu_used;            // 实际使用的采样点
        int nu_count = 0;                               // 采样点数量
        
        // 收敛状态
        bool converged = false;
        std::string stop_reason;                        // 停止原因
    };
    
    // 构造函数：传入配置
    explicit AdaptiveEquationBuilder(const AdaptiveSamplingConfig& config = {});
    
    // 主构建流程
    // regimes: 所有 regime 数据（已 prepare）
    // nimax_lists: 每个 regime 对应的 nimax 索引列表
    // ne: 维度（用于采样器）
    BuildResult build(const std::vector<RegimeData<T>>& regimes,
                     const std::vector<std::vector<int>>& nimax_lists,
                     int ne);

private:
    AdaptiveSamplingConfig config_;
};

// ==========================================
// 独立工具函数：打印变量信息
// ==========================================
template<typename T>
void printVariables(const std::vector<std::vector<int>>& alphas,
                   const std::vector<std::vector<int>>& betas,
                   std::ostream& os = std::cout) {
    os << "=== Unknown Variables (alpha, beta) ===" << std::endl;
    size_t total = alphas.size() * betas.size();
    os << "Total number of variables: " << total << std::endl;

    for (size_t a_idx = 0; a_idx < alphas.size(); ++a_idx) {
        for (size_t b_idx = 0; b_idx < betas.size(); ++b_idx) {
            size_t col = a_idx * betas.size() + b_idx;
            os << "Col " << col << ":\talpha = [";
            for (size_t i = 0; i < alphas[a_idx].size(); ++i) {
                if (i > 0) os << ", ";
                os << alphas[a_idx][i];
            }
            os << "],\tbeta = [";
            for (size_t i = 0; i < betas[b_idx].size(); ++i) {
                if (i > 0) os << ", ";
                os << betas[b_idx][i];
            }
            os << "]" << std::endl;
        }
    }
    os << "========================================" << std::endl;
}

// ==========================================
// RegimeEvaluator 实现
// ==========================================

template<typename T>
void RegimeEvaluator<T>::init(const RegimeData<T>& reg, 
                              int k_max,
                              int lev,
                              int deg,
                              const std::vector<std::vector<int>>& alphas,
                              const std::vector<std::vector<int>>& betas)
{
    // ===== 参数验证 =====
    if (k_max < 0) {
        throw std::invalid_argument("init: k_max must be non-negative");
    }
    if (lev < 0) {
        throw std::invalid_argument("init: lev must be non-negative");
    }
    if (deg < 0) {
        throw std::invalid_argument("init: deg must be non-negative");
    }
    if (alphas.empty() || betas.empty()) {
        throw std::invalid_argument("init: alphas and betas must not be empty");
    }
    if (reg.C == nullptr) {
        throw std::invalid_argument("init: RegimeData.C (coefficient container) is null");
    }
    if (reg.nb <= 0) {
        throw std::invalid_argument("init: RegimeData.nb must be positive");
    }
    
    // ===== 初始化成员变量 =====
    reg_ = &reg;
    C_ = reg.C;
    theta_ = reg.theta;
    nb_ = reg.nb;
    ne_ = theta_.size();
    k_max_ = k_max;
    lev_ = lev;
    deg_ = deg;
    alphas_ = alphas;
    betas_ = betas;
    
    // ===== 创建临时存储器 =====
    g_store_ = std::make_unique<IndexStorage<T>>(ne_, nb_ * (k_max_ + 1));
    f1_store_ = std::make_unique<DualIndexStorage<T>>(ne_, k_max_ + 1);
    f2_store_ = std::make_unique<DualIndexStorage<T>>(ne_, nb_ * (k_max_ + 1));
}

template<typename T>
void RegimeEvaluator<T>::clearTemporaryStorage() {
    if (g_store_) g_store_->clear();
    if (f1_store_) f1_store_->clear();
    if (f2_store_) f2_store_->clear();
}

template<typename T>
std::vector<std::vector<T>> RegimeEvaluator<T>::evaluate(const std::vector<T>& nu, int nimax_idx) {
    // ===== 错误检查 =====
    if (!isInitialized()) {
        throw std::runtime_error("evaluate: RegimeEvaluator not initialized");
    }
    if (nu.empty()) {
        throw std::invalid_argument("evaluate: sampling point nu is empty");
    }
    if (static_cast<int>(nu.size()) != ne_) {
        throw std::invalid_argument("evaluate: nu size mismatch with dimension ne");
    }
    if (alphas_.empty() || betas_.empty()) {
        throw std::invalid_argument("evaluate: alphas or betas not initialized");
    }
    
    // 清理之前的临时数据
    clearTemporaryStorage();
    
    // ===== Step2-3: 计算所需的中间量 =====
    // Step2: 对所有 alpha 计算 g(alpha) = p(alpha) * h(alpha, nu, nimax_idx)
    // 注意：step2_computeG 内部遍历所有 alpha，故只调用一次
    step2_computeG(nu, nimax_idx);
    
    // Step3: 对所有 (alpha, beta) 组合计算 f1(alpha, beta) = (theta + nu)^beta 的多项式展开系数
    // 注意：step3_computeF1 内部遍历所有 (alpha, beta)，故只调用一次
    step3_computeF1(nu);
    
    // ===== Step4: 对所有 (alpha, beta) 组合计算 f2 =====
    // f2(alpha, beta) = f1 * g 的卷积
    for (const auto& alpha : alphas_) {
        for (const auto& beta : betas_) {
            step4_computeF2(alpha, beta);
        }
    }
    
    // ===== Step5: 组装最终矩阵 =====
    return buildFinalMatrix();
}

template<typename T>
void RegimeEvaluator<T>::step2_computeG(const std::vector<T>& nu, int nimax_idx) {
    // 防御性检查：验证存储器已初始化
    if (!g_store_ || !reg_ || !C_) {
        throw std::runtime_error("step2_computeG: prerequisite storage or data not initialized");
    }
    
    // 对所有 alpha 计算 g(alpha)
    for (const auto& alpha : alphas_) {
        // 检查是否已计算
        T* cached = g_store_->retrieve(alpha);
        if (cached != nullptr) continue;
        
        // 获取 p(alpha)（如果不存在则跳过）
        T* p_ptr = reg_->getP(alpha);
        if (!p_ptr) continue;
        
        std::vector<T> g_val(nb_ * (k_max_ + 1), T(0));
        
        // 计算 h_jk
        for (int k = 0; k <= k_max_; ++k) {
            std::vector<T> h_jk(nb_, T(0));
            
            for (int l = 0; l <= k; ++l) {
                long long states = BINOM[l + ne_ - 1][ne_ - 1];
                for (long long cid = 0; cid < states; ++cid) {
                    std::vector<int> gamma = readIndex(static_cast<int>(cid), l, ne_);
                    
                    // 计算 weight = (nu - alpha)^gamma
                    T weight = T(1);
                    for (int m = 0; m < ne_; ++m) {
                        T diff = nu[m] - static_cast<T>(alpha[m]);
                        weight *= detail::power(diff, gamma[m]);
                    }
                    
                    // 累加 h_jk
                    for (int j = 0; j < nb_; ++j) {
                        h_jk[j] += (*C_)(k, l, static_cast<int>(cid), j, nimax_idx) * weight;
                    }
                }
            }
            
            // 矩阵乘：g_i = sum_j p_ij * h_j
            for (int i = 0; i < nb_; ++i) {
                for (int j = 0; j < nb_; ++j) {
                    g_val[i * (k_max_ + 1) + k] += p_ptr[i * nb_ + j] * h_jk[j];
                }
            }
        }
        
        g_store_->insert(alpha, g_val);
    }
}

template<typename T>
void RegimeEvaluator<T>::step3_computeF1(const std::vector<T>& nu) {
    // 防御性检查：验证存储器已初始化
    if (!f1_store_) {
        throw std::runtime_error("step3_computeF1: f1_store not initialized");
    }
    
    // 预分配临时向量以优化多项式乘法（避免重复创建）
    // 主多项式和新多项式共享同一缓冲区
    std::vector<T> work_buffer(k_max_ + 1, T(0));
    
    // 对所有 (alpha, beta) 组合计算 f1
    for (const auto& alpha : alphas_) {
        for (const auto& beta : betas_) {
            // 检查缓存是否已有结果
            if (f1_store_->retrieve(alpha, beta) != nullptr) continue;
            
            std::vector<T> poly(k_max_ + 1, T(0));
            poly[0] = T(1);
            
            for (int i = 0; i < ne_; ++i) {
                int b = beta[i];
                if (b == 0) continue;
                
                std::vector<T> term_poly(b + 1, T(0));
                T th = static_cast<T>(theta_[i]);
                T n_val = nu[i];
                
                // 二项式展开系数
                for (int m = 0; m <= b; ++m) {
                    T coef = static_cast<T>(BINOM[b][m]) *
                             detail::power(th, b - m) *
                             detail::power(n_val, m);
                    term_poly[m] = coef;
                }
                
                // 优化：多项式乘法使用预分配的缓冲区，避免频繁动态分配
                // 清空工作缓冲区
                std::fill(work_buffer.begin(), work_buffer.end(), T(0));
                
                // 进行多项式乘法：work_buffer = poly * term_poly
                for (int p1 = 0; p1 <= k_max_; ++p1) {
                    if (poly[p1] == T(0)) continue;
                    for (int p2 = 0; p2 <= b; ++p2) {
                        int target_idx = p1 + p2;
                        if (target_idx <= k_max_) {
                            work_buffer[target_idx] += poly[p1] * term_poly[p2];
                        }
                    }
                }
                
                // 用计算结果替换 poly
                poly = work_buffer;
            }
            
            // 验证结果后插入存储器
            if (!poly.empty()) {
                f1_store_->insert(alpha, poly, beta);
            }
        }
    }
}

template<typename T>
void RegimeEvaluator<T>::step4_computeF2(const std::vector<int>& alpha, const std::vector<int>& beta) {
    // 防御性检查：验证必需的存储器存在
    if (!g_store_ || !f1_store_ || !f2_store_) {
        throw std::runtime_error("step4_computeF2: required storage not initialized");
    }
    
    // 获取指针并验证
    T* g_ptr = g_store_->retrieve(alpha);
    T* f1_ptr = f1_store_->retrieve(alpha, beta);
    
    // 如果任一数据缺失，说明前序步骤有问题，记录且返回
    if (!g_ptr || !f1_ptr) {
        // 可选：在调试模式下输出警告
        // std::cerr << "Warning: step4_computeF2 missing required data for alpha/beta\n";
        return;
    }
    
    int f2_len = k_max_ + 1;
    std::vector<T> f2_val(nb_ * f2_len, T(0));
    
    // 计算 f2 = f1 * g 的卷积
    for (int i = 0; i < nb_; ++i) {
        for (int l = 0; l <= k_max_; ++l) {
            T f1_val = f1_ptr[l];
            if (f1_val == T(0)) continue;
            
            for (int k = 0; k <= k_max_ - l; ++k) {
                T g_val = g_ptr[i * (k_max_ + 1) + k];
                if (g_val == T(0)) continue;
                
                int target_pow = l + k;
                if (target_pow <= k_max_) {  // 额外的边界检查
                    f2_val[i * f2_len + target_pow] += f1_val * g_val;
                }
            }
        }
    }
    
    // 验证结果非空后再插入
    f2_store_->insert(alpha, f2_val, beta);
}

template<typename T>
std::vector<std::vector<T>> RegimeEvaluator<T>::buildFinalMatrix() {
    // 防御性检查：验证必需的参数和存储器
    if (!f2_store_ || alphas_.empty() || betas_.empty() || nb_ <= 0 || k_max_ < 0) {
        throw std::runtime_error("buildFinalMatrix: invalid state (empty alphas/betas or invalid dimensions)");
    }
    
    int total_k = k_max_ + 1;
    int rows = nb_ * total_k;
    int cols = static_cast<int>(alphas_.size() * betas_.size());
    
    // 验证计算的矩阵维度是否合理
    if (rows <= 0 || cols <= 0) {
        throw std::runtime_error("buildFinalMatrix: invalid matrix dimensions (rows=" + 
                                 std::to_string(rows) + ", cols=" + std::to_string(cols) + ")");
    }
    
    std::vector<std::vector<T>> mat(rows, std::vector<T>(cols, T(0)));
    
    // 填充矩阵：每列对应一个 (alpha, beta) 对
    for (size_t a_idx = 0; a_idx < alphas_.size(); ++a_idx) {
        for (size_t b_idx = 0; b_idx < betas_.size(); ++b_idx) {
            // 获取该 (alpha, beta) 对应的 f2 数据
            T* f2_ptr = f2_store_->retrieve(alphas_[a_idx], betas_[b_idx]);
            if (!f2_ptr) continue;  // 如果数据不存在，跳过该列
            
            int col_idx = static_cast<int>(a_idx * betas_.size() + b_idx);
            
            // 将 f2 数据填充到矩阵对应列
            for (int i = 0; i < nb_; ++i) {
                for (int kt = 0; kt < total_k; ++kt) {
                    int row_idx = i * total_k + kt;
                    // 边界检查（防御性编程）
                    if (row_idx < rows && col_idx < cols) {
                        mat[row_idx][col_idx] = f2_ptr[i * total_k + kt];
                    }
                }
            }
        }
    }
    
    return mat;
}

// ==========================================
// AdaptiveSampler 实现
// ==========================================

template<typename T>
AdaptiveSampler<T>::AdaptiveSampler(int ne, const AdaptiveSamplingConfig& config)
    : ne_(ne), config_(config), position_(0), rng_(std::random_device{}())
{
    // 初始化RNG仅一次，每次generateRandomPoint()不再创建
    generateSequence();
}

template<typename T>
void AdaptiveSampler<T>::generateSequence() {
    sequence_.clear();
    sequence_.reserve(config_.max_nu);
    
    if (config_.use_special_points) {
        // 阶段1: 单位向量 e_i = (0,...,1,...,0)
        for (int i = 0; i < ne_ && (int)sequence_.size() < config_.max_nu; ++i) {
            std::vector<T> pt(ne_, T(0));
            pt[i] = T(1);
            sequence_.push_back(std::move(pt));
        }
        
        // 阶段2: 全1向量 (1,1,...,1)
        if ((int)sequence_.size() < config_.max_nu) {
            sequence_.push_back(std::vector<T>(ne_, T(1)));
        }
        
        // 阶段3: 非对称组合 (1,2,0,...), (2,1,0,...)
        if (ne_ >= 2) {
            for (int i = 0; i < ne_ && (int)sequence_.size() < config_.max_nu; ++i) {
                std::vector<T> pt(ne_, T(0));
                pt[i] = T(1);
                pt[(i+1) % ne_] = T(2);
                sequence_.push_back(std::move(pt));
            }
        }
        
        // 阶段4: 小整数网格 (2,0,0,...), (0,2,0,...)
        for (int i = 0; i < ne_ && (int)sequence_.size() < config_.max_nu; ++i) {
            std::vector<T> pt(ne_, T(0));
            pt[i] = T(2);
            sequence_.push_back(std::move(pt));
        }
    }
    
    // 阶段5: 随机点填充剩余
    int remaining = config_.max_nu - sequence_.size();
    for (int i = 0; i < remaining; ++i) {
        sequence_.push_back(generateRandomPoint());
    }
}

template<typename T>
std::vector<T> AdaptiveSampler<T>::generateRandomPoint() {
    // ✅ 改进1：使用成员变量 rng_ 而不是每次创建新的 random_device 和 mt19937
    // 这样避免了频繁的初始化开销，显著提升性能
    std::vector<T> pt(ne_);
    
    if constexpr (std::is_floating_point_v<T>) {
        std::uniform_real_distribution<double> dis(config_.random_min, config_.random_max);
        for (int i = 0; i < ne_; ++i) {
            pt[i] = static_cast<T>(dis(rng_));
        }
    } else if constexpr (std::is_same_v<T, firefly::FFInt>) {
        // 对于 FFInt，使用模数范围内的随机数
        uint64_t p = firefly::FFInt::p;
        uint64_t low = static_cast<uint64_t>(3);
        uint64_t high = p - 1;
        std::uniform_int_distribution<uint64_t> dis(low, high);
        for (int i = 0; i < ne_; ++i) {
            pt[i] = firefly::FFInt(dis(rng_));
        }
    } else {
        // 其他整数类型
        std::uniform_int_distribution<int64_t> dis(
            static_cast<int64_t>(config_.random_min),
            static_cast<int64_t>(config_.random_max)
        );
        for (int i = 0; i < ne_; ++i) {
            pt[i] = static_cast<T>(dis(rng_));
        }
    }
    
    return pt;
}

template<typename T>
std::vector<T> AdaptiveSampler<T>::next() {
    if (position_ >= sequence_.size()) {
        // 超出预生成范围，动态生成随机点
        return generateRandomPoint();
    }
    return sequence_[position_++];
}

template<typename T>
std::vector<T> AdaptiveSampler<T>::peek() const {
    if (position_ >= sequence_.size()) {
        // 无法预查看动态生成的点
        return {};
    }
    return sequence_[position_];
}

template<typename T>
void AdaptiveSampler<T>::update(int nullity) {
    state_.nu_count++;
    
    if (state_.current_nullity == -1) {
        // 第一次更新
        state_.current_nullity = nullity;
        state_.stable_count = 0;
        return;
    }
    
    if (nullity == state_.current_nullity) {
        state_.stable_count++;
    } else {
        state_.current_nullity = nullity;
        state_.stable_count = 0;
    }
    
    // 检查收敛条件
    if (state_.stable_count >= config_.nullity_stable_threshold && 
        state_.nu_count >= config_.min_nu) {
        state_.converged = true;
    }
}

template<typename T>
bool AdaptiveSampler<T>::shouldCheck() const {
    return state_.nu_count >= config_.min_nu && 
           (state_.nu_count % config_.check_interval == 0);
}

template<typename T>
std::vector<std::vector<T>> AdaptiveSampler<T>::getVerificationPoints(int count) {
    std::vector<std::vector<T>> points;
    points.reserve(count);
    for (int i = 0; i < count; ++i) {
        points.push_back(generateRandomPoint());
    }
    return points;
}

template<typename T>
void AdaptiveSampler<T>::reset() {
    state_ = State{};
    position_ = 0;
    // 可选：重新打乱随机部分或重新生成
}

// ==========================================
// IncrementalNullspaceSolver 实现
// ==========================================

template<typename T>
IncrementalNullspaceSolver<T>::IncrementalNullspaceSolver(size_t num_vars)
    : num_vars_(num_vars) {}

template<typename T>
void IncrementalNullspaceSolver<T>::addRow(const std::vector<T>& row) {
    if (row.size() != num_vars_) {
        throw std::invalid_argument("Row size does not match number of variables");
    }
    matrix_.push_back(row);
    cache_valid_ = false;
}

template<typename T>
void IncrementalNullspaceSolver<T>::addRows(const std::vector<std::vector<T>>& rows) {
    for (const auto& row : rows) {
        addRow(row);
    }
}

template<typename T>
int IncrementalNullspaceSolver<T>::getNullity() {
    if (!cache_valid_) {
        cached_info_ = computeNullspace();
        cache_valid_ = true;
    }
    return cached_info_.nullity;
}

template<typename T>
typename IncrementalNullspaceSolver<T>::NullspaceInfo 
IncrementalNullspaceSolver<T>::getNullspace() {
    if (!cache_valid_) {
        cached_info_ = computeNullspace();
        cache_valid_ = true;
    }
    return cached_info_;
}

template<typename T>
void IncrementalNullspaceSolver<T>::clear() {
    matrix_.clear();
    cached_info_ = NullspaceInfo{};
    cache_valid_ = false;
}

// 辅助函数：判断是否为0（特化版本）
namespace detail {
    template<typename T>
    inline bool isZero(T val) {
        if constexpr (std::is_floating_point_v<T>) {
            return std::abs(val) < 1e-14;
        } else {
            return val == T(0);
        }
    }
    
    template<typename T>
    inline T absValue(T val) {
        if constexpr (std::is_floating_point_v<T>) {
            return std::abs(val);
        } else {
            // 对于有限域，直接返回原值（我们只关心是否为零）
            return val;
        }
    }
}

template<typename T>
int IncrementalNullspaceSolver<T>::gaussianElimination(std::vector<std::vector<T>>& A) {
    if (A.empty() || A[0].empty()) return 0;
    
    int m = A.size();
    int n = A[0].size();
    int rank = 0;
    
    for (int col = 0, row = 0; col < n && row < m; ++col) {
        // 寻找主元
        int pivot = -1;
        for (int i = row; i < m; ++i) {
            if (!detail::isZero(A[i][col])) {
                pivot = i;
                break;
            }
        }
        
        if (pivot == -1) continue;
        
        // 交换行
        std::swap(A[row], A[pivot]);
        
        // 归一化主元行
        T piv_val = A[row][col];
        for (int j = col; j < n; ++j) {
            A[row][j] /= piv_val;
        }
        
        // 消去其他行
        for (int i = 0; i < m; ++i) {
            if (i != row && !detail::isZero(A[i][col])) {
                T factor = A[i][col];
                for (int j = col; j < n; ++j) {
                    A[i][j] -= factor * A[row][j];
                }
            }
        }
        
        rank++;
        row++;
    }
    
    return rank;
}

template<typename T>
typename IncrementalNullspaceSolver<T>::NullspaceInfo 
IncrementalNullspaceSolver<T>::computeNullspace() {
    NullspaceInfo info;
    info.is_valid = false;
    
    if (matrix_.empty()) {
        info.nullity = static_cast<int>(num_vars_);
        info.rank = 0;
        // 零空间是全体空间，基为单位矩阵
        info.basis.resize(num_vars_, std::vector<T>(num_vars_));
        for (size_t i = 0; i < num_vars_; ++i) {
            info.basis[i][i] = T(1);
        }
        info.is_valid = true;
        return info;
    }
    
    // 复制矩阵用于消元
    std::vector<std::vector<T>> A = matrix_;
    int m = A.size();
    int n = static_cast<int>(num_vars_);
    
    // 记录主元列位置
    std::vector<int> pivot_col;
    pivot_col.reserve(std::min(m, n));
    
    int rank = 0;
    for (int col = 0, row = 0; col < n && row < m; ++col) {
        // 寻找主元
        int pivot = -1;
        for (int i = row; i < m; ++i) {
            if (!detail::isZero(A[i][col])) {
                pivot = i;
                break;
            }
        }
        
        if (pivot == -1) continue;
        
        std::swap(A[row], A[pivot]);
        pivot_col.push_back(col);
        
        // 归一化
        T piv_val = A[row][col];
        for (int j = col; j < n; ++j) {
            A[row][j] /= piv_val;
        }
        
        // 消去
        for (int i = 0; i < m; ++i) {
            if (i != row && !detail::isZero(A[i][col])) {
                T factor = A[i][col];
                for (int j = col; j < n; ++j) {
                    A[i][j] -= factor * A[row][j];
                }
            }
        }
        
        rank++;
        row++;
    }
    
    // 构建零空间基
    info.rank = rank;
    info.nullity = n - rank;
    info.basis.clear();
    
    if (info.nullity > 0) {
        // 确定自由变量
        std::vector<bool> is_pivot_col(n, false);
        for (int pc : pivot_col) {
            is_pivot_col[pc] = true;
        }
        
        // 为每个自由变量构建一个基向量
        int free_idx = 0;
        for (int j = 0; j < n; ++j) {
            if (!is_pivot_col[j]) {
                std::vector<T> vec(n, T(0));
                vec[j] = T(1);  // 自由变量设为1
                
                // 主元变量由行最简形确定
                for (int r = 0; r < rank; ++r) {
                    int pc = pivot_col[r];
                    vec[pc] = -A[r][j];
                }
                
                info.basis.push_back(std::move(vec));
                free_idx++;
            }
        }
    }
    
    info.is_valid = true;
    return info;
}

template<typename T>
T IncrementalNullspaceSolver<T>::verifyBasis(
    const std::vector<std::vector<T>>& test_rows,
    const NullspaceInfo& info,
    T tolerance)
{
    if (!info.is_valid || info.basis.empty()) {
        return T(0);
    }
    
    T max_residual = T(0);
    
    for (const auto& row : test_rows) {
        if (row.size() != num_vars_) continue;
        
        // 计算该基向量集是否使行向量为零
        // 实际上我们需要验证：对每个基向量 b，row · b = 0
        for (const auto& basis_vec : info.basis) {
            T dot = T(0);
            for (size_t i = 0; i < num_vars_; ++i) {
                dot += row[i] * basis_vec[i];
            }
            // 使用类型安全的绝对值
            T abs_dot = detail::absValue(dot);
            // 对于浮点数，比较是否大于容差；对于FFInt，直接比较是否非零
            if constexpr (std::is_floating_point_v<T>) {
                if (abs_dot > max_residual) {
                    max_residual = abs_dot;
                }
            } else {
                if (abs_dot != T(0)) {
                    max_residual = T(1);  // FFInt非零即可标记
                }
            }
        }
    }
    
    return max_residual;
}

// ==========================================
// GlobalEquationAssembler 实现
// ==========================================

template<typename T>
void GlobalEquationAssembler<T>::init(int k_max, int lev, int deg,
                                       const std::vector<std::vector<int>>& alphas,
                                       const std::vector<std::vector<int>>& betas)
{
    k_max_ = k_max;
    lev_ = lev;
    deg_ = deg;
    alphas_ = alphas;
    betas_ = betas;
}

template<typename T>
void GlobalEquationAssembler<T>::addRegime(const RegimeData<T>& reg, 
                                            const std::vector<int>& nimax_indices)
{
    RegimeConfig config;
    config.reg = &reg;
    config.nimax_indices = nimax_indices;
    regimes_.push_back(config);
    
    // 创建并初始化 evaluator
    RegimeEvaluator<T> evaluator;
    evaluator.init(reg, k_max_, lev_, deg_, alphas_, betas_);
    evaluators_.push_back(std::move(evaluator));
}

template<typename T>
std::vector<std::vector<T>> GlobalEquationAssembler<T>::evaluateAtNu(const std::vector<T>& nu)
{
    std::vector<std::vector<T>> all_rows;
    
    // 遍历所有 regime
    for (size_t r = 0; r < regimes_.size(); ++r) {
        // 遍历该 regime 的所有 nimax
        for (int nimax_idx : regimes_[r].nimax_indices) {
            // 计算该 (regime, nimax) 在当前 nu 下的贡献
            auto rows = evaluators_[r].evaluate(nu, nimax_idx);
            
            // 累加到全局行集合
            all_rows.insert(all_rows.end(), rows.begin(), rows.end());
        }
    }
    
    return all_rows;
}

template<typename T>
size_t GlobalEquationAssembler<T>::totalRowsPerNu() const
{
    size_t total = 0;
    for (size_t r = 0; r < regimes_.size(); ++r) {
        // 每个 nimax 贡献 evaluator.rowsPerNu() 行
        total += evaluators_[r].rowsPerNu() * regimes_[r].nimax_indices.size();
    }
    return total;
}

// ==========================================
// AdaptiveEquationBuilder 实现
// ==========================================

template<typename T>
AdaptiveEquationBuilder<T>::AdaptiveEquationBuilder(const AdaptiveSamplingConfig& config)
    : config_(config) {}

template<typename T>
typename AdaptiveEquationBuilder<T>::BuildResult 
AdaptiveEquationBuilder<T>::build(
    const std::vector<RegimeData<T>>& regimes,
    const std::vector<std::vector<int>>& nimax_lists,
    int ne)
{
    BuildResult result;
    
    // ===== 改进3：配置参数验证 =====
    // 检查采样配置的合理性
    if (config_.min_nu <= 0 || config_.max_nu <= 0) {
        result.stop_reason = "Invalid config: min_nu and max_nu must be positive";
        return result;
    }
    if (config_.min_nu > config_.max_nu) {
        result.stop_reason = "Invalid config: min_nu must be <= max_nu";
        return result;
    }
    if (config_.nullity_stable_threshold <= 0) {
        result.stop_reason = "Invalid config: nullity_stable_threshold must be positive";
        return result;
    }
    if (config_.check_interval <= 0) {
        result.stop_reason = "Invalid config: check_interval must be positive";
        return result;
    }
    
    // 参数检查
    if (regimes.empty() || nimax_lists.empty() || regimes.size() != nimax_lists.size()) {
        result.stop_reason = "Invalid input: empty regimes or mismatched nimax_lists";
        return result;
    }
    
    // 从第一个 regime 获取 alphas/betas（所有 regime 共享）
    // 注意：这里假设所有 regime 使用相同的 alphas/betas
    // 实际应该从外部传入或从配置获取
    int k_max = regimes[0].C->getKmax();
    
    // 生成 alphas, betas
    std::vector<std::vector<int>> alphas, betas;
    std::vector<int> temp;
    generateAllIndices(ne, config_.lev_hint, temp, alphas, false);
    temp.clear();
    generateAllIndices(ne, config_.deg_hint, temp, betas, false);
    
    size_t num_vars = alphas.size() * betas.size();
    
    // 初始化组件
    AdaptiveSampler<T> sampler(ne, config_);
    GlobalEquationAssembler<T> assembler;
    IncrementalNullspaceSolver<T> solver(num_vars);
    
    // 初始化 assembler
    assembler.init(k_max, config_.lev_hint, config_.deg_hint, alphas, betas);
    for (size_t r = 0; r < regimes.size(); ++r) {
        assembler.addRegime(regimes[r], nimax_lists[r]);
    }
    
    // ✅ 改进2：内存预分配
    // 计算单次采样的总行数（所有 regime 和 nimax 的贡献）
    size_t rows_per_nu = assembler.totalRowsPerNu();
    // 预分配方程矩阵的存储空间，避免多次重新分配
    // 理论最大行数为 rows_per_nu * config_.max_nu
    result.equations.reserve(rows_per_nu * config_.max_nu);
    result.nu_used.reserve(config_.max_nu);
    
    // 主循环：自适应采样
    while (result.nu_count < config_.max_nu) {
        // 获取下一个采样点
        std::vector<T> nu = sampler.next();
        
        // 在该点评估所有 regime 和 nimax
        auto rows = assembler.evaluateAtNu(nu);
        
        // 添加到求解器
        solver.addRows(rows);
        
        // 记录使用的采样点
        result.nu_used.push_back(nu);
        result.nu_count++;
        result.equations.insert(result.equations.end(), rows.begin(), rows.end());
        
        // 检查是否需要进行收敛检测
        if (sampler.shouldCheck()) {
            int current_nullity = solver.getNullity();
            sampler.update(current_nullity);
            
            // 检查是否收敛
            if (sampler.hasConverged()) {
                // 额外验证
                auto verify_points = sampler.getVerificationPoints(config_.verification_points);
                bool verification_passed = true;
                
                auto current_info = solver.getNullspace();
                for (const auto& vnu : verify_points) {
                    auto vrows = assembler.evaluateAtNu(vnu);
                    T residual = solver.verifyBasis(vrows, current_info, T(config_.tolerance));
                    
                    if constexpr (std::is_floating_point_v<T>) {
                        if (residual > T(config_.tolerance)) {
                            verification_passed = false;
                            break;
                        }
                    } else {
                        if (residual != T(0)) {
                            verification_passed = false;
                            break;
                        }
                    }
                }
                
                if (verification_passed) {
                    result.converged = true;
                    result.nullspace = current_info;
                    result.stop_reason = "Converged after " + std::to_string(result.nu_count) + " samples";
                    return result;
                }
                // 验证失败，继续采样
            }
        }
    }
    
    // 达到最大采样数
    result.nullspace = solver.getNullspace();
    result.stop_reason = "Max samples (" + std::to_string(config_.max_nu) + ") reached";
    return result;
}

// ==========================================
// 辅助函数：构建所有区域
// ==========================================
/**
 * 从原始矩阵数据直接构建所有计算区域。
 * 
 * 这是 低级API。大多数情况下，较为优先
 * 选择使用第二个 buildAllRegimes() 重载，他会自动爱抽
 * 这些情报。
 * 
 * @tparam T 数据类型 (通常是 double 或 firefly::FFInt)
 * @param CTable 层递归生成的系数容器列表
 * @param sector 每个区域的 θ 右探上限向量
 * @param A_list 每个区域的 IBP 算数子列表
 * @param Ainv_list 每个区域的 IBP 算数子逆矩阵列表
 * @param ne 数据维度
 * @param nu_per_regime 每个区域的采样点数（对这个函数未昨）
 * @param min_val 随机点最小值（默认 0）
 * @param max_val 随机点最大值（默认 1）
 * @param modulus 有限域模数（仅对 FFInt 有效默认 0）
 * @return 所有区域的 RegimeData 列表
 * 
 * @see buildAllRegimes(vector, vector<RingCell>) 推荐用此重载
 */
template<typename T>
std::vector<RegimeData<T>> buildAllRegimes(
    std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne, int nu_per_regime,
    T min_val = T(0), T max_val = T(1), uint64_t modulus = 0)
{
    std::vector<RegimeData<T>> all_regimes;

    for (size_t r = 0; r < CTable.size(); ++r) {
        for (size_t s = 0; s < CTable[r].size(); ++s) {
            RegimeData<T> reg;
            reg.C = &CTable[r][s];
            reg.theta = sector[r];
            reg.A_ops = A_list[r];
            reg.A_inv_ops = Ainv_list[r];
            reg.nb = reg.C->basis_size();
            all_regimes.push_back(reg);
        }
    }
    return all_regimes;
}

/**
 * 从 RingDataLoader 加载的 RingCell 数据直接构建所有计算区域。
 * @tparam T 数据类型
 * @param CTable 层递归生成的系数容器列表
 * @param ringData RingDataLoader 加载的数据（std::vector<RingCell<T>>）
 * @param ne 多重指标维度（应与 ringData 中的 limitSector 长度一致）
 * @param nu_per_regime 每个区域的采样点数
 * @param min_val 随机采样最小值（默认 0）
 * @param max_val 随机采样最大值（默认 1）
 * @param modulus 有限域模数（仅对 FFInt 有效）
 * @return 所有区域的 RegimeData 列表
 */
template<typename T>
std::vector<RegimeData<T>> buildAllRegimes(
    std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<AlgebraData::RingCell<T>>& ringData,
    int ne, int nu_per_regime,
    T min_val = T(0), T max_val = T(1), uint64_t modulus = 0)
{
    int nreg = static_cast<int>(ringData.size());
    std::vector<std::vector<int>> sector(nreg);
    std::vector<std::vector<std::vector<T>>> A_list(nreg);
    std::vector<std::vector<std::vector<T>>> Ainv_list(nreg);

    for (int i = 0; i < nreg; ++i) {
        sector[i] = ringData[i].limitSector;
        A_list[i] = ringData[i].A_list;
        Ainv_list[i] = ringData[i].Ainv_list;
    }

    // 调用原有实现
    return buildAllRegimes(CTable, sector, A_list, Ainv_list, ne, nu_per_regime,
                           min_val, max_val, modulus);
}

// ==========================================
// RelationCoefficient 类
// 存储每个 (alpha, beta) 对应的系数向量（基础解系）
// ==========================================
template<typename T>
class RelationCoefficient {
private:
    std::vector<std::vector<int>> alphas_;
    std::vector<std::vector<int>> betas_;
    std::vector<std::vector<T>> coeffs_; // 每个元素是一个向量，长度 = 1 + nsol
    // 或者用map快速查找，但这里用顺序存储，通过二分查找
    mutable std::map<std::pair<std::vector<int>, std::vector<int>>, std::vector<T>*> cache_;

    void buildCache() {
        cache_.clear();
        for (size_t a = 0; a < alphas_.size(); ++a) {
            for (size_t b = 0; b < betas_.size(); ++b) {
                size_t idx = a * betas_.size() + b;
                cache_[{alphas_[a], betas_[b]}] = &coeffs_[idx];
            }
        }
    }

public:
    RelationCoefficient() = default;

    // 从alphas, betas和线性求解结果构造
    RelationCoefficient(const std::vector<std::vector<int>>& alphas,
                        const std::vector<std::vector<int>>& betas,
                        const LinearSystemResult<T>& res)
        : alphas_(alphas), betas_(betas)
    {
        // res.Mext 是 cols x (1+nsol) 矩阵，其中 cols = alphas.size() * betas.size()
        size_t cols = alphas.size() * betas.size();
        if (res.Mext.size() != cols) {
            throw std::runtime_error("Mext size mismatch with alpha/beta count");
        }
        // 转置存储：coeffs_ 为 cols 行，每行长度为 1+nsol
        coeffs_.resize(cols);
        for (size_t i = 0; i < cols; ++i) {
            coeffs_[i] = res.Mext[i]; // 直接复制
        }
        buildCache();
    }

    // 访问指定 (alpha, beta) 的系数向量（特解+基础解系）
    const std::vector<T>& operator()(const std::vector<int>& alpha, const std::vector<int>& beta) const {
        auto it = cache_.find({alpha, beta});
        if (it == cache_.end()) {
            throw std::out_of_range("(alpha, beta) pair not found");
        }
        return *it->second;
    }

    // 获取所有 alpha 列表
    const std::vector<std::vector<int>>& getAlphas() const { return alphas_; }
    // 获取所有 beta 列表
    const std::vector<std::vector<int>>& getBetas() const { return betas_; }

    // 获取解空间的维数（不包括特解，即基础解系个数）
    size_t getNumSolutions() const {
        if (coeffs_.empty()) return 0;
        return coeffs_[0].size() - 1; // 第一列是特解
    }

    // 获取特解部分（每个变量的常数项）
    std::vector<T> getParticular() const {
        std::vector<T> part(coeffs_.size());
        for (size_t i = 0; i < coeffs_.size(); ++i) {
            part[i] = coeffs_[i][0];
        }
        return part;
    }

    // 获取第 k 个基础解（k 从 0 开始）
    std::vector<T> getBasisSolution(size_t k) const {
        if (k >= getNumSolutions()) throw std::out_of_range("Basis index out of range");
        std::vector<T> basis(coeffs_.size());
        for (size_t i = 0; i < coeffs_.size(); ++i) {
            basis[i] = coeffs_[i][k + 1];
        }
        return basis;
    }
};

// ==========================================
// 新重构函数：使用自适应构建器（推荐）
// ==========================================

/**
 * 重构约化关系（新版本，使用自适应采样）
 * 
 * @param CTable 级数系数表
 * @param sector sector 列表
 * @param A_list A 矩阵列表
 * @param Ainv_list A 逆矩阵列表
 * @param ne 维度
 * @param lev max |alpha|
 * @param deg max |beta|
 * @param config 自适应采样配置（可选）
 * @return 线性系统结果和关系系数
 */
template<typename T>
static std::pair<LinearSystemResult<T>, RelationCoefficient<T>> reconstructReductionRelation(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne, int lev, int deg,
    const AdaptiveSamplingConfig& config = {})
{
    // 构建所有 regimes（不预生成采样点）
    std::vector<RegimeData<T>> regimes;
    for (size_t r = 0; r < CTable.size(); ++r) {
        for (size_t s = 0; s < CTable[r].size(); ++s) {
            RegimeData<T> reg;
            reg.C = &CTable[r][s];
            reg.theta = sector[r];
            reg.A_ops = A_list[r];
            reg.A_inv_ops = Ainv_list[r];
            reg.nb = reg.C->basis_size();
            regimes.push_back(std::move(reg));  // 使用移动语义避免复制
        }
    }
    
    // 生成 alphas/betas
    std::vector<std::vector<int>> alphas, betas;
    std::vector<int> temp;
    generateAllIndices(ne, lev, temp, alphas, false);
    temp.clear();
    generateAllIndices(ne, deg, temp, betas, false);
    
    if (regimes.empty()) {
        LinearSystemResult<T> empty_result{false, {}, {}};
        RelationCoefficient<T> empty_coeff;
        return {empty_result, empty_coeff};
    }
    
    // prepare 所有 regimes
    for (auto& reg : regimes) {
        reg.prepare(lev, alphas);
    }
    
    // 构建 nimax 列表（从 CTable 结构推断）
    std::vector<std::vector<int>> nimax_lists;
    for (const auto& reg : regimes) {
        std::vector<int> nimax_list;
        int nimax = reg.C->getNimax();
        for (int i = 0; i <= nimax; ++i) {
            nimax_list.push_back(i);
        }
        nimax_lists.push_back(nimax_list);
    }
    
    // 使用自适应构建器
    AdaptiveSamplingConfig cfg = config;
    cfg.lev_hint = lev;
    cfg.deg_hint = deg;
    
    AdaptiveEquationBuilder<T> builder(cfg);
    auto result = builder.build(regimes, nimax_lists, ne);
    
    // 将 BuildResult 转换为 LinearSystemResult
    LinearSystemResult<T> lsr_result;
    lsr_result.hasSolution = result.converged || result.nullspace.is_valid;
    
    // 构建 Mext 格式 [cols] x [1 + nullity]，第一列为特解（齐次问题为0），后续为零空间基
    size_t num_vars = alphas.size() * betas.size();
    int nullity = result.nullspace.nullity;
    
    lsr_result.Mext.resize(num_vars);
    for (size_t i = 0; i < num_vars; ++i) {
        lsr_result.Mext[i].resize(1 + nullity, T(0));
        for (int j = 0; j < nullity; ++j) {
            lsr_result.Mext[i][j + 1] = result.nullspace.basis[j][i];
        }
    }
    
    // 自由变量索引：全部（齐次方程）
    lsr_result.S.resize(num_vars);
    std::iota(lsr_result.S.begin(), lsr_result.S.end(), 0);
    
    // 构造 RelationCoefficient
    RelationCoefficient<T> coeff(alphas, betas, lsr_result);
    
    return {lsr_result, coeff};
}

// ==========================================
// RegimeData 实现
// ==========================================

template<typename T>
void RegimeData<T>::prepare(int lev, const std::vector<std::vector<int>>& alphas) {
    if (is_prepared) return;
    
    int ne = theta.size();
    p_store = std::make_unique<IndexStorage<T>>(ne, nb * nb);
    
    // 递归计算所有 p(alpha)
    for (const auto& alpha : alphas) {
        computePRecursive(alpha, alphas, A_ops, A_inv_ops);
    }
    
    is_prepared = true;
}

template<typename T>
void RegimeData<T>::computePRecursive(
    const std::vector<int>& alpha,
    const std::vector<std::vector<int>>& all_alphas,
    const std::vector<std::vector<T>>& A_ops,
    const std::vector<std::vector<T>>& A_inv_ops)
{
    // 检查是否已计算
    if (p_store->retrieve(alpha) != nullptr) return;
    
    // 基础情况：alpha = 0，单位矩阵
    bool is_zero = true;
    for (int x : alpha) if (x != 0) { is_zero = false; break; }
    if (is_zero) {
        std::vector<T> identity(nb * nb, T(0));
        for (int i = 0; i < nb; ++i) identity[i * nb + i] = T(1);
        p_store->insert(alpha, identity);
        return;
    }
    
    // 找到第一个非零维度
    int idx = -1;
    for (int i = 0; i < (int)alpha.size(); ++i) {
        if (alpha[i] != 0) { idx = i; break; }
    }
    if (idx == -1) return;  // 不应该发生
    
    std::vector<int> prev_alpha = alpha;
    const std::vector<T>* op_ptr = nullptr;
    
    if (alpha[idx] > 0) {
        prev_alpha[idx]--;
        op_ptr = &A_inv_ops[idx];
    } else {
        prev_alpha[idx]++;
        op_ptr = &A_ops[idx];
    }
    
    // 递归计算前一项
    computePRecursive(prev_alpha, all_alphas, A_ops, A_inv_ops);
    T* prev_data = p_store->retrieve(prev_alpha);
    if (!prev_data) return;
    
    // 矩阵乘法：res = Op * prev_data
    std::vector<T> res(nb * nb, T(0));
    const std::vector<T>& op = *op_ptr;
    for (int r = 0; r < nb; ++r) {
        for (int c = 0; c < nb; ++c) {
            T sum = T(0);
            for (int k = 0; k < nb; ++k) {
                sum += op[r * nb + k] * prev_data[k * nb + c];
            }
            res[r * nb + c] = sum;
        }
    }
    p_store->insert(alpha, res);
}

} // namespace RelationSolver

#endif // RELATION_SOLVER_HPP