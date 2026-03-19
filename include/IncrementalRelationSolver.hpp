#ifndef INCREMENTAL_RELATION_SOLVER_HPP
#define INCREMENTAL_RELATION_SOLVER_HPP

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "RelationSolver.hpp"

namespace IncrementalRelation {

// ==========================================
// 变量标识：b[alpha, beta]
// ==========================================
struct VariableId {
    std::vector<int> alpha;
    std::vector<int> beta;
    
    bool operator==(const VariableId& other) const {
        return alpha == other.alpha && beta == other.beta;
    }
    
    bool operator<(const VariableId& other) const {
        if (alpha != other.alpha) return alpha < other.alpha;
        return beta < other.beta;
    }
    
    // 检查是否被另一个变量"覆盖"（用于消去）
    // this 被 cover 覆盖：alpha >= cover.alpha && beta >= cover.beta
    bool isCoveredBy(const VariableId& cover) const {
        if (alpha.size() != cover.alpha.size() || beta.size() != cover.beta.size())
            return false;
        for (size_t i = 0; i < alpha.size(); ++i)
            if (alpha[i] < cover.alpha[i]) return false;
        for (size_t i = 0; i < beta.size(); ++i)
            if (beta[i] < cover.beta[i]) return false;
        return true;
    }
    
    // 检查是否严格大于（用于新解判断）
    bool isStrictlyGreaterThan(const VariableId& other) const {
        return isCoveredBy(other) && (alpha != other.alpha || beta != other.beta);
    }
};

// ==========================================
// 单个 (lev, deg) 的解
// ==========================================
template<typename T>
struct LevelDegreeSolution {
    int lev;
    int deg;
    
    // 新增的变量列表
    std::vector<VariableId> newVariables;
    
    // 解：每个基向量是一个从 VariableId 到系数的映射
    // 实际上存储的是 b[alpha, beta] = sum(coeff * free_var)
    std::vector<std::map<VariableId, std::vector<T>>> basisSolutions;
    
    // 自由变量列表（解空间的基）
    std::vector<VariableId> freeVariables;
    
    // 统计
    int nullity = 0;
    bool converged = false;
    int stableOrder = -1;  // 达到稳定的阶数
};

// ==========================================
// 解管理器：维护所有 level/deg 的解，提供消去功能
// ==========================================
template<typename T>
class SolutionManager {
public:
    // 添加新的解
    void addSolution(int lev, int deg, const LevelDegreeSolution<T>& sol);
    
    // 检查变量是否已被消去（被之前的解覆盖）
    bool isEliminated(const VariableId& var, int currentLev, int currentDeg) const;
    
    // 获取当前 level/deg 下活跃的变量列表（未被消去的）
    std::vector<VariableId> getActiveVariables(
        const std::vector<VariableId>& allVars,
        int lev, int deg) const;
    
    // 获取累积的 bindep（所有已求解的变量）
    std::vector<VariableId> getAllSolvedVariables() const;
    
    // 获取特定 level 的 bindep
    std::vector<VariableId> getLevelSolvedVariables(int lev) const;
    
    // 打印状态
    void printStatus() const;

private:
    // 存储结构：table_[lev][deg] -> solution
    std::map<int, std::map<int, LevelDegreeSolution<T>>> table_;
    
    // 快速查找：所有已求解的变量（bindep）
    std::set<VariableId> allSolvedVars_;
    
    // 每 level 的已求解变量
    std::map<int, std::set<VariableId>> levelSolvedVars_;
};

// ==========================================
// 增量关系求解器（二维表版本）
// ==========================================
template<typename T>
class IncrementalRelationSolver {
public:
    struct Config {
        int maxLev = 3;           // 最大 lev
        int maxDeg = 3;           // 最大 deg
        int maxOrder = 5;         // 展开最大阶数
        AdaptiveSamplingConfig adaptiveConfig;
    };
    
    struct Result {
        // 二维表：result[lev][deg]
        std::map<int, std::map<int, LevelDegreeSolution<T>>> table;
        
        // 收敛信息
        int convergedLev = -1;
        int convergedDeg = -1;
        
        // 总解数
        int totalNullity = 0;
    };
    
    // 构造
    explicit IncrementalRelationSolver(const Config& config = {});
    
    // 初始化 regimes
    void initRegimes(const std::vector<RegimeData<T>>& regimes,
                     const std::vector<std::vector<int>>& nimaxLists,
                     int ne);
    
    // 求解指定 (lev, deg)
    LevelDegreeSolution<T> solveAt(int lev, int deg);
    
    // 完整求解（遍历所有 lev/deg）
    Result solveAll();
    
    // 逐步求解（从当前位置继续）
    LevelDegreeSolution<T> solveNext();  // 按 (1,0)->(1,1)->...->(1,maxDeg)->(2,0)->... 顺序
    
    // 获取当前位置
    std::pair<int, int> getCurrentPosition() const { return {currentLev_, currentDeg_}; }
    
    // 获取所有结果
    const Result& getResult() const { return result_; }

private:
    Config config_;
    SolutionManager<T> solutionManager_;
    Result result_;
    
    // regime 数据
    std::vector<RegimeData<T>> regimes_;
    std::vector<std::vector<int>> nimaxLists_;
    int ne_ = 0;
    
    // 当前位置
    int currentLev_ = 1;
    int currentDeg_ = 0;
    
    // 生成 (lev, deg) 对应的变量列表
    std::vector<VariableId> generateVariables(int lev, int deg);
    
    // 消去已求解变量
    std::vector<VariableId> removeSolvedVariables(
        const std::vector<VariableId>& vars, int lev, int deg);
    
    // 构建方程并求解
    LevelDegreeSolution<T> solveInternal(
        int lev, int deg, const std::vector<VariableId>& activeVars);
};

// ==========================================
// 实现
// ==========================================

template<typename T>
void SolutionManager<T>::addSolution(int lev, int deg, const LevelDegreeSolution<T>& sol) {
    table_[lev][deg] = sol;
    
    // 更新已求解变量集合
    for (const auto& basisSol : sol.basisSolutions) {
        for (const auto& [var, _] : basisSol) {
            allSolvedVars_.insert(var);
            levelSolvedVars_[lev].insert(var);
        }
    }
}

template<typename T>
bool SolutionManager<T>::isEliminated(const VariableId& var, int currentLev, int currentDeg) const {
    // 检查是否被同 level 低 degree 的解消去
    auto levelIt = levelSolvedVars_.find(currentLev);
    if (levelIt != levelSolvedVars_.end()) {
        for (const auto& solvedVar : levelIt->second) {
            // 同 level：alpha 必须相等，beta >= solved.beta
            if (var.alpha == solvedVar.alpha) {
                bool betaGe = true;
                for (size_t i = 0; i < var.beta.size(); ++i) {
                    if (var.beta[i] < solvedVar.beta[i]) {
                        betaGe = false;
                        break;
                    }
                }
                if (betaGe) return true;
            }
        }
    }
    
    // 检查是否被低 level 的解消去（alpha >= solved.alpha && beta >= solved.beta）
    for (const auto& [lev, vars] : levelSolvedVars_) {
        if (lev >= currentLev) continue;
        for (const auto& solvedVar : vars) {
            if (var.isCoveredBy(solvedVar)) return true;
        }
    }
    
    return false;
}

template<typename T>
std::vector<VariableId> SolutionManager<T>::getActiveVariables(
    const std::vector<VariableId>& allVars,
    int lev, int deg) const {
    std::vector<VariableId> active;
    for (const auto& var : allVars) {
        if (!isEliminated(var, lev, deg)) {
            active.push_back(var);
        }
    }
    return active;
}

template<typename T>
std::vector<VariableId> SolutionManager<T>::getAllSolvedVariables() const {
    return std::vector<VariableId>(allSolvedVars_.begin(), allSolvedVars_.end());
}

template<typename T>
std::vector<VariableId> SolutionManager<T>::getLevelSolvedVariables(int lev) const {
    auto it = levelSolvedVars_.find(lev);
    if (it != levelSolvedVars_.end()) {
        return std::vector<VariableId>(it->second.begin(), it->second.end());
    }
    return {};
}

template<typename T>
void SolutionManager<T>::printStatus() const {
    std::cout << "=== Solution Manager Status ===" << std::endl;
    std::cout << "Total solved variables: " << allSolvedVars_.size() << std::endl;
    for (const auto& [lev, vars] : levelSolvedVars_) {
        std::cout << "  Level " << lev << ": " << vars.size() << " variables" << std::endl;
    }
    std::cout << "================================" << std::endl;
}

// IncrementalRelationSolver 实现...

template<typename T>
IncrementalRelationSolver<T>::IncrementalRelationSolver(const Config& config)
    : config_(config) {
    currentLev_ = 1;
    currentDeg_ = 0;
}

template<typename T>
void IncrementalRelationSolver<T>::initRegimes(
    const std::vector<RegimeData<T>>& regimes,
    const std::vector<std::vector<int>>& nimaxLists,
    int ne) {
    regimes_ = regimes;
    nimaxLists_ = nimaxLists;
    ne_ = ne;
}

template<typename T>
std::vector<VariableId> IncrementalRelationSolver<T>::generateVariables(int lev, int deg) {
    std::vector<VariableId> vars;
    
    // 生成所有 |alpha| <= lev, |beta| <= deg 的组合
    std::vector<std::vector<int>> alphas, betas, temp;
    generateAllIndices(ne_, lev, temp, alphas, false);
    temp.clear();
    generateAllIndices(ne_, deg, temp, betas, false);
    
    for (const auto& alpha : alphas) {
        for (const auto& beta : betas) {
            vars.push_back({alpha, beta});
        }
    }
    
    return vars;
}

template<typename T>
std::vector<VariableId> IncrementalRelationSolver<T>::removeSolvedVariables(
    const std::vector<VariableId>& vars, int lev, int deg) {
    return solutionManager_.getActiveVariables(vars, lev, deg);
}

template<typename T>
LevelDegreeSolution<T> IncrementalRelationSolver<T>::solveAt(int lev, int deg) {
    // 1. 生成所有可能的变量
    auto allVars = generateVariables(lev, deg);
    
    // 2. 消去已求解的变量
    auto activeVars = removeSolvedVariables(allVars, lev, deg);
    
    std::cout << "Solving (" << lev << ", " << deg << "): "
              << allVars.size() << " vars, " << activeVars.size() << " active"
              << std::endl;
    
    if (activeVars.empty()) {
        LevelDegreeSolution<T> emptySol;
        emptySol.lev = lev;
        emptySol.deg = deg;
        emptySol.nullity = 0;
        emptySol.converged = true;
        return emptySol;
    }
    
    // 3. 构建方程并求解
    return solveInternal(lev, deg, activeVars);
}

template<typename T>
LevelDegreeSolution<T> IncrementalRelationSolver<T>::solveInternal(
    int lev, int deg, const std::vector<VariableId>& activeVars) {
    
    LevelDegreeSolution<T> sol;
    sol.lev = lev;
    sol.deg = deg;
    sol.newVariables = activeVars;
    
    // 使用 AdaptiveEquationBuilder 求解
    // TODO: 需要修改 AdaptiveEquationBuilder 支持 VariableId 变量列表
    // 这里暂时使用完整求解
    
    AdaptiveEquationBuilder<T> builder(config_.adaptiveConfig);
    auto buildResult = builder.build(regimes_, nimaxLists_, ne_);
    
    sol.nullity = buildResult.nullspace.nullity;
    sol.converged = buildResult.converged;
    
    // 存储解...（需要映射回 VariableId）
    
    // 添加到管理器
    solutionManager_.addSolution(lev, deg, sol);
    
    return sol;
}

template<typename T>
LevelDegreeSolution<T> IncrementalRelationSolver<T>::solveNext() {
    auto sol = solveAt(currentLev_, currentDeg_);
    result_.table[currentLev_][currentDeg_] = sol;
    
    // 更新位置
    if (currentDeg_ < config_.maxDeg) {
        currentDeg_++;
    } else {
        currentLev_++;
        currentDeg_ = 0;
    }
    
    return sol;
}

template<typename T>
typename IncrementalRelationSolver<T>::Result IncrementalRelationSolver<T>::solveAll() {
    for (int lev = 1; lev <= config_.maxLev; ++lev) {
        for (int deg = 0; deg <= config_.maxDeg; ++deg) {
            auto sol = solveAt(lev, deg);
            result_.table[lev][deg] = sol;
            result_.totalNullity += sol.nullity;
            
            // 可以添加收敛判断，如果 nullity 达到预期则提前退出
        }
    }
    return result_;
}

} // namespace IncrementalRelation

#endif // INCREMENTAL_RELATION_SOLVER_HPP
