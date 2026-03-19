#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>  // 新增
#include <map>      // 用于存储 (lev, deg) -> 解

// 新模块头文件
#include "IBPMatrixLoader_Binary.hpp"   // 二进制 IBP 矩阵加载器
#include "LayerRecursion.hpp"            // 层递归入口
#include "RingDataLoader.hpp"            // 环数据加载器
#include "RelationSolver.hpp"            // 关系求解器
#include "Combinatorics.hpp"             // 组合数工具
#include "LinearSolver.hpp"              // 线性求解器
#include "firefly/FFInt.hpp"              // 有限域类型
#include "SeriesCoefficientIO.hpp"


using namespace std;
using namespace firefly;

// ==========================================
// 辅助函数：动态计算采样配置
// ==========================================

/**
 * 根据 (lev, deg) 和系统参数，动态调整採样点数
 * 
 * 原理：
 * - 每个採样点生成 num_regimes × nb × (k_max+1) 个方程
 * - 变量数 = |alphas(lev)| × |betas(deg)|
 * - 确保方程数 ≥ 变量数（通常使用 1.5-2.0 倍的安全系数）
 */
RelationSolver::AdaptiveSamplingConfig 
adjustSamplingConfig(
    const RelationSolver::AdaptiveSamplingConfig& base_config,
    int lev, int deg, int ne, int nb, int k_max, int num_regimes)
{
    auto config = base_config;
    config.lev_hint = lev;
    config.deg_hint = deg;
    
    // 计算该 (lev, deg) 级别的变量总数
    std::vector<int> temp;
    std::vector<std::vector<int>> alphas, betas;
    RelationSolver::generateAllIndices(ne, lev, temp, alphas, false);
    temp.clear();
    RelationSolver::generateAllIndices(ne, deg, temp, betas, false);
    
    int num_variables = alphas.size() * betas.size();
    
    // 计算每个採样点产生的方程数
    int equations_per_sample = num_regimes * nb * (k_max + 1);
    
    // 计算所需的採样点数（包括 1.2x 的安全系数，确保超定系统但不过度）
    static constexpr double SAFETY_FACTOR = 1.2;  // 减小系数从 1.5 到 1.2，平衡准确性和速度
    int min_samples_required = static_cast<int>(
        std::ceil(num_variables * SAFETY_FACTOR / equations_per_sample)
    );
    
    // 确定最终採样数：至少 min_samples_required，但不超过硬件限制
    static constexpr int MAX_SAMPLES_LIMIT = 150;  // 硬件上限：最多 150 个採样点（加快演示）
    config.max_nu = std::min(
        MAX_SAMPLES_LIMIT,
        std::max(
            base_config.max_nu,                    // 基础配置的下界
            min_samples_required                   // 根据变量数计算的下界
        )
    );
    
    // 最小採样数也应该相应增加
    config.min_nu = std::min(
        config.min_nu,
        std::max(3, min_samples_required / 3)  // 至少为 3，但可随变量数调整
    );
    
    // 输出调整信息
    std::cout << "    Adaptive sampling config:" << std::endl;
    std::cout << "      Variables: " << num_variables << std::endl;
    std::cout << "      Equations per sample: " << equations_per_sample << std::endl;
    std::cout << "      Min samples required (" << std::fixed << std::setprecision(1) << (SAFETY_FACTOR*100) << "% safety): " << min_samples_required << std::endl;
    std::cout << "      Adjusted max_nu: " << config.max_nu << std::endl;
    
    return config;
}

// ==========================================
// 结构体：存储指定 (lev, deg) 级别的解
// ==========================================
struct SolutionAtLevel {
    int lev;
    int deg;
    int num_variables;                    // 变量总数 = |alphas| * |betas|
    int solution_dimension;               // 解空间维数
    double computation_time;              // 计算耗时（秒）
    int system_rows;                      // 方程行数
    int system_cols;                      // 方程列数
    int sampling_points_used;             // 实际採样点数
    
    SolutionAtLevel() : sampling_points_used(0) {}
    
    // formatted 输出
    void print(std::ostream& os = std::cout) const {
        double ratio = system_cols > 0 ? static_cast<double>(system_rows) / system_cols : 0.0;
        std::string status = ratio >= 1.0 ? "✓" : "✗";
        os << "  " << status << " (lev=" << lev << ", deg=" << deg << "): "
           << "vars=" << num_variables 
           << " sol_dim=" << solution_dimension
           << " sys=" << system_rows << "x" << system_cols
           << " ratio=" << std::fixed << std::setprecision(2) << ratio
           << " nu=" << sampling_points_used
           << " time=" << std::fixed << std::setprecision(2) << computation_time << "s" << std::endl;
    }
};

int main() {
    // 1. 初始化组合数表
    initBinomial();

    // 2. 参数设置
    const int order = 4;          // 展开阶数
    const int incre = 2;          // 层数增量
    
    // ========== 迭代参数配置 ==========
    // 定义要遍历的 (lev, deg) 范围
    const int lev_min = 1, lev_max = 2;    // lev: 从 1 到 2（缩小范围以加快演示）
    const int deg_min = 0, deg_max = 2;    // deg: 从 0 到 2（先增加 deg，再增加 lev）
    
    const int nu_per_regime = 100; // 每个区间的采样点数

    // 文件路径（假设使用二进制格式）
    const string family = "bub"; // "DBtop"
    const string binIBPFile = "IBPMat_"+family+".bin";   // IBP矩阵二进制文件
    const string binRingFile = "RingData_"+family+".bin"; // 环数据二进制文件
    const string coeffCacheFile = "ExpansionCache"+family+".bin";  // 缓存文件

    try {
        cout << "=== Test 2: Reconstruction of Linear Relations (Finite Field) ===" << endl;

        // ---- 第1a步：加载 IBP 矩阵 ----
        cout << "Loading IBP matrices from " << binIBPFile << " ..." << endl;
        auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(binIBPFile);
        if (ibpmatlist.empty()) {
            cerr << "Error: No IBP matrices loaded." << endl;
            return 1;
        }
        int ne = ibpmatlist[0].ne;
        int nb = ibpmatlist[0].nb;
        cout << "First matrix: ne = " << ne << ", nb = " << nb << endl;
        cout << "Prime modulus: " << FFInt::p << endl;

        // ---- 第1b步：计算/读取展开系数 ----
        vector<vector<seriesCoefficient<FFInt>>> allResults;
        try {
            allResults = SeriesIO::loadAllResults<FFInt>(coeffCacheFile);
            cout << "Loaded allResults from cache file." << endl;
        } catch (const std::exception& e) {
            cout << "Cache not found or invalid, running recursion..." << endl;
            allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
            SeriesIO::saveAllResults(allResults, coeffCacheFile);
            cout << "Saved allResults to cache." << endl;
        }
        cout << "Total solution branches: " << allResults.size() << " matrices, each with "
             << allResults[0].size() << " solutions (first matrix)." << endl;

        // ---- 第2步：加载环数据 ----
        cout << "\nLoading ring data from " << binRingFile << " ..." << endl;
        // 环数据文件应存储整数（模数下的值），传入模数以确保 FFInt 正确初始化
        auto ringData = AlgebraData::LoadBinary<FFInt>(binRingFile, FFInt::p);
        if (ringData.empty()) {
            cerr << "Error: No ring data loaded." << endl;
            return 1;
        }
        cout << "Loaded " << ringData.size() << " regimes." << endl;

        // ---- 数据提取：提前准备，避免重复处理 ----
        std::vector<std::vector<int>> sector_list;
        std::vector<std::vector<std::vector<FFInt>>> A_list, Ainv_list;
        for (const auto& ring : ringData) {
            sector_list.push_back(ring.limitSector);
            A_list.push_back(ring.A_list);
            Ainv_list.push_back(ring.Ainv_list);
        }
        cout << "Extracted " << sector_list.size() << " sectors from ring data." << endl;
        
        auto start = chrono::high_resolution_clock::now();
        cout << "\n=== Iterative Relation Solving ===" << endl;
        cout << "lev range: [" << lev_min << ", " << lev_max << "]" << endl;
        cout << "deg range: [" << deg_min << ", " << deg_max << "]" << endl;
        cout << "Iteration order: deg increases first, then lev" << endl;
        cout << endl;
        
        // 存储所有 (lev, deg) 级别的解
        vector<SolutionAtLevel> solutions;
        double total_computation_time = 0.0;
        
        // 配置采样参数基础版本
        RelationSolver::AdaptiveSamplingConfig base_config;
        base_config.min_nu = 3;
        base_config.max_nu = 50;              // 增加基础上限到 50（会被动态调整覆盖）
        base_config.nullity_stable_threshold = 2;
        base_config.check_interval = 1;
        base_config.verification_points = 1;
        
        // ===== 嵌套循环：外层 lev，内层 deg =====
        for (int current_lev = lev_min; current_lev <= lev_max; ++current_lev) {
            cout << "--- Processing lev = " << current_lev << " ---" << endl;
            
            for (int current_deg = deg_min; current_deg <= deg_max; ++current_deg) {
                cout << "  Solving (lev=" << current_lev << ", deg=" << current_deg << ")..." << endl;
                
                auto iter_start = chrono::high_resolution_clock::now();
                
                // 动态调整採样配置：根据 (lev, deg) 级别自动计算 max_nu
                int k_max = order - 1;  // 展开幂次上限
                int num_regimes = ringData.size();
                auto config = adjustSamplingConfig(
                    base_config,
                    current_lev, current_deg,
                    ne, nb, k_max, num_regimes
                );
                
                try {
                    // 调用高层 API 求解
                    auto [res, rel_coeff] = RelationSolver::reconstructReductionRelation<FFInt>(
                        allResults,
                        sector_list,
                        A_list,
                        Ainv_list,
                        ne, current_lev, current_deg,
                        config
                    );
                    
                    auto iter_end = chrono::high_resolution_clock::now();
                    double elapsed = chrono::duration<double>(iter_end - iter_start).count();
                    
                    // 记录解信息
                    SolutionAtLevel sol;
                    sol.lev = current_lev;
                    sol.deg = current_deg;
                    sol.num_variables = res.Mext.empty() ? 0 : res.Mext.size();
                    sol.solution_dimension = rel_coeff.getNumSolutions();
                    sol.computation_time = elapsed;
                    sol.system_rows = res.Mext.empty() ? 0 : res.Mext[0].size();
                    sol.system_cols = sol.num_variables;
                    sol.sampling_points_used = config.max_nu;  // 记录实际採样点数
                    
                    solutions.push_back(sol);
                    total_computation_time += elapsed;
                    
                    double ratio = sol.system_cols > 0 ? static_cast<double>(sol.system_rows) / sol.system_cols : 0.0;
                    cout << "OK (dim=" << sol.solution_dimension << ", ratio=" 
                         << fixed << setprecision(2) << ratio 
                         << ", nu=" << config.max_nu << ", time=" 
                         << fixed << setprecision(2) << elapsed << "s)" << endl;
                    
                } catch (const exception& e) {
                    cout << "FAILED (" << e.what() << ")" << endl;
                    // 继续处理下一个 (lev, deg)，不中断循环
                }
            }
        }
        
        // ---- 显示求解结果汇总 ----
        cout << "\n=== Solution Summary ===" << endl;
        cout << "Total configurations: " << solutions.size() << endl;
        cout << "Total computation time: " << fixed << setprecision(2) << total_computation_time << "s" << endl;
        cout << endl;
        
        for (const auto& sol : solutions) {
            sol.print();
        }
        
        // ---- 分析渐进性 ----
        if (!solutions.empty()) {
            cout << "\n=== Progressive Analysis ===" << endl;
            
            // 统计解维数的变化趋势
            vector<int> dims;
            for (const auto& sol : solutions) {
                dims.push_back(sol.solution_dimension);
            }
            
            cout << "Solution dimensions: [";
            for (size_t i = 0; i < dims.size(); ++i) {
                if (i > 0) cout << ", ";
                cout << dims[i];
            }
            cout << "]" << endl;
            
            // 显示第一个和最后一个解的详细信息
            cout << "\nFirst solution (lev=" << solutions.front().lev 
                 << ", deg=" << solutions.front().deg << "):" << endl;
            cout << "  Variables: " << solutions.front().num_variables << endl;
            cout << "  Equations: " << solutions.front().system_rows << endl;
            cout << "  Sampling points: " << solutions.front().sampling_points_used << endl;
            cout << "  Equation/Variable ratio: " << std::fixed << std::setprecision(3) 
                 << (static_cast<double>(solutions.front().system_rows) / solutions.front().system_cols) << endl;
            cout << "  Solution dim: " << solutions.front().solution_dimension << endl;
            
            cout << "\nLast solution (lev=" << solutions.back().lev 
                 << ", deg=" << solutions.back().deg << "):" << endl;
            cout << "  Variables: " << solutions.back().num_variables << endl;
            cout << "  Equations: " << solutions.back().system_rows << endl;
            cout << "  Sampling points: " << solutions.back().sampling_points_used << endl;
            cout << "  Equation/Variable ratio: " << std::fixed << std::setprecision(3) 
                 << (static_cast<double>(solutions.back().system_rows) / solutions.back().system_cols) << endl;
            cout << "  Solution dim: " << solutions.back().solution_dimension << endl;
            
            // 添加關鍵觀察
            cout << "\n=== Key Observations ===" << endl;
            int max_nu = 0;
            for (const auto& sol : solutions) {
                max_nu = std::max(max_nu, sol.sampling_points_used);
            }
            cout << "Maximum sampling points used: " << max_nu << endl;
            
            // 计算理论的方程数 vs 实际方程数
            int k_max = order - 1;
            int num_regimes = ringData.size();
            int eq_per_sample = num_regimes * nb * (k_max + 1);
            cout << "Theoretical equations from " << max_nu << " samples: " 
                 << (max_nu * eq_per_sample) << " (=" << max_nu 
                 << " × " << eq_per_sample << ")" << endl;
            cout << "Actual system size (last): " << solutions.back().system_rows << " × " 
                 << solutions.back().system_cols << endl;
            cout << "Status: " << (solutions.back().system_rows >= solutions.back().system_cols ? "✓ Overdetermined" : "✗ Underdetermined") << endl;
        }
        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> total = end - start;
        cout << "\nTotal execution time (including setup): " << total.count() << " seconds." << endl;
        cout << "Test completed successfully." << endl;
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
}