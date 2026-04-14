#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <map>

#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "RingDataLoader.hpp"
#include "RelationSolver.hpp"
#include "Combinatorics.hpp"
#include "LinearSolver.hpp"
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"

using namespace std;
using namespace firefly;

RelationSolver::AdaptiveSamplingConfig 
adjustSamplingConfig(
    const RelationSolver::AdaptiveSamplingConfig& base_config,
    int lev, int deg, int ne, int nb, int k_max, int num_regimes)
{
    auto config = base_config;
    config.lev_hint = lev;
    config.deg_hint = deg;
    
    std::vector<int> temp;
    std::vector<std::vector<int>> alphas, betas;
    RelationSolver::generateAllIndices(ne, lev, temp, alphas, false);
    temp.clear();
    RelationSolver::generateAllIndices(ne, deg, temp, betas, false);
    
    int num_variables = alphas.size() * betas.size();
    int equations_per_sample = num_regimes * nb * (k_max + 1);
    
    static constexpr double SAFETY_FACTOR = 1.2;
    int min_samples_required = static_cast<int>(
        std::ceil(num_variables * SAFETY_FACTOR / equations_per_sample)
    );
    
    static constexpr int MAX_SAMPLES_LIMIT = 150;
    config.max_nu = std::min(
        MAX_SAMPLES_LIMIT,
        std::max(base_config.max_nu, min_samples_required)
    );
    
    config.min_nu = std::min(
        config.min_nu,
        std::max(3, min_samples_required / 3)
    );
    
    std::cout << "    Adaptive sampling config:" << std::endl;
    std::cout << "      Variables: " << num_variables << std::endl;
    std::cout << "      Equations per sample: " << equations_per_sample << std::endl;
    std::cout << "      Min samples required (" << std::fixed << std::setprecision(1) << (SAFETY_FACTOR*100) << "% safety): " << min_samples_required << std::endl;
    std::cout << "      Adjusted max_nu: " << config.max_nu << std::endl;
    
    return config;
}

struct SolutionAtLevel {
    int lev;
    int deg;
    int num_variables;
    int solution_dimension;
    double computation_time;
    int system_rows;
    int system_cols;
    int sampling_points_used;
    
    SolutionAtLevel() : sampling_points_used(0) {}
    
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <family_name> [order] [lev_max] [deg_max]" << endl;
        cerr << "  family_name: e.g. bub, bub0, DBtop" << endl;
        cerr << "  order      : expansion order (default: 4)" << endl;
        cerr << "  lev_max    : max |alpha| (default: 2)" << endl;
        cerr << "  deg_max    : max |beta| (default: 2)" << endl;
        cerr << endl;
        cerr << "Example: " << argv[0] << " bub" << endl;
        cerr << "         " << argv[0] << " DBtop 4 2 1" << endl;
        return 1;
    }

    string family = argv[1];
    int order = (argc > 2) ? stoi(argv[2]) : 4;
    int lev_min = 1;
    int lev_max = (argc > 3) ? stoi(argv[3]) : 2;
    int deg_min = 0;
    int deg_max = (argc > 4) ? stoi(argv[4]) : 2;
    const int incre = 2;
    
    initBinomial();

    const string binIBPFile = "IBPMat_"+family+".bin";
    const string binRingFile = "RingData_"+family+".bin";
    const string coeffCacheFile = "ExpansionCache_"+family+".bin";

    cout << "=== Configuration ===" << endl;
    cout << "Family: " << family << endl;
    cout << "Order: " << order << endl;
    cout << "lev range: [" << lev_min << ", " << lev_max << "]" << endl;
    cout << "deg range: [" << deg_min << ", " << deg_max << "]" << endl;
    cout << "Files: " << binIBPFile << ", " << binRingFile << endl;
    cout << endl;

    try {
        cout << "=== Test: Reconstruction of Linear Relations (Finite Field) ===" << endl;

        // ---- 第1a步：加载 IBP 矩阵 ----
        cout << "Loading IBP matrices from " << binIBPFile << " ..." << endl;
        auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(binIBPFile);
        if (ibpmatlist.empty()) {
            cerr << "Error: No IBP matrices loaded." << endl;
            return 1;
        }
        int ne = ibpmatlist[0].ne;
        int nb = ibpmatlist[0].nb;
        cout << "Loaded: ne=" << ne << ", nb=" << nb << ", mod=" << FFInt::p << endl;

        // ---- 第1b步：计算/读取展开系数 ----
        vector<vector<seriesCoefficient<FFInt>>> allResults;
        try {
            allResults = SeriesIO::loadAllResults<FFInt>(coeffCacheFile);
            cout << "Loaded coefficients from cache." << endl;
        } catch (const std::exception& e) {
            cout << "Computing expansion coefficients (order=" << order << ", incre=" << incre << ")..." << endl;
            auto t0 = chrono::high_resolution_clock::now();
            allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
            auto t1 = chrono::high_resolution_clock::now();
            double dt = chrono::duration<double>(t1 - t0).count();
            cout << "Recursion completed in " << fixed << setprecision(3) << dt << "s" << endl;
            SeriesIO::saveAllResults(allResults, coeffCacheFile);
            cout << "Saved to cache: " << coeffCacheFile << endl;
        }
        cout << "Total branches: " << allResults.size() << " matrices" << endl;

        // ---- 第2步：加载环数据 ----
        cout << "\nLoading ring data from " << binRingFile << " ..." << endl;
        auto ringData = AlgebraData::LoadBinary<FFInt>(binRingFile, FFInt::p);
        if (ringData.empty()) {
            cerr << "Error: No ring data loaded." << endl;
            return 1;
        }
        cout << "Loaded " << ringData.size() << " regimes." << endl;

        // ---- 数据提取 ----
        std::vector<std::vector<int>> sector_list;
        std::vector<std::vector<std::vector<FFInt>>> A_list, Ainv_list;
        for (const auto& ring : ringData) {
            sector_list.push_back(ring.limitSector);
            A_list.push_back(ring.A_list);
            Ainv_list.push_back(ring.Ainv_list);
        }
        cout << "Extracted " << sector_list.size() << " sectors." << endl;
        
        auto start = chrono::high_resolution_clock::now();
        cout << "\n=== Iterative Relation Solving ===" << endl;
        
        vector<SolutionAtLevel> solutions;
        double total_computation_time = 0.0;
        
        RelationSolver::AdaptiveSamplingConfig base_config;
        base_config.min_nu = 3;
        base_config.max_nu = 50;
        base_config.nullity_stable_threshold = 2;
        base_config.check_interval = 1;
        base_config.verification_points = 1;
        
        for (int current_lev = lev_min; current_lev <= lev_max; ++current_lev) {
            cout << "--- Processing lev = " << current_lev << " ---" << endl;
            
            for (int current_deg = deg_min; current_deg <= deg_max; ++current_deg) {
                cout << "  Solving (lev=" << current_lev << ", deg=" << current_deg << ")..." << endl;
                
                auto iter_start = chrono::high_resolution_clock::now();
                
                int k_max = order - 1;
                int num_regimes = ringData.size();
                auto config = adjustSamplingConfig(
                    base_config, current_lev, current_deg,
                    ne, nb, k_max, num_regimes
                );
                
                try {
                    auto [res, rel_coeff] = RelationSolver::reconstructReductionRelation<FFInt>(
                        allResults, sector_list, A_list, Ainv_list,
                        ne, current_lev, current_deg, config
                    );
                    
                    auto iter_end = chrono::high_resolution_clock::now();
                    double elapsed = chrono::duration<double>(iter_end - iter_start).count();
                    
                    SolutionAtLevel sol;
                    sol.lev = current_lev;
                    sol.deg = current_deg;
                    sol.num_variables = res.Mext.empty() ? 0 : res.Mext.size();
                    sol.solution_dimension = rel_coeff.getNumSolutions();
                    sol.computation_time = elapsed;
                    sol.system_rows = res.Mext.empty() ? 0 : res.Mext[0].size();
                    sol.system_cols = sol.num_variables;
                    sol.sampling_points_used = config.max_nu;
                    
                    solutions.push_back(sol);
                    total_computation_time += elapsed;
                    
                    double ratio = sol.system_cols > 0 ? static_cast<double>(sol.system_rows) / sol.system_cols : 0.0;
                    cout << "OK (dim=" << sol.solution_dimension << ", ratio=" 
                         << fixed << setprecision(2) << ratio 
                         << ", nu=" << config.max_nu << ", time=" 
                         << fixed << setprecision(2) << elapsed << "s)" << endl;
                    
                } catch (const exception& e) {
                    cout << "FAILED (" << e.what() << ")" << endl;
                }
            }
        }
        
        cout << "\n=== Solution Summary ===" << endl;
        cout << "Total configurations: " << solutions.size() << endl;
        cout << "Total computation time: " << fixed << setprecision(2) << total_computation_time << "s" << endl;
        for (const auto& sol : solutions) {
            sol.print();
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
