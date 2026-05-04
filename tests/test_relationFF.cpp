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
    int equations_per_sample = num_regimes * nb * (deg + k_max + 1);
    
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

// 辅助函数：将关系系数导出为 MMA 多项式格式 Sum c(a,b) v^b j[a]
template<typename T>
void exportRelationToMMA_Polynomial(
    const RelationSolver::RelationCoefficient<T>& rel_coeff,
    const ::LinearSystemResult<T>& linear_result,
    int lev, int deg, int ne,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);

    const auto& alphas = rel_coeff.getAlphas();
    const auto& betas = rel_coeff.getBetas();
    size_t nAlpha = alphas.size();
    size_t nBeta = betas.size();
    size_t nSols = rel_coeff.getNumSolutions();

    auto isZero = [](const T& val) -> bool {
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            return val.n == 0;
        } else {
            return val == T(0);
        }
    };

    auto writeValue = [&out](const T& val) {
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            out << val.n;
        } else {
            out << val;
        }
    };

    out << "(* C++ Relation Export: Polynomial Form *)\n";
    out << "$RelationResult = <|\n";
    out << "  \"Lev\" -> " << lev << ", \"Deg\" -> " << deg << ",\n";
    out << "  \"HasSolution\" -> " << (linear_result.hasSolution ? "True" : "False") << ",\n";
    out << "  \"NumVariables\" -> " << (nAlpha * nBeta) << ",\n";
    out << "  \"NumSolutions\" -> " << nSols << ",\n";
    out << "  \"Relations\" -> {\n";
    for (size_t sol = 0; sol <= nSols; ++sol) {
        out << "    ";
        bool first_term = true;
        for (size_t a = 0; a < nAlpha; ++a) {
            for (size_t b = 0; b < nBeta; ++b) {
                size_t idx = a * nBeta + b;
                const T& val = linear_result.Mext[idx][sol];
                if (isZero(val)) continue;

                if (!first_term) out << " + ";
                first_term = false;

                writeValue(val);

                for (int i = 0; i < ne; ++i) {
                    int exp = betas[b][i];
                    if (exp == 1) out << "*v" << (i + 1);
                    else if (exp > 1) out << "*v" << (i + 1) << "^" << exp;
                }

                out << "*g[";
                for (int i = 0; i < ne; ++i) {
                    if (i > 0) out << ",";
                    int alpha_i = alphas[a][i];
                    if (alpha_i == 0) {
                        out << "v" << (i + 1);
                    } else if (alpha_i == 1) {
                        out << "v" << (i + 1) << "-1";
                    } else {
                        out << "v" << (i + 1) << "-" << alpha_i;
                    }
                }
                out << "]";
            }
        }
        if (first_term) out << "0";
        if (sol < nSols) out << ",";
        out << "\n";
    }
    out << "  }\n";
    out << "|>;\n";
}

// 统一格式导出函数：包含完整的系数矩阵和元数据
template<typename T>
void exportRelationToMMA_Unified(
    const RelationSolver::RelationCoefficient<T>& rel_coeff,
    const ::LinearSystemResult<T>& linear_result,
    int lev, int deg, int ne, int modulus,
    const std::string& family,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);

    const auto& alphas = rel_coeff.getAlphas();
    const auto& betas = rel_coeff.getBetas();
    size_t nAlpha = alphas.size();
    size_t nBeta = betas.size();
    size_t nSols = rel_coeff.getNumSolutions();
    size_t numVars = nAlpha * nBeta;

    auto isZero = [](const T& val) -> bool {
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            return val.n == 0;
        } else {
            return val == T(0);
        }
    };

    auto writeValue = [&out](const T& val) {
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            out << val.n;
        } else {
            out << val;
        }
    };

    out << "(* C++ Relation Export: Unified Format *)\n";
    out << "(* Family: " << family << ", Lev=" << lev << ", Deg=" << deg << " *)\n";
    out << "$RelationResult = <|\n";
    out << "  \"Family\" -> \"" << family << "\",\n";
    out << "  \"Lev\" -> " << lev << ",\n";
    out << "  \"Deg\" -> " << deg << ",\n";
    out << "  \"NE\" -> " << ne << ",\n";
    out << "  \"Modulus\" -> " << modulus << ",\n";

    // Alphas list
    out << "  \"Alphas\" -> {";
    for (size_t a = 0; a < nAlpha; ++a) {
        out << "{";
        for (int i = 0; i < ne; ++i) {
            if (i > 0) out << ",";
            out << alphas[a][i];
        }
        out << "}";
        if (a < nAlpha - 1) out << ", ";
    }
    out << "},\n";

    // Betas list
    out << "  \"Betas\" -> {";
    for (size_t b = 0; b < nBeta; ++b) {
        out << "{";
        for (int i = 0; i < ne; ++i) {
            if (i > 0) out << ",";
            out << betas[b][i];
        }
        out << "}";
        if (b < nBeta - 1) out << ", ";
    }
    out << "},\n";

    // Coefficients matrix [numVars] x [1 + numSolutions]
    out << "  \"Coefficients\" -> {\n";
    for (size_t i = 0; i < numVars; ++i) {
        out << "    {";
        // Column 0: particular solution (or 0 for homogeneous)
        if (linear_result.hasSolution && i < linear_result.Mext.size() && 0 < linear_result.Mext[i].size()) {
            writeValue(linear_result.Mext[i][0]);
        } else {
            out << "0";
        }
        // Columns 1..nSols: nullspace basis vectors
        for (size_t s = 1; s <= nSols; ++s) {
            out << ", ";
            if (i < linear_result.Mext.size() && s < linear_result.Mext[i].size()) {
                writeValue(linear_result.Mext[i][s]);
            } else {
                out << "0";
            }
        }
        out << "}";
        if (i < numVars - 1) out << ",";
        out << "\n";
    }
    out << "  },\n";

    // Relations (polynomial form for human readability)
    out << "  \"Relations\" -> {\n";
    for (size_t sol = 0; sol <= nSols; ++sol) {
        out << "    ";
        if (sol == 0) {
            out << "\"Particular\" -> ";
        } else {
            out << "\"Basis" << sol << "\" -> ";
        }
        out << "\"";
        bool first_term = true;
        for (size_t a = 0; a < nAlpha; ++a) {
            for (size_t b = 0; b < nBeta; ++b) {
                size_t idx = a * nBeta + b;
                const T& val = (idx < linear_result.Mext.size() && sol < linear_result.Mext[idx].size())
                    ? linear_result.Mext[idx][sol] : T(0);
                if (isZero(val)) continue;

                if (!first_term) out << " + ";
                first_term = false;

                writeValue(val);

                // beta exponents (v variables)
                for (int i = 0; i < ne; ++i) {
                    int exp = betas[b][i];
                    if (exp == 1) out << "*v" << (i + 1);
                    else if (exp > 1) out << "*v" << (i + 1) << "^" << exp;
                }

                // alpha index: g[v1-α1, v2-α2, ...] notation
                out << "*g[";
                for (int i = 0; i < ne; ++i) {
                    if (i > 0) out << ",";
                    int alpha_i = alphas[a][i];
                    if (alpha_i == 0) {
                        out << "v" << (i + 1);
                    } else if (alpha_i == 1) {
                        out << "v" << (i + 1) << "-1";
                    } else {
                        out << "v" << (i + 1) << "-" << alpha_i;
                    }
                }
                out << "]";
            }
        }
        if (first_term) out << "0";
        out << "\"";
        if (sol < nSols) out << ",";
        out << "\n";
    }
    out << "  },\n";

    // Solution info
    out << "  \"HasSolution\" -> " << (linear_result.hasSolution ? "True" : "False") << ",\n";
    out << "  \"NumVariables\" -> " << numVars << ",\n";
    out << "  \"NumSolutions\" -> " << nSols << "\n";
    out << "|>;\n";
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
        bool has_nontrivial = solution_dimension > 0;
        std::string status = has_nontrivial ? "✓" : "✗";
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
    int lev_min = 0;
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
            string mmaCacheFile = "ExpansionMMA_" + family + ".m";
            SeriesIO::exportAllResultsToMMA(allResults, mmaCacheFile);
            cout << "Exported to MMA format: " << mmaCacheFile << endl;
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
        for (size_t i = 0; i < sector_list.size(); ++i) {
            cout << "  Sector " << i << ": [";
            for (size_t j = 0; j < sector_list[i].size(); ++j) {
                if (j > 0) cout << ", ";
                cout << sector_list[i][j];
            }
            cout << "]" << endl;
        }

        auto start = chrono::high_resolution_clock::now();
        cout << "\n=== Iterative Relation Solving ===" << endl;
        
        vector<SolutionAtLevel> solutions;
        double total_computation_time = 0.0;
        
        RelationSolver::AdaptiveSamplingConfig base_config;
        base_config.min_nu = 3;
        base_config.max_nu = 50;
        base_config.nullity_stable_threshold = 5;
        base_config.check_interval = 5;
        base_config.verification_points = 5;
        base_config.random_min = 0;
        base_config.random_max = static_cast<double>(FFInt::p - 1);
        
        for (int current_lev = lev_min; current_lev <= lev_max; ++current_lev) {
            cout << "--- Processing lev = " << current_lev << " ---" << endl;
            
            for (int current_deg = deg_min; current_deg <= deg_max; ++current_deg) {
                cout << "  Solving (lev=" << current_lev << ", deg=" << current_deg << ")..." << endl;
                
                auto iter_start = chrono::high_resolution_clock::now();
                
                int k_max = order;
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

                    if (res.hasSolution) {
                        // Export polynomial format (legacy)
                        string relFilename = "Relations_" + family + "_lev" + to_string(current_lev) + "_deg" + to_string(current_deg) + ".m";
                        exportRelationToMMA_Polynomial(rel_coeff, res, current_lev, current_deg, ne, relFilename);
                        cout << "  Exported relation to: " << relFilename << endl;

                        // Export unified format (for comparison with MMA)
                        string unifiedFilename = "Compare-CPPRelation-" + family + "_lev" + to_string(current_lev) + "_deg" + to_string(current_deg) + ".m";
                        exportRelationToMMA_Unified(rel_coeff, res, current_lev, current_deg, ne,
                            static_cast<int>(FFInt::p), family, unifiedFilename);
                        cout << "  Exported unified format to: " << unifiedFilename << endl;

                        // DIAGNOSTIC: system rank info for deg=0
                        if (current_deg == 0 && !res.Mext.empty()) {
                            cout << "    [DIAG] variables=" << res.Mext.size() 
                                 << ", equations=" << res.Mext[0].size()
                                 << ", solutions=" << rel_coeff.getNumSolutions()
                                 << ", pivots=" << res.pivot_cols.size() << endl;
                        }

                    }

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
