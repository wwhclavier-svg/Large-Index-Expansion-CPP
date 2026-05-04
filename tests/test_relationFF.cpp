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

// 将所有 (lev, deg) 结果导出到一个统一的 .m 文件（对应 MMA 的一次性输出）
template<typename T>
void exportAllResultsToMMA_SingleFile(
    const std::vector<RelationSolver::LevDegResult<T>>& all_results,
    int ne, int order, int modulus,
    const std::string& family,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);

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

    int count_with_solution = 0;
    for (const auto& ld : all_results) {
        if (ld.linear_result.hasSolution) ++count_with_solution;
    }

    out << "(* C++ Relation Export: All Results *)\n";
    out << "(* Family: " << family << ", Order: " << order
        << ", NE: " << ne << ", Modulus: " << modulus << " *)\n";
    out << "(* Total configurations: " << all_results.size()
        << ", with solution: " << count_with_solution << " *)\n";
    out << "$AllRelations = {\n";

    for (size_t idx = 0; idx < all_results.size(); ++idx) {
        const auto& ld = all_results[idx];
        const auto& alphas = ld.coeff.getAlphas();
        const auto& betas = ld.coeff.getBetas();
        size_t nAlpha = alphas.size();
        size_t nBeta = betas.size();
        size_t nSols = ld.coeff.getNumSolutions();
        size_t numVars = nAlpha * nBeta;

        out << "  <|\n";
        out << "    \"Family\" -> \"" << family << "\",\n";
        out << "    \"Lev\" -> " << ld.lev << ",\n";
        out << "    \"Deg\" -> " << ld.deg << ",\n";
        out << "    \"NE\" -> " << ne << ",\n";
        out << "    \"Modulus\" -> " << modulus << ",\n";
        out << "    \"Order\" -> " << order << ",\n";
        out << "    \"StableOrder\" -> " << ld.stable_order << ",\n";
        out << "    \"NumRelations\" -> " << ld.num_relations << ",\n";
        out << "    \"ActiveVars\" -> " << ld.active_vars << ",\n";
        out << "    \"TotalVars\" -> " << ld.total_vars << ",\n";

        // Alphas
        out << "    \"Alphas\" -> {";
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

        // Betas
        out << "    \"Betas\" -> {";
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

        // Coefficients matrix
        out << "    \"Coefficients\" -> {\n";
        for (size_t i = 0; i < numVars; ++i) {
            out << "      {";
            if (ld.linear_result.hasSolution && i < ld.linear_result.Mext.size() && 0 < ld.linear_result.Mext[i].size()) {
                writeValue(ld.linear_result.Mext[i][0]);
            } else {
                out << "0";
            }
            for (size_t s = 1; s <= nSols; ++s) {
                out << ", ";
                if (i < ld.linear_result.Mext.size() && s < ld.linear_result.Mext[i].size()) {
                    writeValue(ld.linear_result.Mext[i][s]);
                } else {
                    out << "0";
                }
            }
            out << "}";
            if (i < numVars - 1) out << ",";
            out << "\n";
        }
        out << "    },\n";

        // Independent pairs
        out << "    \"IndependentPairs\" -> {";
        for (size_t p = 0; p < ld.independent_pairs.size(); ++p) {
            const auto& [alpha, beta] = ld.independent_pairs[p];
            out << "{{";
            for (int i = 0; i < ne; ++i) {
                if (i > 0) out << ",";
                out << alpha[i];
            }
            out << "},{";
            for (int i = 0; i < ne; ++i) {
                if (i > 0) out << ",";
                out << beta[i];
            }
            out << "}}";
            if (p < ld.independent_pairs.size() - 1) out << ", ";
        }
        out << "},\n";

        out << "    \"HasSolution\" -> " << (ld.linear_result.hasSolution ? "True" : "False") << ",\n";
        out << "    \"NumVariables\" -> " << numVars << ",\n";
        out << "    \"NumSolutions\" -> " << nSols << "\n";
        out << "  |>";

        if (idx < all_results.size() - 1) out << ",";
        out << "\n";
    }

    out << "};\n";
    out.close();
    std::cout << "\nExported all " << all_results.size() << " results to: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <family_name> [order] [lev_min] [lev_max] [deg_max]" << endl;
        cerr << "  family_name: e.g. bub, bub0, DBtop" << endl;
        cerr << "  order      : expansion order (default: 4)" << endl;
        cerr << "  lev_min    : min |alpha| (default: 1)" << endl;
        cerr << "  lev_max    : max |alpha| (default: 2)" << endl;
        cerr << "  deg_max    : max |beta| (default: 2)" << endl;
        cerr << endl;
        cerr << "Example: " << argv[0] << " bub" << endl;
        cerr << "         " << argv[0] << " DBtop 5 1 2 1" << endl;
        return 1;
    }

    string family = argv[1];
    int order = (argc > 2) ? stoi(argv[2]) : 4;
    int lev_min = (argc > 3) ? stoi(argv[3]) : 1;
    int lev_max = (argc > 4) ? stoi(argv[4]) : 2;
    int deg_min = 0;
    int deg_max = (argc > 5) ? stoi(argv[5]) : 2;
    const int incre = 2;
    
    initBinomial();

    const string binIBPFile = "IBPMat_"+family+".bin";
    const string binRingFile = "RingData_"+family+".bin";
    const string coeffCacheFile = "ExpansionCache_"+family+"_k"+to_string(order)+".bin";

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
        cout << "\n=== Relation Solving (RemoveSolvedVariables strategy) ===" << endl;

        RelationSolver::AdaptiveSamplingConfig base_config;
        base_config.min_nu = 3;
        base_config.max_nu = 50;
        base_config.nullity_stable_threshold = 5;
        base_config.check_interval = 5;
        base_config.verification_points = 5;
        base_config.plateau_size = 1;       // MMA 默认: 1 个额外阶数确认稳定性
        base_config.random_min = 0;
        base_config.random_max = static_cast<double>(FFInt::p - 1);

        try {
            auto all_results = RelationSolver::reconstructAllRelations<FFInt>(
                allResults, sector_list, A_list, Ainv_list,
                ne, lev_min, lev_max, deg_max, base_config
            );

            // Collect solution summaries and export as single unified file
            vector<SolutionAtLevel> solutions;

            for (const auto& ld : all_results) {
                SolutionAtLevel sol;
                sol.lev = ld.lev;
                sol.deg = ld.deg;
                sol.num_variables = ld.active_vars;
                sol.solution_dimension = ld.coeff.getNumSolutions();
                sol.computation_time = 0;
                sol.system_rows = 0;
                sol.system_cols = ld.active_vars;
                sol.sampling_points_used = 0;
                solutions.push_back(sol);

                cout << "  (lev=" << ld.lev << ", deg=" << ld.deg << "): "
                     << "active=" << ld.active_vars << "/" << ld.total_vars
                     << " sol_dim=" << sol.solution_dimension
                     << " independent=" << ld.independent_pairs.size()
                     << " stable_order=" << ld.stable_order << endl;
            }

            // Export single unified file containing ALL results
            string unifiedFilename = "AllRelations_" + family + "_k" + to_string(order) + ".m";
            exportAllResultsToMMA_SingleFile(all_results, ne, order,
                static_cast<int>(FFInt::p), family, unifiedFilename);

            cout << "\n=== Solution Summary ===" << endl;
            cout << "Total configurations: " << solutions.size() << endl;
            for (const auto& sol : solutions) {
                sol.print();
            }

            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> total = end - start;
            cout << "\nTotal execution time: " << total.count() << " seconds." << endl;
            cout << "Test completed successfully." << endl;
            return 0;

        } catch (const exception& e) {
            cout << "FAILED (" << e.what() << ")" << endl;
        }
    }
    catch (const std::exception& e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
}
