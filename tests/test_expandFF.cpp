#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <filesystem>

#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "Combinatorics.hpp"
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"


using namespace std;
using namespace firefly;

/**
 * 导出 C++ 元数据 + k=0,1 系数到 Compare-CPPMeta-<famname>.m
 * 供 Compare-VerifyLog.wl 读取生成统一验证日志
 */
template <typename T>
void exportMetaToMMA(const vector<vector<seriesCoefficient<T>>>& allResults,
                     const vector<IBPMatrixE<T>>& ibpmatlist,
                     double totalTime, const string& famname) {
    string verifyDir = "verify/" + famname + "/";
    std::filesystem::create_directories(verifyDir);
    string metafile = verifyDir + "Compare-CPPMeta-" + famname + ".m";
    ofstream out(metafile);
    if (!out) {
        cerr << "Warning: Cannot open metadata file for writing: " << metafile << endl;
        return;
    }

    out << "(* C++ Verify Metadata Export *)\n";
    out << "$CPPMeta = <|\n";
    out << "  \"Family\" -> \"" << famname << "\",\n";
    out << "  \"TotalTime\" -> " << fixed << setprecision(6) << totalTime << ",\n";
    out << "  \"Modulus\" -> " << FFInt::p << ",\n";

    // Per-region metadata
    out << "  \"Regions\" -> {\n";
    for (size_t i = 0; i < ibpmatlist.size(); ++i) {
        const auto& mat = ibpmatlist[i];
        int nimax = (i < allResults.size() && !allResults[i].empty())
                    ? allResults[i][0].getNimax() : 0;
        out << "    <| \"RegionIndex\" -> " << (i + 1)
            << ", \"NE\" -> " << mat.ne
            << ", \"NB\" -> " << mat.nb
            << ", \"NIBP\" -> " << mat.nibp
            << ", \"Incre\" -> " << mat.incre
            << ", \"Nimax\" -> " << nimax << " |>";
        if (i < ibpmatlist.size() - 1) out << ",";
        out << "\n";
    }
    out << "  },\n";

    // k=0,1 coefficients from the first region's first solution
    out << "  \"Coefficients\" -> {\n";
    if (!allResults.empty() && !allResults[0].empty()) {
        const auto& coeff = allResults[0][0];
        int kmax = coeff.getKmax();
        int ne = coeff.getNe();
        int nb = coeff.getNb();
        int incre = coeff.getIncre();
        extern long long BINOM[MAX_VAL][MAX_VAL];

        for (int k = 0; k <= min(1, kmax); ++k) {
            out << "    { (* k=" << k << " *)\n";
            int lmax = incre * k;
            bool firstTerm = true;
            for (int l = 0; l <= lmax; ++l) {
                long long nSeeds = getCapacity(ne, l);
                for (long long cid = 0; cid < nSeeds; ++cid) {
                    auto seed = readIndex(cid, l, ne);
                    // 检查该 seed 是否有非零分量
                    bool anyNonZero = false;
                    for (int j = 0; j < nb; ++j) {
                        if (coeff(k, l, cid, j, 0) != T(0)) {
                            anyNonZero = true;
                            break;
                        }
                    }
                    if (!anyNonZero) continue;
                    if (!firstTerm) out << " + ";
                    firstTerm = false;
                    // 输出 nb 维向量表示扩域元素
                    if (nb == 1) {
                        if constexpr (std::is_same_v<T, firefly::FFInt>) {
                            out << coeff(k, l, cid, 0, 0).n;
                        } else {
                            out << coeff(k, l, cid, 0, 0);
                        }
                    } else {
                        out << "{";
                        for (int j = 0; j < nb; ++j) {
                            if (j > 0) out << ", ";
                            const T& val = coeff(k, l, cid, j, 0);
                            if constexpr (std::is_same_v<T, firefly::FFInt>) {
                                out << val.n;
                            } else {
                                out << val;
                            }
                        }
                        out << "}";
                    }
                    for (int v = 0; v < ne; ++v) {
                        if (seed[v] > 0) {
                            out << "*v" << (v + 1);
                            if (seed[v] > 1) out << "^" << seed[v];
                        }
                    }
                }
            }
            if (firstTerm) out << "0";
            out << "\n    }";
            if (k < min(1, kmax)) out << ",";
            out << "\n";
        }
    }
    out << "  }\n";

    out << "|>\n";
    out.close();
    cout << "\nMetadata exported to " << metafile << endl;
}

template <typename T>
void printM1(std::vector<IBPMatrixE<T>> ibpmatlist) {
    if (!ibpmatlist.empty()) {
        const auto& firstMat = ibpmatlist[0];
        std::cout << "\nFirst M1 matrix: dimensions = ["
                << firstMat.M1.size() << ", ";
        if (!firstMat.M1.empty() && !firstMat.M1[0].empty()) {
            std::cout << firstMat.M1[0].size() << ", "
                    << firstMat.M1[0][0].size() << "]" << std::endl;

            // 打印前 10 个非零元素
            int count = 0;
            for (size_t i = 0; i < firstMat.M1.size() && count < 10; ++i) {
                for (size_t j = 0; j < firstMat.M1[i].size() && count < 10; ++j) {
                    for (size_t k = 0; k < firstMat.M1[i][j].size() && count < 10; ++k) {
                        const auto& val = firstMat.M1[i][j][k];
                        if (val != 0 || true) {   // 假设 FFInt 可比较零
                            std::cout << "  M1[" << i << "][" << j << "][" << k << "] = "
                                    << val << std::endl;
                            ++count;
                        }
                    }
                }
            }
            if (count == 0) std::cout << "  M1 has no non-zero elements." << std::endl;
        } else {
            std::cout << "M1 is empty or invalid." << std::endl;
        }

        // Print K1s matrix
        std::cout << "\nK1s matrix: dimensions = ["
                << firstMat.K1s.size() << ", " << firstMat.K1s[0].size() << ", " << firstMat.K1s[0][0].size() << "]" << std::endl;
        for (size_t i = 0; i < firstMat.K1s.size(); ++i) {
            for (size_t j = 0; j < firstMat.K1s[i].size(); ++j) {
                for (size_t k = 0; k < firstMat.K1s[i][j].size(); ++k) {
                    std::cout << "  K1s[" << i << "][" << j << "][" << k << "] = " << firstMat.K1s[i][j][k] << std::endl;
                }
            }
        }

        // Print K2s matrix
        std::cout << "\nK2s matrix: dimensions = ["
                << firstMat.K2s.size() << ", " << firstMat.K2s[0].size() << ", " << firstMat.K2s[0][0].size() << "]" << std::endl;
        for (size_t i = 0; i < firstMat.K2s.size(); ++i) {
            for (size_t j = 0; j < firstMat.K2s[i].size(); ++j) {
                for (size_t k = 0; k < firstMat.K2s[i][j].size(); ++k) {
                    std::cout << "  K2s[" << i << "][" << j << "][" << k << "] = " << firstMat.K2s[i][j][k] << std::endl;
                }
            }
        }
    } else {
        std::cout << "No matrices loaded." << std::endl;
    }
}

template <typename T>
void printResult(const std::vector<std::vector<seriesCoefficient<T>>>& allResults) {
    if (allResults.empty() || allResults[0].empty()) {
        std::cout << "No results to print." << std::endl;
        return;
    }

    // 取第一个矩阵的第一个解
    const auto& sol = allResults[0][0];

    // 获取维度信息（需要上述 getter）
    int kmax = sol.getKmax();
    int incre = sol.getIncre();
    int ne = sol.getNe();
    int nb = sol.getNb();
    int nimax = sol.getNimax();   // 通解个数

    // 组合数表（假设已在某处定义并初始化）
    extern long long BINOM[MAX_VAL][MAX_VAL];

    std::cout << "\nFirst solution coefficients (first 20):" << std::endl;
    int count = 0, max_count = 40;

    // 遍历所有阶、层、种子、底空间分量，只打印第 0 个独立元（nimax 索引 0）
    for (int k = 0; k <= kmax && count < max_count; ++k) {
        int lmax = incre * k;
        for (int l = 0; l <= lmax && count < max_count; ++l) {
            long long num_cid = BINOM[l + ne - 1][ne - 1];  // 该层的种子数
            for (int cid = 0; cid < num_cid && count < max_count; ++cid) {
                std::vector<int> seed = readIndex(cid, l, ne);
                for (int j = 0; j < nb && count < max_count; ++j) {
                    // 获取系数值（只取独立元索引 0）
                    T val = sol(k, l, cid, j, 0);

                    std::ostringstream oss;
                    oss << "[";
                    for (int d = 0; d < ne; ++d) {
                        if (d > 0) oss << ",";
                        oss << seed[d];
                    }
                    oss << "]";
                    std::cout << "  coeff[" << count << "] (k=" << k
                              << ", l=" << l << ", seed=" << oss.str()
                              << ", j=" << j << ") = " << val << std::endl;
                    ++count;
                }
            }
        }
    }

    if (count == 20) {
        std::cout << "  ... (more coefficients not shown)" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // 1. 设置有限域素数
    FFInt::set_new_prime(179424673);

    // 2. 初始化组合数表
    initBinomial();

    // 3. 从命令行参数获取 famname 和 order（默认 bub00, order=4）
    string famname = "bub00";
    if (argc > 1) {
        famname = argv[1];
    }
    int order = 4;
    if (argc > 2) {
        order = atoi(argv[2]);
    }
    cout << "Using family: " << famname << ", order: " << order << endl;
    const int incre = 2;
    string verifyDir = "verify/" + famname + "/";
    std::filesystem::create_directories(verifyDir);
    const string filename = string("data/IBPMat_") + famname + ".bin";
    const string coeffCacheFile = string("data/ExpansionCache_") + famname + "_k" + to_string(order) + ".bin";

    FFInt test_div = FFInt(1) / FFInt(2);
    std::cout << "[SANITY] 1/2 mod p = " << test_div.n << " (expected 89712337 for correct mod inverse)" << std::endl;

    try {
        cout << "=== Test: Expand Coefficients over Finite Field (FFInt) ===" << endl;
        cout << "Using famname: " << famname << endl;
        cout << "Filename: " << filename << endl;
        cout << "CoeffCacheFile: " << coeffCacheFile << endl;

        // 4. 加载 IBP 矩阵数据，模板参数为 FFInt
        //auto ibpmatlist = loadAllIBPMatrices<FFInt>(filename);
        auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(filename);
        if (ibpmatlist.empty()) {
            cerr << "Error: No matrices loaded." << endl;
            return 1;
        }
        cout << "prime: " << FFInt::p << endl;

        printM1(ibpmatlist);

        cout << "Loaded " << ibpmatlist.size() << " matrices." << endl;
        cout << "First matrix: ne = " << ibpmatlist[0].ne
             << ", nb = " << ibpmatlist[0].nb
             << ", nibp = " << ibpmatlist[0].nibp << endl;

        // 5. 执行批处理递归（自动使用 FFInt 版本的求解器）
        auto start = chrono::high_resolution_clock::now();
        auto allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end - start;

        // 6. 输出统计信息
        cout << "\n=== Results ===" << endl;
        cout << "Total time: " << fixed << setprecision(4) << diff.count() << " seconds." << endl;
        cout << "Number of solution branches per matrix:" << endl;
        for (size_t i = 0; i < allResults.size(); ++i) {
            cout << "  Matrix " << i+1 << ": " << allResults[i].size() << " solutions" << endl;
        }

        printResult(allResults);
        SeriesIO::saveAllResults(allResults, coeffCacheFile);
        SeriesIO::exportAllResultsToMMA(allResults, verifyDir + "Compare-CPPResult-" + famname + ".m");
        exportMetaToMMA(allResults, ibpmatlist, diff.count(), famname);

        cout << "\nTest completed successfully." << endl;
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
}
