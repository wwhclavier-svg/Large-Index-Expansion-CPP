/**
 * test_expand_family.cpp
 * 通用多族展开系数验证测试
 * 
 * 用法: ./test_expand_family <family_name> [order]
 * 
 * 支持的族:
 *   - bub, bub00, Bub_Numeric (已有数据)
 *   - DB (已有数据，但计算量大，建议 order<=2)
 *   - tri, box (需要先生成数据)
 * 
 * 验证内容:
 *   1. 矩阵加载成功
 *   2. 展开完成无异常
 *   3. order=0 的初始条件 C(0,0,seed={0..},j)=1
 *   4. 系数可导出为 MMA 格式
 */

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "Combinatorics.hpp"
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"

using namespace std;
using namespace firefly;

struct TestResult {
    string family;
    int order;
    bool loadOk;
    bool expandOk;
    bool initCondOk;
    double timeMs;
    string errorMsg;
    size_t numMatrices;
    int ne, nb, nibp;
};

TestResult runFamilyTest(const string& family, int order) {
    TestResult res;
    res.family = family;
    res.order = order;
    res.loadOk = false;
    res.expandOk = false;
    res.initCondOk = false;
    res.timeMs = 0;
    res.numMatrices = 0;
    res.ne = res.nb = res.nibp = -1;

    try {
        string binFile = "data/IBPMat_" + family + ".bin";
        cout << "\n========================================\n";
        cout << "Testing family: " << family << "\n";
        cout << "  Matrix file: " << binFile << "\n";
        cout << "  Order: " << order << "\n";
        cout << "========================================\n";

        // 1. Load matrices
        auto mats = loadAllIBPMatricesBinary<FFInt>(binFile);
        res.numMatrices = mats.size();
        res.loadOk = true;
        cout << "[PASS] Loaded " << mats.size() << " matrices\n";

        if (!mats.empty()) {
            res.ne = mats[0].ne;
            res.nb = mats[0].nb;
            res.nibp = mats[0].nibp;
            cout << "  Dimensions: ne=" << res.ne 
                 << ", nb=" << res.nb 
                 << ", nibp=" << res.nibp << "\n";
        }

        // 2. Run expansion
        auto start = chrono::high_resolution_clock::now();
        auto allResults = batchProcessRecursion<FFInt>(mats, order, 2);
        auto end = chrono::high_resolution_clock::now();
        res.timeMs = chrono::duration<double, milli>(end - start).count();
        res.expandOk = true;
        cout << "[PASS] Expansion completed in " << fixed << setprecision(2) 
             << res.timeMs << " ms\n";

        // 3. Verify initial condition: C(0,0,seed={0..},j) = 1
        if (!allResults.empty() && !allResults[0].empty()) {
            bool initOk = true;
            for (size_t b = 0; b < allResults[0].size(); ++b) {
                const auto& sol = allResults[0][b];
                int ne = sol.getNe();
                int nb = sol.getNb();
                for (int j = 0; j < nb; ++j) {
                    FFInt c00 = sol(0, 0, 0, j, 0);
                    if (c00.n != 1) {
                        cout << "[FAIL] C(0,0,0," << j << ") = " << c00.n 
                             << " (expected 1)\n";
                        initOk = false;
                    }
                }
            }
            res.initCondOk = initOk;
            if (initOk) {
                cout << "[PASS] Initial condition C(0,0,*,*) = 1 verified\n";
            }
        }

        // 4. Export to MMA format
        string mmaFile = "ExpansionMMA_" + family + "_test.m";
        SeriesIO::exportAllResultsToMMA(allResults, mmaFile);
        cout << "[PASS] Exported to " << mmaFile << "\n";

    } catch (const exception& e) {
        res.errorMsg = e.what();
        cout << "[FAIL] Exception: " << e.what() << "\n";
    }

    return res;
}

void printSummary(const vector<TestResult>& results) {
    cout << "\n\n========================================\n";
    cout << "           TEST SUMMARY\n";
    cout << "========================================\n";
    
    int passed = 0, failed = 0;
    for (const auto& r : results) {
        bool ok = r.loadOk && r.expandOk && r.initCondOk;
        cout << (ok ? "[PASS]" : "[FAIL]") << " " << left << setw(15) << r.family
             << " order=" << r.order
             << " mats=" << r.numMatrices
             << " dims=[" << r.ne << "," << r.nb << "," << r.nibp << "]"
             << " time=" << fixed << setprecision(1) << r.timeMs << "ms";
        if (!ok && !r.errorMsg.empty()) {
            cout << " | Error: " << r.errorMsg;
        }
        cout << "\n";
        if (ok) passed++; else failed++;
    }
    
    cout << "----------------------------------------\n";
    cout << "Total: " << results.size() << " | Passed: " << passed 
         << " | Failed: " << failed << "\n";
    cout << "========================================\n";
}

int main(int argc, char* argv[]) {
    FFInt::set_new_prime(179424673);
    initBinomial();

    vector<pair<string, int>> testCases;
    
    if (argc >= 2) {
        // Single family mode
        string family = argv[1];
        int order = (argc > 2) ? stoi(argv[2]) : 4;
        testCases.push_back({family, order});
    } else {
        // Batch mode: test all available families
        cout << "Running batch test for all available families...\n";
        testCases = {
            {"bub", 4},
            {"bub00", 4},
            {"Bub_Numeric", 4},
        };
        
        // Add larger families if data available
        if (ifstream("data/IBPMat_DB.bin").good()) {
            testCases.push_back({"DB", 1});  // Very large: 63 mats
        }
        if (ifstream("data/IBPMat_Tri.bin").good()) {
            testCases.push_back({"Tri", 4});  // Medium: 2 mats, ne=3
        }
        if (ifstream("data/IBPMat_Box.bin").good()) {
            testCases.push_back({"Box", 4});  // Large: 7 mats, ne=4
        }
    }

    vector<TestResult> results;
    for (const auto& [family, order] : testCases) {
        // Check if binary exists
        if (!ifstream("data/IBPMat_" + family + ".bin").good()) {
            cout << "\n[SKIP] " << family << ": data/IBPMat_" << family 
                 << ".bin not found\n";
            continue;
        }
        results.push_back(runFamilyTest(family, order));
    }

    printSummary(results);

    // Return non-zero if any test failed
    for (const auto& r : results) {
        if (!r.loadOk || !r.expandOk || !r.initCondOk) {
            return 1;
        }
    }
    return 0;
}
