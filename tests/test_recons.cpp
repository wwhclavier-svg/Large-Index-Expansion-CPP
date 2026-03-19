#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

#include "IBPMatrixLoader.hpp"
#include "layerRecursion.hpp"
#include "RingDataLoader.hpp"
#include "RelationRecon.hpp"

using namespace std;

int main() {
    // 1. 初始化组合数表
    initBinomial();

    // 2. 参数设置
    const int order = 5;          // 展开阶数
    const int incre = 2;           // 层数增量
    const int lev = 3;             // |alpha| 上限
    const int deg = 1;             // |beta| 上限
    const int nu_per_regime = 80;  // 每个区间的采样点数

    const string jsonFile = "IBPMatAll_DB.json";
    const string binFile = "IBPRingData_DB.bin";

    try {
        cout << "=== Test 2: Reconstruction of Linear Relations ===" << endl;

        // ---- 第1步：加载 IBP 矩阵并计算展开系数 ----
        cout << "Loading IBP matrices from " << jsonFile << " ..." << endl;
        auto ibpmatlist = loadAllIBPMatrices<double>(jsonFile);
        if (ibpmatlist.empty()) {
            cerr << "Error: No IBP matrices loaded." << endl;
            return 1;
        }
        int ne = ibpmatlist[0].ne;
        int nb = ibpmatlist[0].nb;
        cout << "First matrix: ne = " << ne << ", nb = " << nb << endl;

        cout << "Running batchProcessRecursion (order=" << order << ", incre=" << incre << ") ..." << endl;
        auto start = chrono::high_resolution_clock::now();
        auto allResults = batchProcessRecursion<double>(ibpmatlist, order, incre);
        auto mid = chrono::high_resolution_clock::now();
        chrono::duration<double> recTime = mid - start;
        cout << "Recursion time: " << fixed << setprecision(4) << recTime.count() << " seconds." << endl;
        cout << "Total solution branches: " << allResults.size() << " matrices, each with "
             << allResults[0].size() << " solutions (first matrix)." << endl;

        // ---- 第2步：加载环数据 ----
        cout << "\nLoading ring data from " << binFile << " ..." << endl;
        auto ringData = AlgebraData::RingDataLoader::LoadBinary<double>(binFile);
        if (ringData.empty()) {
            cerr << "Error: No ring data loaded." << endl;
            return 1;
        }
        int nreg = ringData.size();
        cout << "Loaded " << nreg << " regimes." << endl;

        vector<vector<int>> sector(nreg);
        vector<vector<vector<double>>> A_list(nreg), Ainv_list(nreg);
        for (int i = 0; i < nreg; ++i) {
            sector[i] = ringData[i].limitSector;
            A_list[i]   = ringData[i].A_list;
            Ainv_list[i]= ringData[i].Ainv_list;
        }

        // ---- 第3步：构建所有计算区间（regimes） ----
        auto all_regimes = buildAllRegimes<double>(allResults, sector, A_list, Ainv_list, ne, nu_per_regime);
        cout << "Number of regimes built: " << all_regimes.size() << endl;

        // ---- 第4步：初始化求解器并构建全局线性系统 ----
        // 使用第一个矩阵的第一个解作为 Solver 的初始参考（nb 会根据每个 regime 自动调整）
        Solver<double> solver(nb, ne, order, lev, deg, allResults[0][0]);
        cout << "Solver initialized." << endl;

        cout << "Building global linear system (size may be large) ..." << endl;
        auto global_matrix = solver.buildAllEquation(all_regimes);
        cout << "Global system size: " << global_matrix.size() << " x " << global_matrix[0].size() << endl;

        // ---- 第5步：求解齐次线性方程组，得到线性关系 ----
        cout << "Solving homogeneous system ..." << endl;
        auto res = solveLinearSystem(global_matrix);   // 无右侧向量，自动视为零向量
        cout << "Number of independent relations found: " << res.S.size() << endl;

        // 可选：输出部分结果
        if (!res.Mext.empty()) {
            cout << "First row of solution basis (if any):" << endl;
            for (size_t j = 0; j < min<size_t>(res.Mext[0].size(), 10); ++j)
                cout << res.Mext[0][j] << " ";
            cout << "..." << endl;
        }

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> total = end - start;
        cout << "\nTotal execution time: " << total.count() << " seconds." << endl;
        cout << "Test completed successfully." << endl;
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
}