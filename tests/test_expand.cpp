#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

#include "IBPMatrixLoader.hpp"
#include "layerRecursion.hpp"
#include "RelationRecon.hpp"   // 包含 initBinomial() 和 BINOM

using namespace std;

int main() {
    // 1. 初始化组合数表（所有组合数函数依赖此全局表）
    initBinomial();

    // 2. 参数设置
    const int order = 5;        // 最大展开阶数
    const int incre = 2;         // 每阶对应的层数增量
    const string filename = "IBPMatAll_DB.json";

    try {
        cout << "=== Test 1: Expand Coefficients (batchProcessRecursion) ===" << endl;

        // 3. 加载 IBP 矩阵数据
        cout << "Loading IBP matrices from " << filename << " ..." << endl;
        auto ibpmatlist = loadAllIBPMatrices<double>(filename);

        if (ibpmatlist.empty()) {
            cerr << "Error: No matrices loaded." << endl;
            return 1;
        }

        // 4. 打印基本信息
        cout << "Loaded " << ibpmatlist.size() << " matrices." << endl;
        cout << "First matrix: ne = " << ibpmatlist[0].ne
             << ", nb = " << ibpmatlist[0].nb
             << ", nibp = " << ibpmatlist[0].nibp << endl;

        // 5. 执行批处理递归，生成所有展开系数
        auto start = chrono::high_resolution_clock::now();
        auto allResults = batchProcessRecursion<double>(ibpmatlist, order, incre);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end - start;

        // 6. 输出统计信息
        cout << "\n=== Results ===" << endl;
        cout << "Total time: " << fixed << setprecision(4) << diff.count() << " seconds." << endl;
        cout << "Number of solution branches per matrix:" << endl;
        for (size_t i = 0; i < allResults.size(); ++i) {
            cout << "  Matrix " << i+1 << ": " << allResults[i].size() << " solutions" << endl;
            if (!allResults[i].empty()) {
                // 可输出第一个解的维度信息
                auto& C0 = allResults[i][0];
                cout << "     first solution: basis size = " << C0.basis_size()
                     << ", total coefficients = " << C0.total_size() << endl;
            }
        }

        cout << "\nTest completed successfully." << endl;
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
}