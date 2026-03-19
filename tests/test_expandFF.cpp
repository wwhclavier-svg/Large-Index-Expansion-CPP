#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "Combinatorics.hpp"      // 新增
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"


using namespace std;
using namespace firefly;       // 引入 FFInt 命名空间

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

int main() {
    // 1. 设置有限域素数（例如 1000003）
    FFInt::set_new_prime(1000003);

    // 2. 初始化组合数表（与类型无关，仍用 long long）
    initBinomial();

    // 3. 参数设置
    const int order = 4;
    const int incre = 2;
    //const string filename = "IBPMat_DBall.bin";
    //const string coeffCacheFile = "resCache_Expansion_DBall.bin";  // 缓存文件
    //const string filename = "IBPMat_DPpart_TriScale.bin";
    //const string coeffCacheFile = "resCache_Expansion_DPpart_TriScale.bin";  // 缓存文件
    const string filename = "IBPMat_DPpart_QuadriScale.bin";
    const string coeffCacheFile = "resCache_Expansion_DPpart_QuadriScale.bin";  // 缓存文件




    try {
        cout << "=== Test: Expand Coefficients over Finite Field (FFInt) ===" << endl;

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
        

        cout << "\nTest completed successfully." << endl;
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
}