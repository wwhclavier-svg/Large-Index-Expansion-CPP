#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include "layerRecursion.hpp"
#include "RingDataLoader.hpp"
#include "RelationRecon.hpp"
#include "IBPMatrixLoader.hpp"

using namespace std;

int main() {
    // 1. 初始化组合数表
    initBinomial();

    // 2. 参数设置
    int incre = 2;
    // int nibp = 8; // 从文件读取后会更新
    // int ne = 7;      
    // int nb = 1;      
    const int order = 6;

    // 3. 计时
    auto start = chrono::high_resolution_clock::now();

    string filename = "IBPMatAll_DB.json";

    try {
        cout << "--- 正在读取 IBPMatrixE from " + filename + "---" << endl;
        // 注意：loadAllIBPMatrices 需要在某处定义，假设在 IBPMatrixLoader.hpp 中
        auto ibpmatlist = loadAllIBPMatrices<double>(filename);
        
        if(ibpmatlist.empty()) {
            cerr << "Error: No matrices loaded." << endl;
            return 1;
        }

        cout << "Parameters: ne=" << ibpmatlist[0].ne << ", nibp=" << ibpmatlist[0].nibp << ", order=" << order << endl;
        cout << "--- 数据加载成功，开始递归计算 ---" << endl;
        
        int ne = ibpmatlist[0].ne;
        int nb = ibpmatlist[0].nb;
        int nibp = ibpmatlist[0].nibp;

        auto allResults = batchProcessRecursion(ibpmatlist, order, incre);

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end - start;
        cout << "Execution Time: " << fixed << setprecision(4) << diff.count() << " seconds." << endl;

        /*myMat = std::move(ibpmatlist[2]);

        cout << "--- IBP Layer Recursion Efficiency Test ---" << endl;
        cout << "Parameters: ne=" << myMat.ne << ", nb=" << myMat.nb << ", order=" << order << ", nibp=" << myMat.nibp << endl;

        cout << "M1:  " << endl;
        printMatrix(myMat.M1);

        cout << "N1:  " << endl;
        printMatrix(myMat.N1);

        auto result = layerRecursion<Scalar>(myMat, myMat.ne, myMat.nb, myMat.nibp, order, incre);*/

        // --- Ring Data Processing ---
        auto ringData = AlgebraData::RingDataLoader::LoadBinary("IBPRingData_DB.bin");
        int nreg = ringData.size(); 
        vector<vector<int>> sector(nreg); 
        vector<vector<vector<double>>> A_list(nreg), Ainv_list(nreg);
        
        for(int i=0; i < ringData.size(); ++i) {
            sector[i] = ringData[i].limitSector;
            A_list[i] = ringData[i].A_list;
            Ainv_list[i] = ringData[i].Ainv_list;
        }

        int nu_per_regime = 80;
        int lev=3, deg=1;

        // buildAllRegimes 假设在 RelationRecon.hpp 中定义
        auto all_regimes = buildAllRegimes(allResults, sector, A_list, Ainv_list, ne, nu_per_regime); // vector<RegimeData>
        cout << "(1) Regimes built." << endl;

        /*for(int i = 0; i < 3; ++i) {
            cout << "reg("<<i<<")  nb = "<<all_regimes[i].nb << endl;
            cout << "B: " <<endl; printMatrix(all_regimes[i].A_inv_ops);
            cout << "A: " <<endl; printMatrix(all_regimes[i].A_ops);
        }*/

        Solver<double> solver(nb, ne, order, lev, deg, allResults[0][0]);
        cout << "(2) Solver initialized." << endl;
        
        auto global_matrix = solver.buildAllEquation(all_regimes); 
        cout << "size:  " <<global_matrix.size() << "x" << global_matrix[0].size() << endl;
        cout << "(3) System built." << endl;

        auto res = solveLinearSystem(global_matrix);
        cout << "#rel = " << res.S.size() << endl;

    } catch (const std::exception& e) {
        cerr << "Error during recursion: " << e.what() << endl;
        return 1;
    }

    return 0;
}