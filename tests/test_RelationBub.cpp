#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

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

int main() {
    initBinomial();

    const int order = 4;
    const int incre = 2;
    const int lev = 2;
    const int deg = 1;
    const int nu_per_regime = 50;

    const string binIBPFile = "IBPMat_bub.bin";
    const string binRingFile = "RingData_bub.bin";
    const string coeffCacheFile = "ExpansionCache_bub.bin";

    try {
        cout << "=== Test Bub: Reconstruction of Linear Relations ===" << endl;

        // 1. 加载 IBP 矩阵
        cout << "Loading IBP matrices from " << binIBPFile << " ..." << endl;
        auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(binIBPFile);
        if (ibpmatlist.empty()) {
            cerr << "Error: No IBP matrices loaded." << endl;
            return 1;
        }
        int ne = ibpmatlist[0].ne;
        int nb = ibpmatlist[0].nb;
        cout << "Loaded: ne=" << ne << ", nb=" << nb << ", mod=" << FFInt::p << endl;

        // 2. 计算/加载展开系数
        vector<vector<seriesCoefficient<FFInt>>> allResults;
        try {
            allResults = SeriesIO::loadAllResults<FFInt>(coeffCacheFile);
            cout << "Loaded coefficients from cache." << endl;
        } catch (...) {
            cout << "Computing expansion coefficients..." << endl;
            allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
            SeriesIO::saveAllResults(allResults, coeffCacheFile);
        }
        cout << "Total branches: " << allResults.size() << endl;

        // 3. 加载环数据
        cout << "\nLoading ring data from " << binRingFile << " ..." << endl;
        auto ringData = AlgebraData::LoadBinary<FFInt>(binRingFile, FFInt::p);
        if (ringData.empty()) {
            cerr << "Error: No ring data loaded." << endl;
            return 1;
        }
        cout << "Loaded " << ringData.size() << " regimes." << endl;

        // 4. 构建 regimes
        auto all_regimes = RelationSolver::buildAllRegimes<FFInt>(
            allResults, ringData, ne, nu_per_regime);
        cout << "Built " << all_regimes.size() << " regimes for computation." << endl;

        // 5. 配置并求解
        RelationSolver::AdaptiveSamplingConfig config;
        config.min_nu = 3;
        config.max_nu = 50;
        config.nullity_stable_threshold = 2;
        config.lev_hint = lev;
        config.deg_hint = deg;

        cout << "\nSolving (lev=" << lev << ", deg=" << deg << ")..." << endl;
        auto start = chrono::high_resolution_clock::now();

        // 提取 sector, A_list, Ainv_list
        std::vector<std::vector<int>> sector_list;
        std::vector<std::vector<std::vector<FFInt>>> A_list, Ainv_list;
        for (const auto& ring : ringData) {
            sector_list.push_back(ring.limitSector);
            A_list.push_back(ring.A_list);
            Ainv_list.push_back(ring.Ainv_list);
        }

        auto [res, rel_coeff] = RelationSolver::reconstructReductionRelation<FFInt>(
            allResults, sector_list, A_list, Ainv_list,
            ne, lev, deg, config);

        auto end = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(end - start).count();

        cout << "Solution dimension: " << rel_coeff.getNumSolutions() << endl;
        cout << "Computation time: " << fixed << setprecision(2) << elapsed << "s" << endl;

        if (!res.Mext.empty()) {
            cout << "System size: " << res.Mext.size() << " x " << res.Mext[0].size() << endl;
        }

        cout << "\nTest completed successfully." << endl;
        return 0;
    }
    catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }
}
