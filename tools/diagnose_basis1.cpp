// diagnose_basis1.cpp
// 诊断工具：检查 Basis1 是否真的在 M(ν) 的零空间中
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

int main(int argc, char* argv[]) {
    string family = (argc > 1) ? argv[1] : "bub00";
    string verifyDir = "verify/" + family + "/";
    
    initBinomial();
    FFInt::set_new_prime(179424673);
    
    // Load data
    string binIBPFile = verifyDir + "IBPMat_" + family + ".bin";
    string binRingFile = verifyDir + "RingData_" + family + ".bin";
    string coeffCacheFile = verifyDir + "ExpansionCache_" + family + ".bin";
    
    auto ibpMat = loadIBPMatrixBinary<FFInt>(binIBPFile);
    int ne = ibpMat.ne, nb = ibpMat.nb;
    
    // Load or compute expansion
    vector<vector<seriesCoefficient<FFInt>>> allResults;
    try {
        allResults = SeriesIO::loadAllResults<FFInt>(coeffCacheFile);
        cout << "Loaded expansion cache." << endl;
    } catch (...) {
        cout << "Running expansion..." << endl;
        allResults = batchProcessRecursion<FFInt>({ibpMat}, 4, 2);
        SeriesIO::saveAllResults(allResults, coeffCacheFile);
    }
    
    // Load ring data
    auto ringData = AlgebraData::LoadBinary<FFInt>(binRingFile, FFInt::p);
    
    // Build sector list
    auto sector_list = extractSectors(ringData);
    
    // Build A lists
    int numRegs = ringData.size();
    vector<vector<vector<FFInt>>> A_list(numRegs), Ainv_list(numRegs);
    for (int r = 0; r < numRegs; ++r) {
        for (int i = 0; i < ne; ++i) {
            vector<FFInt> A_row(nb), Ainv_row(nb);
            for (int j = 0; j < nb; ++j) {
                FFInt val = ringData[r].VarRule[to_string(i+1)][j];
                A_row[j] = val;
                Ainv_row[j] = val.inv();
            }
            A_list[r].push_back(A_row);
            Ainv_list[r].push_back(Ainv_row);
        }
    }
    
    // Generate alphas and betas for lev=2, deg=0
    vector<int> temp;
    vector<vector<int>> alphas, betas;
    RelationSolver::generateAllIndices(ne, 2, temp, alphas, false);
    temp.clear();
    RelationSolver::generateAllIndices(ne, 0, temp, betas, false);
    
    cout << "\n=== Basis1 Diagnostic ===" << endl;
    cout << "Variables: " << alphas.size() << " α × " << betas.size() << " β = " 
         << alphas.size() * betas.size() << endl;
    cout << "Alphas: ";
    for (auto& a : alphas) cout << "(" << a[0] << "," << a[1] << ") ";
    cout << endl;
    cout << "Betas: ";
    for (auto& b : betas) cout << "(" << b[0] << "," << b[1] << ") ";
    cout << endl;
    
    // Build known Basis1: 78498295*g[v1,v2] + 89712338*g[v1,v2-1] + 1*g[v1,v2-2]
    // alpha=(0,0) beta=(0,0): 78498295
    // alpha=(0,1) beta=(0,0): 89712338  
    // alpha=(0,2) beta=(0,0): 1
    vector<FFInt> basis1(alphas.size() * betas.size(), FFInt(0));
    for (size_t a = 0; a < alphas.size(); ++a) {
        if (alphas[a][0] == 0 && alphas[a][1] == 0)
            basis1[a * betas.size() + 0] = FFInt(78498295);
        if (alphas[a][0] == 0 && alphas[a][1] == 1)
            basis1[a * betas.size() + 0] = FFInt(89712338);
        if (alphas[a][0] == 0 && alphas[a][1] == 2)
            basis1[a * betas.size() + 0] = FFInt(1);
    }
    
    cout << "\nBasis1 vector:" << endl;
    for (size_t i = 0; i < basis1.size(); ++i) {
        if (basis1[i] != FFInt(0))
            cout << "  [" << i << "] α=(" << alphas[i / betas.size()][0] << "," << alphas[i / betas.size()][1] 
                 << ") β=(" << betas[i % betas.size()][0] << "," << betas[i % betas.size()][1] 
                 << "): " << basis1[i].n << endl;
    }
    
    // Build RegimeEvaluator
    vector<RelationSolver::RegimeData<FFInt>> regimes;
    for (int r = 0; r < numRegs; ++r) {
        RelationSolver::RegimeData<FFInt> reg;
        reg.C = &allResults[r][0];
        reg.theta = ringData[r].Theta;
        AppendToOperatorList(A_list[r], Ainv_list[r], reg);
        reg.nb = nb;
        regimes.push_back(reg);
    }
    
    RelationSolver::RegimeEvaluator<FFInt> evaluator;
    evaluator.init(4, 2, 0, alphas, betas);
    for (auto& reg : regimes) {
        vector<int> nimax_list = {allResults[0][0].getNimax()};
        evaluator.addRegime(reg, nimax_list);
    }
    
    // Test at multiple ν points  
    vector<pair<int,int>> test_points = {
        {1, 1}, {1, 2}, {2, 1}, {2, 3}, {5, 7}, {10, 12},
        {3, 5}, {7, 11}, {100, 200}, {1000, 500}
    };
    
    cout << "\n=== Testing M(ν) × Basis1 ===" << endl;
    int passed = 0, failed = 0;
    int total_rows = 0;
    
    for (auto& [n1, n2] : test_points) {
        vector<FFInt> nu = {FFInt(n1), FFInt(n2)};
        auto rows = evaluator.evaluateAtNu(nu);
        
        // rows is a vector of equations (each equation is a vector of coefficients)
        // Check M(ν) × Basis1 = 0
        bool all_zero = true;
        int bad_row = -1;
        for (size_t r = 0; r < rows.size(); ++r) {
            FFInt dot = FFInt(0);
            for (size_t c = 0; c < rows[r].size() && c < basis1.size(); ++c) {
                dot = dot + rows[r][c] * basis1[c];
            }
            if (dot != FFInt(0)) {
                all_zero = false;
                bad_row = r;
                break;
            }
        }
        
        if (all_zero) {
            cout << "  ν=(" << setw(4) << n1 << "," << setw(4) << n2 << "): rows=" << setw(2) << rows.size() << " ✅ M·Basis1=0" << endl;
            passed++;
        } else {
            cout << "  ν=(" << setw(4) << n1 << "," << setw(4) << n2 << "): rows=" << setw(2) << rows.size() << " ❌ M·Basis1≠0 at row " << bad_row << endl;
            // Print the dot product for the failing row
            FFInt dot = FFInt(0);
            for (size_t c = 0; c < rows[bad_row].size() && c < basis1.size(); ++c) {
                dot = dot + rows[bad_row][c] * basis1[c];
            }
            cout << "     Dot product on row " << bad_row << " = " << dot.n << endl;
            // Print the row
            cout << "     Row " << bad_row << ": ";
            for (size_t c = 0; c < rows[bad_row].size(); ++c) {
                if (rows[bad_row][c] != FFInt(0))
                    cout << "[" << c << "]=" << rows[bad_row][c].n << " ";
            }
            cout << endl;
            failed++;
        }
        total_rows += rows.size();
    }
    
    cout << "\n=== Summary ===" << endl;
    cout << "Tested " << test_points.size() << " ν points, " << total_rows << " equations total" << endl;
    cout << "PASS: " << passed << ", FAIL: " << failed << endl;
    
    return failed > 0 ? 1 : 0;
}
