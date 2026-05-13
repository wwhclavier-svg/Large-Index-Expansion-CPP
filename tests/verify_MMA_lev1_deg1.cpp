#include <iostream>
#include <vector>
#include <string>
#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "Combinatorics.hpp"
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"
#include "RingDataLoader.hpp"
#include "RelationSolver.hpp"

using namespace std;
using namespace firefly;

int main() {
    FFInt::set_new_prime(179424673);
    initBinomial();

    string family = "bub00";
    int order = 4;
    int lev = 1, deg = 1;

    // Load IBP and ring data
    string verifyDir = "verify/" + family + "/";
    string ibpFile = verifyDir + "IBPMat_" + family + ".bin";
    string coeffCacheFile = verifyDir + "resCache_Expansion_" + family + ".bin";
    string ringFile = "data/RingData_" + family + ".bin";

    auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(ibpFile);
    if (ibpmatlist.empty()) { cerr << "No IBP matrices" << endl; return 1; }

    // Load or compute expansion
    auto allResults = SeriesIO::loadAllResults<FFInt>(coeffCacheFile);
    if (allResults.empty()) {
        allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, 2);
    }

    auto ringData = AlgebraData::LoadBinary<FFInt>(ringFile, FFInt::p);
    if (ringData.empty()) { cerr << "No ring data" << endl; return 1; }

    // Theta vector from ring data
    cout << "theta = [";
    for (size_t i = 0; i < ringData[0].limitSector.size(); ++i) {
        if (i > 0) cout << ", ";
        cout << ringData[0].limitSector[i];
    }
    cout << "]" << endl;

    int k_max = allResults[0][0].getKmax();
    int nimax = allResults[0][0].getNimax();
    cout << "k_max = " << k_max << endl;
    cout << "nimax = " << nimax << endl;
    cout << "ne = " << allResults[0][0].getNe() << endl;
    cout << "nb = " << allResults[0][0].getNb() << endl;

    // Diagnostic: check how many nimax indices have non-zero coefficients
    cout << "\n--- nimax diagnostic ---" << endl;
    const auto& C = allResults[0][0];
    for (int i = 0; i <= nimax; ++i) {
        bool allZero = true;
        int nonZeroCount = 0;
        for (int k = 0; k <= min(1, k_max); ++k) {
            for (int l = 0; l <= 2*k; ++l) {
                long long nSeeds = getCapacity(2, l);
                for (long long cid = 0; cid < nSeeds; ++cid) {
                    for (int j = 0; j < 1; ++j) {
                        if (C(k, l, cid, j, i) != FFInt(0)) {
                            allZero = false;
                            nonZeroCount++;
                            if (nonZeroCount <= 2)
                                cout << "  i=" << i << " k=" << k << " l=" << l << " cid=" << cid << " val=" << C(k,l,cid,j,i).n << endl;
                        }
                    }
                }
            }
        }
        cout << "i=" << i << ": " << (allZero ? "ALL ZERO" : "non-zero entries=" + to_string(nonZeroCount)) << endl;
    }

    // Generate alphas and betas for lev=1, deg=1
    vector<vector<int>> alphas, betas;
    vector<int> temp;
    RelationSolver::generateAllIndices(2, lev, temp, alphas, false);
    temp.clear();
    RelationSolver::generateAllIndices(2, deg, temp, betas, false);

    cout << "\nAlphas (" << alphas.size() << "): ";
    for (auto& a : alphas) cout << "(" << a[0] << "," << a[1] << ") ";
    cout << "\nBetas (" << betas.size() << "): ";
    for (auto& b : betas) cout << "(" << b[0] << "," << b[1] << ") ";
    cout << endl;

    // Build RegimeData
    RelationSolver::RegimeData<FFInt> reg;
    reg.C = &allResults[0][0];
    reg.theta = ringData[0].limitSector;
    reg.A_ops = ringData[0].A_list;
    reg.A_inv_ops = ringData[0].Ainv_list;
    reg.nb = reg.C->basis_size();
    reg.prepare(lev, alphas);

    // Build evaluator
    RelationSolver::RegimeEvaluator<FFInt> eval;
    eval.init(reg, k_max, lev, deg, alphas, betas);

    // MMA b coefficients for Relation 1 and Relation 2
    // Column order: for each alpha in alphas, for each beta in betas
    // alphas: (0,0), (1,0), (0,1)
    // betas:  (0,0), (1,0), (0,1)
    // Column index = a_idx * 3 + b_idx:
    //   col 0: α=(0,0),β=(0,0)   col 1: α=(0,0),β=(1,0)   col 2: α=(0,0),β=(0,1)
    //   col 3: α=(1,0),β=(0,0)   col 4: α=(1,0),β=(1,0)   col 5: α=(1,0),β=(0,1)
    //   col 6: α=(0,1),β=(0,0)   col 7: α=(0,1),β=(1,0)   col 8: α=(0,1),β=(0,1)

    // MMA Relation 1:
    // (-13801897 + v1) j[(1,0)] + (27603797 + 179424670*v1 + 179424671*v2) j[(0,1)]
    //   + (179424672 + 179424672*v1 + 2*v2) j[(0,0)]
    vector<int64_t> b1_raw = {
        179424672, 179424672, 2,        // α=(0,0): 179424672 + 179424672*v1 + 2*v2
        -13801897, 1, 0,                 // α=(1,0): -13801897 + 1*v1 + 0*v2
        27603797, 179424670, 179424671   // α=(0,1): 27603797 + 179424670*v1 + 179424671*v2
    };

    // MMA Relation 2:
    // (-1 + v2) j[(1,0)] + (165622774 + 2*v1 + v2) j[(0,1)] + (1 + 179424672*v2) j[(0,0)]
    vector<int64_t> b2_raw = {
        1, 0, 179424672,                 // α=(0,0): 1 + 0*v1 + 179424672*v2
        -1, 0, 1,                        // α=(1,0): -1 + 0*v1 + 1*v2
        165622774, 2, 1                  // α=(0,1): 165622774 + 2*v1 + 1*v2
    };

    vector<FFInt> b1(9), b2(9);
    for (int i = 0; i < 9; ++i) {
        b1[i] = FFInt(b1_raw[i]);
        b2[i] = FFInt(b2_raw[i]);
    }

    // Diagnostic: print full M for ν=(1,0) and check nullspace
    {
        vector<FFInt> nu_diag = {FFInt(1), FFInt(0)};
        auto M = eval.evaluate(nu_diag, 0);
        int nrows = M.size(), ncols = M[0].size();
        cout << "\n=== Diagnostic: M matrix at ν=(1,0), nimax=0 ===" << endl;
        cout << "Dimensions: " << nrows << " rows × " << ncols << " cols" << endl;

        for (size_t r = 0; r < M.size(); ++r) {
            cout << "row[" << r << "] ";
            for (size_t c = 0; c < M[r].size(); ++c)
                cout << (M[r][c].n < 1000000 ? " " : "") << M[r][c].n << " ";
            cout << "\n";
        }

        // Simple Gaussian elimination to find nullspace dimension
        // Work on a copy, augmented with identity for RREF
        int R = nrows, C = ncols;
        vector<vector<FFInt>> A = M;  // copy
        vector<int> pivot_col(C, -1);
        int rank = 0;
        int rr = 0;
        for (int cc = 0; cc < C && rr < R; ++cc) {
            // Find pivot row
            int pivot_row = -1;
            for (int r = rr; r < R; ++r) {
                if (A[r][cc] != FFInt(0)) { pivot_row = r; break; }
            }
            if (pivot_row < 0) continue;
            // Swap
            if (pivot_row != rr) swap(A[rr], A[pivot_row]);
            // Normalize
            FFInt inv = FFInt(1) / A[rr][cc];
            for (int c = cc; c < C; ++c) A[rr][c] = A[rr][c] * inv;
            // Eliminate other rows
            for (int r = 0; r < R; ++r) {
                if (r == rr) continue;
                FFInt factor = A[r][cc];
                if (factor == FFInt(0)) continue;
                for (int c = cc; c < C; ++c)
                    A[r][c] = A[r][c] - factor * A[rr][c];
            }
            pivot_col[cc] = rr;
            rank++;
            rr++;
        }
        cout << "\nRank of M(ν=(1,0)): " << rank << " (nullspace dim = " << (C - rank) << ")" << endl;

        // Print free variables and nullspace basis
        vector<int> free_vars;
        for (int c = 0; c < C; ++c) if (pivot_col[c] < 0) free_vars.push_back(c);
        cout << "Free variables (indices): ";
        for (int fv : free_vars) cout << fv << " ";
        cout << endl;

        // For each free variable, extract nullspace vector from RREF
        for (size_t fi = 0; fi < free_vars.size(); ++fi) {
            int fv = free_vars[fi];
            cout << "Nullspace[" << fi << "] (free var = col " << fv << "): ";
            vector<int64_t> ns(C, 0);
            ns[fv] = 1;
            for (int c = 0; c < C; ++c) {
                if (pivot_col[c] >= 0) {
                    ns[c] = (-A[pivot_col[c]][fv]).n;  // solution: x_pivot = -RREF[pivot_row][free]
                }
            }
            for (int c = 0; c < C; ++c) cout << ns[c] << " ";
            cout << endl;

            // Verify this nullspace vector: M · ns = 0?
            bool ok = true;
            for (int r = 0; r < nrows; ++r) {
                FFInt dot(0);
                for (int c = 0; c < C; ++c)
                    dot = dot + M[r][c] * FFInt(ns[c]);
                if (dot.n != 0) { ok = false; cout << "  FAIL row " << r << " dot=" << dot.n; }
            }
            cout << "  Verify: " << (ok ? "PASS" : "FAIL") << endl;

            // Now test this nullspace vector at other ν
            if (ok) {
                vector<vector<FFInt>> other_nus = {
                    {FFInt(0), FFInt(1)}, {FFInt(1), FFInt(1)}, {FFInt(2), FFInt(1)}
                };
                for (auto& nu2 : other_nus) {
                    auto M2 = eval.evaluate(nu2, 0);
                    bool ok2 = true;
                    for (int r = 0; r < (int)M2.size(); ++r) {
                        FFInt dot(0);
                        for (int c = 0; c < C; ++c)
                            dot = dot + M2[r][c] * FFInt(ns[c]);
                        if (dot.n != 0) { ok2 = false; break; }
                    }
                    cout << "    ν=(" << nu2[0].n << "," << nu2[1].n << "): " << (ok2 ? "PASS" : "FAIL") << endl;
                }
            }
        }
    }

    cout << "\n=== Verification: M(nu) · b = 0 ===" << endl;

    // Test at manually chosen nu points and also the random ones from the sampler
    vector<vector<FFInt>> test_nus = {
        {FFInt(1), FFInt(0)},
        {FFInt(0), FFInt(1)},
        {FFInt(1), FFInt(1)},
        {FFInt(2), FFInt(1)},
        {FFInt(1), FFInt(2)},
        {FFInt(3), FFInt(7)},
        {FFInt(42), FFInt(99)},
    };

    bool all_pass_1 = true, all_pass_2 = true;

    for (int nimax_idx = 0; nimax_idx <= allResults[0][0].getNimax(); ++nimax_idx) {
        cout << "\n--- nimax_idx = " << nimax_idx << " ---" << endl;

        for (size_t nui = 0; nui < test_nus.size(); ++nui) {
            auto& nu = test_nus[nui];
            auto M = eval.evaluate(nu, nimax_idx);

            // Print nu
            cout << "ν=(" << nu[0].n << "," << nu[1].n << "): ";

            // M * b1
            bool pass1 = true;
            for (size_t row = 0; row < M.size(); ++row) {
                FFInt dot(0);
                for (size_t col = 0; col < M[row].size(); ++col) {
                    dot = dot + M[row][col] * b1[col];
                }
                if (dot.n != 0) {
                    pass1 = false;
                    cout << "\n  rel1 FAIL row=" << row << " dot=" << dot.n;
                }
            }
            if (pass1) cout << "rel1=0 ";

            // M * b2
            bool pass2 = true;
            for (size_t row = 0; row < M.size(); ++row) {
                FFInt dot(0);
                for (size_t col = 0; col < M[row].size(); ++col) {
                    dot = dot + M[row][col] * b2[col];
                }
                if (dot.n != 0) {
                    pass2 = false;
                    cout << "\n  rel2 FAIL row=" << row << " dot=" << dot.n;
                }
            }
            if (pass2) cout << "rel2=0";
            if (!pass1 || !pass2) cout << " ";

            if (!pass1) all_pass_1 = false;
            if (!pass2) all_pass_2 = false;
        }
    }

    cout << "\n\n=== Summary ===" << endl;
    cout << "Relation 1: " << (all_pass_1 ? "PASS (all zeros)" : "FAIL") << endl;
    cout << "Relation 2: " << (all_pass_2 ? "PASS (all zeros)" : "FAIL") << endl;

    return 0;
}
