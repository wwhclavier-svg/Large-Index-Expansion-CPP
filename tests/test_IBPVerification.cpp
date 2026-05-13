#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "Combinatorics.hpp"
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"
#include "IBPVerification.hpp"

using namespace std;
using namespace firefly;

int main(int argc, char* argv[]) {
    // 1. Setup finite field prime
    FFInt::set_new_prime(179424673);

    // 2. Initialize binomial table
    initBinomial();

    // 3. Get family name from command line (default bub00)
    string famname = "bub00";
    if (argc > 1) {
        famname = argv[1];
    }
    cout << "=== IBP Verification Test for Family: " << famname << " ===" << endl;

    const int order = 4;
    const int incre = 2;
    const string matrixFile = "data/IBPMat_" + famname + ".bin";
    const string coeffCacheFile = "data/resCache_Expansion_" + famname + ".bin";

    try {
        // 4. Load IBP matrix data
        auto ibpmatlist = loadAllIBPMatricesBinary<FFInt>(matrixFile);
        if (ibpmatlist.empty()) {
            cerr << "[FAIL] No matrices loaded." << endl;
            return 1;
        }
        cout << "[OK] Loaded " << ibpmatlist.size() << " matrices" << endl;
        cout << "     First matrix: ne=" << ibpmatlist[0].ne
             << " nb=" << ibpmatlist[0].nb
             << " nibp=" << ibpmatlist[0].nibp << endl;

        // 5. Try to load cached expansion results
        cout << "\n[Step 1] Loading cached expansion results..." << endl;
        vector<vector<seriesCoefficient<FFInt>>> allResults;

        bool haveCache = false;
        {
            ifstream test(coeffCacheFile, ios::binary);
            if (test) {
                test.close();
                try {
                    allResults = SeriesIO::loadAllResults<FFInt>(coeffCacheFile);
                    haveCache = true;
                    cout << "[OK] Loaded cache from " << coeffCacheFile << endl;
                } catch (const exception& e) {
                    cout << "[WARN] Cache loading failed: " << e.what() << endl;
                    haveCache = false;
                }
            }
        }

        if (!haveCache) {
            // 6. Run expansion if no cache
            cout << "\n[Step 2] Running expansion (no cache found)..." << endl;
            auto start = chrono::high_resolution_clock::now();
            allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> diff = end - start;
            cout << "[OK] Expansion completed in " << fixed << setprecision(4)
                 << diff.count() << " seconds" << endl;

            // Save cache
            SeriesIO::saveAllResults(allResults, coeffCacheFile);
            cout << "[OK] Saved cache to " << coeffCacheFile << endl;

            // Export to MMA format
            SeriesIO::exportAllResultsToMMA(allResults, "Compare-CPPResult-" + famname + ".m");
            cout << "[OK] Exported to Compare-CPPResult-" + famname + ".m" << endl;
        }

        // 7. Verify IBP equations for each matrix and solution
        cout << "\n[Step 3] Verifying IBP equations..." << endl;
        cout << "========================================" << endl;

        int totalFailures = 0;
        int matrixIdx = 0;

        for (const auto& matrixSolutions : allResults) {
            cout << "\nMatrix " << (matrixIdx + 1) << ": " << matrixSolutions.size() << " solutions" << endl;

            for (size_t solIdx = 0; solIdx < matrixSolutions.size(); ++solIdx) {
                auto& sol = const_cast<seriesCoefficient<FFInt>&>(matrixSolutions[solIdx]);
                int nimax = sol.getNimax();

                cout << "  Solution " << solIdx << " (nimax=" << nimax << "):" << endl;

                // Create verifier for this matrix
                IBPVerifier<FFInt> verifier(ibpmatlist[matrixIdx], order, incre, nimax);

                // We need to know nindep for verification
                // For now, use nimax as upper bound (actual nindep <= nimax)
                int nindep_estimate = nimax;

                int failures = verifier.verifyAll(sol, nindep_estimate);
                if (failures == 0) {
                    cout << "    [PASS] All IBP equations satisfied" << endl;
                } else {
                    cout << "    [FAIL] " << failures << " equation(s) not satisfied" << endl;
                    totalFailures += failures;
                }
            }
            matrixIdx++;
        }

        cout << "\n========================================" << endl;
        if (totalFailures == 0) {
            cout << "[PASS] All IBP equations verified successfully!" << endl;
            cout << "       The expansion coefficients satisfy the IBP equations." << endl;
        } else {
            cout << "[FAIL] " << totalFailures << " verification(s) failed." << endl;
            cout << "       Some expansion coefficients do NOT satisfy the IBP equations." << endl;
        }
        cout << "========================================" << endl;

        return totalFailures > 0 ? 1 : 0;

    } catch (const std::exception& e) {
        cerr << "[FAIL] Exception: " << e.what() << endl;
        return 1;
    }
}