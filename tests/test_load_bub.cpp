#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursion.hpp"
#include "Combinatorics.hpp"
#include "firefly/FFInt.hpp"
#include "SeriesCoefficientIO.hpp"

using namespace std;
using namespace firefly;

int main(int argc, char* argv[]) {
    try {
        FFInt::set_new_prime(179424673);
        initBinomial();

        std::string family = "bub00";
        if (argc > 1) {
            family = argv[1];
        }

        std::string filename = "data/IBPMat_" + family + ".bin";
        std::cout << "=== Test Loading " << filename << " ===" << std::endl;

        auto ibpmatlist = loadAllIBPMatricesBinary<firefly::FFInt>(filename);
        std::cout << "Loaded " << ibpmatlist.size() << " matrices" << std::endl;
        
        if (!ibpmatlist.empty()) {
            const auto& mat = ibpmatlist[0];
            cout << "First matrix: ne=" << mat.ne 
                      << ", nb=" << mat.nb 
                      << ", nibp=" << mat.nibp << endl;
        }

        // Run expansion
        const int order = 4;
        const int incre = 2;
        cout << "\n=== Running expansion (order=" << order << ", incre=" << incre << ") ===" << endl;
        auto allResults = batchProcessRecursion<FFInt>(ibpmatlist, order, incre);
        
        cout << "Branches: " << allResults[0].size() << endl;

        // Export to MMA
        SeriesIO::exportAllResultsToMMA(allResults, "ExpansionMMA_bub00_debug.m");
        cout << "Exported to ExpansionMMA_bub00_debug.m" << endl;

        // Print raw coefficients for k=0,1,2
        if (!allResults.empty() && !allResults[0].empty()) {
            const auto& sol = allResults[0][0];
            int ne = sol.getNe();
            int nb = sol.getNb();
            cout << "\n=== Raw coefficients ===" << endl;
            for (int k = 0; k <= 2; ++k) {
                int lmax = incre * k;
                for (int l = 0; l <= lmax; ++l) {
                    long long num_cid = BINOM[l + ne - 1][ne - 1];
                    for (int cid = 0; cid < num_cid; ++cid) {
                        vector<int> seed = readIndex(cid, l, ne);
                        for (int j = 0; j < nb; ++j) {
                            FFInt val = sol(k, l, cid, j, 0);
                            cout << "C(k=" << k << ", l=" << l 
                                 << ", seed={" << seed[0] << "," << seed[1] << "}"
                                 << ", j=" << j << ") = " << val.n << endl;
                        }
                    }
                }
            }
        }
        
        cout << "\n=== Success ===" << endl;
        return 0;
    } catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }
}
