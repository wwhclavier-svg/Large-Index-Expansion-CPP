#ifndef IBP_VERIFICATION_HPP
#define IBP_VERIFICATION_HPP

#include <vector>
#include <array>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

#include "Combinatorics.hpp"
#include "SeriesCoefficient.hpp"
#include "IBPMatrixLoader_Binary.hpp"
#include "LayerRecursionCore.hpp"
#include "SeriesCoefficientIO.hpp"
#include "firefly/FFInt.hpp"

using namespace std;
using namespace firefly;

// Forward declaration of Binomial table
extern long long BINOM[MAX_VAL][MAX_VAL];

/**
 * IBP Verification Test
 * ====================
 * Verifies that the expansion coefficients satisfy the IBP equations.
 *
 * The verification checks:
 *   For each (k, l, seed), the inhomogTerms Total should be ZERO
 *   (which means the IBP equation is satisfied)
 *
 * Algorithm:
 *   For each order k = 0..Kmax:
 *     For each level l = 0..min(lmax_k, lmax_{k-1}):
 *       For each seed at level l:
 *         Build inhomogTerms for this (k, l, seed)
 *         Check if Total[m][j][0..nindep] are all ZERO
 *
 * If all checks pass, the expansion is consistent with IBP equations.
 */
template<typename T>
class IBPVerifier {
private:
    const IBPMatrixE<T>& ibpmat;
    int ne, nb, nibp, nimax, incre, kmax;
    LayerRecursionCore::inhomogTerms<T> terms;

public:
    IBPVerifier(const IBPMatrixE<T>& ibpmat, int kmax, int incre, int nimax)
        : ibpmat(ibpmat), kmax(kmax), incre(incre), nimax(nimax),
          terms(kmax, incre, ibpmat.M1.size(), ibpmat.nb, nimax, ibpmat.ne) {}

    /**
     * Verify IBP equation for a specific (k, l, seed) combination
     * @param C The series coefficient container with expansion results
     * @param k Order k
     * @param l Level l
     * @param seed Seed vector (length ne)
     * @param nindep Number of independent solutions
     * @return true if IBP equation is satisfied (Total == 0)
     */
    bool verifySingle(seriesCoefficient<T> C, int k, int l,
                      const vector<int>& seed, int nindep) {
        // Determine ncurr (last non-zero seed index + 1)
        int ncurr = 0;
        for (int i = 0; i < ne; ++i) {
            if (seed[i] > 0) ncurr = i + 1;
        }

        // Build all inhomog terms for this (k, l, seed)
        terms.buildAll(ibpmat, C, k, l, const_cast<vector<int>&>(seed), nindep, ncurr, BINOM);

        // Check if Total is zero for all m, j, and all solution indices
        int nimax_p1 = nimax + 1;
        for (int m = 0; m < nibp; ++m) {
            for (int j = 0; j < nb; ++j) {
                T* total_row = terms.get_Total_row(m, j);
                for (int i = 0; i < nimax_p1; ++i) {
                    if (total_row[i] != T(0)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    int verifyAll(const seriesCoefficient<T>& C_, int nindep) {
        // Create a non-const copy for internal use (inhomogTerms requires non-const)
        seriesCoefficient<T> C = C_;
        int failures = 0;

        for (int k = 0; k <= kmax; ++k) {
            int lmax = incre * k;
            for (int l = 0; l <= lmax; ++l) {
                // Skip levels that don't exist in k-1
                if (k > 0 && l > incre * (k - 1)) continue;

                int nSeeds = (int)BINOM[l + ne - 1][ne - 1];
                for (int cid = 0; cid < nSeeds; ++cid) {
                    vector<int> seed = readIndex(cid, l, ne);

                    bool ok = verifySingle(C, k, l, seed, nindep);
                    if (!ok) {
                        failures++;
                        cout << "  [FAIL] k=" << k << " l=" << l << " seed={";
                        for (int i = 0; i < ne; ++i) {
                            if (i > 0) cout << ",";
                            cout << seed[i];
                        }
                        cout << "}" << endl;
                    }
                }
            }
        }
        return failures;
    }

    /**
     * Print detailed debug info for a specific (k, l, seed)
     */
    void debugSingle(seriesCoefficient<T> C, int k, int l,
                     const vector<int>& seed, int nindep) {
        int ncurr = 0;
        for (int i = 0; i < ne; ++i) {
            if (seed[i] > 0) ncurr = i + 1;
        }

        terms.buildAll(ibpmat, C, k, l, const_cast<vector<int>&>(seed), nindep, ncurr, BINOM);

        cout << "[DEBUG] Verification for k=" << k << " l=" << l << " seed={";
        for (int i = 0; i < ne; ++i) {
            if (i > 0) cout << ",";
            cout << seed[i];
        }
        cout << "} nindep=" << nindep << " ncurr=" << ncurr << endl;

        int nimax_p1 = nimax + 1;
        for (int m = 0; m < min(nibp, 3); ++m) {  // Limit output
            for (int j = 0; j < nb; ++j) {
                T* total_row = terms.get_Total_row(m, j);
                cout << "  Total[" << m << "][" << j << "]=[";
                for (int i = 0; i < min(nimax_p1, 5); ++i) {
                    if (i > 0) cout << ",";
                    cout << total_row[i];
                }
                if (nimax_p1 > 5) cout << ",...";
                cout << "]" << endl;
            }
        }
    }
};

#endif // IBP_VERIFICATION_HPP