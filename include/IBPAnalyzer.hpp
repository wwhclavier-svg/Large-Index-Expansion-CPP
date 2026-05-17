#ifndef IBP_ANALYZER_HPP
#define IBP_ANALYZER_HPP

// ============================================================
// IBPAnalyzer — Large-index conversion + A/B equation generation
//
// Bridges IBPEqGenerator (g-operator form) to RegionSolver
// (A/B polynomial form) and RecursionBuilder (F-tables).
//
// Key transformations:
//   1. Large-index: select leading-order in n per sector
//   2. g→A/B:       g[α] → ∏ A_i^{α_i - v_i}  (v_i ∈ {0,1})
//   3. F-tables:    extract specific A/B monomial coefficients
// ============================================================

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <stdexcept>

#include "IBPEqGenerator.hpp"
#include "PolyArith.hpp"

namespace IBPAnalyzer {

// ============================================================
// F-table output: coefficients extracted from IBP equations
//
// Naming follows MMA's IBPCoefficientMatrix output:
//   F0   — zero-monomial coefficient (all A/B exponents = 0)
//   F1   — coefficient of A_i (N1 in MMA)
//   f2   — coefficient of B_j * A_i (off-diagonal F2)
//   F2D  — v_i-derivative of the zero-monomial (diagonal F2)
//   final F2 = f2 + DiagonalMatrix[F2D]
// ============================================================

struct FTable {
    int nibp, ne;
    std::vector<int64_t>                   F0;    // [nibp]
    std::vector<std::vector<int64_t>>      F2D;   // [nibp][ne]
    std::vector<std::vector<int64_t>>      F1;    // [nibp][ne]
    std::vector<std::vector<std::vector<int64_t>>> f2; // [nibp][ne][ne]
};

// ============================================================
// g-term → A/B monomial string (for Singular)
//
// At v_i = 0: g[α] → ∏_i A_i^{α[i]}
// The coeff is the term coefficient.
// ============================================================

inline std::string gShiftToABString(const std::vector<int>& gShift,
                                     int64_t coeff, int ne, int64_t modulus) {
    if (coeff == 0) return "0";

    // Reduce coeff mod p
    int64_t c = coeff % modulus;
    if (c < 0) c += modulus;

    std::stringstream ss;

    // Check if this is a constant term (all exponents zero)
    bool allZero = true;
    for (int e : gShift) if (e != 0) { allZero = false; break; }

    if (allZero) {
        ss << c;
        return ss.str();
    }

    // Write coefficient
    if (c != 1) {
        ss << c;
    }

    // Write variable factors: A1^e1 * A2^e2 * ...
    for (int i = 0; i < ne; ++i) {
        if (gShift[i] == 0) continue;
        if (c != 1 || i > 0 || (i == 0 && gShift[i] > 0 && ss.tellp() > 0))
            ss << "*";
        else if (c == 1 && ss.tellp() == 0)
            ; // first factor with coeff=1, no prefix needed
        ss << "A" << (i + 1);
        if (gShift[i] > 1) ss << "^" << gShift[i];
    }

    // Handle coeff=-1 case: if ss is empty (shouldn't happen), put -1
    if (c == -1 || c == modulus - 1) {
        std::string s = ss.str();
        if (s.empty()) return std::to_string(modulus - 1);
        return "-" + s;
    }

    return ss.str();
}

// ============================================================
// Compute A/B monomial for a g-term at a specific v_i evaluation.
//
// g[α] at v = (v_1, ..., v_ne):
//   For each i:
//     If α_i ≥ v_i: factor = A_i^{α_i - v_i}
//     Else:          factor = B_i^{v_i - α_i}
//
// Returns exponent vectors in the combined basis:
//   indices 0..ne-1     = B_1..B_ne
//   indices ne..2*ne-1  = A_1..A_ne
// ============================================================

inline void gTermToABExps(const std::vector<int>& gShift,
                           const std::vector<int>& viVals,
                           int ne,
                           std::vector<int>& outBExps,
                           std::vector<int>& outAExps) {
    outBExps.assign(ne, 0);
    outAExps.assign(ne, 0);
    for (int i = 0; i < ne; ++i) {
        int a = gShift[i];
        int v = viVals[i];
        if (a >= v) {
            outAExps[i] = a - v;
        } else {
            outBExps[i] = v - a;
        }
    }
}

// ============================================================
// buildABEquations — generate A/B polynomial strings for Singular
//
// For the leading-order A/B system (used in Groebner basis):
//   - Select d-terms (hasD=true) and active n_i terms
//   - Evaluate g[α] at v = 0 → ∏ A_i^{α_i}
//   - Group like monomials using PolyArith, format as Singular strings
//   - Return one polynomial string per IBP equation
// ============================================================

inline std::vector<std::string> buildABEquations(
    const IBPEqGenerator::IBPEquations& ibp,
    const std::vector<int>& sector,
    int64_t modulus)
{
    int ne = ibp.ne;
    int nibp = ibp.nibp;
    std::vector<std::string> result;
    result.reserve(nibp);

    // Combined variable names: B1..Bne, A1..Ane (matching Singular ordering)
    std::vector<std::string> varNames(2 * ne);
    for (int i = 0; i < ne; ++i) {
        varNames[i]       = "B" + std::to_string(i + 1);
        varNames[ne + i]  = "A" + std::to_string(i + 1);
    }

    // v = all-1: standard IBP shift convention for A/B extraction
    std::vector<int> viVals(ne, 1);

    for (int m = 0; m < nibp; ++m) {
        const auto& eq = ibp.equations[m];
        PolyArith::Polynomial poly;

        for (const auto& term : eq.terms) {
            // d-terms (nIdx=0) are skipped: MMA's expRegSolve2 uses
            // sectorLimitIBP to restrict to the leading order,
            // after sectorLimitIBP, Coefficient[..., "n"] filters them out.
            if (term.nIdx == 0)
                continue;
            // Skip inactive n_i (sector element is 0)
            if (term.nIdx > 0 && sector[term.nIdx - 1] == 0)
                continue;
            int64_t coeff = term.coeff % modulus;
            if (coeff < 0) coeff += modulus;
            // Active n_i: coeff stays as-is (×1, matching MMA's n→0 path)

            // g[α] at v=all-1 → A/B exponents.
            // C++ gShift = absolute ν-index; standard shift s_i = 2*v_i - gShift_i
            std::vector<int> bExps, aExps;
            std::vector<int> stdShift(ne);
            for (int i = 0; i < ne; ++i)
                stdShift[i] = 2 * viVals[i] - term.gShift[i];
            gTermToABExps(stdShift, viVals, ne, bExps, aExps);

            // Combined exponent vector: [B1..Bne, A1..Ane]
            PolyArith::Monomial mon;
            mon.exps.resize(2 * ne, 0);
            for (int i = 0; i < ne; ++i) {
                mon.exps[i]      = bExps[i];
                mon.exps[ne + i] = aExps[i];
            }
            mon.coeff = coeff % modulus;
            if (mon.coeff < 0) mon.coeff += modulus;
            if (mon.coeff != 0) poly.push_back(mon);
        }

        PolyArith::canonicalize(poly, modulus);
        result.push_back(PolyArith::polyToSingularString(poly, varNames));
    }

    return result;
}

// ============================================================
// extractFTable — extract F0, F1, f2, F2D from IBP equations
//
// MMA's IBPCoefficientMatrix computes FTable values from IBP equations.
//
// Key: C++ gShift = z-exponent ∈ {0,1,2,...} (stored by mon2F).
// MMA's g-shift convention is DIFFERENT — the reflection
//   stdShift[k] = 2*v_k - gShift[k]  (at v=all-1) = 2 - gShift[k]
// converts C++ gShift to the standard IBP g-shift convention used by MMA.
//
// This is the SAME reflection used by buildABEquations (verified correct
// by RingData being byte-identical).
//
// At v=1 with reflected shift:
//   F0:    constant     → gShift=[1,...,1]  (= d-term + n_i constant parts)
//   F2D[i]: B_i          → gShift[i]=0, gShift[j≠i]=1
//   F1[i]:  A_i          → gShift = [1,0] → A_2; [0,1] → A_1
//   f2[i][j]: B_j*A_i    → gShift = [0,2] → A_2*B_1; [2,0] → A_1*B_2
// ============================================================

inline FTable extractFTable(
    const IBPEqGenerator::IBPEquations& ibp,
    int64_t modulus)
{
    int ne = ibp.ne;
    int nibp = ibp.nibp;
    FTable ft;
    ft.nibp = nibp;
    ft.ne = ne;
    ft.F0.assign(nibp, 0);
    ft.F2D.assign(nibp, std::vector<int64_t>(ne, 0));
    ft.F1.assign(nibp, std::vector<int64_t>(ne, 0));
    ft.f2.assign(nibp, std::vector<std::vector<int64_t>>(ne, std::vector<int64_t>(ne, 0)));

    for (int m = 0; m < nibp; ++m) {
        const auto& eq = ibp.equations[m];

        for (const auto& term : eq.terms) {
            int64_t c = term.coeff % modulus;
            if (c < 0) c += modulus;
            if (c == 0) continue;

            // g-shift reflection: convert C++ z-exponents to MMA standard IBP convention.
            // Same reflection formula as buildABEquations (verified correct).
            std::vector<int> stdShift(ne);
            for (int k = 0; k < ne; ++k)
                stdShift[k] = 2 - term.gShift[k];  // v_i = 1 → stdShift = 2 - gShift

            // A/B exponents at v=1 with reflected shift:
            //   aExps[k] = max(0, stdShift[k] - 1)
            //   bExps[k] = max(0, 1 - stdShift[k])
            std::vector<int> aExps(ne, 0), bExps(ne, 0);
            for (int k = 0; k < ne; ++k) {
                if (stdShift[k] >= 1) {
                    aExps[k] = stdShift[k] - 1;
                } else {
                    bExps[k] = 1 - stdShift[k];
                }
            }

            // F0: d-term constant (hasD=true, all exps = 0)
            // In MMA, only d-terms contribute to F0.
            {
                bool isConst = true;
                for (int k = 0; k < ne; ++k)
                    if (aExps[k] != 0 || bExps[k] != 0) { isConst = false; break; }
                if (isConst && term.hasD) {
                    ft.F0[m] = (ft.F0[m] + c) % modulus;
                    if (ft.F0[m] < 0) ft.F0[m] += modulus;
                }
                // F2D[i]: n_i constant term (hasD=false, nIdx=i+1, all exps = 0)
                if (isConst && !term.hasD && term.nIdx > 0 && term.nIdx <= ne) {
                    int idx = term.nIdx - 1; // 0-based
                    ft.F2D[m][idx] = (ft.F2D[m][idx] + c) % modulus;
                    if (ft.F2D[m][idx] < 0) ft.F2D[m][idx] += modulus;
                }
            }

            // F1[i]: A_i term (aExps = e_i, bExps all 0) = terms with stdShift[i]=2
            for (int i = 0; i < ne; ++i) {
                if (aExps[i] != 1) continue;
                bool ok = true;
                for (int k = 0; k < ne; ++k) {
                    if (bExps[k] != 0) { ok = false; break; }
                    if (k != i && aExps[k] != 0) { ok = false; break; }
                }
                if (ok) {
                    ft.F1[m][i] = (ft.F1[m][i] + c) % modulus;
                    if (ft.F1[m][i] < 0) ft.F1[m][i] += modulus;
                }
            }

            // f2[i][j]: B_j * A_i term (aExps = e_i, bExps = e_j)
            // = stdShift[i]=2, stdShift[j]=0 (gShift[i]=0, gShift[j]=2)
            for (int i = 0; i < ne; ++i) {
                if (aExps[i] != 1) continue;
                for (int j = 0; j < ne; ++j) {
                    if (bExps[j] != 1) continue;
                    bool ok = true;
                    for (int k = 0; k < ne; ++k) {
                        if (k != i && aExps[k] != 0) { ok = false; break; }
                        if (k != j && bExps[k] != 0) { ok = false; break; }
                    }
                    if (ok) {
                        ft.f2[m][i][j] = (ft.f2[m][i][j] + c) % modulus;
                        if (ft.f2[m][i][j] < 0) ft.f2[m][i][j] += modulus;
                    }
                }
            }
        }
    }

    return ft;
}

// ============================================================
// Debug: print FTable summary
// ============================================================

inline void printFTableSummary(const FTable& ft) {
    std::cout << "[IBPAnalyzer] FTable: nibp=" << ft.nibp << " ne=" << ft.ne << std::endl;

    // Count non-zero entries
    int nzF0 = 0, nzF2D = 0, nzF1 = 0, nzf2 = 0;
    for (int m = 0; m < ft.nibp; ++m) {
        if (ft.F0[m] != 0) ++nzF0;
        for (int i = 0; i < ft.ne; ++i) {
            if (ft.F2D[m][i] != 0) ++nzF2D;
            if (ft.F1[m][i] != 0) ++nzF1;
            for (int j = 0; j < ft.ne; ++j)
                if (ft.f2[m][i][j] != 0) ++nzf2;
        }
    }
    std::cout << "  F0 non-zero:  " << nzF0 << "/" << ft.nibp << std::endl;
    std::cout << "  F2D non-zero: " << nzF2D << "/" << (ft.nibp * ft.ne) << std::endl;
    std::cout << "  F1 non-zero:  " << nzF1 << "/" << (ft.nibp * ft.ne) << std::endl;
    std::cout << "  f2 non-zero:  " << nzf2 << "/" << (ft.nibp * ft.ne * ft.ne) << std::endl;
}

} // namespace IBPAnalyzer

#endif
