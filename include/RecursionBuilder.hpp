#ifndef RECURSION_BUILDER_HPP
#define RECURSION_BUILDER_HPP

// ============================================================
// RecursionBuilder — build recursion matrices from FTable + RegionData
//
// Mirrors MMA's RecursionCoefficientMatrix + recursionMatrixCompanion:
//   1. Compute polynomial entries from FTable + VarRule/FractionRule
//   2. Lift each polynomial to nb×nb matrix via MonomialBasisMatrix
//
// Formulas (v=all-1 basis):
//   F2w[i,j] = f2[i,j] * B_j              (B_j from VarRule/FractionRule)
//   K1[i]    = Σ_j F2w[j,i] + F1[i] * A_i (A_i from VarRule)
//   K2[i]    = Σ_j F2w[i,j]
//   K1s[i]   = K1[i] if sector[i]==1 else 0
//   K2s[i]   = Σ_{j: sector[j]==1} F2w[i,j]
//   M1[i]    = K1s[i] - K2s[i]
//   N1[i]    = F2D[i] + K1[i]
// ============================================================

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

#include "IBPAnalyzer.hpp"
#include "RegionSolver.hpp"
#include "PolyArith.hpp"

namespace RecursionBuilder {

// ============================================================
// Flattened nb×nb matrix (row-major)
// ============================================================

using FlatMatrix = std::vector<int64_t>;  // length nb*nb

inline int64_t& matAt(FlatMatrix& m, int nb, int row, int col) {
    return m[row * nb + col];
}
inline int64_t matAt(const FlatMatrix& m, int nb, int row, int col) {
    return m[row * nb + col];
}

// ============================================================
// FractionRule key builder: "B"[j,i] = A_j * B_i
// ============================================================

inline std::string fractionKey(int j, int i) {
    // MMA convention: FR[j,i] = A_j * B_i, stored as "B"[j,i]
    return "\"B\"[" + std::to_string(j) + "," + std::to_string(i) + "]";
}

// ============================================================
// Look up FractionRule: "B"[j,i] → polynomial in VarIndep
// For nb=1/empty VarIndep: returns constant
// ============================================================

inline PolyArith::Polynomial getFractionRulePoly(
    int j, int i, int ne,
    const std::vector<std::string>& varIndep,
    const std::map<std::string, std::string>& fractionRule,
    const std::map<std::string, std::string>& varRule,
    int64_t modulus)
{
    std::string key = fractionKey(j, i);
    auto it = fractionRule.find(key);
    if (it != fractionRule.end()) {
        try {
            return PolyArith::parseSingularPolynomial(it->second, varIndep, modulus);
        } catch (const std::exception& e) {
            std::cerr << "[getFractionRulePoly] key=" << key << " value='" << it->second << "' varIndep=["; 
            for (auto& v: varIndep) std::cerr << v << ",";
            std::cerr << "] modulus=" << modulus << " err=" << e.what() << std::endl;
            throw;
        }
    }
    // For i==j: A_i * B_i = 1 in the quotient ring
    if (j == i) {
        PolyArith::Polynomial p;
        p.push_back({std::vector<int>(varIndep.size(), 0), int64_t(1)});
        return p;
    }
    // Fallback: try A_j * B_i via VarRule. B_i = 1/A_i.
    // For nb=1: both are constants.
    auto bit = varRule.find("B" + std::to_string(i));
    auto ait = varRule.find("A" + std::to_string(j));
    if (bit != varRule.end() && ait != varRule.end()) {
        // Parse both as constants (empty VarIndep)
        PolyArith::Polynomial a_poly, b_poly;
        try {
            a_poly = PolyArith::parseSingularPolynomial(ait->second, varIndep, modulus);
            b_poly = PolyArith::parseSingularPolynomial(bit->second, varIndep, modulus);
        } catch (const std::exception& e) {
            std::cerr << "[getVarPoly fallback] A" << j << " val='" << ait->second
                      << "' B" << i << " val='" << bit->second
                      << "' varIndep=[] modulus=" << modulus
                      << " err=" << e.what() << std::endl;
            throw;
        }
        PolyArith::Polynomial result;
        for (auto& am : a_poly) {
            for (auto& bm : b_poly) {
                std::vector<int> exps(varIndep.size(), 0);
                for (size_t k = 0; k < varIndep.size(); ++k)
                    exps[k] = am.exps[k] + bm.exps[k];
                int64_t c = (am.coeff * bm.coeff) % modulus;
                if (c < 0) c += modulus;
                if (c != 0) result.push_back({exps, c});
            }
        }
        PolyArith::canonicalize(result, modulus);
        return result;
    }
    return PolyArith::Polynomial(); // zero
}

// ============================================================
// RecursionMatrices — output of RecursionBuilder
//
// Matches IBPMatrixE structure: each entry is a flattened nb×nb matrix.
// Dimensions:
//   M1, N1, K1, K1s, K2s: [nibp][ne]
//   F0:                   [nibp]
//   F2, F2s:              [nibp][ne][ne]
// ============================================================

struct RecursionMatrices {
    int nibp, ne, nb;
    std::vector<std::vector<FlatMatrix>>      M1;   // [nibp][ne]
    std::vector<std::vector<FlatMatrix>>      N1;   // [nibp][ne]
    std::vector<std::vector<FlatMatrix>>      K1;   // [nibp][ne]
    std::vector<std::vector<FlatMatrix>>      K2;   // [nibp][ne]
    std::vector<std::vector<FlatMatrix>>      K1s;  // [nibp][ne]
    std::vector<std::vector<FlatMatrix>>      K2s;  // [nibp][ne]
    std::vector<FlatMatrix>                   F0;   // [nibp]
    std::vector<std::vector<std::vector<FlatMatrix>>> F2;  // [nibp][ne][ne]
    std::vector<std::vector<std::vector<FlatMatrix>>> F2s; // [nibp][ne][ne]
};

// ============================================================
// Get polynomial representation of a variable from RegionData
//
// A_i or B_j is either:
//   - In VarDep:  polynomial(VarIndep) from VarRule
//   - In VarIndep: the generator itself (monomial {e_i})
//   - Otherwise:   constant (should not happen for solved systems)
// ============================================================

inline PolyArith::Polynomial getVarPoly(
    const std::string& varName,
    const std::vector<std::string>& varIndep,
    const std::map<std::string, std::string>& varRule,
    int64_t modulus)
{
    // Check VarRule first
    auto it = varRule.find(varName);
    if (it != varRule.end()) {
        try {
            return PolyArith::parseSingularPolynomial(it->second, varIndep, modulus);
        } catch (const std::exception& e) {
            std::cerr << "[getVarPoly] var=" << varName << " val='" << it->second
                      << "' varIndep=["; for(auto& v:varIndep)std::cerr<<v<<","; 
                      std::cerr << "] modulus=" << modulus << " err=" << e.what() << std::endl;
            throw;
        }
    }

    // Check if it's a generator (in VarIndep)
    for (size_t i = 0; i < varIndep.size(); ++i) {
        if (varIndep[i] == varName) {
            PolyArith::Polynomial p;
            std::vector<int> exps(varIndep.size(), 0);
            exps[i] = 1;
            p.push_back({exps, int64_t(1)});
            return p;
        }
    }

    // Variable not found — return 0 (shouldn't happen)
    return PolyArith::Polynomial();
}

// ============================================================
// Look up a monomial's index in the MonomialBasisIndex
// ============================================================

inline int findBasisIndex(
    const std::vector<int>& exps,
    const std::vector<std::vector<int>>& basisIndex)
{
    for (size_t i = 0; i < basisIndex.size(); ++i) {
        if (basisIndex[i].size() != exps.size()) continue;
        bool match = true;
        for (size_t j = 0; j < exps.size(); ++j) {
            if (basisIndex[i][j] != exps[j]) { match = false; break; }
        }
        if (match) return static_cast<int>(i);
    }
    return -1;
}

// ============================================================
// Lift a polynomial to an nb×nb matrix via MonomialBasisMatrix
//
// For polynomial P = Σ c_k · m_k:
//   Matrix(P) = Σ c_k · MonomialBasisMatrix[idx(m_k)]
// ============================================================

inline FlatMatrix liftPolynomial(
    const PolyArith::Polynomial& poly,
    const std::vector<std::vector<std::vector<int64_t>>>& monomialBasisMatrix,
    const std::vector<std::vector<int>>& monomialBasisIndex,
    int nb,
    int64_t modulus)
{
    FlatMatrix result(nb * nb, 0);

    if (nb == 0) return result;

    // If VarIndep is empty: poly is a constant, result = c * I
    if (monomialBasisMatrix.empty()) {
        if (poly.empty()) return result; // zero
        // poly should have exactly one monomial (constant)
        int64_t c = poly[0].coeff % modulus;
        if (c < 0) c += modulus;
        for (int i = 0; i < nb; ++i)
            matAt(result, nb, i, i) = c;
        return result;
    }

    for (const auto& mon : poly) {
        int64_t c = mon.coeff % modulus;
        if (c < 0) c += modulus;
        if (c == 0) continue;

        // Pad exponent vector to match basis index dimension.
        // Constant terms have empty exps (size 0) but basis indices
        // have size = |VarIndep|. Need zero-padding to find match.
        std::vector<int> lookupExps = mon.exps;
        if (!monomialBasisIndex.empty() && lookupExps.size() < monomialBasisIndex[0].size()) {
            lookupExps.resize(monomialBasisIndex[0].size(), 0);
        }

        int k = findBasisIndex(lookupExps, monomialBasisIndex);
        if (k < 0 || k >= (int)monomialBasisMatrix.size()) {
            // Monomial not in basis — shouldn't happen for reduced polynomials
            continue;
        }

        const auto& matK = monomialBasisMatrix[k];
        for (int i = 0; i < nb && i < (int)matK.size(); ++i) {
            for (int j = 0; j < nb && j < (int)matK[i].size(); ++j) {
                int64_t val = (matK[i][j] * c) % modulus;
                if (val < 0) val += modulus;
                matAt(result, nb, i, j) = (matAt(result, nb, i, j) + val) % modulus;
            }
        }
    }
    return result;
}

// ============================================================
// Zero matrix helper
// ============================================================

inline FlatMatrix zeroMatrix(int nb) {
    return FlatMatrix(nb * nb, 0);
}

// ============================================================
// identityMatrix — nb×nb identity, flattened
// ============================================================

inline FlatMatrix identityMatrix(int nb) {
    FlatMatrix m(nb * nb, 0);
    for (int i = 0; i < nb; ++i)
        matAt(m, nb, i, i) = 1;
    return m;
}

// ============================================================
// addMatrices — element-wise addition mod modulus
// ============================================================

inline void addMatrices(FlatMatrix& dst, const FlatMatrix& src, int nb, int64_t modulus) {
    for (size_t i = 0; i < dst.size() && i < src.size(); ++i) {
        dst[i] = (dst[i] + src[i]) % modulus;
        if (dst[i] < 0) dst[i] += modulus;
    }
}

// ============================================================
// scaleMatrix — scalar multiplication
// ============================================================

inline FlatMatrix scaleMatrix(const FlatMatrix& m, int64_t scalar, int nb, int64_t modulus) {
    FlatMatrix result(nb * nb, 0);
    int64_t s = scalar % modulus;
    if (s < 0) s += modulus;
    for (size_t i = 0; i < m.size(); ++i) {
        result[i] = (m[i] * s) % modulus;
        if (result[i] < 0) result[i] += modulus;
    }
    return result;
}

// ============================================================
// buildRecursionMatrices — main entry point
//
// Takes FTable + RegionData for a single region, produces
// lifted nb×nb matrices for all recursion operators.
// ============================================================

inline RecursionMatrices buildRecursionMatrices(
    const IBPAnalyzer::FTable& ft,
    const RegionSolver::RegionData& reg,
    int64_t modulus)
{
    int nibp = ft.nibp;
    int ne = ft.ne;
    int nb = reg.nb;

    RecursionMatrices rm;
    rm.nibp = nibp;
    rm.ne = ne;
    rm.nb = nb;

    // Allocate result storage
    auto alloc2D = [&](std::vector<std::vector<FlatMatrix>>& dest) {
        dest.resize(nibp, std::vector<FlatMatrix>(ne, zeroMatrix(nb)));
    };
    auto alloc3D = [&](std::vector<std::vector<std::vector<FlatMatrix>>>& dest) {
        dest.resize(nibp,
            std::vector<std::vector<FlatMatrix>>(ne,
                std::vector<FlatMatrix>(ne, zeroMatrix(nb))));
    };

    alloc2D(rm.M1);
    alloc2D(rm.N1);
    alloc2D(rm.K1);
    alloc2D(rm.K2);
    alloc2D(rm.K1s);
    alloc2D(rm.K2s);
    rm.F0.resize(nibp, zeroMatrix(nb));
    alloc3D(rm.F2);
    alloc3D(rm.F2s);

    // Pre-compute poly representations of A_i (B handled via FractionRule)
    std::vector<PolyArith::Polynomial> A_poly(ne);
    for (int i = 0; i < ne; ++i) {
        A_poly[i] = getVarPoly("A" + std::to_string(i + 1),
                               reg.VarIndep, reg.VarRule, modulus);
    }

    // DEBUG: dump f2 matrix for cross-check with MMA
    if (nb > 0 && ne <= 4) {
        std::cerr << "[DEBUG F2] nibp=" << nibp << " ne=" << ne << " nb=" << nb << std::endl;
        for (int m = 0; m < nibp; ++m) {
            std::cerr << "  Eq[" << m << "] f2 (C++: f2[A][B] = coeff of B*_A):";
            for (int i = 0; i < ne; ++i)
                for (int j = 0; j < ne; ++j)
                    if (ft.f2[m][i][j] != 0)
                        std::cerr << " (A" << (i+1) << ",B" << (j+1) << "=" << ft.f2[m][i][j] << ")";
            std::cerr << std::endl;
            // MMA would have f2[[m,i,j]] = coeff of B_j * A_{mmai}
            // where mmai = m in ne==nibp case, else??
            // C++ f2[m][i][j] = coeff of B_j * A_i
            // So MMA f2[[m,i,j]] should equal C++ f2[m][i][j] for same i,j
            std::cerr << "  Eq[" << m << "] f2 for F2w (loop access f2[j][i]):";
            for (int i = 0; i < ne; ++i)
                for (int j = 0; j < ne; ++j)
                    if (ft.f2[m][j][i] != 0)
                        std::cerr << " (i=" << i << ",j=" << j << " coeff=" << ft.f2[m][j][i] << " Fr=i+" << (i+1) << "j+" << (j+1) << ")";
            std::cerr << std::endl;
        }
        std::cerr << "  FractionRule keys:";
        for (auto& [k,v] : reg.FractionRule)
            std::cerr << " " << k << "=" << v;
        std::cerr << std::endl;
    }
    // For each IBP equation
    for (int m = 0; m < nibp; ++m) {
        // --- F2w[i][j] = f2[j][i] * "B"[i,j]  (compensates MMA's transposed f2) ---
        // MMA stores f2 as f2[eq][A][B] = coeff(B_A * A_B), transposed vs C++.
        // C++ f2[i][j] = coeff(B_j * A_i). Using f2[j][i] * Fr(i+1,j+1) matches MMA:
        //   MMA F2w[i][j] = f2[i][j] * B[j,i] = coeff(B_i*A_j) * (A_i*B_j)
        //   C++ F2w[i][j] = f2[j][i] * Fr(i+1,j+1) = coeff(B_i*A_j) * (A_i*B_j)  ✓
        std::vector<std::vector<PolyArith::Polynomial>> F2w_poly(ne,
            std::vector<PolyArith::Polynomial>(ne));

        for (int i = 0; i < ne; ++i) {
            for (int j = 0; j < ne; ++j) {
                if (ft.f2[m][j][i] != 0) {
                    auto frPoly = getFractionRulePoly(i + 1, j + 1, ne,
                        reg.VarIndep, reg.FractionRule, reg.VarRule, modulus);
                    F2w_poly[i][j] = PolyArith::polyScale(
                        frPoly, ft.f2[m][j][i], modulus);
                }
            }
        }

        // --- K1[i] = Σ_j F2w[j][i] + F1[i] * A_i ---
        for (int i = 0; i < ne; ++i) {
            PolyArith::Polynomial k1_poly;

            // Σ_j F2w[j][i] (F2w already includes FractionRule)
            for (int j = 0; j < ne; ++j) {
                k1_poly = PolyArith::polyAdd(k1_poly, F2w_poly[j][i], modulus);
            }

            // + F1[i] * A_i
            if (ft.F1[m][i] != 0) {
                auto term = PolyArith::polyScale(
                    A_poly[i], ft.F1[m][i], modulus);
                k1_poly = PolyArith::polyAdd(k1_poly, term, modulus);
            }

            PolyArith::canonicalize(k1_poly, modulus);
            rm.K1[m][i] = liftPolynomial(k1_poly, reg.MonomialBasisMatrix,
                                         reg.MonomialBasisIndex, nb, modulus);
        }

        // --- K2[i] = Σ_j F2w[i][j] ---
        for (int i = 0; i < ne; ++i) {
            PolyArith::Polynomial k2_poly;
            for (int j = 0; j < ne; ++j) {
                k2_poly = PolyArith::polyAdd(k2_poly, F2w_poly[i][j], modulus);
            }
            PolyArith::canonicalize(k2_poly, modulus);
            rm.K2[m][i] = liftPolynomial(k2_poly, reg.MonomialBasisMatrix,
                                         reg.MonomialBasisIndex, nb, modulus);
        }

        // --- K1s[i] = K1[i] if sector[i]==1 else 0 ---
        for (int i = 0; i < ne; ++i) {
            if (reg.limitSector[i] == 1) {
                rm.K1s[m][i] = rm.K1[m][i];
            }
            // else: stays zero
        }

        // --- K2s[i] = Σ_{j: sector[j]==1} F2w[i][j] ---
        for (int i = 0; i < ne; ++i) {
            PolyArith::Polynomial k2s_poly;
            for (int j = 0; j < ne; ++j) {
                if (reg.limitSector[j] == 1) {
                    k2s_poly = PolyArith::polyAdd(k2s_poly, F2w_poly[i][j], modulus);
                }
            }
            PolyArith::canonicalize(k2s_poly, modulus);
            rm.K2s[m][i] = liftPolynomial(k2s_poly, reg.MonomialBasisMatrix,
                                          reg.MonomialBasisIndex, nb, modulus);
        }

        // --- M1[i] = K1s[i] - K2s[i] ---
        for (int i = 0; i < ne; ++i) {
            for (int r = 0; r < nb; ++r) {
                for (int c = 0; c < nb; ++c) {
                    int64_t val = (matAt(rm.K1s[m][i], nb, r, c)
                                   - matAt(rm.K2s[m][i], nb, r, c)) % modulus;
                    if (val < 0) val += modulus;
                    matAt(rm.M1[m][i], nb, r, c) = val;
                }
            }
        }

        // --- N1[i] = F2D[i] + K1[i] ---
        for (int i = 0; i < ne; ++i) {
            PolyArith::Polynomial n1_poly;
            // F2D[m][i] is a scalar constant
            if (ft.F2D[m][i] != 0) {
                // Constant polynomial
                std::vector<int> zeroExps(reg.VarIndep.size(), 0);
                n1_poly.push_back({zeroExps, ft.F2D[m][i]});
            }
            // Add K1[i] polynomial
            // Need to reconstruct K1 poly — we already lifted it, so let's
            // compute N1 from the already-lifted K1 instead
            // N1 = F2D*I + K1 (as matrices, since F2D is scalar)
            for (int r = 0; r < nb; ++r) {
                for (int c = 0; c < nb; ++c) {
                    int64_t val = matAt(rm.K1[m][i], nb, r, c);
                    if (r == c) {
                        val = (val + ft.F2D[m][i]) % modulus;
                        if (val < 0) val += modulus;
                    }
                    matAt(rm.N1[m][i], nb, r, c) = val;
                }
            }
        }

        // --- F0[m] ---
        {
            PolyArith::Polynomial f0_poly;
            std::vector<int> zeroExps(reg.VarIndep.size(), 0);
            f0_poly.push_back({zeroExps, ft.F0[m]});
            PolyArith::canonicalize(f0_poly, modulus);
            rm.F0[m] = liftPolynomial(f0_poly, reg.MonomialBasisMatrix,
                                      reg.MonomialBasisIndex, nb, modulus);
        }

        // --- F2[i][j] = f2[i][j] * "B"[j,i]  (MMA's F2w, NO F2D on diagonal) ---
        // F2s = sector-filtered F2: F2s[i][j] = F2[i][j] if sector[j]==1 else 0
        for (int i = 0; i < ne; ++i) {
            for (int j = 0; j < ne; ++j) {
                rm.F2[m][i][j] = liftPolynomial(F2w_poly[i][j],
                    reg.MonomialBasisMatrix, reg.MonomialBasisIndex, nb, modulus);

                // F2s: sector-filtered
                if (reg.limitSector[j] == 1) {
                    rm.F2s[m][i][j] = rm.F2[m][i][j];
                }
                // else stays zero
            }
        }
    }

    return rm;
}

// ============================================================
// Debug: print non-zero matrix stats
// ============================================================

inline void printRecursionStats(const RecursionMatrices& rm) {
    auto countNZ = [](const FlatMatrix& m) -> int {
        int nz = 0;
        for (auto v : m) if (v != 0) ++nz;
        return nz;
    };

    std::cout << "[RecursionBuilder] nibp=" << rm.nibp
              << " ne=" << rm.ne << " nb=" << rm.nb << std::endl;

    auto count2D = [&](const std::vector<std::vector<FlatMatrix>>& v) {
        int total = 0, nzTotal = 0;
        for (auto& row : v)
            for (auto& m : row) { ++total; nzTotal += countNZ(m); }
        std::cout << "  total=" << total << " non-zero=" << nzTotal << std::endl;
    };

    std::cout << "M1:"; count2D(rm.M1);
    std::cout << "N1:"; count2D(rm.N1);
    std::cout << "K1:"; count2D(rm.K1);
    std::cout << "K2:"; count2D(rm.K2);
    std::cout << "K1s:"; count2D(rm.K1s);
    std::cout << "K2s:"; count2D(rm.K2s);

    int f0nz = 0;
    for (auto& m : rm.F0) f0nz += countNZ(m);
    std::cout << "F0:  total=" << rm.F0.size() << " non-zero=" << f0nz << std::endl;

    int f2total = 0, f2nz = 0, f2stotal = 0, f2snz = 0;
    for (auto& a : rm.F2) for (auto& b : a) for (auto& m : b) {
        ++f2total; f2nz += countNZ(m);
    }
    for (auto& a : rm.F2s) for (auto& b : a) for (auto& m : b) {
        ++f2stotal; f2snz += countNZ(m);
    }
    std::cout << "F2:  total=" << f2total << " non-zero=" << f2nz << std::endl;
    std::cout << "F2s: total=" << f2stotal << " non-zero=" << f2snz << std::endl;
}

} // namespace RecursionBuilder

#endif
