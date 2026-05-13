#ifndef RING_BUILDER_HPP
#define RING_BUILDER_HPP

// ============================================================
// RingBuilder — compute A_i and Ainv_i matrices from RegionData
//
// Mirrors MMA's ComputeRingData:
//   1. For each A_i (i=1..ne), expand A_i via VarRule → polynomial
//      in VarIndep, lift to nb×nb matrix via MonomialBasisMatrix.
//   2. Compute Ainv_i = Inverse[A_i] mod modulus.
//
// Output can be fed to BinaryRingWriter or directly consumed by
// RelationSolver (via RingCellRaw → RingCell conversion).
// ============================================================

#include <string>
#include <vector>
#include <map>
#include <cstdint>
#include <stdexcept>
#include <iostream>

#include "RegionSolver.hpp"
#include "RecursionBuilder.hpp"
#include "PolyArith.hpp"

namespace RingBuilder {

using FlatMatrix = RecursionBuilder::FlatMatrix;

// ============================================================
// Extended Euclidean algorithm for modular inverse
// ============================================================

inline int64_t modInverse(int64_t a, int64_t modulus) {
    a = a % modulus;
    if (a < 0) a += modulus;
    if (a == 0)
        throw std::runtime_error("RingBuilder::modInverse: zero has no inverse");

    int64_t t = 0, newt = 1;
    int64_t r = modulus, newr = a;

    while (newr != 0) {
        int64_t q = r / newr;
        int64_t tmp = newt;
        newt = t - q * newt;
        t = tmp;
        tmp = newr;
        newr = r - q * newr;
        r = tmp;
    }

    if (r > 1)
        throw std::runtime_error("RingBuilder::modInverse: element not invertible");

    if (t < 0) t += modulus;
    return t;
}

// ============================================================
// Modular matrix inversion via Gauss-Jordan elimination
// ============================================================

inline FlatMatrix invertMatrix(const FlatMatrix& M, int nb, int64_t modulus) {
    if (nb == 0) return {};

    // Augmented matrix [M | I]
    std::vector<std::vector<int64_t>> aug(nb, std::vector<int64_t>(2 * nb, 0));
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            aug[i][j] = RecursionBuilder::matAt(M, nb, i, j) % modulus;
            if (aug[i][j] < 0) aug[i][j] += modulus;
        }
        aug[i][nb + i] = 1;
    }

    for (int col = 0; col < nb; ++col) {
        // Find pivot
        int pivot = -1;
        for (int row = col; row < nb; ++row) {
            if (aug[row][col] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot == -1)
            throw std::runtime_error("RingBuilder::invertMatrix: singular matrix");

        if (pivot != col)
            std::swap(aug[col], aug[pivot]);

        // Normalize pivot row
        int64_t inv = modInverse(aug[col][col], modulus);
        for (int j = 0; j < 2 * nb; ++j) {
            aug[col][j] = (aug[col][j] * inv) % modulus;
            if (aug[col][j] < 0) aug[col][j] += modulus;
        }

        // Eliminate all other rows
        for (int row = 0; row < nb; ++row) {
            if (row == col) continue;
            int64_t factor = aug[row][col];
            if (factor == 0) continue;
            for (int j = 0; j < 2 * nb; ++j) {
                aug[row][j] = (aug[row][j] - factor * aug[col][j]) % modulus;
                if (aug[row][j] < 0) aug[row][j] += modulus;
            }
        }
    }

    FlatMatrix result(nb * nb, 0);
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < nb; ++j)
            RecursionBuilder::matAt(result, nb, i, j) = aug[i][nb + j];
    return result;
}

// ============================================================
// RingMatrices — output of RingBuilder
// ============================================================

struct RingMatrices {
    std::vector<FlatMatrix> A_list;     // [ne] each nb×nb flattened, row-major
    std::vector<FlatMatrix> Ainv_list;  // [ne] each nb×nb flattened, row-major
};

// ============================================================
// computeRingMatrices — main entry point
// ============================================================

inline RingMatrices computeRingMatrices(
    const RegionSolver::RegionData& reg,
    int ne,
    int64_t modulus)
{
    int nb = reg.nb;

    RingMatrices result;
    result.A_list.resize(ne);
    result.Ainv_list.resize(ne);

    for (int i = 0; i < ne; ++i) {
        auto A_poly = RecursionBuilder::getVarPoly(
            "A" + std::to_string(i + 1),
            reg.VarIndep, reg.VarRule, modulus);

        result.A_list[i] = RecursionBuilder::liftPolynomial(
            A_poly, reg.MonomialBasisMatrix,
            reg.MonomialBasisIndex, nb, modulus);

        try {
            result.Ainv_list[i] = invertMatrix(result.A_list[i], nb, modulus);
        } catch (const std::exception& e) {
            std::cerr << "[RingBuilder] invertMatrix failed for A[" << (i+1) << "]"
                      << " nb=" << nb << " region.sector=[";
            for (size_t si = 0; si < reg.limitSector.size(); ++si) {
                if (si) std::cerr << ",";
                std::cerr << reg.limitSector[si];
            }
            std::cerr << "]" << std::endl;
            std::cerr << "  A matrix (" << nb << "x" << nb << "):" << std::endl;
            for (int r = 0; r < nb; ++r) {
                std::cerr << "  ";
                for (int c = 0; c < nb; ++c) {
                    std::cerr << RecursionBuilder::matAt(result.A_list[i], nb, r, c) << " ";
                }
                std::cerr << std::endl;
            }
            std::cerr << "  MonomialBasis: ";
            for (const auto& mb : reg.MonomialBasis) std::cerr << mb << " ";
            std::cerr << std::endl;
            std::cerr << "  VarIndep: ";
            for (const auto& v : reg.VarIndep) std::cerr << v << " ";
            std::cerr << std::endl;
            std::cerr << "  VarDep: ";
            for (const auto& v : reg.VarDep) std::cerr << v << " ";
            std::cerr << std::endl;
            std::cerr << "  VarRule:" << std::endl;
            for (const auto& [k, v] : reg.VarRule)
                std::cerr << "    " << k << " -> " << v << std::endl;
            std::cerr << "  MinPoly (GB):" << std::endl;
            for (const auto& p : reg.MinPoly)
                std::cerr << "    " << p << std::endl;
            std::cerr << "  MonomialBasisIndex:" << std::endl;
            for (size_t mi = 0; mi < reg.MonomialBasisIndex.size(); ++mi) {
                std::cerr << "    [" << mi << "] [";
                for (size_t ei = 0; ei < reg.MonomialBasisIndex[mi].size(); ++ei) {
                    if (ei) std::cerr << ",";
                    std::cerr << reg.MonomialBasisIndex[mi][ei];
                }
                std::cerr << "]  (mono=" << reg.MonomialBasis[mi] << ")" << std::endl;
            }
            std::cerr << "  MonomialBasisMatrix (first few):" << std::endl;
            for (size_t mi = 0; mi < std::min(reg.MonomialBasisMatrix.size(), size_t(3)); ++mi) {
                std::cerr << "    M[" << mi << "] (" << reg.MonomialBasis[mi] << "):" << std::endl;
                for (int r = 0; r < nb; ++r) {
                    std::cerr << "      ";
                    for (int c = 0; c < nb; ++c)
                        std::cerr << reg.MonomialBasisMatrix[mi][r][c] << " ";
                    std::cerr << std::endl;
                }
            }
            // Verify M[0] should be I
            bool m0_ok = true;
            for (int r = 0; r < nb && m0_ok; ++r)
                for (int c = 0; c < nb && m0_ok; ++c)
                    if (reg.MonomialBasisMatrix[0][r][c] != (r==c ? 1 : 0))
                        m0_ok = false;
            std::cerr << "  M[0] is identity: " << (m0_ok ? "YES" : "NO - BUG!") << std::endl;
            throw;
        }
    }

    return result;
}

// ============================================================
// Debug: print ring matrices and verify A * Ainv = I
// ============================================================

inline void printRingStats(const RingMatrices& rm, int nb) {
    int ne = static_cast<int>(rm.A_list.size());
    std::cout << "[RingBuilder] ne=" << ne << " nb=" << nb << std::endl;

    for (int i = 0; i < ne; ++i) {
        int nz = 0;
        for (auto v : rm.A_list[i]) if (v != 0) ++nz;
        std::cout << "  A[" << (i + 1) << "] non-zero=" << nz
                  << "/" << (nb * nb) << std::flush;

        // Verify A * Ainv = I
        bool ok = true;
        int64_t mod = 179424673; // consistent with the rest of the pipeline
        for (int r = 0; r < nb && ok; ++r) {
            for (int c = 0; c < nb && ok; ++c) {
                int64_t sum = 0;
                for (int k = 0; k < nb; ++k) {
                    sum = (sum + RecursionBuilder::matAt(rm.A_list[i], nb, r, k)
                                 * RecursionBuilder::matAt(rm.Ainv_list[i], nb, k, c)) % mod;
                }
                if (sum < 0) sum += mod;
                int64_t expected = (r == c) ? 1 : 0;
                if (sum != expected) ok = false;
            }
        }
        std::cout << " verify=" << (ok ? "OK" : "FAIL") << std::endl;
    }
}

} // namespace RingBuilder

#endif
