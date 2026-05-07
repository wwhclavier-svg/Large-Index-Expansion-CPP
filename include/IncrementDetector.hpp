#ifndef INCREMENT_DETECTOR_HPP
#define INCREMENT_DETECTOR_HPP

#include <vector>
#include <cmath>
#include "IBPMatrixLoader_Binary.hpp"
#include "firefly/FFInt.hpp"

// Flatten [nibp][ne][nb*nb] tensor into (nibp*nb) × (ne*nb) matrix
template<typename T>
std::vector<std::vector<T>> flattenMatrix(
    const std::vector<std::vector<std::vector<T>>>& tensor,
    int nibp, int ne, int nb)
{
    int rows = nibp * nb;
    int cols = ne * nb;
    std::vector<std::vector<T>> mat(rows, std::vector<T>(cols, T(0)));
    for (int i = 0; i < nibp; ++i) {
        for (int j = 0; j < ne; ++j) {
            for (int a = 0; a < nb; ++a) {
                for (int b = 0; b < nb; ++b) {
                    T val = tensor[i][j][a * nb + b];
                    if constexpr (std::is_same_v<T, firefly::FFInt>) {
                        if (val.n != 0)
                            mat[i * nb + a][j * nb + b] = val;
                    } else {
                        if (std::abs(val) > 1e-15)
                            mat[i * nb + a][j * nb + b] = val;
                    }
                }
            }
        }
    }
    return mat;
}

// Check if all entries in N1 tensor are zero
template<typename T>
bool isN1Zero(const std::vector<std::vector<std::vector<T>>>& N1, int nibp, int ne, int nb)
{
    int nblock = nb * nb;
    for (int i = 0; i < nibp; ++i)
        for (int j = 0; j < ne; ++j)
            for (int k = 0; k < nblock; ++k)
                if constexpr (std::is_same_v<T, firefly::FFInt>) {
                    if (N1[i][j][k].n != 0) return false;
                } else {
                    if (std::abs(N1[i][j][k]) > 1e-15) return false;
                }
    return true;
}

// Matrix rank via Gaussian elimination (row echelon form)
template<typename T>
int matrixRank(std::vector<std::vector<T>> A)
{
    if (A.empty()) return 0;
    int rows = static_cast<int>(A.size());
    int cols = static_cast<int>(A[0].size());

    int rank = 0;
    for (int col = 0; col < cols && rank < rows; ++col) {
        // Find pivot
        int pivot = -1;
        for (int r = rank; r < rows; ++r) {
            bool nonZero = false;
            if constexpr (std::is_same_v<T, firefly::FFInt>) {
                nonZero = (A[r][col].n != 0);
            } else {
                nonZero = (std::abs(A[r][col]) > 1e-10);
            }
            if (nonZero) { pivot = r; break; }
        }
        if (pivot == -1) continue;

        // Swap to pivot row
        if (pivot != rank) std::swap(A[rank], A[pivot]);

        // Normalize pivot row
        T pivotVal = A[rank][col];
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            T inv = T(1) / pivotVal;
            for (int c = col; c < cols; ++c) A[rank][c] = A[rank][c] * inv;
        } else {
            for (int c = col; c < cols; ++c) A[rank][c] /= pivotVal;
        }

        // Eliminate rows below
        for (int r = rank + 1; r < rows; ++r) {
            bool nonZero = false;
            if constexpr (std::is_same_v<T, firefly::FFInt>) {
                nonZero = (A[r][col].n != 0);
            } else {
                nonZero = (std::abs(A[r][col]) > 1e-10);
            }
            if (!nonZero) continue;
            T factor = A[r][col];
            for (int c = col; c < cols; ++c) A[r][c] = A[r][c] - factor * A[rank][c];
        }

        rank++;
    }
    return rank;
}

// Check if augmented matrix [M | N] has higher rank than M alone
template<typename T>
bool isIncompatiable(const std::vector<std::vector<T>>& M,
                     const std::vector<std::vector<T>>& N)
{
    int rows = static_cast<int>(M.size());
    int mcols = static_cast<int>(M[0].size());
    int ncols = static_cast<int>(N[0].size());

    // Build augmented matrix [M | N]
    std::vector<std::vector<T>> aug(rows, std::vector<T>(mcols + ncols, T(0)));
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < mcols; ++c) aug[r][c] = M[r][c];
        for (int c = 0; c < ncols; ++c) aug[r][mcols + c] = N[r][c];
    }

    int rankM = matrixRank(M);
    int rankAug = matrixRank(aug);
    return rankAug > rankM;
}

// Detect optimal increment from IBP matrix data
// Returns: 1 (tight recursion), 2 (default), 3 (incompatible system)
template<typename T>
int detectIncrement(const IBPMatrixE<T>& ibpmat)
{
    int nibp = ibpmat.nibp;
    int ne = ibpmat.ne;
    int nb = ibpmat.nb;

    auto M1flattened = flattenMatrix(ibpmat.M1, nibp, ne, nb);
    auto N1flattened = flattenMatrix(ibpmat.N1, nibp, ne, nb);

    int rankM = matrixRank(M1flattened);
    int fullRank = ne * nb;
    bool fullRankFlag = (rankM >= fullRank);
    bool incompatFlag = isIncompatiable(M1flattened, N1flattened);
    bool N1zeroFlag = isN1Zero(ibpmat.N1, nibp, ne, nb);

    int incre = 2; // default
    if (incompatFlag) {
        incre = 3;
    } else if (fullRankFlag && N1zeroFlag) {
        incre = 1;
    }

    std::cout << "  IncrementDetector: M1_rank=" << rankM << "/" << fullRank
              << " full=" << (fullRankFlag ? "Y" : "N")
              << " incompat=" << (incompatFlag ? "Y" : "N")
              << " N1zero=" << (N1zeroFlag ? "Y" : "N")
              << " → incre=" << incre << std::endl;

    return incre;
}

#endif // INCREMENT_DETECTOR_HPP
