#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

namespace utility {

// 符号函数：返回 (-1)^l 的泛型版本
template<typename T>
inline T sgn(int l) {
    return static_cast<T>((l % 2 == 0) ? 1 : -1);
}

// 标量加法
template<typename T>
inline void addTo(T& A, const T& B) {
    A += B;
}

// 向量加法（逐元素，只加到较小长度）
template<typename T>
inline void addTo(std::vector<T>& A, const std::vector<T>& B) {
    size_t n = std::min(A.size(), B.size());
    for (size_t i = 0; i < n; ++i) A[i] += B[i];
}

// 标量乘法（就地）
template<typename T>
inline void multiplyTo(T& A, const T& c) {
    A *= c;
}

// 向量乘法（就地，逐元素乘标量）
template<typename T>
inline void multiplyTo(std::vector<T>& A, const T& c) {
    for (size_t i = 0; i < A.size(); ++i) A[i] *= c;
}

// 向量乘法（返回新向量）
template<typename T>
inline std::vector<T> multiplyBy(const std::vector<T>& A, const T& c) {
    std::vector<T> res(A.size());
    for (size_t i = 0; i < A.size(); ++i) res[i] = A[i] * c;
    return res;
}

// 矩阵乘法（就地，逐元素乘标量）
template<typename T>
inline void multiplyTo(std::vector<std::vector<T>>& A, const T& c) {
    for (auto& row : A) multiplyTo(row, c);
}

// 打印一维向量 (替代原 printVec)
template<typename T>
void printVector(const std::vector<T>& v, 
                 const std::string& name = "", 
                 const std::string& separator = ", ", 
                 std::ostream& os = std::cout) {
    if (!name.empty()) os << name << ": ";
    os << "{";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1) os << separator;
    }
    os << "}";
}

// 打印二维矩阵 (vector of vectors)
template<typename T>
void printMatrix(const std::vector<std::vector<T>>& mat, 
                 const std::string& name = "", 
                 int width = 5, 
                 std::ostream& os = std::cout) {
    if (!name.empty()) os << name << ":\n";
    for (const auto& row : mat) {
        os << "[ ";
        for (const auto& val : row) {
            os << std::setw(width) << val << " ";
        }
        os << "]\n";
    }
}

// 打印三维矩阵 (vector of vector of vectors)
template<typename T>
void printMatrix3D(const std::vector<std::vector<std::vector<T>>>& mat,
                   const std::string& name = "",
                   int width = 5,
                   std::ostream& os = std::cout) {
    if (!name.empty()) os << name << ":\n";
    for (size_t i = 0; i < mat.size(); ++i) {
        os << "Layer " << i << ":\n";
        for (size_t j = 0; j < mat[i].size(); ++j) {
            os << "  Row " << j << ": [ ";
            for (const auto& val : mat[i][j]) {
                os << std::setw(width) << val << " ";
            }
            os << "]\n";
        }
    }
}

// 打印 C 风格数组（指针形式），替代原 printMat
template<typename T>
void printArray(const T* data, int rows, int cols, 
                const std::string& prefix = "", 
                int width = 5, 
                std::ostream& os = std::cout) {
    for (int i = 0; i < rows; ++i) {
        os << prefix << "[ ";
        for (int j = 0; j < cols; ++j) {
            os << std::setw(width) << data[i * cols + j] << " ";
        }
        os << "]\n";
    }
}

} // namespace utility

#endif // UTILITIES_HPP