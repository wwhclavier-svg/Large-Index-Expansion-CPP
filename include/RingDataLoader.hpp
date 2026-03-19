#ifndef RING_DATA_LOADER_HPP
#define RING_DATA_LOADER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <type_traits>   // 添加
#include <cstdint>        // 添加
#include "firefly/FFInt.hpp"

namespace AlgebraData {

// 定义核心数据结构，方便外部引用
template <typename T>
struct RingCell {
    std::vector<int> limitSector;
    std::vector<std::vector<T>> A_list;    // [ne][nb*nb]
    std::vector<std::vector<T>> Ainv_list; // [ne][nb*nb]
    int matDim; // nb
};

// 自由函数声明
template <typename T>
std::vector<RingCell<T>> LoadBinary(const std::string& filename);

template <typename T>
std::vector<RingCell<T>> LoadBinary(const std::string& filename, uint64_t modulus = 0);

inline int Idx(int row, int col, int nb) { return row * nb + col; }

template<typename T>
std::vector<RingCell<T>> LoadBinary(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    int numCells = 0, globalNe;
    file.read(reinterpret_cast<char*>(&numCells), sizeof(numCells));
    file.read(reinterpret_cast<char*>(&globalNe), sizeof(int));
    if (!file) throw std::runtime_error("Failed to read number of cells.");

    std::vector<RingCell<T>> cells;
    cells.reserve(numCells);

    for (int i = 0; i < numCells; ++i) {
        RingCell<T> cell;
        
        // 1.读取 limitSector 长度
        int secLen = 0;
        file.read(reinterpret_cast<char*>(&secLen), sizeof(secLen));
        if (!file) throw std::runtime_error("Failed to read limitSector length.");

        // 读取 limitSector
        cell.limitSector.resize(secLen);
        file.read(reinterpret_cast<char*>(cell.limitSector.data()), secLen * sizeof(int));
        if (!file) throw std::runtime_error("Failed to read limitSector.");

        // 2.读取矩阵维度 matDim (nb)
        int nb = 0;
        file.read(reinterpret_cast<char*>(&nb), sizeof(nb));
        if (!file) throw std::runtime_error("Failed to read matDim.");
        cell.matDim = nb;

        int ne = secLen;                // limitSector 的长度即为 ne
        int matSize = nb * nb;           // 每个矩阵的元素个数

        // 准备临时缓冲区读取double数据
        std::vector<double> temp(matSize);

        // 3.读取A_list: ne行, 每行matSize个double
        cell.A_list.resize(ne);
        for (int m = 0; m < ne; ++m) {
            file.read(reinterpret_cast<char*>(temp.data()), matSize * sizeof(double));
            if (!file) throw std::runtime_error("Failed to read A_list element.");
            cell.A_list[m].resize(matSize);
            for (int k = 0; k < matSize; ++k) {
                cell.A_list[m][k] = static_cast<T>(temp[k]);
            }
        }

        // 4.读取Ainv_list
        cell.Ainv_list.resize(ne);
        for (int m = 0; m < ne; ++m) {
            file.read(reinterpret_cast<char*>(temp.data()), matSize * sizeof(double));
            if (!file) throw std::runtime_error("Failed to read Ainv_list element.");
            cell.Ainv_list[m].resize(matSize);
            for (int k = 0; k < matSize; ++k) {
                cell.Ainv_list[m][k] = static_cast<T>(temp[k]);
            }
        }

        cells.push_back(std::move(cell));
    }

    file.close();
    return cells;
}
                                                                            
template<typename T>
std::vector<RingCell<T>> LoadBinary(const std::string& filename, uint64_t modulus) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open file: " + filename);

    // 读取头部：单元数、全局 ne、模数（可选）
    int numCells = 0, globalNe = 0;
    file.read(reinterpret_cast<char*>(&numCells), sizeof(numCells));
    file.read(reinterpret_cast<char*>(&globalNe), sizeof(globalNe));
    if (!file) throw std::runtime_error("Failed to read header.");

    // 如果文件包含模数，可在此读取（假设放在头部）
    // uint64_t fileModulus = 0;
    // file.read(reinterpret_cast<char*>(&fileModulus), sizeof(fileModulus));
    // 若提供了 modulus 参数，可用它覆盖或校验

    std::vector<RingCell<T>> cells;
    cells.reserve(numCells);

    for (int i = 0; i < numCells; ++i) {
        RingCell<T> cell;

        // 读取 limitSector（始终为整数）
        int secLen = 0;
        file.read(reinterpret_cast<char*>(&secLen), sizeof(secLen));
        if (!file) throw std::runtime_error("Failed to read limitSector length.");
        if (secLen != globalNe) throw std::runtime_error("limitSector length mismatch");

        cell.limitSector.resize(secLen);
        file.read(reinterpret_cast<char*>(cell.limitSector.data()), secLen * sizeof(int));
        if (!file) throw std::runtime_error("Failed to read limitSector.");

        // 读取矩阵维度
        int nb = 0;
        file.read(reinterpret_cast<char*>(&nb), sizeof(nb));
        if (!file) throw std::runtime_error("Failed to read matDim.");
        cell.matDim = nb;

        int matSize = nb * nb;

        // 根据数据类型选择读取方式
        if constexpr (std::is_same_v<T, double>) {
            // 读取 double
            std::vector<double> temp(matSize);
            cell.A_list.resize(globalNe);
            for (int m = 0; m < globalNe; ++m) {
                file.read(reinterpret_cast<char*>(temp.data()), matSize * sizeof(double));
                if (!file) throw std::runtime_error("Failed to read A_list.");
                cell.A_list[m].assign(temp.begin(), temp.end());
            }
            cell.Ainv_list.resize(globalNe);
            for (int m = 0; m < globalNe; ++m) {
                file.read(reinterpret_cast<char*>(temp.data()), matSize * sizeof(double));
                if (!file) throw std::runtime_error("Failed to read Ainv_list.");
                cell.Ainv_list[m].assign(temp.begin(), temp.end());
            }
        }
        else if constexpr (std::is_same_v<T, firefly::FFInt>) {
            // 设置模数（如果尚未设置）
            static bool primeSet = false;
            if (modulus != 0 && !primeSet) {
                firefly::FFInt::set_new_prime(modulus);
                primeSet = true;
            }
            // 读取 int64_t
            std::vector<int64_t> temp(matSize);
            cell.A_list.resize(globalNe);
            for (int m = 0; m < globalNe; ++m) {
                file.read(reinterpret_cast<char*>(temp.data()), matSize * sizeof(int64_t));
                if (!file) throw std::runtime_error("Failed to read A_list.");
                cell.A_list[m].resize(matSize);
                for (int k = 0; k < matSize; ++k)
                    cell.A_list[m][k] = firefly::FFInt(temp[k]);
            }
            cell.Ainv_list.resize(globalNe);
            for (int m = 0; m < globalNe; ++m) {
                file.read(reinterpret_cast<char*>(temp.data()), matSize * sizeof(int64_t));
                if (!file) throw std::runtime_error("Failed to read Ainv_list.");
                cell.Ainv_list[m].resize(matSize);
                for (int k = 0; k < matSize; ++k)
                    cell.Ainv_list[m][k] = firefly::FFInt(temp[k]);
            }
        }
        else {
            static_assert(sizeof(T) == 0, "Unsupported type for RingDataLoader");
        }

        cells.push_back(std::move(cell));
    }
    return cells;
}


} // namespace AlgebraData

#endif // RING_DATA_LOADER_HPP