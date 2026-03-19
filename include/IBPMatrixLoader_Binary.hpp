#ifndef IBPMATRIX_LOADER_HPP
#define IBPMATRIX_LOADER_HPP

#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <stdexcept>

#include <fstream>
#include <string>
#include <cstring>
#include <regex>
#include <sstream>

#include <cstdint>
#include <functional>
#include <cmath>

#include <set>
#include <vector>
#include <array>
#include <unordered_map>

#include <algorithm>
#include "firefly/FFInt.hpp"

//namespace firefly { class FFInt; }

using namespace std;


// 辅助函数：读取大端序 32 位整数
inline int32_t readBE32(std::ifstream& in) {
    uint8_t bytes[4];
    in.read(reinterpret_cast<char*>(bytes), 4);
    if (!in) throw std::runtime_error("Failed to read 32-bit integer");
    return (static_cast<int32_t>(bytes[0]) << 24) |
           (static_cast<int32_t>(bytes[1]) << 16) |
           (static_cast<int32_t>(bytes[2]) << 8)  |
           static_cast<int32_t>(bytes[3]);
}

// 辅助函数：读取大端序 64 位整数
inline int64_t readBE64(std::ifstream& in) {
    uint8_t bytes[8];
    in.read(reinterpret_cast<char*>(bytes), 8);
    if (!in) throw std::runtime_error("Failed to read 64-bit integer");
    return (static_cast<int64_t>(bytes[0]) << 56) |
           (static_cast<int64_t>(bytes[1]) << 48) |
           (static_cast<int64_t>(bytes[2]) << 40) |
           (static_cast<int64_t>(bytes[3]) << 32) |
           (static_cast<int64_t>(bytes[4]) << 24) |
           (static_cast<int64_t>(bytes[5]) << 16) |
           (static_cast<int64_t>(bytes[6]) << 8)  |
           static_cast<int64_t>(bytes[7]);
}

// 辅助函数：读取大端序双精度浮点数
inline double readBEDouble(std::ifstream& in) {
    uint8_t bytes[8];
    in.read(reinterpret_cast<char*>(bytes), 8);
    if (!in) throw std::runtime_error("Failed to read double");
    // 反转字节顺序 (大端 -> 小端)
    uint8_t swapped[8];
    for (int i = 0; i < 8; ++i) swapped[i] = bytes[7 - i];
    double result;
    std::memcpy(&result, swapped, 8);
    return result;
}


template<typename T>
struct IBPMatrixE
{
    vector<vector<vector<T>>> N1;
    vector<vector<vector<T>>> K1;
    vector<vector<T>> F0;
    vector<vector<vector<vector<T>>>> F2;
    vector<vector<vector<T>>> M1;
    vector<vector<vector<T>>> K1s;
    vector<vector<vector<T>>> K2s;
    vector<vector<vector<vector<T>>>> F2s;
    int nibp, ne, nb;
};
//  nb*nb matrix structure is stored in 1-dim vector

// ==========================================
// 1. 模板元编程助手：计算 Vector 嵌套深度
// ==========================================
template<typename T>
struct tensor_depth {
    static const int value = 0;
};

template<typename T>
struct tensor_depth<vector<T>> {
    static const int value = 1 + tensor_depth<T>::value;
};

// ==========================================
// 2. 递归调整大小 (Recursive Resizer)
// ==========================================
// 目的：根据 JSON 的 dims 数组，初始化 C++ 的嵌套 vector
// 规则：前 K-1 维一一对应，最后一维容纳 JSON 剩余维度的乘积

// 基础情况：处理最底层的 vector<T> (非 vector 的 T)
template<typename T>
void resizeTensor(vector<T>& vec, const vector<int>& dims, size_t current_dim_idx) {
    // 计算剩余维度的乘积作为线性大小
    // 例如：dims=[10, 9, 2, 2], current=2 (即对应2,2) -> size = 4
    int linear_size = 1;
    for (size_t i = current_dim_idx; i < dims.size(); ++i) {
        linear_size *= dims[i];
    }
    vec.assign(linear_size, static_cast<T>(0));
}

// 递归情况：处理嵌套 vector<vector<U>>
template<typename U>
void resizeTensor(vector<vector<U>>& vec, const vector<int>& dims, size_t current_dim_idx) {
    // 维度检查
    if (current_dim_idx >= dims.size()) return;

    // 当前层 resize
    int current_size = dims[current_dim_idx];
    vec.resize(current_size);

    // 递归对每一个元素进行 resize
    for (auto& sub_vec : vec) {
        resizeTensor(sub_vec, dims, current_dim_idx + 1);
    }
}

// ==========================================
// 3. 递归插入数值 (Recursive Inserter)
// ==========================================
// 目的：根据坐标 vector<int> coords，找到对应位置填入 value

// 基础情况：最底层 vector<T>，执行线性索引计算并赋值
template<typename T, typename V>
void insertValue(vector<T>& vec, const vector<int>& coords, V value, const vector<int>& dims, size_t dim_offset) {
    // 此时 coords 和 dims 包含了剩余维度的信息
    // 需要计算 flattened index
    // 公式: idx = c_0 * (d_1*...*d_n) + c_1 * (...) + c_n
    
    int flat_idx = 0;
    int stride = 1;
    
    // 倒序遍历剩余维度来计算 stride
    // 假设 coords 是 [c_k, c_k+1, ..., c_n]
    // 注意：MMA 导出通常是 1-based，需要 -1
    for (int k = (int)dims.size() - 1; k >= (int)dim_offset; --k) {
        int coord_val = coords[k - 1]; // coords 索引通常比 dims 少 1 (因为行指针隐含了第0维)
                                       // 但在此逻辑中，coords 是完整的剩余坐标
                                       // 让我们修正一下传入的 coords 定义：
                                       // coords 是 columnIndices 中的一行，例如 [ne_idx, nb_row, nb_col]
        
        int current_idx = coords[k - dim_offset]; 
        flat_idx += current_idx * stride;
        stride *= dims[k];
    }
    
    if (flat_idx >= 0 && flat_idx < vec.size()) {
        vec[flat_idx] = static_cast<T>(value);
    }
}

// 递归情况：剥离一个坐标，深入下一层
template<typename U, typename V>
void insertValue(vector<vector<U>>& vec, const vector<int>& coords, V value, const vector<int>& dims, size_t dim_offset) {
    // 获取当前层对应的坐标
    // coords 对应的是 columnIndices 中的内容，对应 dims[1], dims[2]...
    // dim_offset 是当前 C++ vector 在 dims 中的层级索引 (从 1 开始，因为 0 是行指针)
    
    // 这里的逻辑主要处理 C++ 的中间层级
    // 比如 F2 (4D): vec[i][j][k][l]
    // i 由 rowPointer 控制，这里的 vec 已经是 vec[i]
    // coords[0] 对应 j, coords[1] 对应 k...
    
    int coord_idx = coords[0]; 

    // 构造剩余坐标 (传递子切片)
    // 效率优化：实际生产中应传递指针或 iterator 避免拷贝，这里为了清晰用 vector
    vector<int> next_coords(coords.begin() + 1, coords.end());
    
    if (coord_idx >= 0 && coord_idx < vec.size()) {
        insertValue(vec[coord_idx], next_coords, value, dims, dim_offset + 1);
    }
}



template<typename T>
std::vector<IBPMatrixE<T>> loadAllIBPMatricesBinary0(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    // 读取魔数
    char magic[4];
    in.read(magic, 4);
    if (std::string(magic, 4) != "IBP1")
        throw std::runtime_error("Invalid magic number");

    // 读取组数
    int32_t numRegs;
    in.read(reinterpret_cast<char*>(&numRegs), sizeof(numRegs));

    std::vector<IBPMatrixE<T>> results;
    results.reserve(numRegs);

    // 辅助 lambda：读取一个 CSR 稀疏矩阵并填充到目标向量
    auto readSparse = [&](auto& targetVec) {
        uint8_t exists;
        in.read(reinterpret_cast<char*>(&exists), sizeof(exists));
        if (!exists) return;

        // 读取 dims
        int32_t dims_len;
        in.read(reinterpret_cast<char*>(&dims_len), sizeof(dims_len));
        std::vector<int32_t> dims(dims_len);
        in.read(reinterpret_cast<char*>(dims.data()), dims_len * sizeof(int32_t));

        // 读取 rowPtr
        int32_t rowPtr_len;
        in.read(reinterpret_cast<char*>(&rowPtr_len), sizeof(rowPtr_len));
        std::vector<int32_t> rowPtr(rowPtr_len);
        in.read(reinterpret_cast<char*>(rowPtr.data()), rowPtr_len * sizeof(int32_t));

        // 读取 colIdx
        int32_t colIdx_len;
        in.read(reinterpret_cast<char*>(&colIdx_len), sizeof(colIdx_len));
        std::vector<int32_t> colIdx(colIdx_len);
        in.read(reinterpret_cast<char*>(colIdx.data()), colIdx_len * sizeof(int32_t));

        // 读取 values
        int32_t values_len;
        in.read(reinterpret_cast<char*>(&values_len), sizeof(values_len));
        std::vector<T> values(values_len);

        if constexpr (std::is_same_v<T, double>) {
            std::vector<double> tmp(values_len);
            in.read(reinterpret_cast<char*>(tmp.data()), values_len * sizeof(double));
            std::copy(tmp.begin(), tmp.end(), values.begin());
        } else {
            std::vector<int64_t> tmp(values_len);
            in.read(reinterpret_cast<char*>(tmp.data()), values_len * sizeof(int64_t));
            std::copy(tmp.begin(), tmp.end(), values.begin());
        }

        // 分配目标向量空间
        resizeTensor(targetVec, dims, 0);

        // 根据 CSR 数据填充
        int coord_len = dims_len - 1;  // 每个坐标的元素个数
        for (int i = 0; i < dims[0]; ++i) {
            int start = rowPtr[i];
            int end = rowPtr[i + 1];
            for (int p = start; p < end; ++p) {
                // 提取坐标（从 1-based 转为 0-based）
                std::vector<int> coords(coord_len);
                for (int d = 0; d < coord_len; ++d) {
                    coords[d] = colIdx[p * coord_len + d] - 1;
                }
                T val = values[p];
                // 插入值
                insertValue(targetVec[i], coords, val, dims, 1);
            }
        }
    };

    // 算子顺序必须与写入顺序一致
    const std::vector<std::string> opOrder = {"M1", "N1", "K1", "F0", "F2", "K1s", "K2s", "F2s"};

    for (int r = 0; r < numRegs; ++r) {
        IBPMatrixE<T> mat;

        // 读取元数据
        int32_t nibp, ne, nb, incre;
        int64_t modulus;
        in.read(reinterpret_cast<char*>(&nibp), sizeof(nibp));
        in.read(reinterpret_cast<char*>(&ne), sizeof(ne));
        in.read(reinterpret_cast<char*>(&nb), sizeof(nb));
        in.read(reinterpret_cast<char*>(&incre), sizeof(incre));
        in.read(reinterpret_cast<char*>(&modulus), sizeof(modulus));

        mat.nibp = nibp;
        mat.ne = ne;
        mat.nb = nb;
        // incre 可存储或忽略

        // 如果是 FFInt 类型，设置素数
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            static bool prime_set = false;
            static uint64_t prime = 0;
            if (!prime_set) {
                firefly::FFInt::set_new_prime(static_cast<uint64_t>(modulus));
                prime_set = true;
                prime = static_cast<uint64_t>(modulus);
            } else if (prime != static_cast<uint64_t>(modulus)) {
                throw std::runtime_error("Inconsistent prime modulus across matrices");
            }
        }

        // 建立算子名称到成员变量的映射
        std::unordered_map<std::string, std::function<void()>> readers;
        readers["M1"]   = [&]() { readSparse(mat.M1); };
        readers["N1"]   = [&]() { readSparse(mat.N1); };
        readers["K1"]   = [&]() { readSparse(mat.K1); };
        readers["F0"]   = [&]() { readSparse(mat.F0); };
        readers["F2"]   = [&]() { readSparse(mat.F2); };
        readers["K1s"]  = [&]() { readSparse(mat.K1s); };
        readers["K2s"]  = [&]() { readSparse(mat.K2s); };
        readers["F2s"]  = [&]() { readSparse(mat.F2s); };

        // 按顺序读取
        for (const auto& op : opOrder) {
            readers[op]();
        }

        results.push_back(std::move(mat));
    }

    return results;
}


// 主函数：加载所有 IBP 矩阵
template<typename T>
std::vector<IBPMatrixE<T>> loadAllIBPMatricesBinary(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);


    // 读取魔数 (4 字节字符，无需转换)
    char magic[4];
    in.read(magic, 4);
    if (std::string(magic, 4) != "IBP1")
        throw std::runtime_error("Invalid magic number");

    // 读取组数 (大端序 32 位)
    int32_t numRegs = readBE32(in);
    if (numRegs <= 0 || numRegs > 10000)  // 合理性检查
        throw std::runtime_error("Invalid numRegs: " + std::to_string(numRegs));
    cout << "valid numRegs = " << numRegs << endl; 

    std::vector<IBPMatrixE<T>> results;
    results.reserve(numRegs);

    // 算子顺序必须与写入顺序一致
    const std::vector<std::string> opOrder = {"M1", "N1", "K1", "F0", "F2", "K1s", "K2s", "F2s"};

    for (int r = 0; r < numRegs; ++r) {
        IBPMatrixE<T> mat;

        // 读取元数据 (均为大端序)
        int32_t nibp = readBE32(in);
            cout << "valid nibp = " << nibp << endl; 
        int32_t ne   = readBE32(in);
            cout << "valid ne = " << ne << endl; 
        int32_t nb   = readBE32(in);
            cout << "valid nb = " << nb << endl; 
        int32_t incre = readBE32(in);
            cout << "valid incre = " << incre << endl; 
        int64_t modulus = readBE64(in);
            cout << "valid modolus = " << modulus << endl; 

        mat.nibp = nibp;
        mat.ne = ne;
        mat.nb = nb;
        // incre 可存储或忽略

        // 如果是 FFInt 类型，设置素数
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            static bool prime_set = false;
            static uint64_t prime = 0;
            if (!prime_set) {
                firefly::FFInt::set_new_prime(static_cast<uint64_t>(modulus));
                prime_set = true;
                prime = static_cast<uint64_t>(modulus);
            } else if (prime != static_cast<uint64_t>(modulus)) {
                throw std::runtime_error("Inconsistent prime modulus across matrices");
            }
        }
        cout << "new prime set"<<endl;

        // 定义读取稀疏矩阵的 lambda
        auto readSparse = [&](auto& targetVec, int nb) {
            // 存在标志 (1 字节，无需转换)
            uint8_t exists;
            in.read(reinterpret_cast<char*>(&exists), sizeof(exists));
            if (!exists) return;

            // 读取 dims 数组
            int32_t dims_len = readBE32(in);
            std::vector<int32_t> dims(dims_len), dims_flattened(dims_len-1);
            for (int i = 0; i < dims_len; ++i) dims[i] = readBE32(in);
            cout << "dims = [";
            for (int i = 0; i < dims_len - 1; ++i) {
                if(i < dims_len - 2) { 
                    dims_flattened[i] = dims[i];
                    cout << dims_flattened[i] << " ";
                } else {
                    dims_flattened[i] = dims[i]*dims[i+1];
                    cout << dims_flattened[i];
                }
            }
            cout << "]" << endl;

            // 读取 rowPtr 数组
            int32_t rowPtr_len = readBE32(in);
            std::vector<int32_t> rowPtr(rowPtr_len);
            for (int i = 0; i < rowPtr_len; ++i) rowPtr[i] = readBE32(in);

            // 读取 colIdx 数组
            int32_t colIdx_len = readBE32(in);
            std::vector<int32_t> colIdx(colIdx_len);
            for (int i = 0; i < colIdx_len; ++i) colIdx[i] = readBE32(in);

            // 读取 values 数组
            int32_t values_len = readBE32(in);
            cout << "values_len: "<< values_len << endl;
            std::vector<T> values(values_len);

            if constexpr (std::is_same_v<T, double>) {
                for (int i = 0; i < values_len; ++i) values[i] = readBEDouble(in);
            } else if constexpr (std::is_same_v<T, firefly::FFInt>) {
                // 对于 FFInt，文件存储的是 int64_t 整数值（模数下）
                for (int i = 0; i < values_len; ++i) {
                    int64_t val = readBE64(in);
                    values[i] = firefly::FFInt(val);  // 假设 FFInt 可从 int64_t 构造
                }
            } else {
                // 其他类型可类似处理，这里假设只有 double 和 FFInt
                static_assert(std::is_same_v<T, double> || std::is_same_v<T, firefly::FFInt>,
                              "Unsupported type");
            }

            // 分配目标向量空间
            resizeTensor(targetVec, dims_flattened, 0);

            // 根据 CSR 数据填充
            int coord_len = dims_len - 1;  // 每个坐标的分量个数
            for (int i = 0; i < dims[0]; ++i) {
                int start = rowPtr[i];
                int end = rowPtr[i + 1];
                for (int p = start; p < end; ++p) {
                    // 提取坐标（从 1-based 转为 0-based）
                    std::vector<int> coords(coord_len), coords_flattened(coord_len - 1);
                    for (int d = 0; d < coord_len; ++d) coords[d] = colIdx[p * coord_len + d] - 1;
                    for (int d = 0; d < coord_len - 1; ++d) {
                        if(d < coord_len - 2) {
                            coords_flattened[d] = coords[d];
                        } else {
                            coords_flattened[d] = coords[d]*nb + coords[d+1];
                        }
                    }
                    T val = values[p];
                    // 插入值
                    insertValue(targetVec[i], coords_flattened, val, dims_flattened, 1);
                }
            }
        };

        // 建立算子名称到成员变量的映射
        std::unordered_map<std::string, std::function<void()>> readers;
        readers["M1"]   = [&]() { readSparse(mat.M1, nb); };
        readers["N1"]   = [&]() { readSparse(mat.N1, nb); };
        readers["K1"]   = [&]() { readSparse(mat.K1, nb); };
        readers["F0"]   = [&]() { readSparse(mat.F0, nb); };
        readers["F2"]   = [&]() { readSparse(mat.F2, nb); };
        readers["K1s"]  = [&]() { readSparse(mat.K1s, nb); };
        readers["K2s"]  = [&]() { readSparse(mat.K2s, nb); };
        readers["F2s"]  = [&]() { readSparse(mat.F2s, nb); };

        // 按顺序读取
        for (const auto& op : opOrder) {
            readers[op]();
            cout << op << " <- read" << endl;
        }
        cout << "all read" << endl;
        //cout << "first: " << mat.M1[0][0][0] << endl;

        results.push_back(std::move(mat));
    }

    return results;
}

#endif