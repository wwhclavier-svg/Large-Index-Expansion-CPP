#ifndef IBPMATRIX_LOADER_HPP
#define IBPMATRIX_LOADER_HPP

#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

#include <fstream>
#include <string>
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

#include "json.hpp"

//namespace firefly { class FFInt; }

using json = nlohmann::json;
using namespace std;

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
        
        int current_idx = coords[k - dim_offset] - 1; // 1-based -> 0-based
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
    
    int coord_idx = coords[0] - 1; // 1-based -> 0-based
    
    // 构造剩余坐标 (传递子切片)
    // 效率优化：实际生产中应传递指针或 iterator 避免拷贝，这里为了清晰用 vector
    vector<int> next_coords(coords.begin() + 1, coords.end());
    
    if (coord_idx >= 0 && coord_idx < vec.size()) {
        insertValue(vec[coord_idx], next_coords, value, dims, dim_offset + 1);
    }
}

// ==========================================
// 4. 通用加载入口
// ==========================================
template<typename T, typename VecType>
void loadGenericMatrix(const json& jRoot, const string& key, VecType& targetVec, IBPMatrixE<T>& info) {
    if (!jRoot.contains(key)) return;
    const auto& data = jRoot.at(key);

    // 1. 读取维度
    vector<int> dims = data.at("dims").get<vector<int>>();
    
    // 更新全局信息 (以 M1 为准，或者每次都更新)
    if (key == "M1" && dims.size() >= 4) {
        info.nibp = dims[0];
        info.ne = dims[1];
        info.nb = dims[2]; // 假设最后两维是 nb x nb
    }

    // 2. 初始化目标结构
    // C++ vector 的第 0 维对应 JSON dims[0]
    resizeTensor(targetVec, dims, 0);

    // 3. 读取稀疏数据
    const auto& rowPtr = data.at("rowPointers").get<vector<int>>();
    const auto& colIdxList = data.at("columnIndices").get<vector<vector<int>>>();

    vector<T> values;
    if constexpr (std::is_same_v<T, firefly::FFInt>) {
        auto int_vals = data.at("values").get<vector<int64_t>>();
        values.assign(int_vals.begin(), int_vals.end());
    } else {
        values = data.at("values").get<vector<T>>();
    }

    // 4. 遍历行 (第一维)
    int dim0 = dims[0];
    for (int i = 0; i < dim0; ++i) {
        int start = rowPtr[i];
        int end = rowPtr[i + 1];

        for (int ptr = start; ptr < end; ++ptr) {
            const auto& coords = colIdxList[ptr];
            T val = values[ptr];

            // 调用递归插入器
            // targetVec[i] 是入口，coords 对应 dims[1]...dims[end]
            // dim_offset = 1，因为我们已经进入了第 0 维
            insertValue(targetVec[i], coords, val, dims, 1);
        }
    }
    
    //cout << "Loaded " << key << " (Depth: " << tensor_depth<decltype(targetVec)>::value << ")" << endl;
}

// ==========================================
// 5. 顶层加载函数
// ==========================================
template<typename T>
void loadIBPMatrices(const string& filename, IBPMatrixE<T>& mat) {
    ifstream f(filename);
    json j;
    f >> j;

    // 自动处理所有成员，无论深度如何
    loadGenericMatrix(j, "N1", mat.N1, mat);
    loadGenericMatrix(j, "K1", mat.K1, mat);
    loadGenericMatrix(j, "F0", mat.F0, mat);   // 2D 自动处理
    loadGenericMatrix(j, "F2", mat.F2, mat);   // 4D 自动处理
    loadGenericMatrix(j, "M1", mat.M1, mat);
    loadGenericMatrix(j, "K1s", mat.K1s, mat);
    loadGenericMatrix(j, "K2s", mat.K2s, mat);
    loadGenericMatrix(j, "F2s", mat.F2s, mat);
}

template<typename T>
vector<IBPMatrixE<T>> loadAllIBPMatrices(const string& filename) {
    ifstream f(filename);
    if (!f.is_open()) throw runtime_error("Could not open file: " + filename);

    json jRoot;
    f >> jRoot;

    if constexpr (std::is_same_v<T, firefly::FFInt>) {
    uint64_t prime = 0;
    bool prime_set = false;
    for (const auto& jItem : jRoot.is_array() ? jRoot : json::array({jRoot})) {
        // 从第一个矩阵中获取 modulus（假设每个矩阵都有 M1 算子）
        uint64_t mat_prime = 0;
        if (jItem.contains("M1")) {
            mat_prime = jItem["M1"]["modulus"].get<uint64_t>();
        } else if (jItem.contains("N1")) {
            mat_prime = jItem["N1"]["modulus"].get<uint64_t>();
        } // 可补充其他算子，但通常至少有一个
        else {
            throw runtime_error("No modulus found in matrix");
        }
        if (!prime_set) {
            firefly::FFInt::set_new_prime(mat_prime);
            prime_set = true;
            prime = mat_prime;
        } else if (prime != mat_prime) {
            throw runtime_error("Inconsistent prime modulus across matrices");
        }
    }
}

    vector<IBPMatrixE<T>> results;
    int regCount = 0;

    // 检查根节点是否为列表
    if (jRoot.is_array()) {
        for (const auto& jItem : jRoot) {
            IBPMatrixE<T> mat;
            
            // 为每个矩阵对象调用之前的通用读取逻辑
            // 注意：loadGenericMatrix 内部会根据 dims 自动 resize mat 里的各个 vector
            loadGenericMatrix(jItem, "N1", mat.N1, mat);
            loadGenericMatrix(jItem, "K1", mat.K1, mat);
            loadGenericMatrix(jItem, "F0", mat.F0, mat);
            loadGenericMatrix(jItem, "F2", mat.F2, mat);
            loadGenericMatrix(jItem, "M1", mat.M1, mat);
            loadGenericMatrix(jItem, "K1s", mat.K1s, mat);
            loadGenericMatrix(jItem, "K2s", mat.K2s, mat);
            loadGenericMatrix(jItem, "F2s", mat.F2s, mat);
            
            // 基础元数据解析（如果 JSON 顶层有这些字段）
            if (jItem.contains("nibp")) mat.nibp = jItem["nibp"];
            if (jItem.contains("ne")) mat.ne = jItem["ne"];
            if (jItem.contains("nb")) mat.nb = jItem["nb"];

            regCount++;
            cout << "reg " << regCount << "  nb = " << mat.nb << endl;

            results.push_back(std::move(mat)); // 使用 move 提高效率
        }
    } else {
        // 如果只有一个矩阵对象而非列表，也兼容处理
        IBPMatrixE<T> mat;

        loadGenericMatrix(jRoot, "N1", mat.N1, mat);
        loadGenericMatrix(jRoot, "K1", mat.K1, mat);
        loadGenericMatrix(jRoot, "F0", mat.F0, mat);
        loadGenericMatrix(jRoot, "F2", mat.F2, mat);
        loadGenericMatrix(jRoot, "M1", mat.M1, mat);
        loadGenericMatrix(jRoot, "K1s", mat.K1s, mat);
        loadGenericMatrix(jRoot, "K2s", mat.K2s, mat);
        loadGenericMatrix(jRoot, "F2s", mat.F2s, mat);

        results.push_back(std::move(mat));
    }

    return results;
}




template<typename T>
std::vector<IBPMatrixE<T>> loadAllIBPMatricesBinary(const std::string& filename) {
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



#endif