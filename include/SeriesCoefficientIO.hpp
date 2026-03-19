#ifndef SERIES_COEFFICIENT_IO_HPP
#define SERIES_COEFFICIENT_IO_HPP

#include <fstream>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include "SeriesCoefficient.hpp"
#include "firefly/FFInt.hpp"

namespace SeriesIO {

// 魔数用于文件格式识别
constexpr char MAGIC[] = "SERCOEF";
constexpr uint32_t VERSION = 1;

// 写入单个 seriesCoefficient 对象（内部使用）
template<typename T>
void writeCoefficient(std::ofstream& out, const seriesCoefficient<T>& coeff) {
    // 写入元数据
    int32_t kmax = coeff.getKmax();
    int32_t incre = coeff.getIncre();
    int32_t ne = coeff.getNe();
    int32_t nb = coeff.getNb();
    int32_t nimax = coeff.getNimax();
    out.write(reinterpret_cast<const char*>(&kmax), sizeof(kmax));
    out.write(reinterpret_cast<const char*>(&incre), sizeof(incre));
    out.write(reinterpret_cast<const char*>(&ne), sizeof(ne));
    out.write(reinterpret_cast<const char*>(&nb), sizeof(nb));
    out.write(reinterpret_cast<const char*>(&nimax), sizeof(nimax));

    // 写入数据长度
    uint64_t dataSize = coeff.total_size();
    out.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));

    // 根据类型写入数据
    if constexpr (std::is_same_v<T, double>) {
        // double 直接写二进制
        out.write(reinterpret_cast<const char*>(coeff.raw_data().data()),
                  dataSize * sizeof(double));
    }
    else if constexpr (std::is_same_v<T, firefly::FFInt>) {
        // FFInt 转换为 uint64_t 再写入
        std::vector<uint64_t> buffer(dataSize);
        for (size_t i = 0; i < dataSize; ++i) {
            buffer[i] = coeff.raw_data()[i].n;  // 假设 FFInt 有公共成员 n
        }
        out.write(reinterpret_cast<const char*>(buffer.data()),
                  dataSize * sizeof(uint64_t));
    }
    else {
        static_assert(sizeof(T) == 0, "Unsupported type for serialization");
    }
}

// 读取单个 seriesCoefficient 对象（内部使用）
template<typename T>
seriesCoefficient<T> readCoefficient(std::ifstream& in) {
    // 读取元数据
    int32_t kmax, incre, ne, nb, nimax;
    in.read(reinterpret_cast<char*>(&kmax), sizeof(kmax));
    in.read(reinterpret_cast<char*>(&incre), sizeof(incre));
    in.read(reinterpret_cast<char*>(&ne), sizeof(ne));
    in.read(reinterpret_cast<char*>(&nb), sizeof(nb));
    in.read(reinterpret_cast<char*>(&nimax), sizeof(nimax));
    if (!in) throw std::runtime_error("Failed to read coefficient metadata");

    // 读取数据长度
    uint64_t dataSize;
    in.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));
    if (!in) throw std::runtime_error("Failed to read data size");

    // 构造 seriesCoefficient 对象（自动计算 offsets 并分配零初始化内存）
    seriesCoefficient<T> coeff(kmax, incre, ne, nb, nimax, BINOM);

    // 确保数据大小一致
    if (coeff.total_size() != dataSize)
        throw std::runtime_error("Data size mismatch while loading");

    // 根据类型读取数据
    if constexpr (std::is_same_v<T, double>) {
        in.read(reinterpret_cast<char*>(coeff.raw_data().data()),
                dataSize * sizeof(double));
    }
    else if constexpr (std::is_same_v<T, firefly::FFInt>) {
        std::vector<uint64_t> buffer(dataSize);
        in.read(reinterpret_cast<char*>(buffer.data()),
                dataSize * sizeof(uint64_t));
        for (size_t i = 0; i < dataSize; ++i) {
            coeff.raw_data()[i] = firefly::FFInt(buffer[i]);
        }
    }
    else {
        static_assert(sizeof(T) == 0, "Unsupported type for deserialization");
    }

    if (!in) throw std::runtime_error("Failed to read coefficient data");
    return coeff;
}

/**
 * 将 allResults（vector<vector<seriesCoefficient<T>>>）保存到二进制文件
 */
template<typename T>
void saveAllResults(const std::vector<std::vector<seriesCoefficient<T>>>& allResults,
                    const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);

    // 写入魔数和版本
    out.write(MAGIC, sizeof(MAGIC) - 1);  // 不包含末尾 '\0'
    uint32_t version = VERSION;
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));

    // 写入外层大小
    uint32_t nOuter = static_cast<uint32_t>(allResults.size());
    out.write(reinterpret_cast<const char*>(&nOuter), sizeof(nOuter));

    for (const auto& innerVec : allResults) {
        uint32_t nInner = static_cast<uint32_t>(innerVec.size());
        out.write(reinterpret_cast<const char*>(&nInner), sizeof(nInner));

        for (const auto& coeff : innerVec) {
            writeCoefficient(out, coeff);
        }
    }
}

/**
 * 从二进制文件加载 allResults
 */
template<typename T>
std::vector<std::vector<seriesCoefficient<T>>> loadAllResults(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file for reading: " + filename);

    // 验证魔数
    char magic[sizeof(MAGIC) - 1];
    in.read(magic, sizeof(magic));
    if (std::string(magic, sizeof(magic)) != std::string(MAGIC, sizeof(MAGIC) - 1))
        throw std::runtime_error("Invalid magic number in file");

    uint32_t version;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != VERSION)
        throw std::runtime_error("Unsupported file version");

    uint32_t nOuter;
    in.read(reinterpret_cast<char*>(&nOuter), sizeof(nOuter));
    std::vector<std::vector<seriesCoefficient<T>>> result(nOuter);

    for (uint32_t i = 0; i < nOuter; ++i) {
        uint32_t nInner;
        in.read(reinterpret_cast<char*>(&nInner), sizeof(nInner));
        result[i].reserve(nInner);
        for (uint32_t j = 0; j < nInner; ++j) {
            result[i].push_back(readCoefficient<T>(in));
        }
    }

    return result;
}

} // namespace SeriesIO

#endif // SERIES_COEFFICIENT_IO_HPP