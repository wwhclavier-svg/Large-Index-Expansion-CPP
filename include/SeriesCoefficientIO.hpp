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
constexpr uint32_t VERSION = 2;

// 写入单个 seriesCoefficient 对象（内部使用）
template<typename T>
void writeCoefficient(std::ofstream& out, const seriesCoefficient<T>& coeff) {
    // 写入元数据
    int32_t kmax = coeff.getKmax();
    int32_t incre = coeff.getIncre();
    int32_t ne = coeff.getNe();
    int32_t nb = coeff.getNb();
    int32_t nimax = coeff.getNimax();
    int32_t active_i_max = coeff.getActiveIMax();
    out.write(reinterpret_cast<const char*>(&kmax), sizeof(kmax));
    out.write(reinterpret_cast<const char*>(&incre), sizeof(incre));
    out.write(reinterpret_cast<const char*>(&ne), sizeof(ne));
    out.write(reinterpret_cast<const char*>(&nb), sizeof(nb));
    out.write(reinterpret_cast<const char*>(&nimax), sizeof(nimax));
    out.write(reinterpret_cast<const char*>(&active_i_max), sizeof(active_i_max));

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
seriesCoefficient<T> readCoefficient(std::ifstream& in, uint32_t version) {
    // 读取元数据
    int32_t kmax, incre, ne, nb, nimax, active_i_max;
    in.read(reinterpret_cast<char*>(&kmax), sizeof(kmax));
    in.read(reinterpret_cast<char*>(&incre), sizeof(incre));
    in.read(reinterpret_cast<char*>(&ne), sizeof(ne));
    in.read(reinterpret_cast<char*>(&nb), sizeof(nb));
    in.read(reinterpret_cast<char*>(&nimax), sizeof(nimax));
    if (version >= 2) {
        in.read(reinterpret_cast<char*>(&active_i_max), sizeof(active_i_max));
    } else {
        active_i_max = nimax; // 旧格式：假设全部 i 都活跃
    }
    if (!in) throw std::runtime_error("Failed to read coefficient metadata");

    // 读取数据长度
    uint64_t dataSize;
    in.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));
    if (!in) throw std::runtime_error("Failed to read data size");

    // 构造 seriesCoefficient 对象（自动计算 offsets 并分配零初始化内存）
    seriesCoefficient<T> coeff(kmax, incre, ne, nb, nimax, BINOM);
    coeff.setActiveIMax(active_i_max);

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
    if (version < 1 || version > VERSION)
        throw std::runtime_error("Unsupported file version: " + std::to_string(version));

    uint32_t nOuter;
    in.read(reinterpret_cast<char*>(&nOuter), sizeof(nOuter));
    std::vector<std::vector<seriesCoefficient<T>>> result(nOuter);

    for (uint32_t i = 0; i < nOuter; ++i) {
        uint32_t nInner;
        in.read(reinterpret_cast<char*>(&nInner), sizeof(nInner));
        result[i].reserve(nInner);
        for (uint32_t j = 0; j < nInner; ++j) {
            result[i].push_back(readCoefficient<T>(in, version));
        }
    }

    return result;
}

/**
 * 将 allResults 导出为 Mathematica 可读的 .m 文件格式。
 * 输出格式：
 *   (1) 每个 seed={a1,a2,...} 处的系数转为单项式 c*v1^a1*v2^a2*...；
 *   (2) 同一 k 下不同 l（即不同 seed 长度）的单项式归并相加为多项式；
 *   (3) 解索引 i 翻到最外层，完全为零的解被剔除。
 */
template<typename T>
void exportAllResultsToMMA(const std::vector<std::vector<seriesCoefficient<T>>>& allResults,
                           const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);

    out << "(* C++ Expansion Cache Export *)\n";
    out << "$ExpansionResults = {\n";

    for (size_t mi = 0; mi < allResults.size(); ++mi) {
        out << "  { (* Matrix " << (mi + 1) << " *)\n";
        for (size_t si = 0; si < allResults[mi].size(); ++si) {
            const auto& coeff = allResults[mi][si];
            int kmax = coeff.getKmax();
            int incre = coeff.getIncre();
            int ne = coeff.getNe();
            int nb = coeff.getNb();
            int nimax = coeff.getNimax();

            out << "    <| \"Kmax\" -> " << kmax << ", \"Incre\" -> " << incre
                << ", \"NE\" -> " << ne << ", \"NB\" -> " << nb
                << ", \"Nimax\" -> " << nimax << ",\n";
            out << "       \"Solutions\" -> {\n";

            // 先收集哪些解 (i) 非全零
            std::vector<int> active_i;
            for (int i = 0; i <= nimax; ++i) {
                bool allZero = true;
                for (int k = 0; k <= kmax && allZero; ++k) {
                    int lmax = incre * k;
                    for (int l = 0; l <= lmax && allZero; ++l) {
                        long long nSeeds = getCapacity(ne, l);
                        for (long long cid = 0; cid < nSeeds && allZero; ++cid) {
                            for (int j = 0; j < nb; ++j) {
                                if (coeff(k, l, cid, j, i) != T(0)) {
                                    allZero = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (!allZero) active_i.push_back(i);
            }

            for (size_t ii = 0; ii < active_i.size(); ++ii) {
                int i = active_i[ii];
                out << "         <| \"i\" -> " << i << ", \"H\" -> {\n";
                for (int k = 0; k <= kmax; ++k) {
                    out << "           ";
                    int lmax = incre * k;
                    bool firstTerm = true;
                    for (int l = 0; l <= lmax; ++l) {
                        long long nSeeds = getCapacity(ne, l);
                        for (long long cid = 0; cid < nSeeds; ++cid) {
                            auto seed = readIndex(cid, l, ne);
                            // 检查该 seed 是否有非零分量
                            bool anyNonZero = false;
                            for (int j = 0; j < nb; ++j) {
                                if (coeff(k, l, cid, j, i) != T(0)) {
                                    anyNonZero = true;
                                    break;
                                }
                            }
                            if (!anyNonZero) continue;
                            if (!firstTerm) out << " + ";
                            firstTerm = false;
                            // 输出 nb 维向量表示扩域元素
                            if (nb == 1) {
                                if constexpr (std::is_same_v<T, firefly::FFInt>) {
                                    out << coeff(k, l, cid, 0, i).n;
                                } else {
                                    out << coeff(k, l, cid, 0, i);
                                }
                            } else {
                                out << "{";
                                for (int j = 0; j < nb; ++j) {
                                    if (j > 0) out << ", ";
                                    if constexpr (std::is_same_v<T, firefly::FFInt>) {
                                        out << coeff(k, l, cid, j, i).n;
                                    } else {
                                        out << coeff(k, l, cid, j, i);
                                    }
                                }
                                out << "}";
                            }
                            for (int v = 0; v < ne; ++v) {
                                if (seed[v] > 0) {
                                    out << "*v" << (v + 1);
                                    if (seed[v] > 1) out << "^" << seed[v];
                                }
                            }
                        }
                    }
                    if (firstTerm) out << "0";
                    out << " (* k=" << k << " *)";
                    if (k < kmax) out << ",";
                    out << "\n";
                }
                out << "         }|>";
                if (ii < active_i.size() - 1) out << ",";
                out << "\n";
            }

            out << "       }\n";
            out << "    |>";
            if (si < allResults[mi].size() - 1) out << ",";
            out << "\n";
        }
        out << "  }";
        if (mi < allResults.size() - 1) out << ",";
        out << "\n";
    }

    out << "};\n";
    out.close();
}

} // namespace SeriesIO

#endif // SERIES_COEFFICIENT_IO_HPP