#ifndef BINARY_RING_WRITER_HPP
#define BINARY_RING_WRITER_HPP

#include <fstream>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <type_traits>

#include "firefly/FFInt.hpp"
#include "RingDataLoader.hpp"

// ==========================================
// Raw ring cell data (original un-flattened format)
// ==========================================
template<typename T>
struct RingCellRaw {
    std::vector<int>   limitSector;
    int                nb = 0;
    std::vector<T>     A_flat;     // ne * nb*nb elements, row-major [e][nb*nb]
    std::vector<T>     Ainv_flat;  // ne * nb*nb elements, row-major [e][nb*nb]
};

// ==========================================
// Write RingData binary
// Format (native endian, no magic):
//   [numCells: int32]
//   [globalNe: int32]
//   Per cell:
//     [secLen: int32] [limitSector: int32[secLen]]
//     [nb: int32]
//     [A_list: ne * nb*nb × elem]
//     [Ainv_list: ne * nb*nb × elem]
//   elem = double (if modulus==0) or int64_t (for FFInt)
// ==========================================
template<typename T>
void writeRingDataBinary(const std::string& filename,
                         const std::vector<RingCellRaw<T>>& cells,
                         int globalNe) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);

    int32_t count = static_cast<int32_t>(cells.size());
    int32_t ne = static_cast<int32_t>(globalNe);
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));
    out.write(reinterpret_cast<const char*>(&ne), sizeof(ne));

    for (const auto& cell : cells) {
        // limitSector
        int32_t secLen = static_cast<int32_t>(cell.limitSector.size());
        out.write(reinterpret_cast<const char*>(&secLen), sizeof(secLen));
        out.write(reinterpret_cast<const char*>(cell.limitSector.data()),
                  secLen * sizeof(int));

        // nb
        int32_t nb = static_cast<int32_t>(cell.nb);
        out.write(reinterpret_cast<const char*>(&nb), sizeof(nb));

        // A_list and Ainv_list are stored as flat vectors
        // Write element-by-element to handle type conversion
        for (auto& v : cell.A_flat) {
            if constexpr (std::is_same_v<T, firefly::FFInt>) {
                int64_t val = static_cast<int64_t>(v.n);
                out.write(reinterpret_cast<const char*>(&val), sizeof(val));
            } else {
                double val = static_cast<double>(v);
                out.write(reinterpret_cast<const char*>(&val), sizeof(val));
            }
        }
        for (auto& v : cell.Ainv_flat) {
            if constexpr (std::is_same_v<T, firefly::FFInt>) {
                int64_t val = static_cast<int64_t>(v.n);
                out.write(reinterpret_cast<const char*>(&val), sizeof(val));
            } else {
                double val = static_cast<double>(v);
                out.write(reinterpret_cast<const char*>(&val), sizeof(val));
            }
        }
    }

    out.close();
    std::cout << "Written " << cells.size() << " ring cell(s) to " << filename
              << std::endl;
}

// ==========================================
// Read RingData into raw format
// ==========================================
template<typename T>
std::vector<RingCellRaw<T>> readRingDataRaw(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    int32_t numCells, globalNe;
    in.read(reinterpret_cast<char*>(&numCells), sizeof(numCells));
    in.read(reinterpret_cast<char*>(&globalNe), sizeof(globalNe));
    if (!in) throw std::runtime_error("Failed to read ring data header");

    std::vector<RingCellRaw<T>> cells;
    cells.reserve(numCells);

    for (int c = 0; c < numCells; ++c) {
        RingCellRaw<T> cell;

        // limitSector
        int32_t secLen;
        in.read(reinterpret_cast<char*>(&secLen), sizeof(secLen));
        cell.limitSector.resize(secLen);
        in.read(reinterpret_cast<char*>(cell.limitSector.data()),
                secLen * sizeof(int));

        // nb
        int32_t nb;
        in.read(reinterpret_cast<char*>(&nb), sizeof(nb));
        cell.nb = nb;

        int32_t matSize = globalNe * nb * nb;

        // A_list
        cell.A_flat.resize(matSize);
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            for (int i = 0; i < matSize; ++i) {
                int64_t val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                cell.A_flat[i] = firefly::FFInt(val);
            }
        } else {
            for (int i = 0; i < matSize; ++i) {
                double val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                cell.A_flat[i] = static_cast<T>(val);
            }
        }

        // Ainv_list
        cell.Ainv_flat.resize(matSize);
        if constexpr (std::is_same_v<T, firefly::FFInt>) {
            for (int i = 0; i < matSize; ++i) {
                int64_t val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                cell.Ainv_flat[i] = firefly::FFInt(val);
            }
        } else {
            for (int i = 0; i < matSize; ++i) {
                double val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                cell.Ainv_flat[i] = static_cast<T>(val);
            }
        }

        cells.push_back(std::move(cell));
    }

    return cells;
}

// ==========================================
// Round-trip verification for RingData
// ==========================================
template<typename T>
bool verifyRingRoundTrip(const std::string& filename) {
    auto raw = readRingDataRaw<T>(filename);
    std::string tmpFile = filename + ".roundtrip_test";
    writeRingDataBinary<T>(tmpFile, raw, raw.empty() ? 0 : static_cast<int>(raw[0].limitSector.size()));

    std::ifstream a(filename, std::ios::binary | std::ios::ate);
    std::ifstream b(tmpFile, std::ios::binary | std::ios::ate);
    if (a.tellg() != b.tellg()) {
        std::cerr << "FAIL: size mismatch (" << a.tellg() << " vs " << b.tellg()
                  << ")" << std::endl;
        std::remove(tmpFile.c_str());
        return false;
    }
    a.seekg(0); b.seekg(0);
    size_t pos = 0;
    char ca, cb;
    while (a.get(ca) && b.get(cb)) {
        if (ca != cb) {
            std::cerr << "FAIL: byte mismatch at offset " << pos << std::endl;
            std::remove(tmpFile.c_str());
            return false;
        }
        ++pos;
    }
    std::remove(tmpFile.c_str());
    std::cout << "PASS: " << filename << " ring round-trip byte-identical" << std::endl;
    return true;
}

#endif
