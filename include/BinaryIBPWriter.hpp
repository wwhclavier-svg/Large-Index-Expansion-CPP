#ifndef BINARY_IBP_WRITER_HPP
#define BINARY_IBP_WRITER_HPP

#include <fstream>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <type_traits>

#include "firefly/FFInt.hpp"
#include "IBPMatrixLoader_Binary.hpp"

// ==========================================
// Big-endian write helpers
// ==========================================
inline void writeBE32(std::ofstream& out, int32_t val) {
    uint8_t b[4];
    b[0] = static_cast<uint8_t>((val >> 24) & 0xFF);
    b[1] = static_cast<uint8_t>((val >> 16) & 0xFF);
    b[2] = static_cast<uint8_t>((val >> 8)  & 0xFF);
    b[3] = static_cast<uint8_t>( val        & 0xFF);
    out.write(reinterpret_cast<const char*>(b), 4);
}

inline void writeBE64(std::ofstream& out, int64_t val) {
    uint8_t b[8];
    b[0] = static_cast<uint8_t>((val >> 56) & 0xFF);
    b[1] = static_cast<uint8_t>((val >> 48) & 0xFF);
    b[2] = static_cast<uint8_t>((val >> 40) & 0xFF);
    b[3] = static_cast<uint8_t>((val >> 32) & 0xFF);
    b[4] = static_cast<uint8_t>((val >> 24) & 0xFF);
    b[5] = static_cast<uint8_t>((val >> 16) & 0xFF);
    b[6] = static_cast<uint8_t>((val >> 8)  & 0xFF);
    b[7] = static_cast<uint8_t>( val        & 0xFF);
    out.write(reinterpret_cast<const char*>(b), 8);
}

inline void writeBEDouble(std::ofstream& out, double val) {
    uint8_t raw[8], swapped[8];
    std::memcpy(raw, &val, 8);
    for (int i = 0; i < 8; ++i) swapped[i] = raw[7 - i];
    out.write(reinterpret_cast<const char*>(swapped), 8);
}

// ==========================================
// Sparse tensor in original (un-flattened) format
// Matches Mathematica's SparseArray CSR representation
// ==========================================
template<typename T>
struct SparseTensorRaw {
    std::vector<int32_t> dims;          // original dims (not C++ flattened)
    std::vector<int32_t> rowPtr;        // dims[0] + 1 entries
    std::vector<int32_t> colIdx;        // flat: nz * (dims.size()-1) entries, 1-based
    std::vector<T>        values;       // non-zero values

    bool empty() const { return dims.empty(); }

    void writeBinary(std::ofstream& out) const {
        // dims_len + dims[]
        writeBE32(out, static_cast<int32_t>(dims.size()));
        for (auto d : dims) writeBE32(out, d);

        // rowPtr
        writeBE32(out, static_cast<int32_t>(rowPtr.size()));
        for (auto r : rowPtr) writeBE32(out, r);

        // colIdx
        writeBE32(out, static_cast<int32_t>(colIdx.size()));
        for (auto c : colIdx) writeBE32(out, c);

        // values
        writeBE32(out, static_cast<int32_t>(values.size()));
        for (auto& v : values) {
            if constexpr (std::is_same_v<T, firefly::FFInt>)
                writeBE64(out, static_cast<int64_t>(v.n));
            else
                writeBEDouble(out, static_cast<double>(v));
        }
    }
};

// ==========================================
// Raw regime data (one-per-regime, per-operator sparse tensors)
// ==========================================
template<typename T>
struct RegimeSparseData {
    int32_t nibp, ne, nb, incre;
    int64_t modulus;

    // Operators in canonical order
    SparseTensorRaw<T> M1, N1, K1, F0, F2, K1s, K2s, F2s;

    bool empty() const { return M1.empty() && N1.empty(); }
};

// ==========================================
// Main write function
// ==========================================
template<typename T>
void writeIBPMatrixBinary(const std::string& filename,
                          const std::vector<RegimeSparseData<T>>& regimes) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);

    // Magic "IBP1"
    out.write("IBP1", 4);

    // Number of regimes
    writeBE32(out, static_cast<int32_t>(regimes.size()));

    static const std::vector<std::string> opOrder =
        {"M1", "N1", "K1", "F0", "F2", "K1s", "K2s", "F2s"};

    for (const auto& reg : regimes) {
        // Metadata
        writeBE32(out, reg.nibp);
        writeBE32(out, reg.ne);
        writeBE32(out, reg.nb);
        writeBE32(out, reg.incre);
        writeBE64(out, reg.modulus);

        // Operators in order
        for (const auto& op : opOrder) {
            const SparseTensorRaw<T>* tensor = nullptr;
            if (op == "M1")  tensor = &reg.M1;
            if (op == "N1")  tensor = &reg.N1;
            if (op == "K1")  tensor = &reg.K1;
            if (op == "F0")  tensor = &reg.F0;
            if (op == "F2")  tensor = &reg.F2;
            if (op == "K1s") tensor = &reg.K1s;
            if (op == "K2s") tensor = &reg.K2s;
            if (op == "F2s") tensor = &reg.F2s;

            uint8_t exists = (tensor && !tensor->empty()) ? 1 : 0;
            out.write(reinterpret_cast<const char*>(&exists), 1);
            if (exists)
                tensor->writeBinary(out);
        }
    }

    out.close();
    std::cout << "Written " << regimes.size() << " regime(s) to " << filename
              << std::endl;
}

// ==========================================
// Read .bin into raw sparse format (preserves original dims/coords)
// Used for verification: read → write → compare byte-for-byte
// ==========================================
template<typename T>
std::vector<RegimeSparseData<T>> readIBPMatrixRaw(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    char magic[4];
    in.read(magic, 4);
    if (std::string(magic, 4) != "IBP1")
        throw std::runtime_error("Invalid magic number");

    int32_t numRegs = readBE32(in);
    std::vector<RegimeSparseData<T>> results;
    results.reserve(numRegs);

    static const std::vector<std::string> opOrder =
        {"M1", "N1", "K1", "F0", "F2", "K1s", "K2s", "F2s"};

    for (int r = 0; r < numRegs; ++r) {
        RegimeSparseData<T> reg;
        reg.nibp   = readBE32(in);
        reg.ne     = readBE32(in);
        reg.nb     = readBE32(in);
        reg.incre  = readBE32(in);
        reg.modulus = readBE64(in);

        for (const auto& op : opOrder) {
            uint8_t exists;
            in.read(reinterpret_cast<char*>(&exists), 1);
            if (!exists) continue;

            SparseTensorRaw<T> tensor;

            // Read dims
            int32_t dims_len = readBE32(in);
            tensor.dims.resize(dims_len);
            for (int i = 0; i < dims_len; ++i) tensor.dims[i] = readBE32(in);

            // Read rowPtr
            int32_t rowPtr_len = readBE32(in);
            tensor.rowPtr.resize(rowPtr_len);
            for (int i = 0; i < rowPtr_len; ++i) tensor.rowPtr[i] = readBE32(in);

            // Read colIdx
            int32_t colIdx_len = readBE32(in);
            tensor.colIdx.resize(colIdx_len);
            for (int i = 0; i < colIdx_len; ++i) tensor.colIdx[i] = readBE32(in);

            // Read values
            int32_t values_len = readBE32(in);
            tensor.values.resize(values_len);
            if constexpr (std::is_same_v<T, firefly::FFInt>) {
                for (int i = 0; i < values_len; ++i) {
                    int64_t val = readBE64(in);
                    tensor.values[i] = firefly::FFInt(val);
                }
            } else {
                for (int i = 0; i < values_len; ++i)
                    tensor.values[i] = readBEDouble(in);
            }

            // Assign to struct
            if (op == "M1")  reg.M1  = std::move(tensor);
            if (op == "N1")  reg.N1  = std::move(tensor);
            if (op == "K1")  reg.K1  = std::move(tensor);
            if (op == "F0")  reg.F0  = std::move(tensor);
            if (op == "F2")  reg.F2  = std::move(tensor);
            if (op == "K1s") reg.K1s = std::move(tensor);
            if (op == "K2s") reg.K2s = std::move(tensor);
            if (op == "F2s") reg.F2s = std::move(tensor);
        }

        results.push_back(std::move(reg));
    }

    return results;
}

// ==========================================
// Test utility: verify round-trip byte-for-byte
// ==========================================
template<typename T>
bool verifyRoundTrip(const std::string& filename) {
    auto raw = readIBPMatrixRaw<T>(filename);
    std::string tmpFile = filename + ".roundtrip_test";
    writeIBPMatrixBinary<T>(tmpFile, raw);

    // Compare files byte-for-byte
    std::ifstream a(filename, std::ios::binary | std::ios::ate);
    std::ifstream b(tmpFile, std::ios::binary | std::ios::ate);
    if (a.tellg() != b.tellg()) {
        std::cerr << "FAIL: size mismatch" << std::endl;
        std::remove(tmpFile.c_str());
        return false;
    }
    a.seekg(0); b.seekg(0);
    size_t pos = 0;
    char ca, cb;
    while (a.get(ca) && b.get(cb)) {
        if (ca != cb) {
            std::cerr << "FAIL: byte mismatch at offset " << pos
                      << " (0x" << std::hex << (int)(unsigned char)ca
                      << " vs 0x" << (int)(unsigned char)cb << std::dec << ")" << std::endl;
            std::remove(tmpFile.c_str());
            return false;
        }
        ++pos;
    }
    std::remove(tmpFile.c_str());
    std::cout << "PASS: " << filename << " round-trip byte-identical" << std::endl;
    return true;
}

#endif
