// family_generate.cpp — CLI: family.json → .bin files
//
// Full pipeline:
//   FamilyConfig → IBPEqGenerator → IBPAnalyzer → RegionSolver
//   → RecursionBuilder → RingBuilder → BinaryWriter
//
// Usage:
//   ./family_generate families/bub00.json
//   ./family_generate families/bub00.json --output /tmp/out --diff --source-dir .

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <sys/stat.h>

#include "FamilyConfig.hpp"
#include "IBPEqGenerator.hpp"
#include "IBPAnalyzer.hpp"
#include "RegionSolver.hpp"
#include "RecursionBuilder.hpp"
#include "RingBuilder.hpp"
#include "BinaryIBPWriter.hpp"
#include "BinaryRingWriter.hpp"

using namespace std;
using namespace firefly;

// ============================================================
// Dense → sparse CSR conversion
//
// RecursionBuilder produces dense FlatMatrix (nb×nb) per operator entry.
// BinaryIBPWriter expects SparseTensorRaw (CSR format matching MMA).
// ============================================================

// Helper: build CSR from a 4D array [nibp][ne] of FlatMatrix (nb×nb)
template<typename T>
SparseTensorRaw<T> denseToSparse4D(
    const vector<vector<RecursionBuilder::FlatMatrix>>& data,  // [nibp][ne]
    int nibp, int ne, int nb)
{
    SparseTensorRaw<T> result;
    result.dims = {int32_t(nibp), int32_t(ne), int32_t(nb), int32_t(nb)};

    result.rowPtr.resize(nibp + 1, 0);
    int32_t nz = 0;

    for (int m = 0; m < nibp; ++m) {
        for (int i = 0; i < ne; ++i) {
            for (int r = 0; r < nb; ++r) {
                for (int c = 0; c < nb; ++c) {
                    int64_t val = RecursionBuilder::matAt(data[m][i], nb, r, c);
                    if (val != 0) {
                        // MMA IBP1 format: colIdx stores multi-dimensional coords (1-based)
                        result.colIdx.push_back(i + 1);
                        result.colIdx.push_back(r + 1);
                        result.colIdx.push_back(c + 1);
                        if constexpr (is_same_v<T, FFInt>)
                            result.values.push_back(FFInt(static_cast<uint64_t>(val)));
                        else
                            result.values.push_back(static_cast<T>(val));
                        ++nz;
                    }
                }
            }
        }
        result.rowPtr[m + 1] = nz;
    }
    return result;
}

// 3D: [nibp] of FlatMatrix (for F0)
template<typename T>
SparseTensorRaw<T> denseToSparse3D(
    const vector<RecursionBuilder::FlatMatrix>& data,  // [nibp]
    int nibp, int nb)
{
    SparseTensorRaw<T> result;
    result.dims = {int32_t(nibp), int32_t(nb), int32_t(nb)};

    result.rowPtr.resize(nibp + 1, 0);
    int32_t nz = 0;

    for (int m = 0; m < nibp; ++m) {
        for (int r = 0; r < nb; ++r) {
            for (int c = 0; c < nb; ++c) {
                int64_t val = RecursionBuilder::matAt(data[m], nb, r, c);
                if (val != 0) {
                    // MMA IBP1 format: colIdx stores multi-dimensional coords (1-based)
                    result.colIdx.push_back(r + 1);
                    result.colIdx.push_back(c + 1);
                    if constexpr (is_same_v<T, FFInt>)
                        result.values.push_back(FFInt(static_cast<uint64_t>(val)));
                    else
                        result.values.push_back(static_cast<T>(val));
                    ++nz;
                }
            }
        }
        result.rowPtr[m + 1] = nz;
    }
    return result;
}

// 5D: [nibp][ne][ne] of FlatMatrix (for F2, F2s)
template<typename T>
SparseTensorRaw<T> denseToSparse5D(
    const vector<vector<vector<RecursionBuilder::FlatMatrix>>>& data,  // [nibp][ne][ne]
    int nibp, int ne, int nb)
{
    SparseTensorRaw<T> result;
    result.dims = {int32_t(nibp), int32_t(ne), int32_t(ne), int32_t(nb), int32_t(nb)};

    result.rowPtr.resize(nibp + 1, 0);
    int32_t nz = 0;

    for (int m = 0; m < nibp; ++m) {
        for (int i = 0; i < ne; ++i) {
            for (int j = 0; j < ne; ++j) {
                for (int r = 0; r < nb; ++r) {
                    for (int c = 0; c < nb; ++c) {
                        int64_t val = RecursionBuilder::matAt(data[m][i][j], nb, r, c);
                        if (val != 0) {
                            // MMA IBP1 format: colIdx stores multi-dimensional coords (1-based)
                            result.colIdx.push_back(i + 1);
                            result.colIdx.push_back(j + 1);
                            result.colIdx.push_back(r + 1);
                            result.colIdx.push_back(c + 1);
                            if constexpr (is_same_v<T, FFInt>)
                                result.values.push_back(FFInt(static_cast<uint64_t>(val)));
                            else
                                result.values.push_back(static_cast<T>(val));
                            ++nz;
                        }
                    }
                }
            }
        }
        result.rowPtr[m + 1] = nz;
    }
    return result;
}

// Convert one region's RecursionMatrices → RegimeSparseData
template<typename T>
RegimeSparseData<T> convertToRegimeSparse(
    const RecursionBuilder::RecursionMatrices& rm,
    int64_t modulus)
{
    RegimeSparseData<T> reg;
    reg.nibp    = static_cast<int32_t>(rm.nibp);
    reg.ne      = static_cast<int32_t>(rm.ne);
    reg.nb      = static_cast<int32_t>(rm.nb);
    reg.incre   = 2;  // default in MMA
    reg.modulus = modulus;

    reg.M1  = denseToSparse4D<T>(rm.M1,  rm.nibp, rm.ne, rm.nb);
    reg.N1  = denseToSparse4D<T>(rm.N1,  rm.nibp, rm.ne, rm.nb);
    reg.K1  = denseToSparse4D<T>(rm.K1,  rm.nibp, rm.ne, rm.nb);
    reg.K1s = denseToSparse4D<T>(rm.K1s, rm.nibp, rm.ne, rm.nb);
    reg.K2s = denseToSparse4D<T>(rm.K2s, rm.nibp, rm.ne, rm.nb);
    reg.F0  = denseToSparse3D<T>(rm.F0,  rm.nibp, rm.nb);
    reg.F2  = denseToSparse5D<T>(rm.F2,  rm.nibp, rm.ne, rm.nb);
    reg.F2s = denseToSparse5D<T>(rm.F2s, rm.nibp, rm.ne, rm.nb);

    return reg;
}

// Convert RingMatrices → RingCellRaw
template<typename T>
RingCellRaw<T> convertToRingCell(
    const RingBuilder::RingMatrices& rm,
    const vector<int>& limitSector)
{
    RingCellRaw<T> cell;
    cell.limitSector = limitSector;

    int ne = static_cast<int>(rm.A_list.size());
    int nb = 0;
    if (ne > 0 && !rm.A_list.empty()) {
        nb = static_cast<int>(sqrt(rm.A_list[0].size()));
    }
    cell.nb = nb;

    // Flatten: concatenate A_list[0..ne-1] → A_flat
    for (int i = 0; i < ne; ++i) {
        for (auto v : rm.A_list[i]) {
            if constexpr (is_same_v<T, FFInt>)
                cell.A_flat.push_back(FFInt(static_cast<uint64_t>(v)));
            else
                cell.A_flat.push_back(static_cast<T>(v));
        }
    }
    for (int i = 0; i < ne; ++i) {
        for (auto v : rm.Ainv_list[i]) {
            if constexpr (is_same_v<T, FFInt>)
                cell.Ainv_flat.push_back(FFInt(static_cast<uint64_t>(v)));
            else
                cell.Ainv_flat.push_back(static_cast<T>(v));
        }
    }

    return cell;
}

// ============================================================
// File utilities
// ============================================================

inline bool mkdirP(const string& path) {
    string cmd = "mkdir -p \"" + path + "\"";
    return system(cmd.c_str()) == 0;
}

inline bool filesIdentical(const string& a, const string& b) {
    ifstream fa(a, ios::binary | ios::ate);
    ifstream fb(b, ios::binary | ios::ate);
    if (!fa || !fb) return false;
    if (fa.tellg() != fb.tellg()) return false;
    fa.seekg(0); fb.seekg(0);
    return equal(istreambuf_iterator<char>(fa), istreambuf_iterator<char>(),
                 istreambuf_iterator<char>(fb));
}

// ============================================================
// Main
// ============================================================
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << R"(Usage: family_generate <family.json> [options]

  Reads a family.json config, generates IBPMat_*.bin and RingData_*.bin.

Options:
  --source-dir DIR   Comparison source directory for --diff (default: .)
  --output DIR       Output directory (default: .)
  --diff             Compare generated .bin files against source-dir originals
)";
        return 1;
    }

    string jsonFile  = argv[1];
    string sourceDir = ".";
    string outputDir = ".";
    bool   doDiff    = false;

    for (int i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--source-dir" && i + 1 < argc) {
            sourceDir = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            outputDir = argv[++i];
        } else if (arg == "--diff") {
            doDiff = true;
        } else {
            cerr << "Unknown flag: " << arg << endl;
            return 1;
        }
    }

    using Clock = std::chrono::high_resolution_clock;
    using Ms    = std::chrono::milliseconds;
    auto t_total_0 = Clock::now();

    if (outputDir != ".") mkdirP(outputDir);

    // ========================================
    // 1. Parse family config
    // ========================================
    cout << "=== FamilyGenerate ===" << endl;
    FamilyDef fam;
    try {
        fam = parseFamilyConfig(jsonFile);
    } catch (const exception& e) {
        cerr << "Failed to parse " << jsonFile << ": " << e.what() << endl;
        return 1;
    }

    cout << "Family: " << fam.name
         << " (L=" << fam.nLoop() << " E=" << fam.nExt()
         << " N=" << fam.nProp() << ")" << endl;
    cout << "Modulus: " << fam.modulus << endl;

    // Set FFInt modulus
    FFInt::set_new_prime(static_cast<uint64_t>(fam.modulus));

    // ========================================
    // 2. Generate IBP equations
    // ========================================
    cout << "\n[Step 1] Generating IBP equations..." << endl;
    auto t1_0 = Clock::now();
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    auto t1_1 = Clock::now();
    int ne = ibp.ne;
    cout << "  ne=" << ne << " nibp=" << ibp.nibp
         << " subsectors=" << ibp.sectorlist.size() << endl;
    cout << "  [TIMING] Step 1 (IBP eq gen): "
         << std::chrono::duration_cast<Ms>(t1_1 - t1_0).count() / 1000.0 << " s" << endl;

    // ========================================
    // 3. Extract FTable
    // ========================================
    cout << "\n[Step 2] Extracting FTable..." << endl;
    auto t2_0 = Clock::now();
    auto ft = IBPAnalyzer::extractFTable(ibp, fam.modulus);
    auto t2_1 = Clock::now();
    cout << "  nibp=" << ft.nibp << " ne=" << ft.ne << endl;
    cout << "  [TIMING] Step 2 (FTable): "
         << std::chrono::duration_cast<Ms>(t2_1 - t2_0).count() / 1000.0 << " s" << endl;

    // ========================================
    // 4. Solve all sectors → regions
    // ========================================
    cout << "\n[Step 3] Solving all sectors..." << endl;
    auto t3_0 = Clock::now();
    auto topSector = fam.topSector;  // top sector from family config (all active propagators)
    auto regions = RegionSolver::solveAllSectors(ibp, topSector, fam.modulus);
    auto t3_1 = Clock::now();
    cout << "  Found " << regions.size() << " region(s)" << endl;
    cout << "  [TIMING] Step 3 (Region solve): "
         << std::chrono::duration_cast<Ms>(t3_1 - t3_0).count() / 1000.0 << " s" << endl;

    if (regions.empty()) {
        cerr << "No zero-dimensional regions found. Aborting." << endl;
        return 1;
    }

    // Sort regions by sector to match MMA's output order
    // MMA order: {1,0}, {0,1}, then {1,1} groups (ascending by sum, then descending by first element)
    std::sort(regions.begin(), regions.end(), [](const RegionSolver::RegionData& a, const RegionSolver::RegionData& b) {
        const auto& sa = a.limitSector;
        const auto& sb = b.limitSector;
        int suma = 0, sumb = 0;
        for (int v : sa) suma += v;
        for (int v : sb) sumb += v;
        if (suma != sumb) return suma < sumb;
        // For equal sums, descending by first element (so {1,0} before {0,1})
        return sa[0] > sb[0];
    });

    // ========================================
    // 5. Build recursion + ring matrices per region
    // ========================================
    cout << "\n[Step 4] Building recursion + ring matrices..." << endl;
    auto t4_0 = Clock::now();

    vector<RegimeSparseData<FFInt>> ibpRegimes;
    vector<RingCellRaw<FFInt>>      ringCells;

    for (size_t ridx = 0; ridx < regions.size(); ++ridx) {
        auto& reg = regions[ridx];
        cout << "  Region " << ridx << ": nb=" << reg.nb << endl;

        // Build recursion matrices
        auto recMats = RecursionBuilder::buildRecursionMatrices(ft, reg, fam.modulus);

        // Build ring matrices
        auto ringMats = RingBuilder::computeRingMatrices(reg, ne, fam.modulus);

        // Convert to sparse + raw formats
        ibpRegimes.push_back(convertToRegimeSparse<FFInt>(recMats, fam.modulus));
        ringCells.push_back(convertToRingCell<FFInt>(ringMats, reg.limitSector));
    }
    auto t4_1 = Clock::now();
    cout << "  [TIMING] Step 4 (Recursion+Ring build): "
         << std::chrono::duration_cast<Ms>(t4_1 - t4_0).count() / 1000.0 << " s" << endl;

    // ========================================
    // 6. Write output .bin files
    // ========================================
    string ibpOut  = outputDir + "/IBPMat_"  + fam.name + ".bin";
    string ringOut = outputDir + "/RingData_" + fam.name + ".bin";

    cout << "\n[Step 5] Writing .bin files..." << endl;
    auto t5_0 = Clock::now();
    writeIBPMatrixBinary<FFInt>(ibpOut, ibpRegimes);
    writeRingDataBinary<FFInt>(ringOut, ringCells, ne);
    auto t5_1 = Clock::now();
    cout << "  [TIMING] Step 5 (Write .bin): "
         << std::chrono::duration_cast<Ms>(t5_1 - t5_0).count() / 1000.0 << " s" << endl;

    auto t_total_1 = Clock::now();
    double total_s = std::chrono::duration_cast<Ms>(t_total_1 - t_total_0).count() / 1000.0;
    cout << "\n[TIMING] Total wall-clock: " << total_s << " s" << endl;

    // ========================================
    // 7. Optional diff
    // ========================================
    if (doDiff) {
        string ibpSrc  = sourceDir + "/IBPMat_"  + fam.name + ".bin";
        string ringSrc = sourceDir + "/RingData_" + fam.name + ".bin";
        bool ibpOk  = filesIdentical(ibpSrc, ibpOut);
        bool ringOk = filesIdentical(ringSrc, ringOut);

        cout << "\n--- Diff ---" << endl;
        cout << "IBP:  " << (ibpOk  ? "IDENTICAL" : "DIFFER") << endl;
        cout << "Ring: " << (ringOk ? "IDENTICAL" : "DIFFER") << endl;
        cout << "=== " << (ibpOk && ringOk ? "DIFF PASS" : "DIFF FAIL") << " ===" << endl;
        return ibpOk && ringOk ? 0 : 1;
    }

    cout << "\nDone. Output:" << endl;
    cout << "  " << ibpOut << endl;
    cout << "  " << ringOut << endl;
    return 0;
}
