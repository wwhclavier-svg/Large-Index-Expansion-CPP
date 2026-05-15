// family_generate.cpp — CLI: family.json → .bin files
//
// Full pipeline:
//   FamilyConfig → IBPEqGenerator → IBPAnalyzer → RegionSolver
//   → RecursionBuilder → RingBuilder → BinaryWriter
//
// Modes:
//   Default:               Full pipeline, all sectors → IBPMat + RingData
//   --sector <bin> [--cache-dir <dir>]:  Only solve one sector
//   --merge-sectors <dir> [--output <dir>]:  Merge partial outputs from --sector runs
//
// Usage:
//   ./family_generate families/bub00.json
//   ./family_generate families/bub00.json --sector 11 --cache-dir cache/sectors
//   ./family_generate families/bub00.json --merge-sectors cache/sectors

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <cstdlib>
#include <sys/stat.h>
#include <glob.h>
#include <filesystem>

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
// Dense → sparse CSR conversion (same as before)
// ============================================================

template<typename T>
SparseTensorRaw<T> denseToSparse4D(
    const vector<vector<RecursionBuilder::FlatMatrix>>& data,
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

template<typename T>
SparseTensorRaw<T> denseToSparse3D(
    const vector<RecursionBuilder::FlatMatrix>& data,
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

template<typename T>
SparseTensorRaw<T> denseToSparse5D(
    const vector<vector<vector<RecursionBuilder::FlatMatrix>>>& data,
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

template<typename T>
RegimeSparseData<T> convertToRegimeSparse(
    const RecursionBuilder::RecursionMatrices& rm,
    int64_t modulus)
{
    RegimeSparseData<T> reg;
    reg.nibp    = static_cast<int32_t>(rm.nibp);
    reg.ne      = static_cast<int32_t>(rm.ne);
    reg.nb      = static_cast<int32_t>(rm.nb);
    reg.incre   = 2;
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

template<typename T>
RingCellRaw<T> convertToRingCell(
    const RingBuilder::RingMatrices& rm,
    const vector<int>& limitSector)
{
    RingCellRaw<T> cell;
    cell.limitSector = limitSector;
    int ne = static_cast<int>(rm.A_list.size());
    int nb = 0;
    if (ne > 0 && !rm.A_list.empty())
        nb = static_cast<int>(sqrt(rm.A_list[0].size()));
    cell.nb = nb;
    for (int i = 0; i < ne; ++i)
        for (auto v : rm.A_list[i])
            cell.A_flat.push_back(FFInt(static_cast<uint64_t>(v)));
    for (int i = 0; i < ne; ++i)
        for (auto v : rm.Ainv_list[i])
            cell.Ainv_flat.push_back(FFInt(static_cast<uint64_t>(v)));
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

// Glob: find files matching pattern
inline vector<string> globFiles(const string& pattern) {
    vector<string> results;
    glob_t g;
    if (::glob(pattern.c_str(), GLOB_TILDE, nullptr, &g) == 0) {
        for (size_t i = 0; i < g.gl_pathc; ++i)
            results.push_back(g.gl_pathv[i]);
        globfree(&g);
    }
    return results;
}

// ============================================================
// Merge partial sector files
// ============================================================

// Read raw bytes from a file, starting after an initial header skip
inline vector<char> readBytesAfter(const string& path, int64_t skipBytes) {
    ifstream f(path, ios::binary | ios::ate);
    if (!f) throw runtime_error("Cannot open: " + path);
    int64_t total = f.tellg();
    int64_t remaining = total - skipBytes;
    if (remaining <= 0) return {};
    f.seekg(skipBytes);
    vector<char> buf(remaining);
    f.read(buf.data(), remaining);
    return buf;
}

// Read just the regime count from an IBP1 partial file
inline int32_t readRegimeCount(const string& path) {
    ifstream f(path, ios::binary);
    char magic[4];
    f.read(magic, 4);
    if (string(magic, 4) != "IBP1") return 0;
    // Read BE int32
    uint8_t b[4];
    f.read(reinterpret_cast<char*>(b), 4);
    return (b[0] << 24) | (b[1] << 16) | (b[2] << 8) | b[3];
}

// Read cell count from RingData partial file (native endian int32)
inline int32_t readCellCount(const string& path) {
    ifstream f(path, ios::binary);
    int32_t count;
    f.read(reinterpret_cast<char*>(&count), sizeof(count));
    return count;
}

// Read globalNe from RingData partial file (native endian int32)
inline int32_t readRingNe(const string& path) {
    ifstream f(path, ios::binary);
    f.seekg(4); // skip cell count
    int32_t ne;
    f.read(reinterpret_cast<char*>(&ne), sizeof(ne));
    return ne;
}

void mergeIBPSectors(const string& pattern, const string& output) {
    auto files = globFiles(pattern);
    if (files.empty()) {
        cerr << "[MERGE] No partial IBP files found matching: " << pattern << endl;
        return;
    }
    sort(files.begin(), files.end());

    int32_t totalRegimes = 0;
    for (const auto& f : files)
        totalRegimes += readRegimeCount(f);

    ofstream out(output, ios::binary);
    if (!out) throw runtime_error("Cannot write: " + output);

    out.write("IBP1", 4);
    // Write BE int32
    uint8_t hdr[4] = { uint8_t((totalRegimes >> 24) & 0xFF),
                       uint8_t((totalRegimes >> 16) & 0xFF),
                       uint8_t((totalRegimes >> 8) & 0xFF),
                       uint8_t(totalRegimes & 0xFF) };
    out.write(reinterpret_cast<const char*>(hdr), 4);

    for (const auto& f : files) {
        auto payload = readBytesAfter(f, 8); // skip magic(4) + count(4)
        out.write(payload.data(), payload.size());
    }
    out.close();
    cout << "[MERGE] Merged " << files.size() << " partial IBP files → " << output
         << " (" << totalRegimes << " regimes)" << endl;
}

void mergeRingSectors(const string& pattern, const string& output) {
    auto files = globFiles(pattern);
    if (files.empty()) {
        cerr << "[MERGE] No partial RingData files found matching: " << pattern << endl;
        return;
    }
    sort(files.begin(), files.end());

    int32_t totalCells = 0;
    int32_t globalNe = 0;
    for (const auto& f : files) {
        totalCells += readCellCount(f);
        if (globalNe == 0) globalNe = readRingNe(f);
    }

    ofstream out(output, ios::binary);
    if (!out) throw runtime_error("Cannot write: " + output);

    out.write(reinterpret_cast<const char*>(&totalCells), sizeof(totalCells));
    out.write(reinterpret_cast<const char*>(&globalNe), sizeof(globalNe));

    for (const auto& f : files) {
        auto payload = readBytesAfter(f, 8); // skip cellCount(4) + ne(4)
        out.write(payload.data(), payload.size());
    }
    out.close();
    cout << "[MERGE] Merged " << files.size() << " partial RingData files → " << output
         << " (" << totalCells << " cells)" << endl;
}

// ============================================================
// Main
// ============================================================
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << R"(Usage: family_generate <family.json> [options]

  Reads a family.json config, generates IBPMat_*.bin and RingData_*.bin.

Options:
  --sector <binary>        Only solve one subsector (e.g. "110")
  --cache-dir DIR          Cache directory for per-sector partial output
  --merge-sectors DIR      Merge partial sector outputs from DIR
  --output DIR             Output directory (default: data/)
  --source-dir DIR         Comparison source for --diff (default: data/)
  --diff                   Compare against source-dir originals
)";
        return 1;
    }

    string jsonFile     = argv[1];
    string sourceDir    = "data";
    string outputDir    = "data";
    string sectorFilter = "";
    string cacheDir     = "";
    string mergeDir     = "";
    bool   doDiff       = false;

    for (int i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--source-dir" && i + 1 < argc) sourceDir = argv[++i];
        else if (arg == "--output" && i + 1 < argc) outputDir = argv[++i];
        else if (arg == "--sector" && i + 1 < argc) sectorFilter = argv[++i];
        else if (arg == "--cache-dir" && i + 1 < argc) cacheDir = argv[++i];
        else if (arg == "--merge-sectors" && i + 1 < argc) mergeDir = argv[++i];
        else if (arg == "--diff") doDiff = true;
        else { cerr << "Unknown flag: " << arg << endl; return 1; }
    }

    // ── Merge mode ────────────────────────────────────────────────
    if (!mergeDir.empty()) {
        if (outputDir == ".") outputDir = mergeDir; // default: write back to merge dir
        string famName = "";
        auto ibpFiles = globFiles(mergeDir + "/IBPMat_*_sector_*.bin");
        if (!ibpFiles.empty()) {
            // Extract family name from first file: IBPMat_<fam>_sector_<key>.bin
            string fn = filesystem::path(ibpFiles[0]).filename().string();
            auto p1 = fn.find("_sector_");
            if (p1 != string::npos)
                famName = fn.substr(7, p1 - 7); // after "IBPMat_" (7 chars)
        }
        string ibpOut  = outputDir + "/IBPMat_"  + famName + ".bin";
        string ringOut = outputDir + "/RingData_" + famName + ".bin";

        mergeIBPSectors(mergeDir + "/IBPMat_*_sector_*.bin", ibpOut);
        mergeRingSectors(mergeDir + "/RingData_*_sector_*.bin", ringOut);
        return 0;
    }

    // ── Normal / sector mode ──────────────────────────────────────
    using Clock = std::chrono::high_resolution_clock;
    using Ms    = std::chrono::milliseconds;
    auto t_total_0 = Clock::now();

    if (outputDir != ".") mkdirP(outputDir);

    // 1. Parse family config
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

    FFInt::set_new_prime(static_cast<uint64_t>(fam.modulus));

    // 2. Generate IBP equations
    cout << "\n[Step 1] Generating IBP equations..." << endl;
    auto t1_0 = Clock::now();
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    auto t1_1 = Clock::now();
    int ne = ibp.ne;
    cout << "  ne=" << ne << " nibp=" << ibp.nibp
         << " subsectors=" << ibp.sectorlist.size() << endl;
    cout << "  [TIMING] Step 1: "
         << chrono::duration_cast<Ms>(t1_1 - t1_0).count() / 1000.0 << " s" << endl;

    // 3. Extract FTable
    cout << "\n[Step 2] Extracting FTable..." << endl;
    auto t2_0 = Clock::now();
    auto ft = IBPAnalyzer::extractFTable(ibp, fam.modulus);
    auto t2_1 = Clock::now();
    cout << "  nibp=" << ft.nibp << " ne=" << ft.ne << endl;
    cout << "  [TIMING] Step 2: "
         << chrono::duration_cast<Ms>(t2_1 - t2_0).count() / 1000.0 << " s" << endl;

    // 4. Solve sectors → regions
    cout << "\n[Step 3] Solving sectors..." << endl;
    auto t3_0 = Clock::now();
    auto topSector = fam.topSector;
    vector<RegionSolver::RegionData> regions;

    if (!sectorFilter.empty()) {
        // Single-sector mode: only solve the specified sector
        cout << "  Filtering: sector = " << sectorFilter << endl;
        auto subsectors = RegionSolver::generateSubsectors(topSector, ibp.nl);
        bool found = false;
        for (const auto& sub : subsectors) {
            string key;
            for (int b : sub) key += to_string(b);
            if (key == sectorFilter) {
                found = true;
                auto abEqs = IBPAnalyzer::buildABEquations(ibp, sub, fam.modulus);
                bool hasContent = false;
                for (const auto& eq : abEqs) {
                    for (char ch : eq)
                        if (ch != '0' && ch != '+' && ch != '-' && ch != '*') { hasContent = true; break; }
                    if (hasContent) break;
                    if (!eq.empty() && eq != "0") hasContent = true;
                }
                if (hasContent) {
                    cerr << "  Solving sector [" << key << "]" << endl;
                    auto t_s0 = Clock::now();
                    auto sregions = RegionSolver::solveRegion(abEqs, sub, ne, fam.modulus);
                    auto t_s1 = Clock::now();
                    double sec = chrono::duration_cast<Ms>(t_s1 - t_s0).count() / 1000.0;
                    cerr << "  -> " << sregions.size() << " region(s) in " << sec << " s" << endl;
                    for (auto& r : sregions)
                        regions.push_back(move(r));
                }
                break;
            }
        }
        if (!found)
            cerr << "  Sector " << sectorFilter << " not found in family" << endl;
    } else {
        // Normal mode: solve all sectors
        regions = RegionSolver::solveAllSectors(ibp, topSector, fam.modulus);
    }
    auto t3_1 = Clock::now();
    cout << "  Found " << regions.size() << " region(s) total" << endl;
    cout << "  [TIMING] Step 3: "
         << chrono::duration_cast<Ms>(t3_1 - t3_0).count() / 1000.0 << " s" << endl;

    vector<RegimeSparseData<FFInt>> ibpRegimes;
    vector<RingCellRaw<FFInt>>      ringCells;

    if (regions.empty()) {
        cout << "  No regions found. Writing empty output." << endl;
    } else {
        // Sort regions by sector (MMA ordering: weight ascending, then binary descending)
        // Within same sector, sort by nb descending then by raw prime content for deterministic order
        sort(regions.begin(), regions.end(),
            [](const RegionSolver::RegionData& a, const RegionSolver::RegionData& b) {
                const auto& sa = a.limitSector;
                const auto& sb = b.limitSector;
                int ne = (int)sa.size();
                int suma = 0, sumb = 0;
                int va = 0, vb = 0;
                for (int i = 0; i < ne; ++i) {
                    suma += sa[i]; sumb += sb[i];
                    va = (va << 1) | sa[i];
                    vb = (vb << 1) | sb[i];
                }
                if (suma != sumb) return suma < sumb;
                if (va != vb) return va > vb;  // descending binary within same weight
                // Same sector: tiebreak by nb (descending) for deterministic order
                if (a.nb != b.nb) return a.nb > b.nb;
                // Same nb: tiebreak by MonomialBasis content
                if (a.MonomialBasis.size() != b.MonomialBasis.size())
                    return a.MonomialBasis.size() < b.MonomialBasis.size();
                for (size_t i = 0; i < a.MonomialBasis.size(); ++i) {
                    if (a.MonomialBasis[i] != b.MonomialBasis[i])
                        return a.MonomialBasis[i] < b.MonomialBasis[i];
                }
                return false;  // truly equal
            });

        // 5. Build recursion + ring matrices per region
        cout << "\n[Step 4] Building recursion + ring matrices..." << endl;
        auto t4_0 = Clock::now();
        for (size_t ridx = 0; ridx < regions.size(); ++ridx) {
            auto& reg = regions[ridx];
            cout << "  Region " << ridx << ": nb=" << reg.nb << endl;
            auto recMats = RecursionBuilder::buildRecursionMatrices(ft, reg, fam.modulus);
            auto ringMats = RingBuilder::computeRingMatrices(reg, ne, fam.modulus);
            ibpRegimes.push_back(convertToRegimeSparse<FFInt>(recMats, fam.modulus));
            ringCells.push_back(convertToRingCell<FFInt>(ringMats, reg.limitSector));
        }
        auto t4_1 = Clock::now();
        cout << "  [TIMING] Step 4: "
             << chrono::duration_cast<Ms>(t4_1 - t4_0).count() / 1000.0 << " s" << endl;
    }

    // 6. Write output
    string ibpOut, ringOut;
    if (!sectorFilter.empty() && !cacheDir.empty()) {
        // Sector mode: write partial output with sector key in filename
        mkdirP(cacheDir);
        ibpOut  = cacheDir + "/IBPMat_"  + fam.name + "_sector_" + sectorFilter + ".bin";
        ringOut = cacheDir + "/RingData_" + fam.name + "_sector_" + sectorFilter + ".bin";
    } else {
        ibpOut  = outputDir + "/IBPMat_"  + fam.name + ".bin";
        ringOut = outputDir + "/RingData_" + fam.name + ".bin";
    }

    cout << "\n[Step 5] Writing .bin files..." << endl;
    auto t5_0 = Clock::now();
    writeIBPMatrixBinary<FFInt>(ibpOut, ibpRegimes);
    writeRingDataBinary<FFInt>(ringOut, ringCells, ne);
    auto t5_1 = Clock::now();
    cout << "  [TIMING] Step 5: "
         << chrono::duration_cast<Ms>(t5_1 - t5_0).count() / 1000.0 << " s" << endl;

    auto t_total_1 = Clock::now();
    double total_s = chrono::duration_cast<Ms>(t_total_1 - t_total_0).count() / 1000.0;
    cout << "\n[TIMING] Total wall-clock: " << total_s << " s" << endl;

    // 7. Optional diff
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
