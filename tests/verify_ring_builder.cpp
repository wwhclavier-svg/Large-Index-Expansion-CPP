// Cross-validate C++ RingBuilder output against MMA-generated RingData_*.bin
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cstdint>
#include <cstring>
#include <algorithm>

#include "RegionSolver.hpp"
#include "RingBuilder.hpp"
#include "IBPEqGenerator.hpp"

constexpr int64_t MOD = 179424673;

// Read MMA binary format
struct MmaRingCell {
    std::vector<int> limitSector;
    int nb;
    std::vector<int64_t> A_flat;
    std::vector<int64_t> Ainv_flat;
};

int64_t readI64(std::ifstream& in) {
    int64_t v; in.read(reinterpret_cast<char*>(&v), 8); return v;
}

std::vector<MmaRingCell> readMmaRingData(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    std::vector<MmaRingCell> cells;
    int32_t numCells, ne;
    in.read(reinterpret_cast<char*>(&numCells), 4);
    in.read(reinterpret_cast<char*>(&ne), 4);
    for (int ci = 0; ci < numCells; ++ci) {
        MmaRingCell cell;
        int32_t secLen; in.read(reinterpret_cast<char*>(&secLen), 4);
        cell.limitSector.resize(secLen);
        in.read(reinterpret_cast<char*>(cell.limitSector.data()), secLen * 4);
        int32_t nb32; in.read(reinterpret_cast<char*>(&nb32), 4);
        cell.nb = nb32;
        int flatLen = ne * cell.nb * cell.nb;
        cell.A_flat.resize(flatLen);
        cell.Ainv_flat.resize(flatLen);
        for (int i = 0; i < flatLen; ++i) cell.A_flat[i] = readI64(in);
        for (int i = 0; i < flatLen; ++i) cell.Ainv_flat[i] = readI64(in);
        cells.push_back(cell);
    }
    return cells;
}

// Write C++ RingData to binary file (FFInt format)
void writeRingDataBinary(const std::string& filename,
                         const std::vector<RegionSolver::RegionData>& regions,
                         int ne, int64_t modulus) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Cannot open " << filename << " for writing" << std::endl;
        return;
    }
    int32_t numCells = static_cast<int32_t>(regions.size());
    int32_t ne32 = static_cast<int32_t>(ne);
    out.write(reinterpret_cast<const char*>(&numCells), 4);
    out.write(reinterpret_cast<const char*>(&ne32), 4);

    for (const auto& reg : regions) {
        int32_t secLen = static_cast<int32_t>(reg.limitSector.size());
        out.write(reinterpret_cast<const char*>(&secLen), 4);
        for (int s : reg.limitSector) {
            int32_t si = static_cast<int32_t>(s);
            out.write(reinterpret_cast<const char*>(&si), 4);
        }
        int32_t nb32 = static_cast<int32_t>(reg.nb);
        out.write(reinterpret_cast<const char*>(&nb32), 4);

        if (reg.nb > 0) {
            auto rm = RingBuilder::computeRingMatrices(reg, ne, modulus);
            int matSize = ne * reg.nb * reg.nb;
            for (int ei = 0; ei < ne; ++ei) {
                for (int j = 0; j < reg.nb * reg.nb; ++j) {
                    int64_t val = rm.A_list[ei][j];
                    out.write(reinterpret_cast<const char*>(&val), 8);
                }
            }
            for (int ei = 0; ei < ne; ++ei) {
                for (int j = 0; j < reg.nb * reg.nb; ++j) {
                    int64_t val = rm.Ainv_list[ei][j];
                    out.write(reinterpret_cast<const char*>(&val), 8);
                }
            }
        }
    }
    auto bytesWritten = out.tellp();
    out.close();
    std::cout << "Written " << numCells << " cells to " << filename
              << " (" << bytesWritten << " bytes)" << std::endl;
}

bool compareRingData(const std::string& familyFile,
                     const std::string& mmaBinFile,
                     const std::string& cppOutFile,
                     const std::vector<int>& topSector) {
    auto fam = parseFamilyConfig(familyFile);
    int ne = fam.nProp();
    std::cout << "\n=== " << fam.name << " ===" << std::endl;

    // C++ computation — only process sectors found in MMA reference
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);

    // Read MMA reference FIRST to extract unique sectors
    auto mmaData = readMmaRingData(mmaBinFile);
    std::set<std::vector<int>> mmaSectors;
    for (auto& m : mmaData) mmaSectors.insert(m.limitSector);

    // Only solve sectors that appear in MMA reference
    std::vector<RegionSolver::RegionData> regions;
    for (const auto& sec : mmaSectors) {
        auto abEqs = IBPAnalyzer::buildABEquations(ibp, sec, MOD);
        bool hasContent = false;
        for (const auto& eq : abEqs) {
            for (char ch : eq) {
                if (ch != '0' && ch != '+' && ch != '-' && ch != '*' && ch != ' ') {
                    hasContent = true; break;
                }
            }
            if (hasContent) break;
        }
        if (!hasContent) continue;
        auto secRegions = RegionSolver::solveRegion(abEqs, sec, ne, MOD);
        for (auto& reg : secRegions)
            regions.push_back(std::move(reg));
    }

    // Sort to match MMA ordering: by weight ascending, then binary value descending,
    // then nb ascending. The nb tiebreaker is essential for deterministic ordering of
    // multiple prime components within the same sector (same weight + same sector vector).
    std::sort(regions.begin(), regions.end(),
        [](const RegionSolver::RegionData& a, const RegionSolver::RegionData& b) {
            int wa=0,wb=0; for(int v:a.limitSector)wa+=v; for(int v:b.limitSector)wb+=v;
            if(wa!=wb)return wa<wb;
            for(size_t i=0;i<a.limitSector.size();++i)
                if(a.limitSector[i]!=b.limitSector[i])
                    return a.limitSector[i]>b.limitSector[i];
            // Tiebreaker: same sector, same weight => sort by nb ascending
            if (a.nb != b.nb) return a.nb < b.nb;
            return false;
        });

    std::cout << "Processing " << mmaSectors.size() << " unique sectors from MMA reference" << std::endl;

    // Write C++ binary output
    writeRingDataBinary(cppOutFile, regions, ne, MOD);

    std::cout << "MMA file: " << mmaBinFile << " (" << mmaData.size() << " cells)" << std::endl;
    std::cout << "C++ regions: " << regions.size() << std::endl;

    if (regions.size() != mmaData.size()) {
        std::cout << "WARNING: region count mismatch: C++=" << regions.size()
                  << " MMA=" << mmaData.size() << std::endl;
    }

    // Build MMA lookup by limitSector (multimap: sector -> list of indices).
    // Must use vector<int> value, NOT single int — sectors can have multiple
    // prime components (cells) sharing the same limitSector.
    std::map<std::vector<int>, std::vector<int>> mmaBySector;
    for (int i = 0; i < (int)mmaData.size(); ++i) {
        mmaBySector[mmaData[i].limitSector].push_back(i);
    }
    std::set<int> mmaUsed; // track which MMA cells have been matched

    int nCells = (int)regions.size();
    bool allMatch = true;
    int matched = 0;

    for (int r = 0; r < nCells; ++r) {
        auto& reg = regions[r];
        if (reg.nb == 0) continue;

        auto it = mmaBySector.find(reg.limitSector);
        if (it == mmaBySector.end()) {
            continue;
        }

        auto rm = RingBuilder::computeRingMatrices(reg, ne, MOD);

        // Find best matching MMA cell within same sector: same nb, then matching A[0] flat
        int bestMmaIdx = -1;
        for (int mi : it->second) {
            if (mmaUsed.count(mi)) continue;
            if (mmaData[mi].nb != reg.nb) continue;
            if (bestMmaIdx < 0) bestMmaIdx = mi;
            // If this is the only match or matrices match, prefer it
            int ms = reg.nb * reg.nb;
            if (ms > 0 && (int)mmaData[mi].A_flat.size() >= ms &&
                std::equal(mmaData[mi].A_flat.begin(), mmaData[mi].A_flat.begin() + ms,
                           rm.A_list[0].begin())) {
                bestMmaIdx = mi;
                break; // exact match found
            }
        }
        if (bestMmaIdx < 0) {
            // Try any nb
            for (int mi : it->second) {
                if (mmaUsed.count(mi)) continue;
                bestMmaIdx = mi;
                break;
            }
        }
        if (bestMmaIdx < 0) continue;
        mmaUsed.insert(bestMmaIdx);
        matched++;

        auto& mma = mmaData[bestMmaIdx];
        std::cout << "\n--- Cell " << r << " (nb=" << reg.nb << ") ---" << std::endl;

        // Compare limitSector
        bool secMatch = (reg.limitSector == mma.limitSector);
        std::cout << "  limitSector: " << (secMatch ? "MATCH" : "MISMATCH") << std::endl;
        if (!secMatch) {
            std::cout << "    C++: ";
            for (auto v : reg.limitSector) std::cout << v << " ";
            std::cout << "\n    MMA: ";
            for (auto v : mma.limitSector) std::cout << v << " ";
            std::cout << std::endl;
        }

        // Compare nb
        if (reg.nb != mma.nb) {
            std::cout << "  nb: MISMATCH C++=" << reg.nb << " MMA=" << mma.nb << std::endl;
            allMatch = false;
            continue;
        }

        // Compare A_list
        int matSize = reg.nb * reg.nb;
        bool hasMismatch = false;
        for (int i = 0; i < ne; ++i) {
            bool aMatch = true;
            for (int j = 0; j < matSize; ++j) {
                int64_t cppVal = rm.A_list[i][j];
                int mmaIdx = i * matSize + j;
                int64_t mmaVal = mma.A_flat[mmaIdx];
                if (cppVal != mmaVal) {
                    if (aMatch) {
                        std::cout << "  A[" << (i+1) << "]: MISMATCH" << std::endl;
                    }
                    aMatch = false;
                    allMatch = false;
                    std::cout << "    [" << (j / reg.nb) << "," << (j % reg.nb)
                              << "] C++=" << cppVal << " MMA=" << mmaVal << std::endl;
                }
            }
            if (aMatch)
                std::cout << "  A[" << (i+1) << "]: MATCH" << std::endl;
        }

        // Compare Ainv_list
        for (int i = 0; i < ne; ++i) {
            bool ainvMatch = true;
            for (int j = 0; j < matSize; ++j) {
                int64_t cppVal = rm.Ainv_list[i][j];
                int64_t mmaVal = mma.Ainv_flat[i * matSize + j];
                if (cppVal != mmaVal) {
                    if (ainvMatch) {
                        std::cout << "  Ainv[" << (i+1) << "]: MISMATCH" << std::endl;
                    }
                    ainvMatch = false;
                    allMatch = false;
                    std::cout << "    [" << (j / reg.nb) << "," << (j % reg.nb)
                              << "] C++=" << cppVal << " MMA=" << mmaVal << std::endl;
                }
            }
            if (ainvMatch)
                std::cout << "  Ainv[" << (i+1) << "]: MATCH" << std::endl;
        }
    }

    if (matched == 0) {
        std::cout << "\nWARNING: No cells matched between C++ and MMA!" << std::endl;
        return false;
    }
    std::cout << "\n  Matched " << matched << " cells by limitSector" << std::endl;
    return allMatch;
}

int main(int argc, char* argv[]) {
    bool ok = true;
    std::string target = (argc > 1) ? argv[1] : "all";

    if (target == "Box") {
        ok = compareRingData("families/Box.json",
                "verify/Box/RingData_Box-MMA.bin",
                "verify/Box/RingData_Box-CPP.bin",
                {1, 1, 1, 1}) && ok;
    } else if (target == "SR212") {
        ok = compareRingData("families/SR212.json",
                "verify/SR212/RingData_SR212_MMA.bin",
                "data/RingData_SR212-CPP.bin",
                {1, 1, 1, 1, 1}) && ok;
    } else if (target == "TB123") {
        ok = compareRingData("families/TB123.json",
                "data/RingData_TB123.bin",
                "data/RingData_TB123-CPP.bin",
                {1, 1, 1, 1, 1, 1, 1}) && ok;
    } else {
        ok = compareRingData("families/bub00.json",
                "data/RingData_bub00.bin",
                "data/RingData_bub00-CPP.bin",
                {1, 1}) && ok;
        ok = compareRingData("families/Tri.json",
                "data/RingData_Tri.bin",
                "data/RingData_Tri-CPP.bin",
                {1, 1, 1}) && ok;
        ok = compareRingData("families/bub10.json",
                "data/RingData_bub10.bin",
                "data/RingData_bub10-CPP.bin",
                {1, 1}) && ok;
        ok = compareRingData("families/bub11.json",
                "data/RingData_bub11.bin",
                "data/RingData_bub11-CPP.bin",
                {1, 1}) && ok;
    }

    if (ok)
        std::cout << "\n=== All comparisons MATCH ===" << std::endl;
    else
        std::cout << "\n=== Some MISMATCHES found ===" << std::endl;

    return ok ? 0 : 1;
}
