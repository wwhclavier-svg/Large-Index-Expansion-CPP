// Test RingBuilder: compute A_i and Ainv_i matrices from RegionData
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include "../include/RegionSolver.hpp"
#include "../include/RecursionBuilder.hpp"
#include "../include/RingBuilder.hpp"
#include <iostream>
#include <cassert>
#include <cstdint>

static const int64_t MOD = 179424673;

void testFamily(const std::string& familyFile, const std::vector<int>& topSector) {
    auto fam = parseFamilyConfig(familyFile);
    int ne = fam.nProp();
    std::cout << "\n=== " << fam.name << " (ne=" << ne << ") ===" << std::endl;

    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, topSector, MOD);
    auto ft = IBPAnalyzer::extractFTable(ibp, MOD);
    auto regions = RegionSolver::solveRegion(abEqs, topSector, ne, MOD);

    std::cout << "Regions: " << regions.size() << std::endl;
    assert(!regions.empty());

    int totalOK = 0, totalFail = 0;

    for (size_t r = 0; r < regions.size(); ++r) {
        auto& reg = regions[r];
        if (reg.nb == 0) continue;

        std::cout << "\n--- Region " << r << " (nb=" << reg.nb
                  << ", VarIndep=" << reg.VarIndep.size() << ") ---" << std::endl;

        // Print VarRule for A_i
        std::cout << "VarRule for A_i:" << std::endl;
        for (int i = 1; i <= ne; ++i) {
            std::string key = "A" + std::to_string(i);
            auto it = reg.VarRule.find(key);
            if (it != reg.VarRule.end())
                std::cout << "  " << key << " = " << it->second << std::endl;
            else {
                // Check if in VarIndep
                bool inIndep = false;
                for (auto& v : reg.VarIndep) if (v == key) { inIndep = true; break; }
                if (inIndep)
                    std::cout << "  " << key << " = <generator>" << std::endl;
                else
                    std::cout << "  " << key << " = <not found>" << std::endl;
            }
        }

        // Compute ring matrices
        auto rm = RingBuilder::computeRingMatrices(reg, ne, MOD);

        // Verify A_i * Ainv_i = I for each propagator
        for (int i = 0; i < ne; ++i) {
            bool ok = true;
            for (int rrow = 0; rrow < reg.nb && ok; ++rrow) {
                for (int c = 0; c < reg.nb && ok; ++c) {
                    int64_t sum = 0;
                    for (int k = 0; k < reg.nb; ++k) {
                        sum = (sum + RecursionBuilder::matAt(rm.A_list[i], reg.nb, rrow, k)
                                   * RecursionBuilder::matAt(rm.Ainv_list[i], reg.nb, k, c)) % MOD;
                    }
                    if (sum < 0) sum += MOD;
                    int64_t expected = (rrow == c) ? 1 : 0;
                    if (sum != expected) {
                        std::cout << "  FAIL: A[" << (i+1) << "] * Ainv[" << (i+1)
                                  << "][" << rrow << "][" << c << "] = "
                                  << sum << " != " << expected << std::endl;
                        ok = false;
                    }
                }
            }
            if (ok) {
                std::cout << "  A[" << (i+1) << "] * Ainv[" << (i+1) << "] = I: PASS" << std::endl;
                ++totalOK;
            } else {
                ++totalFail;
            }
        }

        // Print A_i matrix for debugging
        for (int i = 0; i < ne; ++i) {
            std::cout << "  A[" << (i+1) << "] = [";
            for (int rrow = 0; rrow < reg.nb; ++rrow) {
                if (rrow > 0) std::cout << "; ";
                for (int c = 0; c < reg.nb; ++c) {
                    if (c > 0) std::cout << ",";
                    std::cout << RecursionBuilder::matAt(rm.A_list[i], reg.nb, rrow, c);
                }
            }
            std::cout << "]" << std::endl;
        }
    }

    std::cout << "\n  Summary: " << totalOK << " OK, " << totalFail << " FAIL" << std::endl;
    assert(totalFail == 0);
}

int main() {
    // Test 1: bub00 — nb=1, VarIndep empty, all A_i are constants
    testFamily("families/bub00.json", {1, 1});

    // Test 2: Tri — nb=2, VarIndep={A3}, A1/A2 are constants, A3 is generator
    testFamily("families/Tri.json", {1, 1, 1});

    std::cout << "\n=== All RingBuilder tests PASSED ===" << std::endl;
    return 0;
}
