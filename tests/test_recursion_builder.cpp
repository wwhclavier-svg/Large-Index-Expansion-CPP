// Integration test: IBPAnalyzer → RegionSolver → RecursionBuilder for bub00
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include "../include/RegionSolver.hpp"
#include "../include/RecursionBuilder.hpp"
#include <iostream>
#include <cassert>
#include <cstdint>

int main() {
    int64_t mod = 179424673;

    auto fam = parseFamilyConfig("families/bub00.json");
    std::cout << "Family: " << fam.name << " (ne=" << fam.nProp() << ")" << std::endl;

    // Step 1: Generate IBP equations
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    std::cout << "IBP equations: " << ibp.equations.size() << std::endl;

    // Step 2: Build A/B equations for top sector [1,1]
    std::vector<int> topSector = {1, 1};
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, topSector, mod);
    std::cout << "\nA/B equations:" << std::endl;
    for (size_t i = 0; i < abEqs.size(); ++i)
        std::cout << "  [" << i << "] " << abEqs[i] << std::endl;

    // Step 3: Extract FTable (v=all-1 basis)
    auto ft = IBPAnalyzer::extractFTable(ibp, mod);
    IBPAnalyzer::printFTableSummary(ft);

    // Step 4: Run RegionSolver
    auto regions = RegionSolver::solveRegion(abEqs, topSector, fam.nProp(), mod);
    std::cout << "\nFound " << regions.size() << " regions" << std::endl;

    assert(!regions.empty());

    // Step 5: Run RecursionBuilder for each region
    for (size_t r = 0; r < regions.size(); ++r) {
        auto& reg = regions[r];
        std::cout << "\n=== Region " << r << " ===" << std::endl;
        std::cout << "  nb=" << reg.nb << " VarIndep=" << reg.VarIndep.size() << std::endl;
        std::cout << "  VarDep: ";
        for (auto& v : reg.VarDep) std::cout << v << " ";
        std::cout << std::endl;
        std::cout << "  VarRule:" << std::endl;
        for (auto& kv : reg.VarRule)
            std::cout << "    " << kv.first << " = " << kv.second << std::endl;
        std::cout << "  MonomialBasis: ";
        for (auto& mb : reg.MonomialBasis) std::cout << mb << " ";
        std::cout << std::endl;

        // Build recursion matrices
        auto rm = RecursionBuilder::buildRecursionMatrices(ft, reg, mod);
        RecursionBuilder::printRecursionStats(rm);

        // Verify key properties
        int nb = reg.nb;

        // For nb=1 (VarIndep=∅), each matrix is 1×1 = a scalar.
        // Verify M1 + K2s = K1s (should hold since M1 = K1s - K2s)
        for (int m = 0; m < rm.nibp; ++m) {
            for (int i = 0; i < rm.ne; ++i) {
                int64_t m1 = rm.M1[m][i][0];
                int64_t k1s = rm.K1s[m][i][0];
                int64_t k2s = rm.K2s[m][i][0];
                int64_t check = (m1 + k2s) % mod;
                if (check < 0) check += mod;
                if (check != k1s) {
                    std::cout << "  FAIL: M1+K2s != K1s at m=" << m << " i=" << i
                              << ": " << m1 << "+" << k2s << "=" << check
                              << " vs " << k1s << std::endl;
                    return 1;
                }
            }
        }

        // Verify N1 = F2D + K1 (already built into construction, but cross-check)
        for (int m = 0; m < rm.nibp; ++m) {
            for (int i = 0; i < rm.ne; ++i) {
                int64_t n1 = rm.N1[m][i][0];
                int64_t k1 = rm.K1[m][i][0];
                int64_t f2d_i = ft.F2D[m][i];
                int64_t expected = (f2d_i + k1) % mod;
                if (expected < 0) expected += mod;
                if (n1 != expected) {
                    std::cout << "  FAIL: N1 != F2D+K1 at m=" << m << " i=" << i
                              << ": " << n1 << " vs " << expected << std::endl;
                    return 1;
                }
            }
        }

        std::cout << "  All consistency checks PASSED" << std::endl;
    }

    std::cout << "\n=== All RecursionBuilder tests PASSED ===" << std::endl;
    return 0;
}
