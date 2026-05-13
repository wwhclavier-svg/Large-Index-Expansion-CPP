// Integration test: IBPAnalyzer → RegionSolver for bub00
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include "../include/RegionSolver.hpp"
#include <iostream>
#include <cassert>
#include <cstdint>

int main() {
    int64_t mod = 179424673;
    auto fam = parseFamilyConfig("families/bub00.json");
    std::cout << "Family: " << fam.name << std::endl;

    // Generate IBP equations
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    std::cout << "IBP equations: " << ibp.equations.size() << std::endl;

    // Build A/B equations for top sector [1,1]
    std::vector<int> topSector = {1, 1};
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, topSector, mod);
    std::cout << "\nA/B equations:" << std::endl;
    for (size_t i = 0; i < abEqs.size(); ++i)
        std::cout << "  [" << i << "] " << abEqs[i] << std::endl;

    // Run RegionSolver
    auto regions = RegionSolver::solveRegion(abEqs, topSector, fam.nProp(), mod);

    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Found " << regions.size() << " regions" << std::endl;
    for (size_t r = 0; r < regions.size(); ++r) {
        const auto& reg = regions[r];
        std::cout << "\nRegion " << r << ":" << std::endl;
        std::cout << "  nb = " << reg.nb << std::endl;
        std::cout << "  VarIndep: ";
        for (auto& v : reg.VarIndep) std::cout << v << " ";
        std::cout << std::endl;
        std::cout << "  VarDep: ";
        for (auto& v : reg.VarDep) std::cout << v << " ";
        std::cout << std::endl;
        std::cout << "  VarDeg: ";
        for (int d : reg.VarDeg) std::cout << d << " ";
        std::cout << std::endl;

        // Print VarRule
        if (!reg.VarRule.empty()) {
            std::cout << "  VarRule:" << std::endl;
            for (auto& [k, v] : reg.VarRule)
                std::cout << "    " << k << " -> " << v << std::endl;
        }
        // Print FractionRule
        if (!reg.FractionRule.empty()) {
            std::cout << "  FractionRule (" << reg.FractionRule.size() << " entries):" << std::endl;
            for (auto& [k, v] : reg.FractionRule)
                std::cout << "    " << k << " -> " << v << std::endl;
        }
        // Print MonomialBasis
        if (!reg.MonomialBasis.empty()) {
            std::cout << "  MonomialBasis: ";
            for (auto& b : reg.MonomialBasis) std::cout << b << "  ";
            std::cout << std::endl;
        }
        // Print MonomialBasisMatrix size
        if (!reg.MonomialBasisMatrix.empty()) {
            std::cout << "  MonomialBasisMatrix: " << reg.MonomialBasisMatrix.size()
                      << " x " << (reg.MonomialBasisMatrix.empty() ? 0 : reg.MonomialBasisMatrix[0].size())
                      << " x " << (reg.MonomialBasisMatrix.empty() || reg.MonomialBasisMatrix[0].empty() ? 0 : reg.MonomialBasisMatrix[0][0].size())
                      << std::endl;
        }
    }

    std::cout << "\nDone." << std::endl;
    return 0;
}
