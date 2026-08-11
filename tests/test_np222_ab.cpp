// test_np222_ab.cpp — Dump A/B eigen equations for NP222 sector 1111100
// and run RegionSolver to compare with RelationMeta_NP222.m regime A-values.
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include "../include/RegionSolver.hpp"
#include <iostream>
#include <cstdint>

int main() {
    int64_t mod = 179424673;
    auto fam = parseFamilyConfig("families/NP222.json");
    std::cout << "Family: " << fam.name << "  nProp=" << fam.nProp() << std::endl;

    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    std::cout << "IBP equations: " << ibp.equations.size()
              << "  ne=" << ibp.ne << "  nibp=" << ibp.nibp << std::endl;

    std::vector<int> sector = {1, 1, 1, 1, 1, 0, 0};
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, sector, mod);
    std::cout << "\n=== A/B equations (sector 1111100) ===" << std::endl;
    for (size_t i = 0; i < abEqs.size(); ++i)
        std::cout << "[" << i << "] " << abEqs[i] << std::endl;

    auto regions = RegionSolver::solveRegion(abEqs, sector, fam.nProp(), mod);
    std::cout << "\n=== Regions: " << regions.size() << " ===" << std::endl;
    for (size_t r = 0; r < regions.size(); ++r) {
        const auto& reg = regions[r];
        std::cout << "\nRegion " << r << ": nb=" << reg.nb << std::endl;
        std::cout << "  VarRule:" << std::endl;
        for (auto& [k, v] : reg.VarRule)
            std::cout << "    " << k << " -> " << v << std::endl;
        std::cout << "  MinPoly:" << std::endl;
        for (auto& p : reg.MinPoly)
            std::cout << "    " << p << std::endl;
    }
    return 0;
}
