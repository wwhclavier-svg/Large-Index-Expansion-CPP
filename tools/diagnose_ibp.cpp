// Diagnostic: dump C++ IBP + A/B equations for bub00 and compare with MMA
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <iostream>

int main() {
    auto fam = parseFamilyConfig("families/bub00.json");
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);

    std::cout << "\n=== C++ IBP Equations (g-operator form) ===" << std::endl;
    std::cout << "ne=" << ibp.ne << " nibp=" << ibp.nibp << std::endl;

    for (size_t m = 0; m < ibp.equations.size(); ++m) {
        auto& eq = ibp.equations[m];
        std::cout << "\nEq[" << m << "] (j=" << eq.j << ", k=" << eq.k << "):";
        for (auto& t : eq.terms) {
            std::cout << "\n  g[";
            for (int i = 0; i < ibp.ne; ++i) {
                if (i > 0) std::cout << ",";
                std::cout << t.gShift[i];
            }
            std::cout << "] * " << t.coeff;
            if (t.hasD) std::cout << " * d";
            if (t.nIdx > 0) std::cout << " * n" << t.nIdx;
        }
    }
    std::cout << "\n" << std::endl;

    // A/B equations for top sector
    std::vector<int> topSector = {1, 1};
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, topSector, fam.modulus);

    std::cout << "=== C++ A/B Equations (top sector) ===" << std::endl;
    std::cout << "modulus=" << fam.modulus << std::endl;
    for (size_t i = 0; i < abEqs.size(); ++i) {
        std::cout << "AB[" << i << "]: " << abEqs[i] << std::endl;
    }

    return 0;
}
