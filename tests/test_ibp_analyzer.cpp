// Test IBPAnalyzer: A/B conversion and F-table extraction
// Compares with known output for bub00 family
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <iostream>
#include <cassert>
#include <cstdint>

using namespace IBPEqGenerator;

int main() {
    int64_t mod = 179424673;

    // Load bub00 family
    auto fam = parseFamilyConfig("families/bub00.json");
    std::cout << "Family: " << fam.name << " (ne=" << fam.nProp()
              << ", nl=" << fam.nLoop() << ", nE=" << fam.nExt() << ")" << std::endl;

    // Generate IBP equations
    auto ibp = generateIBPEquations(fam);
    std::cout << "Generated " << ibp.equations.size() << " IBP equations" << std::endl;
    std::cout << "nibp=" << ibp.nibp << " ne=" << ibp.ne << std::endl;

    // Print one equation for inspection
    for (size_t m = 0; m < ibp.equations.size(); ++m) {
        const auto& eq = ibp.equations[m];
        std::cout << "Eq[" << m << "] (j=" << eq.j << ",k=" << eq.k << "): "
                  << eq.terms.size() << " terms" << std::endl;
        for (const auto& t : eq.terms) {
            std::cout << "  g[";
            for (size_t i = 0; i < t.gShift.size(); ++i) {
                if (i > 0) std::cout << ",";
                std::cout << t.gShift[i];
            }
            std::cout << "] * " << t.coeff;
            if (t.nIdx > 0) std::cout << " * n" << t.nIdx;
            if (t.hasD) std::cout << " * d";
            std::cout << std::endl;
        }
    }

    // Test buildABEquations with top sector [1,1]
    std::vector<int> topSector = {1, 1};
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, topSector, mod);
    std::cout << "\n=== A/B Equations (top sector) ===" << std::endl;
    for (size_t i = 0; i < abEqs.size(); ++i) {
        std::cout << "  [" << i << "] " << abEqs[i] << std::endl;
    }

    // Test extractFTable
    auto ft = IBPAnalyzer::extractFTable(ibp, mod);
    IBPAnalyzer::printFTableSummary(ft);

    // Print F0 values
    std::cout << "\nF0: ";
    for (int m = 0; m < ft.nibp; ++m) std::cout << ft.F0[m] << " ";
    std::cout << std::endl;

    // Print F1 values
    std::cout << "F1:" << std::endl;
    for (int m = 0; m < ft.nibp; ++m) {
        std::cout << "  [" << m << "] ";
        for (int i = 0; i < ft.ne; ++i) std::cout << ft.F1[m][i] << " ";
        std::cout << std::endl;
    }

    // Print F2D values
    std::cout << "F2D:" << std::endl;
    for (int m = 0; m < ft.nibp; ++m) {
        std::cout << "  [" << m << "] ";
        for (int i = 0; i < ft.ne; ++i) std::cout << ft.F2D[m][i] << " ";
        std::cout << std::endl;
    }

    // Print f2 values
    std::cout << "f2:" << std::endl;
    for (int m = 0; m < ft.nibp; ++m) {
        for (int i = 0; i < ft.ne; ++i) {
            for (int j = 0; j < ft.ne; ++j) {
                if (ft.f2[m][i][j] != 0)
                    std::cout << "  [" << m << "][" << i << "][" << j << "] = "
                              << ft.f2[m][i][j] << std::endl;
            }
        }
    }

    std::cout << "\nAll IBPAnalyzer tests completed." << std::endl;
    return 0;
}
