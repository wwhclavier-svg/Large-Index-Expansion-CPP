// dump_ibp.cpp — dump IBP equations + AB equations for any family
// Usage: dump_ibp <family.json> [sector like 1111110]
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <iostream>
#include <cstdint>

using namespace IBPEqGenerator;

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "Usage: dump_ibp <family.json> [sector]\n"; return 1; }
    int64_t mod = 179424673;
    auto fam = parseFamilyConfig(argv[1]);
    std::vector<int> sector = fam.topSector;
    if (argc > 2) {
        sector.clear();
        for (char c : std::string(argv[2])) sector.push_back(c - '0');
    }
    auto ibp = generateIBPEquations(fam);
    std::cout << "Family: " << fam.name << " nibp=" << ibp.nibp << " ne=" << ibp.ne << std::endl;
    for (size_t m = 0; m < ibp.equations.size(); ++m) {
        const auto& eq = ibp.equations[m];
        std::cout << "Eq[" << m << "] (j=" << eq.j << ",k=" << eq.k << "): "
                  << eq.terms.size() << " terms" << std::endl;
        for (const auto& t : eq.terms) {
            std::cout << "  g[";
            for (size_t i = 0; i < t.gShift.size(); ++i) {
                if (i) std::cout << ",";
                std::cout << t.gShift[i];
            }
            std::cout << "] * " << t.coeff;
            if (t.nIdx > 0) std::cout << " * n" << t.nIdx;
            if (t.hasD) std::cout << " * d";
            std::cout << std::endl;
        }
    }
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, sector, mod);
    std::cout << "\n=== A/B Equations (sector ";
    for (int b : sector) std::cout << b;
    std::cout << ") ===" << std::endl;
    for (size_t i = 0; i < abEqs.size(); ++i)
        std::cout << "AB[" << i << "] " << abEqs[i] << std::endl;
    return 0;
}
