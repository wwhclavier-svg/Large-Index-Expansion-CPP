// test_sr212_ab.cpp — Dump SR212 top-sector 11111 的 AB 特征系统
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <iostream>
#include <cstdint>

int main() {
    int64_t mod = 179424673;
    auto fam = parseFamilyConfig("families/SR212.json");
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    std::vector<int> sector = {1, 1, 1, 1, 1};
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, sector, mod);
    for (size_t i = 0; i < abEqs.size(); ++i)
        std::cout << (i ? ",\n" : "") << abEqs[i];
    std::cout << std::endl;
    return 0;
}
