// test_np222_true_ab.cpp — dump 匹配族 F1 的 top-sector AB 特征系统
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <iostream>

int main() {
    int64_t mod = 179424673;
    FamilyDef fam;
    fam.name = "NP222true";
    fam.propagators = {"(l1-l2)^2", "(l2+p)^2", "(l2-p)^2", "(l1-l2+p)^2", "(l1-l2-p)^2"};
    fam.loopMomenta = {"l1", "l2"};
    fam.externalMomenta = {"p"};
    fam.kinematicRules = {{"p^2", "s"}};
    fam.topSector = {1,1,1,1,1};
    fam.numeric = {{"s","1"},{"d","1/13"}};
    fam.modulus = mod;
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, fam.topSector, mod);
    for (size_t i = 0; i < abEqs.size(); ++i)
        std::cout << (i ? ",\n" : "") << abEqs[i];
    std::cout << std::endl;
    return 0;
}
