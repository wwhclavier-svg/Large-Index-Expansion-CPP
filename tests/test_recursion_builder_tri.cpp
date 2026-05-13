// Test RecursionBuilder with Tri family (nb=2, VarIndep non-empty)
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

    auto fam = parseFamilyConfig("families/Tri.json");
    std::cout << "Family: " << fam.name << " (ne=" << fam.nProp() << ")" << std::endl;

    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    std::cout << "IBP equations: " << ibp.equations.size() << std::endl;

    std::vector<int> topSector = {1, 1, 1};
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, topSector, mod);

    auto ft = IBPAnalyzer::extractFTable(ibp, mod);

    // Print FTable details for debugging
    std::cout << "\nF0: ";
    for (int m = 0; m < ft.nibp; ++m) std::cout << ft.F0[m] << " ";
    std::cout << "\nF1:";
    for (int m = 0; m < ft.nibp; ++m) {
        std::cout << "\n  m=" << m << ": ";
        for (int i = 0; i < ft.ne; ++i) std::cout << ft.F1[m][i] << " ";
    }
    std::cout << "\nF2D:";
    for (int m = 0; m < ft.nibp; ++m) {
        std::cout << "\n  m=" << m << ": ";
        for (int i = 0; i < ft.ne; ++i) std::cout << ft.F2D[m][i] << " ";
    }
    std::cout << "\nf2 non-zero:";
    for (int m = 0; m < ft.nibp; ++m)
        for (int i = 0; i < ft.ne; ++i)
            for (int j = 0; j < ft.ne; ++j)
                if (ft.f2[m][i][j] != 0)
                    std::cout << "\n  f2[" << m << "][" << i << "][" << j << "]=" << ft.f2[m][i][j];
    std::cout << std::endl;

    auto regions = RegionSolver::solveRegion(abEqs, topSector, fam.nProp(), mod);
    std::cout << "\nFound " << regions.size() << " regions" << std::endl;
    assert(!regions.empty());

    for (size_t r = 0; r < regions.size(); ++r) {
        auto& reg = regions[r];
        std::cout << "\n=== Region " << r << " ===" << std::endl;
        std::cout << "  nb=" << reg.nb << " VarIndep: ";
        for (auto& v : reg.VarIndep) std::cout << v << " ";
        std::cout << "\n  VarDep: ";
        for (auto& v : reg.VarDep) std::cout << v << " ";
        std::cout << "\n  VarRule:";
        for (auto& kv : reg.VarRule)
            std::cout << "\n    " << kv.first << " = " << kv.second;
        std::cout << "\n  MonomialBasisIndex:";
        for (auto& idx : reg.MonomialBasisIndex) {
            std::cout << " [";
            for (size_t i = 0; i < idx.size(); ++i) {
                if (i > 0) std::cout << ",";
                std::cout << idx[i];
            }
            std::cout << "]";
        }
        std::cout << "\n  MonomialBasisMatrix: nb=" << reg.nb
                  << " dims=" << reg.MonomialBasisMatrix.size();
        if (!reg.MonomialBasisMatrix.empty()) {
            std::cout << "x" << reg.MonomialBasisMatrix[0].size()
                      << "x" << reg.MonomialBasisMatrix[0][0].size();
        }
        std::cout << std::endl;

        if (reg.nb == 0) continue;

        // Print MonomialBasisMatrix entries
        for (int k = 0; k < (int)reg.MonomialBasisMatrix.size(); ++k) {
            std::cout << "  MBM[" << k << "]:";
            for (int i = 0; i < reg.nb; ++i) {
                std::cout << " [";
                for (int j = 0; j < reg.nb; ++j) {
                    if (j > 0) std::cout << ",";
                    std::cout << reg.MonomialBasisMatrix[k][i][j];
                }
                std::cout << "]";
            }
            std::cout << std::endl;
        }

        auto rm = RecursionBuilder::buildRecursionMatrices(ft, reg, mod);
        RecursionBuilder::printRecursionStats(rm);

        // Print specific K1 and N1 entries
        std::cout << "  K1[0][0] flat:";
        for (auto v : rm.K1[0][0]) std::cout << " " << v;
        std::cout << "\n  K1[0][1] flat:";
        for (auto v : rm.K1[0][1]) std::cout << " " << v;
        std::cout << "\n  K1[0][2] flat:";
        for (auto v : rm.K1[0][2]) std::cout << " " << v;
        std::cout << "\n  N1[0][0] flat:";
        for (auto v : rm.N1[0][0]) std::cout << " " << v;
        std::cout << "\n  N1[0][1] flat:";
        for (auto v : rm.N1[0][1]) std::cout << " " << v;
        std::cout << "\n  N1[0][2] flat:";
        for (auto v : rm.N1[0][2]) std::cout << " " << v;

        // Verify consistency
        bool ok = true;
        for (int m = 0; m < rm.nibp && ok; ++m) {
            for (int i = 0; i < rm.ne && ok; ++i) {
                for (int rr = 0; rr < reg.nb; ++rr) {
                    for (int c = 0; c < reg.nb; ++c) {
                        int64_t m1 = RecursionBuilder::matAt(rm.M1[m][i], reg.nb, rr, c);
                        int64_t k1s = RecursionBuilder::matAt(rm.K1s[m][i], reg.nb, rr, c);
                        int64_t k2s = RecursionBuilder::matAt(rm.K2s[m][i], reg.nb, rr, c);
                        int64_t check = (m1 + k2s) % mod;
                        if (check < 0) check += mod;
                        if (check != k1s) {
                            std::cout << "\n  FAIL M1+K2s!=K1s at m=" << m << " i=" << i
                                      << "[" << rr << "][" << c << "]" << std::endl;
                            ok = false;
                        }
                    }
                }
            }
        }
        if (ok) std::cout << "\n  M1+K2s==K1s: PASS" << std::endl;
    }

    std::cout << "\n=== DONE ===" << std::endl;
    return 0;
}
