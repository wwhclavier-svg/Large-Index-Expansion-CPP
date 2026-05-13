// debug_region_solver.cpp — Step-by-step RegionSolver tracing for bub00
// Focus: trace solveRegion pipeline to find A_i constant mismatch root cause
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include "../include/RegionSolver.hpp"
#include "../include/RecursionBuilder.hpp"
#include "../include/RingBuilder.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdint>

static const int64_t MOD = 179424673;

// Dump a polynomial vector with label
void dumpPolys(const std::string& label, const std::vector<std::string>& polys) {
    std::cout << "\n--- " << label << " (" << polys.size() << " polynomials) ---" << std::endl;
    for (size_t i = 0; i < polys.size(); ++i) {
        std::cout << "  [" << i << "] " << polys[i] << std::endl;
    }
}

// Dump RegionData
void dumpRegion(int idx, const RegionSolver::RegionData& reg) {
    std::cout << "\n=== Region " << idx << " ===" << std::endl;
    std::cout << "  limitSector: ";
    for (int v : reg.limitSector) std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "  nb=" << reg.nb << std::endl;
    std::cout << "  VarIndep: ";
    for (auto& v : reg.VarIndep) std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "  VarDep: ";
    for (auto& v : reg.VarDep) std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "  VarDeg: ";
    for (int d : reg.VarDeg) std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "  MinPoly (" << reg.MinPoly.size() << "):" << std::endl;
    for (auto& p : reg.MinPoly) std::cout << "    " << p << std::endl;
    std::cout << "  MonomialBasis: ";
    for (auto& b : reg.MonomialBasis) std::cout << b << " ";
    std::cout << std::endl;
    std::cout << "  MonomialBasisIndex:" << std::endl;
    for (size_t i = 0; i < reg.MonomialBasisIndex.size(); ++i) {
        std::cout << "    [" << i << "] ";
        for (int e : reg.MonomialBasisIndex[i]) std::cout << e << ",";
        std::cout << std::endl;
    }
    std::cout << "  VarRule (" << reg.VarRule.size() << " entries):" << std::endl;
    for (auto& [k, v] : reg.VarRule) {
        std::cout << "    " << k << " -> " << v << std::endl;
    }
    std::cout << "  FractionRule (" << reg.FractionRule.size() << " entries):" << std::endl;
    for (auto& [k, v] : reg.FractionRule) {
        std::cout << "    " << k << " -> " << v << std::endl;
    }
    if (!reg.MonomialBasisMatrix.empty()) {
        std::cout << "  MonomialBasisMatrix: dims=[" << reg.MonomialBasisMatrix.size()
                  << "][" << (reg.MonomialBasisMatrix[0].empty() ? 0 : reg.MonomialBasisMatrix[0].size())
                  << "][" << (reg.MonomialBasisMatrix[0].empty() ? 0 : reg.MonomialBasisMatrix[0][0].size())
                  << "]" << std::endl;
        for (size_t k = 0; k < std::min((size_t)3, reg.MonomialBasisMatrix.size()); ++k) {
            std::cout << "    MBM[" << k << "][0][*]: ";
            for (size_t j = 0; j < reg.MonomialBasisMatrix[k][0].size(); ++j)
                std::cout << reg.MonomialBasisMatrix[k][0][j] << " ";
            std::cout << std::endl;
        }
    }
}

// Read MMA RingData for comparison
struct MmaRingData {
    std::vector<int> limitSector;
    int nb;
    std::vector<int64_t> A_flat;
    std::vector<int64_t> Ainv_flat;
};

std::vector<MmaRingData> readMmaRingData(const std::string& filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + filename);
    auto readI32 = [&]() -> int32_t { int32_t v; f.read(reinterpret_cast<char*>(&v), 4); return v; };
    auto readI64 = [&]() -> int64_t { int64_t v; f.read(reinterpret_cast<char*>(&v), 8); return v; };
    int numCells = readI32();
    int globalNe = readI32();
    std::vector<MmaRingData> cells;
    for (int c = 0; c < numCells; ++c) {
        MmaRingData cell;
        int secLen = readI32();
        cell.limitSector.resize(secLen);
        for (int i = 0; i < secLen; ++i) cell.limitSector[i] = readI32();
        cell.nb = readI32();
        int matSize = cell.nb * cell.nb;
        int flatLen = globalNe * matSize;
        cell.A_flat.resize(flatLen);
        for (int i = 0; i < flatLen; ++i) cell.A_flat[i] = readI64();
        cell.Ainv_flat.resize(flatLen);
        for (int i = 0; i < flatLen; ++i) cell.Ainv_flat[i] = readI64();
        cells.push_back(cell);
    }
    return cells;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "RegionSolver Debug Trace — bub00" << std::endl;
    std::cout << "========================================" << std::endl;

    // Step 1: Load family
    auto fam = parseFamilyConfig("families/bub00.json");
    std::cout << "\nFamily: " << fam.name << " ne=" << fam.nProp() << std::endl;

    // Step 2: Generate IBP equations
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    std::vector<int> topSector = {1, 1};

    // Step 3: Build A/B equations
    auto abEqs = IBPAnalyzer::buildABEquations(ibp, topSector, MOD);
    dumpPolys("A/B Equations", abEqs);

    // Step 4: Trace solveRegion step by step
    std::cout << "\n========================================" << std::endl;
    std::cout << "Step-by-step solveRegion" << std::endl;
    std::cout << "========================================" << std::endl;

    int ne = fam.nProp();
    std::vector<std::string> allVars;
    for (int i = 1; i <= ne; ++i) allVars.push_back("B" + std::to_string(i));
    for (int i = 1; i <= ne; ++i) allVars.push_back("A" + std::to_string(i));
    std::cout << "\nallVars: ";
    for (auto& v : allVars) std::cout << v << " ";
    std::cout << std::endl;

    // Add A_i*B_i-1 constraints
    std::vector<std::string> fullIdeal = abEqs;
    for (int i = 1; i <= ne; ++i)
        fullIdeal.push_back("A" + std::to_string(i) + "*B" + std::to_string(i) + "-1");
    dumpPolys("Full Ideal (with A_i*B_i-1)", fullIdeal);

    // Groebner basis
    auto gb = RegionSolver::computeGroebnerBasis(fullIdeal, allVars, MOD);
    dumpPolys("Groebner Basis (full)", gb);

    // Minimal associated primes
    auto [primelist, dims] = RegionSolver::computeMinAssPrimes(gb, allVars, MOD);
    std::cout << "\n--- Minimal Associated Primes ---" << std::endl;
    std::cout << "Number of components: " << primelist.size() << std::endl;
    for (size_t p = 0; p < primelist.size(); ++p) {
        std::cout << "\n  Component " << p << " (dim=" << dims[p] << "):" << std::endl;
        for (size_t i = 0; i < primelist[p].size(); ++i) {
            std::cout << "    [" << i << "] " << primelist[p][i] << std::endl;
        }
    }

    // Per-component analysis
    std::vector<std::string> Avars, Bvars;
    for (int i = 1; i <= ne; ++i) {
        Bvars.push_back("B" + std::to_string(i));
        Avars.push_back("A" + std::to_string(i));
    }

    for (size_t p = 0; p < primelist.size(); ++p) {
        if (dims[p] != 0) {
            std::cout << "\n  Skipping component " << p << " (dim=" << dims[p] << " != 0)" << std::endl;
            continue;
        }
        std::cout << "\n========================================" << std::endl;
        std::cout << "Component " << p << " processing" << std::endl;
        std::cout << "========================================" << std::endl;

        // Recompute GB for this component (what MMA calls compGb)
        auto compGb = RegionSolver::computeGroebnerBasis(primelist[p], allVars, MOD);
        dumpPolys("compGb (GB of prime component)", compGb);

        // classifyVariablesPrimeA
        std::vector<std::string> primeA, vargen, varpar, ximinpoly;
        std::vector<int> varDeg;
        RegionSolver::classifyVariablesPrimeA(primelist[p], Avars, Bvars,
            primeA, vargen, varpar, varDeg, ximinpoly);

        std::cout << "\n--- classifyVariablesPrimeA results ---" << std::endl;
        std::cout << "primeA (" << primeA.size() << "):" << std::endl;
        for (auto& p_str : primeA) std::cout << "  " << p_str << std::endl;
        std::cout << "vargen: ";
        for (auto& v : vargen) std::cout << v << " ";
        std::cout << std::endl;
        std::cout << "varpar: ";
        for (auto& v : varpar) std::cout << v << " ";
        std::cout << std::endl;
        std::cout << "varDeg: ";
        for (int d : varDeg) std::cout << d << " ";
        std::cout << std::endl;
        std::cout << "ximinpoly (" << ximinpoly.size() << "):" << std::endl;
        for (auto& p_str : ximinpoly) std::cout << "  " << p_str << std::endl;

        // MonomialBasisIndex
        auto basisIndex = RegionSolver::computeMonomialBasisIndex(ximinpoly, vargen);
        std::cout << "\nMonomialBasisIndex (" << basisIndex.size() << "):" << std::endl;
        for (size_t i = 0; i < basisIndex.size(); ++i) {
            std::cout << "  [" << i << "] ";
            for (int e : basisIndex[i]) std::cout << e << ",";
            std::cout << std::endl;
        }

        // VarRule solving
        std::map<std::string, std::string> varRule;
        if (!varpar.empty()) {
            std::cout << "\n--- solveVarRule ---" << std::endl;
            std::cout << "Calling solveVarRule with:" << std::endl;
            std::cout << "  varpar: ";
            for (auto& v : varpar) std::cout << v << " ";
            std::cout << std::endl;
            std::cout << "  vargen: ";
            for (auto& v : vargen) std::cout << v << " ";
            std::cout << std::endl;

            RegionSolver::solveVarRule(compGb, varpar, vargen, MOD, varRule);
            std::cout << "VarRule result (" << varRule.size() << " entries):" << std::endl;
            for (auto& [k, v] : varRule) {
                std::cout << "  " << k << " -> " << v << std::endl;
            }

            // Verify: what does the GB say about these variables?
            std::cout << "\n--- GB linear equations for varpar ---" << std::endl;
            for (const auto& poly : compGb) {
                std::string lm = RegionSolver::extractLeadingMonomial(poly);
                auto lmExps = RegionSolver::parseExponents(lm, allVars);
                // Check if this looks like a linear eq for a varpar variable
                bool isLinear = false;
                int varIdx = -1;
                for (size_t vi = 0; vi < varpar.size(); ++vi) {
                    for (int ai = 0; ai < ne; ++ai) {
                        if (varpar[vi] == Avars[ai] && lmExps[ne + ai] == 1) {
                            bool allZero = true;
                            for (int ki = 0; ki < ne * 2 && allZero; ++ki)
                                if (ki != ne + ai && lmExps[ki] != 0) allZero = false;
                            if (allZero) { isLinear = true; varIdx = vi; break; }
                        }
                    }
                    if (isLinear) break;
                }
                if (isLinear) {
                    std::cout << "  LINEAR[" << varpar[varIdx] << "]: " << poly << std::endl;
                    std::cout << "    LM = " << lm << std::endl;
                    // Parse and show non-LM terms
                    auto p = PolyArith::parseSingularPolynomial(poly, allVars, MOD);
                    std::cout << "    Parsed monomials:" << std::endl;
                    for (auto& m : p) {
                        std::cout << "      exps=[";
                        for (int e : m.exps) std::cout << e << ",";
                        std::cout << "] coeff=" << m.coeff << std::endl;
                    }
                }
            }
        }

        // FractionRule
        std::map<std::string, std::string> fractionRule;
        RegionSolver::computeFractionRule(compGb, allVars, allVars, vargen, ne, MOD, fractionRule);
        std::cout << "\nFractionRule (" << fractionRule.size() << " entries):" << std::endl;
        for (auto& [k, v] : fractionRule) {
            std::cout << "  " << k << " -> " << v << std::endl;
        }

        // MonomialBasisMatrix
        int nb = (int)basisIndex.size();
        std::vector<std::vector<std::vector<int64_t>>> mbm;
        if (nb > 0) {
            if (vargen.empty()) {
                mbm.resize(1, std::vector<std::vector<int64_t>>(1, std::vector<int64_t>(1, 1)));
            } else {
                RegionSolver::computeMonomialBasisMatrix(ximinpoly, vargen, basisIndex, nb, MOD, mbm);
            }
        }
        std::cout << "\nMonomialBasisMatrix: nb=" << nb << std::endl;
        if (!mbm.empty() && !mbm[0].empty()) {
            std::cout << "  MBM[0] (dim " << mbm.size() << "x" << mbm[0].size() << "x" << mbm[0][0].size() << "):" << std::endl;
            for (size_t k = 0; k < std::min((size_t)nb, mbm.size()); ++k) {
                std::cout << "    k=" << k << ": ";
                for (size_t i = 0; i < mbm[k].size(); ++i) {
                    for (size_t j = 0; j < mbm[k][i].size(); ++j) {
                        std::cout << mbm[k][i][j];
                        if (j < mbm[k][i].size() - 1) std::cout << ",";
                    }
                    if (i < mbm[k].size() - 1) std::cout << "; ";
                }
                std::cout << std::endl;
            }
        }

        // Ring matrices
        std::cout << "\n--- RingBuilder computeRingMatrices ---" << std::endl;
        RegionSolver::RegionData reg;
        reg.limitSector = topSector;
        reg.nb = nb;
        reg.VarIndep = vargen;
        reg.VarDep = varpar;
        reg.VarDeg = varDeg;
        reg.MinPoly = ximinpoly;
        reg.VarRule = varRule;
        reg.FractionRule = fractionRule;
        reg.MonomialBasisMatrix = mbm;
        // Build MonomialBasis strings
        for (const auto& idx : basisIndex) {
            PolyArith::Polynomial mpoly;
            mpoly.push_back({idx, int64_t(1)});
            reg.MonomialBasis.push_back(PolyArith::polyToSingularString(mpoly, vargen));
        }

        auto rm = RingBuilder::computeRingMatrices(reg, ne, MOD);
        for (int i = 0; i < ne; ++i) {
            std::cout << "  A[" << (i+1) << "] = [";
            for (int j = 0; j < nb * nb; ++j) {
                if (j > 0) std::cout << ",";
                std::cout << rm.A_list[i][j];
            }
            std::cout << "]" << std::endl;
        }
    }

    // Compare with MMA
    std::cout << "\n========================================" << std::endl;
    std::cout << "MMA Reference Comparison" << std::endl;
    std::cout << "========================================" << std::endl;
    try {
        auto mmaData = readMmaRingData("data/RingData_bub00.bin");
        std::cout << "MMA data/RingData_bub00.bin: " << mmaData.size() << " cells" << std::endl;
        for (size_t c = 0; c < mmaData.size(); ++c) {
            auto& mma = mmaData[c];
            std::cout << "\nCell " << c << ":" << std::endl;
            std::cout << "  limitSector: ";
            for (int v : mma.limitSector) std::cout << v << " ";
            std::cout << std::endl;
            std::cout << "  nb=" << mma.nb << std::endl;
            int ne_local = 2;
            for (int i = 0; i < ne_local; ++i) {
                std::cout << "  A[" << (i+1) << "] = [";
                for (int j = 0; j < mma.nb * mma.nb; ++j) {
                    if (j > 0) std::cout << ",";
                    std::cout << mma.A_flat[i * mma.nb * mma.nb + j];
                }
                std::cout << "]" << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cout << "Could not read MMA data: " << e.what() << std::endl;
    }

    std::cout << "\n=== Debug trace complete ===" << std::endl;
    return 0;
}