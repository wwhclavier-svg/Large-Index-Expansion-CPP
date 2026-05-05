#ifndef REGION_SOLVER_SINGULAR_HPP
#define REGION_SOLVER_SINGULAR_HPP

// Module C: Region Solver — algebraic decomposition of IBP equations
//
// Mirrors MMA's LIECoreAlgebra.wl:expRegSolve2 pipeline:
//   1. Convert IBP eqs → A/B variable equations
//   2. Groebner basis (Mathematica or Singular)
//   3. Minimal associated primes (Singular minAssGTZ)
//   4. Per-component Groebner basis + variable classification
//   5. Monomial basis extraction (quotient ring)
//   6. Recursion matrix construction
//
// Dependencies: SingularRunner.hpp, IBPEqGenerator.hpp

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cstdint>
#include <cctype>
#include <stdexcept>

#include "SingularRunner.hpp"

namespace RegionSolver {

// ==========================================
// Region data (mirrors MMA's bivarPrimeInfo output)
// ==========================================
struct RegionData {
    std::vector<int>          limitSector;     // subsector mask
    int                       nb = 0;          // monomial basis dimension
    std::vector<std::string>  VarIndep;        // independent (generator) variables
    std::vector<std::string>  VarDep;          // dependent (parametrized) variables
    std::vector<std::string>  MinPoly;         // minimal polynomials (GB for vargen)
    std::vector<std::string>  MonomialBasis;   // monomial basis for quotient ring
    std::vector<std::vector<int>> MonomialBasisIndex; // exponent vectors for basis
    std::map<std::string, std::string> VarRule;  // varDep → polynomial(varIndep)
    std::vector<int>          VarDeg;          // degrees of each variable
};

// ==========================================
// Utilities
// ==========================================

// Convert "A"[i] ↔ A[i] format (Singular-safe ↔ MMA-readable)
inline std::string makeSingularVar(const std::string& name, int idx) {
    return name + std::to_string(idx);
}

inline std::string makeBracketVar(const std::string& name, int idx) {
    return "\"" + name + "\"[" + std::to_string(idx) + "]";
}

// ==========================================
// Variable ordering for Singular ring declaration
//
// varorder = {{"B"[1],...,"B"[ne]}, {"A"[1],...,"A"[ne]}, param_vars}
// This encodes: A's > B's > params (dp block ordering)
// In dp: variables in LATER blocks are SMALLER (less significant)
// So: (dp(nB), dp(nA), dp(nParam), C) means:
//   params(1) < ... < params(nParam) < A[1] < ... < A[ne] < B[1] < ... < B[ne]
// ==========================================

inline std::string buildVarDeclaration(int ne, const std::vector<std::string>& paramNames,
                                        const std::string& monomialOrder = "dp") {
    std::stringstream ss;

    // Count variables
    int nA = ne, nB = ne;

    // Build variable list: B_1,...,B_nB, A_1,...,A_nA, param_1,...
    std::vector<std::string> vars;
    for (int i = 1; i <= nB; ++i) vars.push_back("B" + std::to_string(i));
    for (int i = 1; i <= nA; ++i) vars.push_back("A" + std::to_string(i));
    for (size_t i = 0; i < paramNames.size(); ++i)
        vars.push_back(paramNames[i]);

    // Block ordering: (dp(nB), dp(nA), dp(nParams), C)
    // NOTE: Singular uses dp(N) with variables in REVERSE order for block ordering
    // But for simplicity, I'll use a flat dp ordering first
    ss << "(";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i > 0) ss << ",";
        ss << vars[i];
    }
    ss << "),(" << monomialOrder << ")";

    return ss.str();
}

// ==========================================
// Build A/B equations from IBP equations
//
// From expRegSolve2:
//   aeqs0 = Coefficient[ibpeqs, n] /. g[α] → ∏ A_i^{α_i - v_i}
//   aeqs  = aeqs0 /. {1/A_i → B_i}  ∪  {A_i B_i - 1}
//
// In C++, we work with the stub pre-computed data for now.
// ==========================================

inline std::string buildABEquationsScript(
    const std::vector<std::string>& ibpPolynomials,
    int ne, int64_t modulus,
    const std::vector<std::string>& paramNames)
{
    std::stringstream ss;

    ss << "// Build A/B variable system from IBP equations\n";
    ss << "LIB \"primdec.lib\";\n\n";

    // Ring: B[1..ne], A[1..ne], params — dp ordering
    ss << "ring r0 = " << modulus << ", (";
    for (int i = 1; i <= ne; ++i) ss << "B" << i << ",";
    for (int i = 1; i <= ne; ++i) { ss << "A" << i; if (i < ne) ss << ","; }
    for (size_t i = 0; i < paramNames.size(); ++i) ss << "," << paramNames[i];
    ss << "), (dp);\n\n";

    // Equations
    ss << "ideal aeqs = ";
    for (size_t i = 0; i < ibpPolynomials.size(); ++i) {
        if (i > 0) ss << ",";
        ss << ibpPolynomials[i];
    }
    ss << ";\n";

    // Add A_i * B_i - 1 = 0 relations
    ss << "// Add A_i*B_i - 1 relations\n";
    ss << "ideal ab = aeqs";
    for (int i = 1; i <= ne; ++i)
        ss << ",A" << i << "*B" << i << "-1";
    ss << ";\n\n";

    return ss.str();
}

// ==========================================
// Groebner basis computation via Singular
// ==========================================

inline std::vector<std::string> computeGroebnerBasis(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& varOrder,
    int64_t modulus)
{
    return SingularRunner::groebnerBasis(polynomials, varOrder, modulus);
}

// ==========================================
// Minimal associated primes via Singular
//
// Mirrors expRegSolve2's call to SingularMinAssPrime.
// Returns zero-dimensional prime components only.
// ==========================================

inline std::pair<std::vector<std::vector<std::string>>, std::vector<int>>
computeMinAssPrimes(
    const std::vector<std::string>& ideal,
    const std::vector<std::string>& varOrder,
    int64_t modulus)
{
    auto [primelist, dims] = SingularRunner::minimalAssPrimes(
        ideal, varOrder, modulus, "lp");

    // Filter to zero-dimensional components
    std::vector<std::vector<std::string>> zdComps;
    std::vector<int> zdDims;
    for (size_t i = 0; i < primelist.size(); ++i) {
        if (dims[i] == 0) {
            zdComps.push_back(primelist[i]);
            zdDims.push_back(dims[i]);
        }
    }

    std::cout << "[RegionSolver] " << primelist.size() << " prime components, "
              << zdComps.size() << " zero-dimensional" << std::endl;

    return {zdComps, zdDims};
}

// ==========================================
// Monomial basis extraction (quotient ring basis)
//
// Mirrors quotientRingBasisPower:
//   - From Groebner basis, extract leading monomials
//   - Compute max degrees for each variable
//   - Generate all monomials below the degree bound
// ==========================================

// Simple parser for leading monomial extraction from Singular output.
// For lex ordering, the LM is just the first term of each polynomial.

inline std::string extractLeadingMonomial(const std::string& poly) {
    // The leading monomial for lex is the first term in the Singular output.
    // Singular writes: x1^2-1 which means LM = x1^2.
    // We need to extract the first term before any + or -.

    std::string lm;
    size_t i = 0;
    // Skip leading sign
    if (i < poly.size() && poly[i] == '-')
        { lm += '-'; ++i; }
    // Read until we hit + or - (but not ^ for exponent)
    while (i < poly.size()) {
        if (poly[i] == '+' || (poly[i] == '-' && i > 0 && poly[i-1] != '^')) break;
        lm += poly[i];
        ++i;
    }
    return lm;
}

// Count variable occurrences: parse x1^2*x3 from a monomial
inline std::vector<int> parseExponents(const std::string& monomial,
                                        const std::vector<std::string>& varOrder) {
    int n = static_cast<int>(varOrder.size());
    std::vector<int> exps(n, 0);

    std::string s = monomial;
    // Remove sign
    if (!s.empty() && s[0] == '-') s.erase(0, 1);

    // Find each variable pattern xi^e or xi
    for (int i = 0; i < n; ++i) {
        const std::string& v = varOrder[i];
        size_t pos = 0;
        while ((pos = s.find(v, pos)) != std::string::npos) {
            size_t after = pos + v.size();
            if (after < s.size() && s[after] == '^') {
                // Read exponent
                ++after; // skip ^
                std::string expStr;
                while (after < s.size() && std::isdigit(s[after]))
                    expStr += s[after++];
                exps[i] += std::stoi(expStr);
            } else {
                exps[i] += 1;
            }
            pos = after;
        }
    }
    return exps;
}

inline std::vector<std::vector<int>> computeMonomialBasisIndex(
    const std::vector<std::string>& gb,
    const std::vector<std::string>& varOrder)
{
    int n = static_cast<int>(varOrder.size());

    // Empty varOrder → single basis element (constant 1)
    if (n == 0)
        return {{}};

    // Extract leading monomials and compute max degree per variable
    std::vector<int> maxDeg(n, 0);
    for (const auto& poly : gb) {
        std::string lm = extractLeadingMonomial(poly);
        auto exps = parseExponents(lm, varOrder);
        for (int i = 0; i < n; ++i)
            maxDeg[i] = std::max(maxDeg[i], exps[i]);
    }

    // Generate all monomial indices:
    // For each var i: exponent ∈ [0, maxDeg[i]-1] (since GB gives degree bound)
    std::vector<std::vector<int>> basis;
    std::vector<int> current(n, 0);

    std::function<void(int)> generate = [&](int dim) {
        if (dim == n) {
            // Check if this monomial is NOT divisible by any leading monomial
            bool divisible = false;
            for (const auto& poly : gb) {
                std::string lm = extractLeadingMonomial(poly);
                auto lmExps = parseExponents(lm, varOrder);
                bool div = true;
                for (int j = 0; j < n; ++j) {
                    if (current[j] < lmExps[j]) { div = false; break; }
                }
                if (div) { divisible = true; break; }
            }
            if (!divisible)
                basis.push_back(current);
            return;
        }
        // Bounds: exponent from 0 to maxDeg[dim]-1, inclusive both ends.
        // If maxDeg[dim]==1, range is [0, 0] → one value.
        for (int e = 0; e < maxDeg[dim]; ++e) {
            current[dim] = e;
            generate(dim + 1);
        }
    };
    generate(0);

    return basis;
}

// ==========================================
// Variable classification into generators and parametrized
//
// Mirrors MMA's logic:
//   - For each polynomial in the GB, check leading term degree
//   - If deg > 1: variable is a "generator" (vargen)
//   - If deg == 1: variable is "parametrized" (varpar, solvable linearly)
// ==========================================

inline void classifyVariables(
    const std::vector<std::string>& gb,
    const std::vector<std::string>& varOrder,
    std::vector<std::string>& vargen,
    std::vector<std::string>& varpar,
    std::vector<int>& varDeg)
{
    // Count degree per variable from leading terms
    int n = static_cast<int>(varOrder.size());
    std::vector<int> maxTermDeg(n, 0);
    std::vector<std::string> maxTermPoly(n); // which poly has the highest-degree LM for this var

    for (const auto& poly : gb) {
        std::string lm = extractLeadingMonomial(poly);
        auto exps = parseExponents(lm, varOrder);
        for (int i = 0; i < n; ++i) {
            if (exps[i] > maxTermDeg[i]) {
                maxTermDeg[i] = exps[i];
                maxTermPoly[i] = poly;
            }
        }
    }

    vargen.clear();
    varpar.clear();
    varDeg.clear();

    for (int i = 0; i < n; ++i) {
        varDeg.push_back(maxTermDeg[i]);
        if (maxTermDeg[i] > 1) {
            vargen.push_back(varOrder[i]);
        } else if (maxTermDeg[i] == 1) {
            varpar.push_back(varOrder[i]);
        }
        // maxTermDeg[i] == 0 → variable doesn't appear in any LM (free parameter)
    }
}

// ==========================================
// Full region solving pipeline
//
// mirrors expRegSolve2 + bivarPrimeInfo
// ==========================================

inline std::vector<RegionData> solveRegion(
    const std::vector<std::string>& abEquations,   // A/B system
    const std::vector<int>& limitSector,
    int ne,
    int64_t modulus)
{
    std::vector<RegionData> results;

    // Build variable ordering for this region
    // The GB variables are only the A_i and B_i where the sector is active
    // For the full system: all A_i and B_i
    std::vector<std::string> allVars;
    for (int i = 1; i <= ne; ++i) allVars.push_back("B" + std::to_string(i));
    for (int i = 1; i <= ne; ++i) allVars.push_back("A" + std::to_string(i));

    std::cout << "[RegionSolver] Computing Groebner basis for sector ";
    for (int s : limitSector) std::cout << s;
    std::cout << " (ne=" << ne << ", neqs=" << abEquations.size() << ")" << std::endl;

    // Step 1: Groebner basis
    auto gb = computeGroebnerBasis(abEquations, allVars, modulus);
    std::cout << "  GB done: " << gb.size() << " polynomials" << std::endl;

    // Step 2: Minimal associated primes (zero-dimensional only)
    auto [primelist, dims] = computeMinAssPrimes(gb, allVars, modulus);

    if (primelist.empty()) {
        std::cout << "  No zero-dimensional prime components found." << std::endl;
        return results;
    }

    // Step 3: Per-component processing
    for (size_t p = 0; p < primelist.size(); ++p) {
        std::cout << "  Component " << (p+1) << "/" << primelist.size()
                  << " (" << primelist[p].size() << " gens, dim=" << dims[p] << ")" << std::endl;

        // Recompute GB for this component (with appropriate ordering)
        auto compGb = computeGroebnerBasis(primelist[p], allVars, modulus);

        // Classify variables
        RegionData reg;
        reg.limitSector = limitSector;

        classifyVariables(compGb, allVars, reg.VarIndep, reg.VarDep, reg.VarDeg);
        reg.MinPoly = compGb;

        // Compute monomial basis
        reg.MonomialBasisIndex = computeMonomialBasisIndex(compGb, reg.VarIndep);
        reg.nb = static_cast<int>(reg.MonomialBasisIndex.size());

        // Build monomial basis strings
        for (const auto& idx : reg.MonomialBasisIndex) {
            std::string term = "1";
            for (size_t i = 0; i < idx.size(); ++i) {
                if (idx[i] > 0) {
                    if (idx[i] == 1)
                        term += "*" + reg.VarIndep[i];
                    else
                        term += "*" + reg.VarIndep[i] + "^" + std::to_string(idx[i]);
                }
            }
            if (term == "1" && !reg.VarIndep.empty()) term = "1"; // the constant 1
            reg.MonomialBasis.push_back(term);
        }

        std::cout << "    VarIndep: ";
        for (auto& v : reg.VarIndep) std::cout << v << " ";
        std::cout << "\n    VarDep: ";
        for (auto& v : reg.VarDep) std::cout << v << " ";
        std::cout << "\n    nb = " << reg.nb << std::endl;

        results.push_back(reg);
    }

    return results;
}

// ==========================================
// Stub: extract region data from existing .bin file
//
// Phase 3 stub: reads IBP + Ring binary data, reconstructs RegionData.
// In full implementation, this will be replaced by Singular computation.
// ==========================================

// For the stub, region data is embedded in the binary files.
// We provide a wrapper that loads from the .bin and converts to RegionData.

} // namespace RegionSolver
#endif
