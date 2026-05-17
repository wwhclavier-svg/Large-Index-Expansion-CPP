#ifndef REGION_SOLVER_SINGULAR_HPP
#define REGION_SOLVER_SINGULAR_HPP

// Module C: Region Solver — algebraic decomposition of IBP equations
//
// Mirrors MMA's LIECoreAlgebra.wl:expRegSolve2 pipeline:
//   1. Convert IBP eqs → A/B variable equations (done by IBPAnalyzer)
//   2. Groebner basis (Singular)
//   3. Minimal associated primes (Singular primdecGTZ)
//   4. Per-component: GB, variable classification, VarRule, FractionRule
//   5. Monomial basis extraction + multiplication matrices (quotient ring)
//   6. Subsector decomposition: enumerate all subsectors, solve each
//
// Dependencies: SingularRunner.hpp, PolyArith.hpp, IBPAnalyzer.hpp

#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cstdint>
#include <cctype>
#include <stdexcept>

#include "SingularRunner.hpp"
#include "PolyArith.hpp"
#include "IBPAnalyzer.hpp"

namespace RegionSolver {

// ==========================================
// ISO 8601 timestamp for correlating Singular calls with system monitoring
// ==========================================
inline std::string nowStamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm;
    localtime_r(&t, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm);
    return buf;
}

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
    std::map<std::string, std::string> FractionRule; // "B"[j,i] → polynomial(varIndep)
    std::vector<std::vector<std::vector<int64_t>>> MonomialBasisMatrix; // [nb][nb][nb]
    std::vector<int>          VarDeg;          // degrees of each variable
};

// ==========================================
// Utilities
// ==========================================

inline std::string makeSingularVar(const std::string& name, int idx) {
    return name + std::to_string(idx);
}

// ==========================================
// Variable ordering for Singular ring declaration
// ==========================================

inline std::string buildVarDeclaration(int ne, const std::vector<std::string>& paramNames,
                                        const std::string& monomialOrder = "dp") {
    std::stringstream ss;
    int nA = ne, nB = ne;
    std::vector<std::string> vars;
    for (int i = 1; i <= nB; ++i) vars.push_back("B" + std::to_string(i));
    for (int i = 1; i <= nA; ++i) vars.push_back("A" + std::to_string(i));
    for (size_t i = 0; i < paramNames.size(); ++i)
        vars.push_back(paramNames[i]);
    ss << "(";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i > 0) ss << ",";
        ss << vars[i];
    }
    ss << "),(" << monomialOrder << ")";
    return ss.str();
}

// ==========================================
// Build A/B equations Singular script
// ==========================================

inline std::string buildABEquationsScript(
    const std::vector<std::string>& ibpPolynomials,
    int ne, int64_t modulus,
    const std::vector<std::string>& paramNames)
{
    std::stringstream ss;
    ss << "LIB \"primdec.lib\";\n\n";
    ss << "ring r0 = " << modulus << ", (";
    for (int i = 1; i <= ne; ++i) ss << "B" << i << ",";
    for (int i = 1; i <= ne; ++i) { ss << "A" << i; if (i < ne) ss << ","; }
    for (size_t i = 0; i < paramNames.size(); ++i) ss << "," << paramNames[i];
    ss << "), (dp);\n\n";

    ss << "ideal aeqs = ";
    for (size_t i = 0; i < ibpPolynomials.size(); ++i) {
        if (i > 0) ss << ",";
        ss << ibpPolynomials[i];
    }
    ss << ";\n";

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
// Minimal associated primes via Singular (combined GB + primdec + per-comp GB)
// ==========================================

inline SingularRunner::PrimdecResult computePrimdec(
    const std::vector<std::string>& ideal,
    const std::vector<std::string>& varOrder,
    int64_t modulus,
    bool useMinAss = true)
{
    const char* env = getenv("SINGULAR_DECOMP");
    if (env && std::string(env) == "primdecGTZ") useMinAss = false;
    std::cerr << "  [" << nowStamp() << "] Singular " << (useMinAss ? "minAssGTZ" : "primdecGTZ")
              << " starting..." << std::endl;
    auto result = SingularRunner::combinedPrimdec(ideal, varOrder, modulus, useMinAss);
    std::cerr << "  [" << nowStamp() << "] Singular " << (useMinAss ? "minAssGTZ" : "primdecGTZ")
              << " done" << std::endl;
    return result;
}

// ==========================================
// Monomial basis extraction (unchanged)
// ==========================================

inline std::string extractLeadingMonomial(const std::string& poly) {
    std::string lm;
    size_t i = 0;
    if (i < poly.size() && poly[i] == '-') { lm += '-'; ++i; }
    while (i < poly.size()) {
        if (poly[i] == '+' || (poly[i] == '-' && i > 0 && poly[i-1] != '^')) break;
        lm += poly[i];
        ++i;
    }
    return lm;
}

inline std::vector<int> parseExponents(const std::string& monomial,
                                        const std::vector<std::string>& varOrder) {
    int n = static_cast<int>(varOrder.size());
    std::vector<int> exps(n, 0);
    std::string s = monomial;
    if (!s.empty() && s[0] == '-') s.erase(0, 1);
    for (int i = 0; i < n; ++i) {
        const std::string& v = varOrder[i];
        size_t pos = 0;
        while ((pos = s.find(v, pos)) != std::string::npos) {
            size_t after = pos + v.size();
            if (after < s.size() && s[after] == '^') {
                ++after;
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
    if (n == 0) return {{}};

    std::vector<int> maxDeg(n, 0);
    for (const auto& poly : gb) {
        std::string lm = extractLeadingMonomial(poly);
        auto exps = parseExponents(lm, varOrder);
        for (int i = 0; i < n; ++i)
            maxDeg[i] = std::max(maxDeg[i], exps[i]);
    }

    std::vector<std::vector<int>> basis;
    std::vector<int> current(n, 0);
    std::function<void(int)> generate = [&](int dim) {
        if (dim == n) {
            bool divisible = false;
            for (const auto& poly : gb) {
                std::string lm = extractLeadingMonomial(poly);
                auto lmExps = parseExponents(lm, varOrder);
                // Skip LMs that don't involve any VarIndep variable
                // (their LM in the subring is 1, which divides everything)
                bool hasVar = false;
                for (int e : lmExps) if (e > 0) { hasVar = true; break; }
                if (!hasVar) continue;
                bool div = true;
                for (int j = 0; j < n; ++j) {
                    if (current[j] < lmExps[j]) { div = false; break; }
                }
                if (div) { divisible = true; break; }
            }
            if (!divisible) basis.push_back(current);
            return;
        }
        for (int e = 0; e < maxDeg[dim]; ++e) {
            current[dim] = e;
            generate(dim + 1);
        }
    };
    generate(0);
    return basis;
}

// ==========================================
// Variable classification
// ==========================================

inline void classifyVariables(
    const std::vector<std::string>& gb,
    const std::vector<std::string>& varOrder,
    std::vector<std::string>& vargen,
    std::vector<std::string>& varpar,
    std::vector<int>& varDeg)
{
    int n = static_cast<int>(varOrder.size());
    std::vector<int> maxTermDeg(n, 0);
    for (const auto& poly : gb) {
        std::string lm = extractLeadingMonomial(poly);
        auto exps = parseExponents(lm, varOrder);
        for (int i = 0; i < n; ++i) {
            if (exps[i] > maxTermDeg[i]) maxTermDeg[i] = exps[i];
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
    }
}

// MMA-compatible: classify variables from A-only polynomials.
inline void classifyVariablesPrimeA(
    const std::vector<std::string>& compGb,
    const std::vector<std::string>& Avars,
    const std::vector<std::string>& Bvars,
    std::vector<std::string>& primeA,
    std::vector<std::string>& vargen,
    std::vector<std::string>& varpar,
    std::vector<int>& varDeg,
    std::vector<std::string>& ximinpoly)
{
    primeA.clear(); vargen.clear(); varpar.clear(); varDeg.clear(); ximinpoly.clear();
    for (const auto& poly : compGb) {
        bool hasB = false;
        for (const auto& bv : Bvars)
            if (poly.find(bv) != std::string::npos) { hasB = true; break; }
        if (!hasB) primeA.push_back(poly);
    }
    if (primeA.empty()) return;
    struct LmInfo { std::string lm; int deg; };
    std::vector<LmInfo> lminfo;
    for (const auto& poly : primeA) {
        std::string lm = extractLeadingMonomial(poly);
        auto exps = parseExponents(lm, Avars);
        int deg = 0;
        for (int e : exps) deg += e;
        lminfo.push_back({lm, deg});
    }
    std::set<std::string> genVars, parVars;
    for (size_t i = 0; i < primeA.size(); ++i) {
        auto exps = parseExponents(lminfo[i].lm, Avars);
        for (size_t vi = 0; vi < Avars.size(); ++vi)
            if (exps[vi] > 0)
                (lminfo[i].deg > 1 ? genVars : parVars).insert(Avars[vi]);
    }
    vargen.assign(genVars.begin(), genVars.end());
    varpar.assign(parVars.begin(), parVars.end());
    // If a variable is classified as both generator AND parameter, keep only generator
    for (const auto& v : genVars)
        varpar.erase(std::remove(varpar.begin(), varpar.end(), v), varpar.end());
    for (size_t i = 0; i < primeA.size(); ++i)
        if (lminfo[i].deg > 1) ximinpoly.push_back(primeA[i]);
    // Ensure all Avars not appearing as generator are treated as parameters.
    // Variables in non-leading positions (body of a polynomial) are missed by
    // the LM-based classification but still need to be in VarDep for VarRule.
    std::set<std::string> classified(genVars.begin(), genVars.end());
    classified.insert(varpar.begin(), varpar.end());
    for (const auto& av : Avars) {
        if (classified.find(av) == classified.end())
            varpar.push_back(av);
    }
}


// ==========================================
// Singular polynomial reduction helper
//
// Reduces a list of polynomials modulo a GB in a given ring.
// Returns the reduced polynomials as strings.
// ==========================================

inline std::vector<std::string> reducePolynomialsSingular(
    const std::vector<std::string>& polys,
    const std::vector<std::string>& gb,
    const std::vector<std::string>& varOrder,
    int64_t modulus)
{
    if (polys.empty()) return {};

    std::stringstream ss;
    ss << "ring r = " << modulus << ", (";
    for (size_t i = 0; i < varOrder.size(); ++i) {
        if (i > 0) ss << ",";
        ss << varOrder[i];
    }
    ss << "), (lp);\n";

    // GB
    ss << "ideal gb = ";
    for (size_t i = 0; i < gb.size(); ++i) {
        if (i > 0) ss << ",";
        ss << gb[i];
    }
    ss << ";\n\n";

    // Reduce each polynomial
    for (size_t i = 0; i < polys.size(); ++i) {
        ss << "poly p" << i << " = " << polys[i] << ";\n";
        ss << "poly r" << i << " = reduce(p" << i << ", gb);\n";
    }
    ss << "string results;\n";

    // Build the output string with pipe separators
    ss << "results = ";
    for (size_t i = 0; i < polys.size(); ++i) {
        if (i > 0) ss << " + \"|\" + ";
        ss << "string(r" << i << ")";
    }
    ss << ";\n";

    auto out = SingularRunner::runSingular(ss.str(), {"results"});
    if (out.count("results") == 0) return {};

    // Parse pipe-delimited results
    std::vector<std::string> result;
    std::string data = out["results"];
    // Remove trailing newline
    while (!data.empty() && (data.back() == '\n' || data.back() == '\r'))
        data.pop_back();

    size_t pos = 0;
    while (pos < data.size()) {
        size_t next = data.find('|', pos);
        if (next == std::string::npos) {
            result.push_back(data.substr(pos));
            break;
        }
        result.push_back(data.substr(pos, next - pos));
        pos = next + 1;
    }
    return result;
}

// ==========================================
// solveVarRule — extract VarRule from Groebner basis
//
// For lex-ordered GB, each varpar variable has a polynomial with LM = varpar^1.
// This polynomial is of the form: varpar_i + P(vargen) = 0
// So varpar_i = -P(vargen) = (mod - P(vargen)) mod modulus.
//
// We parse the Singular polynomial string, extract the non-LM terms,
// and build the VarRule mapping.
// ==========================================

inline void solveVarRule(
    const std::vector<std::string>& compGb,
    const std::vector<std::string>& varpar,
    const std::vector<std::string>& vargen,
    int64_t modulus,
    std::map<std::string, std::string>& varRule)
{
    varRule.clear();

    // Combined variable list: vargen first, then varpar
    std::vector<std::string> allVars = vargen;
    allVars.insert(allVars.end(), varpar.begin(), varpar.end());

    for (const auto& dpVar : varpar) {
        // Find the GB polynomial whose LM (first term) starts with this variable
        for (const auto& poly : compGb) {
            std::string lm = extractLeadingMonomial(poly);
            // Check if LM is exactly dpVar^1
            auto lmExps = parseExponents(lm, allVars);
            int varIdx = -1;
            for (int i = 0; i < (int)allVars.size(); ++i) {
                if (allVars[i] == dpVar) { varIdx = i; break; }
            }
            if (varIdx < 0) continue;

            bool isLinear = (lmExps[varIdx] == 1);
            for (int i = 0; i < (int)allVars.size(); ++i)
                if (i != varIdx && lmExps[i] != 0) { isLinear = false; break; }
            if (!isLinear) continue;

            // Parse full polynomial and extract non-LM part
            // The polynomial is: dpVar + other_terms = 0  →  dpVar = -other_terms
            PolyArith::Polynomial p = PolyArith::parseSingularPolynomial(poly, allVars, modulus);

            // Remove the LM term (dpVar with coeff 1)
            PolyArith::Polynomial rhs;
            for (const auto& m : p) {
                bool isLM = (m.exps[varIdx] == 1);
                for (int i = 0; i < (int)allVars.size(); ++i)
                    if (i != varIdx && m.exps[i] != 0) isLM = false;
                if (isLM) continue; // skip the LM term
                // Negate
                int64_t c = (-m.coeff) % modulus;
                if (c < 0) c += modulus;
                if (c != 0) rhs.push_back({m.exps, c});
            }
            PolyArith::canonicalize(rhs, modulus);

            // Remap exponents to vargen-only space (first nVargen entries of allVars)
            int nVargen = (int)vargen.size();
            PolyArith::Polynomial projected;
            for (const auto& m : rhs) {
                std::vector<int> newExps(nVargen, 0);
                for (int ei = 0; ei < nVargen && ei < (int)m.exps.size(); ++ei)
                    newExps[ei] = m.exps[ei];
                projected.push_back({newExps, m.coeff});
            }
            PolyArith::canonicalize(projected, modulus);

            // Format as string in vargen variables only
            // When vargen is empty, ensure constant polynomials format correctly
            std::string resultStr;
            if (vargen.empty()) {
                // For empty vargen: result should be a constant
                // If projected is non-empty, it must be a constant polynomial
                if (!projected.empty()) {
                    resultStr = std::to_string(projected[0].coeff % modulus);
                    if (resultStr[0] == '-') {
                        // negative constant - convert to mod representation
                        int64_t cv = projected[0].coeff % modulus;
                        if (cv < 0) cv += modulus;
                        resultStr = std::to_string(cv);
                    }
                } else {
                    resultStr = "0";
                }
            } else {
                resultStr = PolyArith::polyToSingularString(projected, vargen);
            }
            varRule[dpVar] = resultStr;
            break;
        }
    }
}

// ==========================================
// computeFractionRule — compute B[j,i] = A_i * B_j reduced mod GB
//
// For each (i,j) with i≠j, reduces A_i * B_j modulo the prime ideal GB.
// Returns polynomial expressions in VarIndep.
// ==========================================

inline void computeFractionRule(
    const std::vector<std::string>& compGb,
    const std::vector<std::string>& varOrder,  // all A/B vars
    const std::vector<std::string>& abVars,    // A1..Ane, B1..Bne in order
    const std::vector<std::string>& vargen,    // VarIndep for formatting output
    int ne,
    int64_t modulus,
    std::map<std::string, std::string>& fractionRule)
{
    fractionRule.clear();

    // Build list of products A_i * B_j for i != j
    std::vector<std::string> products;
    std::vector<std::pair<int,int>> pairs;
    for (int i = 1; i <= ne; ++i) {
        for (int j = 1; j <= ne; ++j) {
            if (i == j) continue;
            products.push_back("A" + std::to_string(i) + "*B" + std::to_string(j));
            pairs.push_back({i, j});
        }
    }

    if (products.empty()) return;

    auto reduced = reducePolynomialsSingular(products, compGb, varOrder, modulus);

    for (size_t k = 0; k < reduced.size() && k < pairs.size(); ++k) {
        int i = pairs[k].first;
        int j = pairs[k].second;
        std::string key = "\"B\"[" + std::to_string(j) + "," + std::to_string(i) + "]";
        // Parse with full varOrder, then project to vargen
        PolyArith::Polynomial poly = PolyArith::parseSingularPolynomial(reduced[k], varOrder, modulus);
        // Remap exponents: extract only vargen components
        // Build index map from varOrder → vargen
        std::map<std::string, int> varOrderIdx;
        for (size_t vi = 0; vi < varOrder.size(); ++vi) varOrderIdx[varOrder[vi]] = (int)vi;
        std::vector<int> vargenIdx(varOrder.size(), -1);
        for (size_t vi = 0; vi < vargen.size(); ++vi) {
            auto it = varOrderIdx.find(vargen[vi]);
            if (it != varOrderIdx.end()) vargenIdx[it->second] = (int)vi;
        }
        PolyArith::Polynomial projected;
        for (const auto& m : poly) {
            std::vector<int> newExps(vargen.empty() ? 0 : vargen.size(), 0);
            for (size_t ei = 0; ei < m.exps.size() && ei < varOrder.size(); ++ei) {
                if (m.exps[ei] != 0 && vargenIdx[ei] >= 0)
                    newExps[vargenIdx[ei]] = m.exps[ei];
            }
            projected.push_back({newExps, m.coeff});
        }
        PolyArith::canonicalize(projected, modulus);
        fractionRule[key] = PolyArith::polyToSingularString(projected, vargen);
    }
}

// ==========================================
// computeMonomialBasisMatrix — compute multiplication matrices
//
// For each basis monomial m_k, computes the nb×nb matrix representing
// multiplication by m_k in the quotient ring.
//
// Matrix[k][i][j] = coefficient of basis[j] in reduce(basis[i] * basis[k], gb)
//
// Uses Singular's reduce() to handle the polynomial reduction.
// ==========================================

// Helper: extract all variable names referenced in GB polynomials
inline std::vector<std::string> extractAllVarsFromGB(
    const std::vector<std::string>& compGb,
    const std::vector<std::string>& knownVars)
{
    std::set<std::string> vars(knownVars.begin(), knownVars.end());
    for (const auto& poly : compGb) {
        std::string cur;
        for (char ch : poly) {
            if (std::isalpha(ch) || ch == '_') {
                cur += ch;
            } else {
                if (!cur.empty() && std::all_of(cur.begin(), cur.end(),
                        [](char c) { return std::isalpha(c) || c == '_'; })) {
                    if (cur != "x" && cur != "y" && cur != "z") // skip Singular keywords
                        vars.insert(cur);
                }
                cur.clear();
            }
        }
        if (!cur.empty() && std::all_of(cur.begin(), cur.end(),
                [](char c) { return std::isalpha(c) || c == '_'; })) {
            if (cur != "x" && cur != "y" && cur != "z")
                vars.insert(cur);
        }
    }
    return std::vector<std::string>(vars.begin(), vars.end());
}

inline void computeMonomialBasisMatrix(
    const std::vector<std::string>& compGb,
    const std::vector<std::string>& vargen,
    const std::vector<std::vector<int>>& basisIndex,
    int nb,
    int64_t modulus,
    std::vector<std::vector<std::vector<int64_t>>>& basisMatrix)
{
    basisMatrix.clear();
    if (nb == 0 || vargen.empty()) return;

    basisMatrix.resize(nb,
        std::vector<std::vector<int64_t>>(nb,
            std::vector<int64_t>(nb, 0)));

    // Build monomial strings from basis index (in vargen space)
    std::vector<std::string> basisPolys;
    for (int i = 0; i < nb; ++i) {
        PolyArith::Polynomial p;
        p.push_back({basisIndex[i], int64_t(1)});
        basisPolys.push_back(PolyArith::polyToSingularString(p, vargen));
    }

    // Build index map: vargen name → position in vargen vector
    std::map<std::string, int> vargenIdx;
    for (size_t vi = 0; vi < vargen.size(); ++vi)
        vargenIdx[vargen[vi]] = (int)vi;

    // Extract ALL variables referenced in GB to build a complete Singular ring
    auto allGbVars = extractAllVarsFromGB(compGb, vargen);

    // For each basis element k, compute products with all basis elements
    for (int k = 0; k < nb; ++k) {
        std::vector<std::string> products;
        for (int i = 0; i < nb; ++i) {
            if (basisPolys[i] == "1") {
                products.push_back(basisPolys[k]);
            } else if (basisPolys[k] == "1") {
                products.push_back(basisPolys[i]);
            } else {
                products.push_back(basisPolys[i] + "*" + basisPolys[k]);
            }
        }

        // Use the full variable list from GB for the Singular ring, so that
        // any variable referenced in GB polynomials is properly defined.
        std::cerr << "    [" << nowStamp() << "] MBM reduction k=" << k << "/" << nb
                  << " (" << nb << " products)..." << std::endl;
        auto reduced = reducePolynomialsSingular(products, compGb, allGbVars, modulus);

        if (reduced.empty()) {
            std::cerr << "[computeMonomialBasisMatrix] WARNING: Singular reduction "
                         "returned empty for k=" << k << std::endl;
            continue;
        }

        // Parse each reduced polynomial and express in the monomial basis
        for (int i = 0; i < nb && i < (int)reduced.size(); ++i) {
            PolyArith::Polynomial poly = PolyArith::parseSingularPolynomial(
                reduced[i], vargen, modulus);

            if (poly.empty()) continue;

            for (const auto& mon : poly) {
                for (int j = 0; j < nb; ++j) {
                    bool match = true;
                    for (size_t vi = 0; vi < basisIndex[j].size(); ++vi) {
                        if (basisIndex[j][vi] != mon.exps[vi]) {
                            match = false;
                            break;
                        }
                    }
                    if (match) {
                        basisMatrix[k][i][j] = (basisMatrix[k][i][j] + mon.coeff) % modulus;
                        if (basisMatrix[k][i][j] < 0) basisMatrix[k][i][j] += modulus;
                        break;
                    }
                }
            }
        }
    }

    // Transpose each k-matrix to match MMA convention:
    // MMA: M[a][c][b] = coeff of basis[c] in (basis[a] * basis[b])
    // C++ built: M[a][b][c] = coeff of basis[c] in (basis[b] * basis[a])
    // After transpose, C++ M_T[a][c][b] matches MMA M[a][c][b]
    // (commutative multiplication, so basis[a]*basis[b] = basis[b]*basis[a])
    for (int k = 0; k < nb; ++k) {
        for (int i = 0; i < nb; ++i) {
            for (int j = i + 1; j < nb; ++j) {
                std::swap(basisMatrix[k][i][j], basisMatrix[k][j][i]);
            }
        }
    }
}

// ==========================================
// Full region solving pipeline
//
// Takes A/B polynomial strings (from IBPAnalyzer::buildABEquations),
// computes Groebner basis, primary decomposition, variable classification,
// VarRule, FractionRule, and monomial basis matrices.
// ==========================================

inline std::vector<RegionData> solveRegion(
    const std::vector<std::string>& abEquations,
    const std::vector<int>& limitSector,
    int ne,
    int64_t modulus)
{
    std::vector<RegionData> results;

    std::vector<std::string> allVars;
    for (int i = 1; i <= ne; ++i) allVars.push_back("B" + std::to_string(i));
    for (int i = 1; i <= ne; ++i) allVars.push_back("A" + std::to_string(i));

    // Add A_i*B_i-1 constraints to make the ideal zero-dimensional
    std::vector<std::string> fullIdeal = abEquations;
    for (int i = 1; i <= ne; ++i)
        fullIdeal.push_back("A" + std::to_string(i) + "*B" + std::to_string(i) + "-1");

    // Single Singular invocation: GB + dim check + primdecGTZ + per-comp GB
    auto pd = computePrimdec(fullIdeal, allVars, modulus);
    auto& primelist = pd.primelist;
    auto& compGbs = pd.compGbs;

    if (primelist.empty())
        return results;

    // Per-component processing (MMA-compatible)
    for (size_t p = 0; p < primelist.size(); ++p) {
        std::cerr << "  [" << nowStamp() << "] Component " << p << "/" << primelist.size()
                  << " (" << primelist[p].size() << " polys)" << std::endl;

        // Debug: print raw prime
        if (ne <= 4) {
            std::cerr << "\n[DEBUG rawPrime] sector=[";
            for (size_t si=0;si<limitSector.size();++si){if(si)std::cerr<<",";std::cerr<<limitSector[si];}
            std::cerr << "] comp=" << p << " nPoly=" << primelist[p].size() << std::endl;
            for (size_t gi=0; gi<primelist[p].size(); ++gi)
                std::cerr << "  [" << gi << "] " << primelist[p][gi] << std::endl;
        }

        const auto& compGb = compGbs[p];

        RegionData reg;
        reg.limitSector = limitSector;

        // A and B variable lists
        std::vector<std::string> Avars, Bvars;
        for (int i = 1; i <= ne; ++i) {
            Bvars.push_back("B" + std::to_string(i));
            Avars.push_back("A" + std::to_string(i));
        }

        // MMA-compatible: classify from RAW component (not compGb)
        std::vector<std::string> primeA, ximinpoly;
        classifyVariablesPrimeA(primelist[p], Avars, Bvars,
            primeA, reg.VarIndep, reg.VarDep, reg.VarDeg, ximinpoly);
        reg.MinPoly = ximinpoly;

        // MMA-compatible: basis index from ximinpoly only
        reg.MonomialBasisIndex = computeMonomialBasisIndex(ximinpoly, reg.VarIndep);
        reg.nb = static_cast<int>(reg.MonomialBasisIndex.size());

        for (const auto& idx : reg.MonomialBasisIndex) {
            PolyArith::Polynomial mpoly;
            mpoly.push_back({idx, int64_t(1)});
            reg.MonomialBasis.push_back(PolyArith::polyToSingularString(mpoly, reg.VarIndep));
        }

        // Solve VarRule
        if (!reg.VarDep.empty()) {
            // Debug: print GB for nb>1 regions to compare with MMA
            if (reg.nb > 1) {
                std::cerr << "\n[DEBUG solveVarRule] sector=[";
                for (size_t si=0;si<limitSector.size();++si){if(si)std::cerr<<",";std::cerr<<limitSector[si];}
                std::cerr << "] nb=" << reg.nb << " VarIndep=[";
                for (size_t vi=0;vi<reg.VarIndep.size();++vi){if(vi)std::cerr<<",";std::cerr<<reg.VarIndep[vi];}
                std::cerr << "] VarDep=[";
                for (size_t vi=0;vi<reg.VarDep.size();++vi){if(vi)std::cerr<<",";std::cerr<<reg.VarDep[vi];}
                std::cerr << "]" << std::endl;
                std::cerr << "  compGb:" << std::endl;
                for (size_t gi=0; gi<compGb.size(); ++gi)
                    std::cerr << "    [" << gi << "] " << compGb[gi] << std::endl;
            }
            solveVarRule(compGb, reg.VarDep, reg.VarIndep, modulus, reg.VarRule);
            if (reg.nb > 1) {
                std::cerr << "  VarRule result:" << std::endl;
                for (auto& [k,v] : reg.VarRule)
                    std::cerr << "    " << k << " -> " << v << std::endl;
            }
        }

        // Substitute VarDep variables out of ximinpoly using VarRule.
        // ximinpoly (A-only minipoly) may reference VarDep variables
        // that are not in VarIndep. Substituting them out ensures
        // computeMonomialBasisMatrix works in a pure VarIndep ring.
        std::vector<std::string> ximinpolySubsted;
        for (const auto& poly : ximinpoly) {
            std::string p = poly;
            if (!reg.VarDep.empty()) {
                for (const auto& [var, val] : reg.VarRule) {
                    size_t pos = 0;
                    while ((pos = p.find(var, pos)) != std::string::npos) {
                        bool wordBoundary = true;
                        if (pos > 0 && (std::isalpha(p[pos-1]) || std::isdigit(p[pos-1])))
                            wordBoundary = false;
                        size_t end = pos + var.size();
                        if (end < p.size() && (std::isalpha(p[end]) || std::isdigit(p[end])))
                            wordBoundary = false;
                        if (wordBoundary) {
                            p.replace(pos, var.size(), "(" + val + ")");
                            pos += val.size() + 2;
                        } else {
                            pos += var.size();
                        }
                    }
                }
            }
            ximinpolySubsted.push_back(p);
        }

        // FractionRule
        computeFractionRule(compGb, allVars, allVars, reg.VarIndep, ne, modulus, reg.FractionRule);

        // MonomialBasisMatrix from ximinpoly (A-only minimal polynomials)
        if (reg.nb > 0) {
            if (reg.VarIndep.empty()) {
                // nb=1 identity matrix for constant ring
                reg.MonomialBasisMatrix.resize(1,
                    std::vector<std::vector<int64_t>>(1,
                        std::vector<int64_t>(1, 1)));
            } else {
                computeMonomialBasisMatrix(ximinpolySubsted, reg.VarIndep,
                    reg.MonomialBasisIndex, reg.nb, modulus, reg.MonomialBasisMatrix);
            }
        }

        results.push_back(reg);
    }
    return results;
}

// ==========================================
// Subsector enumeration
//
// Generates all binary masks m such that m[i] <= topSector[i].
// For a top sector of {1,1,...,1}, this enumerates all 2^ne subsectors.
// ==========================================

// Generate subsectors matching MMA's Subsets[activeIndices, {nl, ne}]
inline std::vector<std::vector<int>> generateSubsectors(const std::vector<int>& topSector, int nl) {
    int ne = static_cast<int>(topSector.size());
    std::vector<int> activeProps;
    for (int i = 0; i < ne; ++i)
        if (topSector[i] == 1) activeProps.push_back(i);

    int nActive = static_cast<int>(activeProps.size());
    int minSize = std::max(1, nl);
    int maxSize = ne;
    std::vector<std::vector<int>> result;
    std::function<void(int, std::vector<int>&)> dfs = [&](int start, std::vector<int>& chosen) {
        if ((int)chosen.size() >= minSize && (int)chosen.size() <= maxSize) {
            std::vector<int> sector(ne, 0);
            for (int idx : chosen) sector[idx] = 1;
            result.push_back(sector);
        }
        if ((int)chosen.size() >= maxSize) return;
        for (int i = start; i < nActive; ++i) {
            chosen.push_back(activeProps[i]);
            dfs(i + 1, chosen);
            chosen.pop_back();
        }
    };
    std::vector<int> chosen;
    dfs(0, chosen);

    // Sort to match MMA ordering: by weight ascending, then binary value descending
    std::sort(result.begin(), result.end(),
        [ne](const std::vector<int>& a, const std::vector<int>& b) {
            int wa = 0, wb = 0;
            int va = 0, vb = 0;
            for (int i = 0; i < ne; ++i) {
                wa += a[i];
                wb += b[i];
                va = (va << 1) | a[i];
                vb = (vb << 1) | b[i];
            }
            if (wa != wb) return wa < wb;
            return va > vb;  // descending binary within same weight
        });

    return result;
}

// ==========================================
// solveAllSectors — full subsector decomposition
//
// Mirrors MMA's expRegSolve2 outer loop over limitSector values.
//
// For each subsector of the top sector:
//   1. Build A/B equations via IBPAnalyzer::buildABEquations
//   2. Run solveRegion (GB + primdec + VarRule + FractionRule + MBM)
//   3. Collect results with correct limitSector labels
//
// Skips subsectors that produce no zero-dimensional components.
// When deduplicate=true (default), sorts subsectors top-first and deduplicates
// identical regions across subsectors by nb+FracRule+MBM
// (matches MMA's LIESolveRegion).
// ==========================================

inline std::vector<RegionData> solveAllSectors(
    const IBPEqGenerator::IBPEquations& ibp,
    const std::vector<int>& topSector,
    int64_t modulus,
    bool deduplicate = false)
{
    int ne = ibp.ne;
    int nl = ibp.nl;
    std::vector<RegionData> allResults;
    auto subsectors = generateSubsectors(topSector, nl);

    if (deduplicate) {
        std::sort(subsectors.begin(), subsectors.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
                int sa = 0, sb = 0;
                for (int v : a) sa += v;
                for (int v : b) sb += v;
                if (sa != sb) return sa > sb;
                return a > b;
            });
    }

    int sectorIdx = 0, totalSectors = (int)subsectors.size();
    for (const auto& sub : subsectors) {
        ++sectorIdx;
        auto abEqs = IBPAnalyzer::buildABEquations(ibp, sub, modulus);

        // DEBUG: print A/B for ne<=4 families
        if (ne <= 4) {
            std::cout << "\n[DEBUG A/B] sector=[";
            for (size_t si=0;si<sub.size();++si){if(si)std::cout<<",";std::cout<<sub[si];}
            std::cout << "] ne=" << ne << " nibp=" << ibp.nibp << std::endl;
            for (size_t ei=0;ei<abEqs.size();++ei) std::cout << "  AB[" << ei << "] = " << abEqs[ei] << std::endl;
        }
        // Skip empty subsectors (all equations are constant 0)
        bool hasContent = false;
        for (const auto& eq : abEqs) {
            for (char ch : eq) {
                if (ch != '0' && ch != '+' && ch != '-' && ch != '*') {
                    hasContent = true;
                    break;
                }
            }
            if (hasContent) break;
            if (!eq.empty() && eq != "0") hasContent = true;
        }
        if (!hasContent) continue;

        std::cerr << "[sector " << sectorIdx << "/" << totalSectors << "] solving [";
        for (size_t si=0;si<sub.size();++si){if(si)std::cerr<<",";std::cerr<<sub[si];}
        std::cerr << "]" << std::endl;

        auto t_solve_0 = std::chrono::steady_clock::now();
        auto regions = solveRegion(abEqs, sub, ne, modulus);
        auto t_solve_1 = std::chrono::steady_clock::now();
        double sec_s = std::chrono::duration<double>(t_solve_1 - t_solve_0).count();
        int nreg = (int)regions.size();
        std::cerr << "  -> " << nreg << " region(s) in " << sec_s << " s" << std::endl;

        if (deduplicate) {
            for (auto& reg : regions) {
                bool isDuplicate = false;
                for (const auto& existing : allResults) {
                    if (reg.nb != existing.nb) continue;
                    if (reg.FractionRule != existing.FractionRule) continue;
                    if (reg.MonomialBasisIndex != existing.MonomialBasisIndex) continue;
                    isDuplicate = true;
                    std::cerr << "  -> skipping duplicate region (nb=" << reg.nb
                              << ", sector=[" << (reg.limitSector.empty()?-1:reg.limitSector[0])
                              << "...])" << std::endl;
                    break;
                }
                if (!isDuplicate)
                    allResults.push_back(std::move(reg));
            }
        } else {
            for (auto& reg : regions)
                allResults.push_back(std::move(reg));
        }
    }

    return allResults;
}

} // namespace RegionSolver
#endif
