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
// Minimal associated primes via Singular
// ==========================================

inline std::pair<std::vector<std::vector<std::string>>, std::vector<int>>
computeMinAssPrimes(
    const std::vector<std::string>& ideal,
    const std::vector<std::string>& varOrder,
    int64_t modulus)
{
    auto [primelist, dims] = SingularRunner::minimalAssPrimes(
        ideal, varOrder, modulus, "lp");
    std::vector<std::vector<std::string>> zdComps;
    std::vector<int> zdDims;
    for (size_t i = 0; i < primelist.size(); ++i) {
        if (dims[i] == 0) {
            zdComps.push_back(primelist[i]);
            zdDims.push_back(dims[i]);
        }
    }
    return {zdComps, zdDims};
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
            varRule[dpVar] = PolyArith::polyToSingularString(projected, vargen);
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

        // Use FULL variable ring (allVars from compGb context) to avoid
        // "undefined variable" errors. The GB polynomials may reference
        // VarDep variables that aren't in vargen alone.
        // We build the full varOrder from the compGb by extracting variable names.
        // For now, use the vargen ring if the GB only involves vargen variables;
        // otherwise use the full ring and project results.
        //
        // Strategy: use vargen ring for reduction. If the GB references non-vargen
        // vars, Singular will error → reduced comes back empty → we fall through
        // and basisMatrix stays zero.
        //
        // The caller (solveRegion) ensures compGb is projected to vargen space
        // by substituting VarRule before calling this function. See the
        // solveRegion code where we build vargenGb.
        auto reduced = reducePolynomialsSingular(products, compGb, vargen, modulus);

        if (reduced.empty()) {
            std::cerr << "[computeMonomialBasisMatrix] WARNING: Singular reduction "
                         "returned empty for k=" << k
                      << ". GB may reference undefined variables." << std::endl;
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

    // Step 1: Groebner basis
    auto gb = computeGroebnerBasis(fullIdeal, allVars, modulus);

    // Step 2: Minimal associated primes (zero-dimensional only)
    auto [primelist, dims] = computeMinAssPrimes(gb, allVars, modulus);

    if (primelist.empty())
        return results;

    // Step 3: Per-component processing
    for (size_t p = 0; p < primelist.size(); ++p) {

        // Recompute GB for this component
        auto compGb = computeGroebnerBasis(primelist[p], allVars, modulus);

        RegionData reg;
        reg.limitSector = limitSector;

        // Classify variables
        classifyVariables(compGb, allVars, reg.VarIndep, reg.VarDep, reg.VarDeg);
        reg.MinPoly = compGb;

        // Compute monomial basis
        reg.MonomialBasisIndex = computeMonomialBasisIndex(compGb, reg.VarIndep);
        reg.nb = static_cast<int>(reg.MonomialBasisIndex.size());

        // Build monomial basis strings
        for (const auto& idx : reg.MonomialBasisIndex) {
            PolyArith::Polynomial mpoly;
            mpoly.push_back({idx, int64_t(1)});
            reg.MonomialBasis.push_back(PolyArith::polyToSingularString(mpoly, reg.VarIndep));
        }

        // Solve VarRule: express VarDep in terms of VarIndep
        if (!reg.VarDep.empty()) {
            solveVarRule(compGb, reg.VarDep, reg.VarIndep, modulus, reg.VarRule);
        }

        // Compute FractionRule: A_i * B_j reduced modulo GB
        computeFractionRule(compGb, allVars, allVars, reg.VarIndep, ne, modulus, reg.FractionRule);

        // Compute MonomialBasisMatrix: multiplication matrices
        if (reg.nb > 0) {
            if (reg.VarIndep.empty()) {
                reg.MonomialBasisMatrix.resize(1,
                    std::vector<std::vector<int64_t>>(1,
                        std::vector<int64_t>(1, 1)));
            } else {
                // Project GB to VarIndep space by substituting VarRule.
                // The compGb contains all A/B variables, but Singular reduction
                // for the monomial basis needs only VarIndep variables in the ring.
                std::vector<std::string> vargenGb;
                // Build replacement list: for each allVars variable, its
                // expression as a polynomial in VarIndep.
                std::vector<PolyArith::Polynomial> replacements(allVars.size());
                for (size_t vi = 0; vi < allVars.size(); ++vi) {
                    const std::string& vname = allVars[vi];
                    auto rit = reg.VarRule.find(vname);
                    if (rit != reg.VarRule.end()) {
                        replacements[vi] = PolyArith::parseSingularPolynomial(
                            rit->second, reg.VarIndep, modulus);
                    } else {
                        // VarIndep variable: it maps to itself in VarIndep space
                        auto vit = std::find(reg.VarIndep.begin(), reg.VarIndep.end(), vname);
                        if (vit != reg.VarIndep.end()) {
                            int pos = (int)(vit - reg.VarIndep.begin());
                            std::vector<int> exps(reg.VarIndep.size(), 0);
                            exps[pos] = 1;
                            replacements[vi].push_back({exps, int64_t(1)});
                        } else {
                            // Variable not in VarRule or VarIndep → constant 0
                        }
                    }
                }
                for (const auto& gbPoly : compGb) {
                    auto poly = PolyArith::parseSingularPolynomial(gbPoly, allVars, modulus);
                    auto projected = PolyArith::polySubstituteMulti(
                        poly, replacements, modulus, (int)reg.VarIndep.size());
                    PolyArith::canonicalize(projected, modulus);
                    if (!projected.empty()) {
                        vargenGb.push_back(
                            PolyArith::polyToSingularString(projected, reg.VarIndep));
                    }
                }
                computeMonomialBasisMatrix(vargenGb, reg.VarIndep,
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
// ==========================================

inline std::vector<RegionData> solveAllSectors(
    const IBPEqGenerator::IBPEquations& ibp,
    const std::vector<int>& topSector,
    int64_t modulus)
{
    int ne = ibp.ne;
    int nl = ibp.nl;
    std::vector<RegionData> allResults;
    auto subsectors = generateSubsectors(topSector, nl);

    for (const auto& sub : subsectors) {
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

        auto regions = solveRegion(abEqs, sub, ne, modulus);
        for (auto& reg : regions) {
            allResults.push_back(std::move(reg));
        }
    }

    return allResults;
}

} // namespace RegionSolver
#endif
