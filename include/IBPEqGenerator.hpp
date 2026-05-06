#ifndef IBP_EQ_GENERATOR_HPP
#define IBP_EQ_GENERATOR_HPP

// ============================================================
// Module B: IBP Equation Generator — C++ orchestration with
// Singular back-end for symbolic algebra (Phase 2).
//
// Pipeline:
//   FamilyDef  ──[parsing]──→  ScalarProductBasis
//              ──[Singular]──→ SP2PD transformation
//              ──[Singular]──→ IBP + LI derivative matrix
//              ──[assembly]──→ IBPEquations (g-operator form)
// ============================================================

#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <cctype>
#include <stdexcept>
#include <functional>

#include "FamilyConfig.hpp"
#include "SingularRunner.hpp"

namespace IBPEqGenerator {

// ============================================================
// Step 0 — Propagator parser
// ============================================================

// Parse a single scalar product token like "l1^2", "-2*l1*p2", "msq"
struct SPTerm {
    int64_t coeff = 0;      // integer coefficient
    int li = 0;             // loop momentum index (1-based), 0 if none
    int lj = 0;             // second loop index for l_i·l_j, 0 if p or no loop
    int pk = 0;             // external momentum index (1-based), 0 if loop or none
    std::string param;      // parameter name like "msq", "s", etc. (empty if none)
};

// Full expansion: propagator = Σ coeff_k × sp_k + const_k
struct PropagatorExpansion {
    std::vector<SPTerm> terms; // all terms in the expansion
};

// Scalar product basis for one family
struct ScalarProductBasis {
    int nl, nE, ne;
    // llIdx[i][j] = index of l_i·l_j in the basis (1-based, for i≤j)
    // lpIdx[i][j] = index of l_i·p_j
    // Total: nll = nl*(nl+1)/2, nlp = nl*nE, M = nll + nlp
    std::vector<std::vector<int>> llIdx; // [1+nl][1+nl], 0 = unused
    std::vector<std::vector<int>> lpIdx; // [1+nl][1+nE]
    int nSP = 0; // total scalar products

    // Propagator expansions
    std::vector<PropagatorExpansion> pdExpansions;

    // Scalar product names for Singular output
    std::string spName(int idx) const {
        return "sp" + std::to_string(idx);
    }
};

// ============================================================
// Propagator auto-expansion parser
// ============================================================

// Find index of a momentum name. Returns {isLoop, idx} where idx is 1-based.
inline std::pair<bool, int> findMomentum(const std::string& name,
                                          const std::vector<std::string>& loopMomenta,
                                          const std::vector<std::string>& externalMomenta) {
    for (int i = 0; i < (int)loopMomenta.size(); ++i)
        if (loopMomenta[i] == name) return {true, i + 1};
    for (int i = 0; i < (int)externalMomenta.size(); ++i)
        if (externalMomenta[i] == name) return {false, i + 1};
    return {false, -1};
}

// Parse a linear combination of momenta: "l1+p", "k1-p1", "l1+l2+p", "k1-p1-p2-p3"
// Returns vector of {coefficient, momentum_name}
inline std::vector<std::pair<int, std::string>> parseMomentumLinearCombination(const std::string& expr) {
    std::vector<std::pair<int, std::string>> result;
    std::string current;
    int sign = 1;

    // The first term has no leading sign (it's positive)
    for (size_t i = 0; i <= expr.size(); ++i) {
        char c = (i < expr.size()) ? expr[i] : '\0';
        if (c == '+' || c == '-' || c == '\0') {
            if (!current.empty()) {
                // Trim whitespace (both ends)
                while (!current.empty() && current.back() == ' ') current.pop_back();
                size_t pos = 0;
                while (pos < current.size() && current[pos] == ' ') ++pos;
                result.push_back({sign, current.substr(pos)});
                current.clear();
            }
            sign = (c == '-') ? -1 : 1;
        } else {
            current += c;
        }
    }
    return result;
}

// Expand (Σ c_i v_i)^2 into scalar products and kinematic invariants
// Returns: {spCoeffs (sp_index → coeff), constExpr (parameter expression string)}
struct MomentumSquareExpansion {
    std::map<int, int64_t> spCoeffs;       // scalar product index → coefficient
    std::map<std::string, int64_t> kinCoeffs; // kinematic key → coefficient (e.g., "p^2" → -2)
};

inline MomentumSquareExpansion expandMomentumSquare(
    const std::vector<std::pair<int, std::string>>& momenta,
    const ScalarProductBasis& spb,
    const std::vector<std::string>& loopMomenta,
    const std::vector<std::string>& externalMomenta)
{
    MomentumSquareExpansion result;
    int n = (int)momenta.size();

    // For each momentum, determine {isLoop, index}
    struct MomInfo { bool isLoop; int idx; int coeff; std::string name; };
    std::vector<MomInfo> info;
    for (auto& [coeff, name] : momenta) {
        auto [isLoop, idx] = findMomentum(name, loopMomenta, externalMomenta);
        info.push_back({isLoop, idx, coeff, name});
    }

    // Self terms: c_i^2 * (v_i·v_i)
    for (int i = 0; i < n; ++i) {
        int64_t c2 = (int64_t)info[i].coeff * info[i].coeff;
        if (info[i].isLoop) {
            int li = info[i].idx;
            int spIdx = spb.llIdx[li][li];
            result.spCoeffs[spIdx] += c2;
        } else {
            // External momentum squared → kinematic invariant
            std::string key = info[i].name + "^2";
            result.kinCoeffs[key] += c2;
        }
    }

    // Cross terms: 2*c_i*c_j * (v_i·v_j) for i < j
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int64_t cij = 2 * (int64_t)info[i].coeff * info[j].coeff;
            if (info[i].isLoop && info[j].isLoop) {
                int li = info[i].idx, lj = info[j].idx;
                if (li > lj) std::swap(li, lj);
                int spIdx = spb.llIdx[li][lj];
                result.spCoeffs[spIdx] += cij;
            } else if (info[i].isLoop && !info[j].isLoop) {
                int spIdx = spb.lpIdx[info[i].idx][info[j].idx];
                result.spCoeffs[spIdx] += cij;
            } else if (!info[i].isLoop && info[j].isLoop) {
                int spIdx = spb.lpIdx[info[j].idx][info[i].idx];
                result.spCoeffs[spIdx] += cij;
            } else {
                // Both external momenta
                std::string key = info[i].name + "*" + info[j].name;
                result.kinCoeffs[key] += cij;
            }
        }
    }
    return result;
}

// Parse a single squared expression (either "name^2" or "(expr)^2")
// sign: overall sign multiplier
inline void parseSquaredTerm(const std::string& expr, int sign,
                              const ScalarProductBasis& spb,
                              const std::vector<std::string>& loopMomenta,
                              const std::vector<std::string>& externalMomenta,
                              std::vector<int64_t>& spCoeffs,
                              std::map<std::string, int64_t>& constCoeffs) {
    std::string inner;
    if (!expr.empty() && expr[0] == '(') {
        // Compound squared: (l1+p)^2
        size_t end = expr.find(')');
        if (end == std::string::npos) {
            throw std::runtime_error("Unmatched parenthesis in propagator: " + expr);
        }
        inner = expr.substr(1, end - 1);
    } else {
        // Simple squared: l1^2, k1^2, p^2
        size_t hat = expr.find('^');
        if (hat != std::string::npos) {
            inner = expr.substr(0, hat);
        } else {
            inner = expr; // just a param name like "msq"
        }
    }

    // Check if it's a compound expression (contains + or - operators)
    bool isSimple = (inner.find('+') == std::string::npos && inner.find('-') == std::string::npos);

    if (!isSimple) {
        // Compound: expand the linear combination squared
        auto momenta = parseMomentumLinearCombination(inner);
        auto expansion = expandMomentumSquare(momenta, spb, loopMomenta, externalMomenta);
        for (auto& [spIdx, coeff] : expansion.spCoeffs)
            spCoeffs[spIdx] += sign * coeff;
        for (auto& [key, coeff] : expansion.kinCoeffs)
            constCoeffs[key] += sign * coeff;
    } else {
        // Simple: could be a single momentum or a parameter
        auto [isLoop, idx] = findMomentum(inner, loopMomenta, externalMomenta);
        if (idx > 0 && isLoop) {
            // Loop momentum squared: l_i^2 → scalar product
            spCoeffs[spb.llIdx[idx][idx]] += sign;
        } else if (idx > 0 && !isLoop) {
            // External momentum squared: p_i^2 → kinematic invariant
            constCoeffs[inner + "^2"] += sign;
        } else {
            // Parameter name like "msq"
            constCoeffs[inner] += sign;
        }
    }
}

// Parse a full propagator string into scalar product coefficients and constant expression.
// The constant expression is a string like "-s-msq" representing Σ constCoeffs.
inline void parsePropagator(const std::string& propStr,
                             const ScalarProductBasis& spb,
                             const std::vector<std::string>& loopMomenta,
                             const std::vector<std::string>& externalMomenta,
                             const std::map<std::string, std::string>& kinematicRules,
                             std::vector<int64_t>& spCoeffs,    // length nSP
                             std::string& constExpr)             // e.g., "-s-msq"
{
    spCoeffs.assign(spb.nSP + 1, 0); // 1-indexed: spCoeffs[1..nSP]
    std::map<std::string, int64_t> rawConstCoeffs; // raw kinematic key → coefficient

    // Tokenize: split the propagator string into signed terms
    // Strategy: walk the string; a '+' or '-' at top-level depth splits terms
    struct Term { int sign; std::string expr; };
    std::vector<Term> terms;

    int depth = 0;
    size_t i = 0;
    int sign = 1; // default sign for the first term

    while (i < propStr.size()) {
        // Skip whitespace
        while (i < propStr.size() && propStr[i] == ' ') ++i;
        if (i >= propStr.size()) break;

        // At depth 0, check for sign
        if (depth == 0 && (propStr[i] == '+' || propStr[i] == '-')) {
            sign = (propStr[i] == '-') ? -1 : 1;
            ++i;
            // Skip whitespace after sign
            while (i < propStr.size() && propStr[i] == ' ') ++i;
        }

        // Read expression until next top-level operator or end
        std::string expr;
        while (i < propStr.size()) {
            char c = propStr[i];
            if (c == '(') {
                ++depth;
                expr += c;
            } else if (c == ')') {
                --depth;
                expr += c;
            } else if (depth == 0 && (c == '+' || c == '-')) {
                break; // next term
            } else {
                expr += c;
            }
            ++i;
        }

        // Trim trailing whitespace
        while (!expr.empty() && expr.back() == ' ') expr.pop_back();

        if (!expr.empty())
            terms.push_back({sign, expr});
    }

    // Process each term
    std::map<std::string, int64_t> termConstCoeffs;
    for (auto& [tSign, expr] : terms) {
        // Check if expr contains ^2 or is a squared expression
        bool isSquared = false;

        // Contains ^2 somewhere (may be inside or outside parens)
        if (expr.find("^2") != std::string::npos) {
            isSquared = true;
        }

        if (isSquared) {
            parseSquaredTerm(expr, tSign, spb, loopMomenta, externalMomenta,
                            spCoeffs, termConstCoeffs);
        } else {
            // It's a parameter name without ^2
            termConstCoeffs[expr] += tSign;
        }
    }

    // Apply kinematic rules to raw constant coefficients
    // e.g., "p^2" → "s", "p1*p2" → "(s3-s1-s2)/2"
    // If no rule found, keep the key as-is (it IS a parameter name like "msq")
    for (auto& [key, coeff] : termConstCoeffs) {
        auto it = kinematicRules.find(key);
        std::string resolvedKey;
        if (it != kinematicRules.end()) {
            resolvedKey = it->second;
            // Normalize for Singular: remove spaces, insert * for implicit
            // multiplication (e.g., "2 t" → "2*t"). Singular does not parse
            // "2t" as 2*t — it treats it as a single identifier.
            std::string norm;
            for (char ch : resolvedKey) if (ch != ' ') norm += ch;
            resolvedKey.clear();
            for (size_t ci = 0; ci < norm.size(); ++ci) {
                if (ci > 0 && isdigit(norm[ci-1]) && isalpha(norm[ci]))
                    resolvedKey += '*';
                resolvedKey += norm[ci];
            }
        } else {
            resolvedKey = key; // parameter name, like "msq"
        }
        // Accumulate into rawConstCoeffs with the resolved key
        rawConstCoeffs[resolvedKey] += coeff;
    }

    // Build constant expression string
    std::stringstream css;
    bool first = true;
    for (auto& [key, coeff] : rawConstCoeffs) {
        if (coeff == 0) continue;
        if (!first) {
            css << (coeff > 0 ? "+" : "");
        }
        if (coeff == 1) {
            if (!first) css << "1*"; // needs explicit coeff if not first
            css << key;
        } else if (coeff == -1) {
            css << "-" << key;
        } else {
            css << coeff << "*" << key;
        }
        first = false;
    }
    constExpr = css.str();
    if (constExpr.empty()) constExpr = "0";
}

// ============================================================
// Tokeniser / parser for propagator strings
// ============================================================

inline ScalarProductBasis buildScalarProductBasis(const FamilyDef& fam) {
    ScalarProductBasis spb;
    spb.nl = fam.nLoop();
    spb.nE = fam.nExt();
    spb.ne = fam.nProp();

    int nl = spb.nl, nE = spb.nE;

    // Build llIdx and lpIdx
    spb.llIdx.resize(nl + 1, std::vector<int>(nl + 1, 0));
    spb.lpIdx.resize(nl + 1, std::vector<int>(nE + 1, 0));

    int idx = 0;
    // l_i·l_j for i ≤ j
    for (int i = 1; i <= nl; ++i)
        for (int j = i; j <= nl; ++j)
            spb.llIdx[i][j] = ++idx;
    // l_i·p_j
    for (int i = 1; i <= nl; ++i)
        for (int j = 1; j <= nE; ++j)
            spb.lpIdx[i][j] = ++idx;
    spb.nSP = idx;

    // Parse each propagator — compute C matrix rows
    spb.pdExpansions.resize(spb.ne);

    return spb;
}

// ============================================================
// Step 1 — Singular ring generation
// ============================================================

// Build the Singular ring variable string
inline std::string buildRingVars(const ScalarProductBasis& spb,
                                  const std::map<std::string, std::string>& numericParams) {
    (void)numericParams;
    std::stringstream ss;

    // Dimension parameter
    ss << "d";

    // Scalar product variables (indexed: sp(1)..sp(nSP))
    if (spb.nSP > 0)
        ss << ",sp(1.." << spb.nSP << ")";

    // z variables (propagator denominators, indexed: z(1)..z(ne))
    if (spb.ne > 0)
        ss << ",z(1.." << spb.ne << ")";

    return ss.str();
}

// ============================================================
// Step 2 — SP2PD computation (Singular script with real coefficients)
// ============================================================

// Build the Singular script that:
//   1. Expands each propagator in the scalar-product basis
//   2. Solves the linear system: sp[*] = f(z[*], params)
//   3. Outputs the SP2PD transformation matrix
//
// The key: D_i = -z_i + const_i (where const_i comes from masses, kinematic invars)
// In the sp basis: D_i = Σ_j C_{ij} × sp_j + const_i
// So: Σ_j C_{ij} × sp_j + const_i = -z_i → Σ_j C_{ij} × sp_j = -z_i - const_i
// Invert: sp_j = Σ_i Inv_{ji} × (-z_i - const_i)

inline std::string buildSP2PDScript(const ScalarProductBasis& spb,
                                     const std::vector<std::string>& propagatorStrings,
                                     int64_t modulus,
                                     const std::map<std::string, std::string>& numericParams,
                                     const std::vector<std::string>& loopMomenta,
                                     const std::vector<std::string>& externalMomenta,
                                     const std::map<std::string, std::string>& kinematicRules) {
    std::stringstream ss;

    int ne = spb.ne, nSP = spb.nSP;

    ss << "// SP2PD: Scalar Product → Propagator Denominator\n";
    ss << "// Build ring with scalar products + z-vars + parameters\n\n";

    ss << "LIB \"linalg.lib\";\n";
    ss << "ring rSP = " << modulus << ", (" << buildRingVars(spb, numericParams) << "), (dp);\n\n";

    // Define kinematic/numeric parameters as ring coefficients
    for (const auto& [k, v] : numericParams) {
        if (k != "d") {
            ss << "number " << k << " = " << v << ";\n";
        }
    }
    ss << "\n";

    ss << "// ========================================\n";
    ss << "// Auto-generated propagator coefficients\n";
    ss << "// ========================================\n\n";

    ss << "int ne = " << ne << ";\n";
    ss << "int nl = " << spb.nl << ";\n";
    ss << "int nE = " << spb.nE << ";\n";
    ss << "int nSP = " << nSP << ";\n\n";

    ss << "// C[i,j] = coefficient of sp[j] in propagator D[i]\n";
    ss << "matrix C[" << ne << "][" << nSP << "];\n";
    ss << "matrix Cconst[" << ne << "][1];\n\n";

    // Auto-generate coefficients by parsing each propagator
    for (int i = 0; i < ne; ++i) {
        std::vector<int64_t> spCoeffs; // 1-indexed
        std::string constExpr;
        parsePropagator(propagatorStrings[i], spb, loopMomenta, externalMomenta,
                       kinematicRules, spCoeffs, constExpr);

        ss << "// D[" << (i+1) << "] = " << propagatorStrings[i] << "\n";
        for (int j = 1; j <= nSP; ++j) {
            if (spCoeffs[j] != 0)
                ss << "C[" << (i+1) << "," << j << "] = " << spCoeffs[j] << ";\n";
        }
        ss << "Cconst[" << (i+1) << ",1] = " << constExpr << ";\n\n";
    }

    ss << "// ========================================\n";
    ss << "// Linear solve: sp[*] = Inv × (-z[*] - const[*])\n";
    ss << "//   C * sp + Cconst = -z\n";
    ss << "//   C * sp = -z - Cconst\n";
    ss << "//   sp = C^{-1} * (-z - Cconst)  [if ne == nSP]\n";
    ss << "// ========================================\n\n";

    ss << "// SP2PD rules: sp[j] expressed in z[*] basis\n";
    ss << "ideal sp2pd;\n";
    ss << "string sp2pd_result = \"\";\n";
    ss << "if (ne == nSP) {\n";
    ss << "  matrix C_inv = inverse(C);\n";
    ss << "  sp2pd_result = string(nSP) + \"|\" + string(ne);\n";
    ss << "  // sp[j] = Σ_i C_inv[j,i] * (-z[i] - Cconst[i,1])\n";
    ss << "  for (int j = 1; j <= nSP; j++) {\n";
    ss << "    poly rule_j;\n";
    ss << "    for (int i = 1; i <= ne; i++) {\n";
    ss << "      rule_j = rule_j + C_inv[j,i] * (-z(i) - Cconst[i,1]);\n";
    ss << "    }\n";
    ss << "    sp2pd[j] = sp(j) - rule_j;\n";
    ss << "    sp2pd_result = sp2pd_result + \"|\" + string(rule_j);\n";
    ss << "  }\n";
    ss << "} else {\n";
    ss << "  sp2pd_result = \"FAIL: ne != nSP\";\n";
    ss << "}\n";

    return ss.str();
}

// ============================================================
// Step 2b — Execute SP2PD via Singular and parse output
// ============================================================

struct SP2PDResult {
    bool success = false;
    int ne, nSP;
    // sp2pd[j] = Σ_i coeff_ji * z_i + const_j
    // coeff[j][i] = coefficient of z[i] in sp[j] (1-based indexing)
    std::vector<std::vector<std::string>> coeffs; // [nSP+1][ne+1] — polynomial strings
    std::vector<std::string> constTerms;           // [nSP+1] — constant term strings
    std::vector<std::string> fullRules;            // [nSP+1] — full expression string
};

inline SP2PDResult runSP2PD(const FamilyDef& fam) {
    auto spb = buildScalarProductBasis(fam);
    std::string script = buildSP2PDScript(spb, fam.propagators, fam.modulus,
                                           fam.numeric, fam.loopMomenta,
                                           fam.externalMomenta, fam.kinematicRules);

    SP2PDResult sp2pd;
    sp2pd.ne = spb.ne;
    sp2pd.nSP = spb.nSP;

    // Run Singular via the shared runner
    std::map<std::string, std::string> result;
    try {
        result = SingularRunner::runSingular(script, {"sp2pd_result"}, "/tmp");
    } catch (const std::exception& e) {
        std::cerr << "[SP2PD] Singular execution failed: " << e.what() << std::endl;
        return sp2pd;
    }

    auto it = result.find("sp2pd_result");
    if (it == result.end()) {
        std::cerr << "[SP2PD] No output from Singular for family " << fam.name << std::endl;
        return sp2pd;
    }

    std::string data = it->second;
    // Trim whitespace
    while (!data.empty() && (data[0] == ' ' || data[0] == '\n')) data.erase(0, 1);
    while (!data.empty() && (data.back() == ' ' || data.back() == '\n')) data.pop_back();

    if (data.find("FAIL") == 0) {
        std::cerr << "[SP2PD] ne != nSP for family " << fam.name << std::endl;
        return sp2pd;
    }

    // Parse: pipe-delimited: nSP|ne|rule1|rule2|...|ruleN
    std::vector<std::string> fields;
    {
        std::stringstream lss(data);
        std::string field;
        while (std::getline(lss, field, '|')) {
            // Trim whitespace
            while (!field.empty() && field[0] == ' ') field.erase(0, 1);
            while (!field.empty() && (field.back() == '\r' || field.back() == '\n' || field.back() == ' '))
                field.pop_back();
            fields.push_back(field);
        }
    }

    if (fields.size() < 3) {
        std::cerr << "[SP2PD] Unexpected output format from Singular (got "
                  << fields.size() << " fields)" << std::endl;
        return sp2pd;
    }

    int parsedNSP = std::stoi(fields[0]);
    int parsedNE  = std::stoi(fields[1]);

    if (parsedNSP != spb.nSP || parsedNE != spb.ne) {
        std::cerr << "[SP2PD] Dimension mismatch: got nSP=" << parsedNSP
                  << " ne=" << parsedNE << std::endl;
        return sp2pd;
    }

    sp2pd.success = true;
    sp2pd.fullRules.resize(spb.nSP + 1, "0");

    for (int j = 1; j <= spb.nSP; ++j) {
        if (j + 1 < (int)fields.size()) {
            sp2pd.fullRules[j] = fields[j + 1]; // fields[0]=nSP, fields[1]=ne
        }
    }

    return sp2pd;
}

// ============================================================
// Step 3 — IBP derivative matrix computation
// ============================================================

// For each IBP identity: ∫ d^d k ∂/∂l_j^μ ( q_k^μ / ∏ D_i^{n_i} ) = 0
//
// The derivative matrix derivL[i][j][k]:
//   i = propagator index (1..ne)
//   j = loop momentum index (1..nl)
//   k = momentum component index (1..nl+nE, loops first then externals)
//
// derivL[i][j][k] = (∂D_i/∂l_j^μ) · q_k^μ  expressed in SP2PD basis
//
// After SP2PD substitution: deriv → linear combination of z[*] + const.

// Build a Singular script that computes the IBP derivative matrix.
// C++ pre-computes derivative contributions in terms of scalar product indices;
// Singular does polynomial substitution via SP2PD rules to get z-basis expressions.
inline std::string buildIBPDerivativeScript(const ScalarProductBasis& spb,
                                              const std::vector<std::string>& propagatorStrings,
                                              int64_t modulus,
                                              const std::map<std::string, std::string>& numericParams,
                                              const std::vector<std::string>& loopMomenta,
                                              const std::vector<std::string>& externalMomenta,
                                              const std::map<std::string, std::string>& kinematicRules) {
    std::stringstream ss;

    int nl = spb.nl, nE = spb.nE, ne = spb.ne, nSP = spb.nSP;
    int nTotal = nl + nE;

    // Build reverse lookup: sp index → (type, a, b)
    // type 0 = ll[a][b], type 1 = lp[a][b]
    struct SPInfo { int type; int a; int b; };
    std::vector<SPInfo> spInfo(nSP + 1, {0, 0, 0});
    for (int a = 1; a <= nl; ++a)
        for (int b = a; b <= nl; ++b)
            spInfo[spb.llIdx[a][b]] = {0, a, b};
    for (int a = 1; a <= nl; ++a)
        for (int b = 1; b <= nE; ++b)
            spInfo[spb.lpIdx[a][b]] = {1, a, b};

    // Parse propagator coefficients
    std::vector<std::vector<int64_t>> C_rows(ne, std::vector<int64_t>(nSP + 1, 0));

    for (int i = 0; i < ne; ++i) {
        std::string constExpr;
        parsePropagator(propagatorStrings[i], spb, loopMomenta, externalMomenta,
                       kinematicRules, C_rows[i], constExpr);
    }

    // For each (i,j,k), pre-compute derivative as linear combination of sp indices
    // deriv_sp[idx][s] = coefficient of sp[s] in the derivative
    // deriv_const[idx] = string of constant term (kinematic invariants)
    int nibp = ne * nl * nTotal;
    std::vector<std::vector<int64_t>> derivSp(nibp);
    std::vector<std::string> derivConst(nibp);

    int idx = 0;
    for (int i = 0; i < ne; ++i) {
        for (int j = 1; j <= nl; ++j) {
            for (int k = 1; k <= nTotal; ++k) {
                std::vector<int64_t> dSp(nSP + 1, 0);
                std::map<std::string, int64_t> constTerms;

                // Σ_s C[i,s] * (∂sp[s]/∂l_j · q_k)
                for (int s = 1; s <= nSP; ++s) {
                    int64_t c = C_rows[i][s];
                    if (c == 0) continue;

                    auto& info = spInfo[s];
                    if (info.type == 0) {
                        // sp[s] = l_a·l_b (a≤b)
                        // ∂(l_a·l_b)/∂l_j · q_k = δ_aj*(l_b·q_k) + δ_bj*(l_a·q_k)
                        // Note: when a==b, both δ terms coincide, giving 2*δ_aj*(l_a·q_k)
                        int a = info.a, b = info.b;
                        int64_t factor = (a == b) ? 2 : 1;
                        if (j == a) {
                            if (k <= nl) {
                                int mn = std::min(b, k), mx = std::max(b, k);
                                dSp[spb.llIdx[mn][mx]] += c * factor;
                            } else {
                                dSp[spb.lpIdx[b][k - nl]] += c * factor;
                            }
                        }
                        if (a != b && j == b) {
                            if (k <= nl) {
                                int mn = std::min(a, k), mx = std::max(a, k);
                                dSp[spb.llIdx[mn][mx]] += c;
                            } else {
                                dSp[spb.lpIdx[a][k - nl]] += c;
                            }
                        }
                    } else {
                        // sp[s] = l_a·p_b
                        int a = info.a, b = info.b;
                        if (j == a) {
                            if (k <= nl) {
                                dSp[spb.lpIdx[k][b]] += c;
                            } else {
                                // p_b·p_{k-nl}: kinematic invariant
                                std::string e1 = externalMomenta[b - 1];
                                std::string e2 = externalMomenta[k - nl - 1];
                                std::string key;
                                if (e1 == e2) key = e1 + "^2";
                                else if (e1 < e2) key = e1 + "*" + e2;
                                else key = e2 + "*" + e1;
                                constTerms[key] += c;
                            }
                        }
                    }
                }

                derivSp[idx] = dSp;
                // Apply kinematicRules to constant term keys
                std::map<std::string, int64_t> resolvedConstTerms;
                for (auto& [key, coeff] : constTerms) {
                    auto it = kinematicRules.find(key);
                    std::string resolvedKey;
                    if (it != kinematicRules.end()) {
                        resolvedKey = it->second;
                        // Normalize for Singular: remove spaces, insert * for
                        // implicit multiplication (e.g., "2 t" → "2*t")
                        std::string norm;
                        for (char ch : resolvedKey) if (ch != ' ') norm += ch;
                        resolvedKey.clear();
                        for (size_t ci = 0; ci < norm.size(); ++ci) {
                            if (ci > 0 && isdigit(norm[ci-1]) && isalpha(norm[ci]))
                                resolvedKey += '*';
                            resolvedKey += norm[ci];
                        }
                    } else {
                        resolvedKey = key;
                    }
                    resolvedConstTerms[resolvedKey] += coeff;
                }
                // Build constant expression string
                std::stringstream css;
                bool first = true;
                for (auto& [key, coeff] : resolvedConstTerms) {
                    if (coeff == 0) continue;
                    if (!first) css << (coeff > 0 ? "+" : "");
                    if (coeff == 1) { if (!first) css << "+"; }
                    else if (coeff == -1) css << "-";
                    else if (coeff > 0) css << coeff << "*";
                    else css << coeff << "*";
                    css << key;
                    first = false;
                }
                derivConst[idx] = css.str();
                if (derivConst[idx].empty()) derivConst[idx] = "0";

                ++idx;
            }
        }
    }


    // Now build the Singular script
    ss << "// ========================================\n";
    ss << "// IBP Derivative Matrix Computation\n";
    ss << "// ========================================\n\n";

    ss << "LIB \"linalg.lib\";\n";
    ss << "ring rDeriv = " << modulus << ", (" << buildRingVars(spb, numericParams) << "), (dp);\n\n";

    // Numeric parameters
    for (const auto& [k, v] : numericParams) {
        if (k != "d") ss << "number " << k << " = " << v << ";\n";
    }
    ss << "\n";

    ss << "int ne = " << ne << ";\n";
    ss << "int nl = " << nl << ";\n";
    ss << "int nE = " << nE << ";\n";
    ss << "int nSP = " << nSP << ";\n";
    ss << "int nTotal = " << nTotal << ";\n";
    ss << "int nibp = ne * nl * nTotal;\n\n";

    // C matrix and Cconst
    ss << "matrix C[" << ne << "][" << nSP << "];\n";
    ss << "matrix Cconst[" << ne << "][1];\n";
    for (int i = 0; i < ne; ++i) {
        for (int s = 1; s <= nSP; ++s)
            if (C_rows[i][s] != 0)
                ss << "C[" << (i+1) << "," << s << "] = " << C_rows[i][s] << ";\n";
        // Re-parse constExpr for the ring context
        std::string ce;
        parsePropagator(propagatorStrings[i], spb, loopMomenta, externalMomenta,
                       kinematicRules, C_rows[i], ce);
        ss << "Cconst[" << (i+1) << ",1] = " << ce << ";\n";
    }
    ss << "\n";

    // SP2PD polynomials (use list to store polynomials)
    ss << "matrix C_inv = inverse(C);\n";
    ss << "list spL;\n";
    ss << "for (int ss = 1; ss <= nSP; ss++) {\n";
    ss << "  poly p = 0;\n";
    ss << "  for (int ii = 1; ii <= ne; ii++) {\n";
    ss << "    p = p + C_inv[ss,ii] * (-z(ii) - Cconst[ii,1]);\n";
    ss << "  }\n";
    ss << "  spL = spL + list(p);\n";
    ss << "}\n";
    ss << "\n";

    // Derivative assembly
    ss << "// Derivative: for each (i,j,k), deriv = Σ_s dSp[s] * poly(spL[s] + const\n";
    ss << "string deriv_result = \"\";\n";
    ss << "int idx = 0;\n";
    ss << "for (int i = 1; i <= ne; i++) {\n";
    ss << "  for (int j = 1; j <= nl; j++) {\n";
    ss << "    for (int k = 1; k <= nTotal; k++) {\n";
    ss << "      idx++;\n";
    ss << "      poly deriv = 0;\n";

    idx = 0;
    for (int i = 0; i < ne; ++i) {
        for (int j = 1; j <= nl; ++j) {
            for (int k = 1; k <= nTotal; ++k) {
                ss << "      if (idx == " << (idx+1) << ") {\n";
                // Add scalar product contributions
                bool hasTerms = false;
                for (int s = 1; s <= nSP; ++s) {
                    if (derivSp[idx][s] != 0) {
                        if (derivSp[idx][s] == 1)
                            ss << "        deriv = deriv + poly(spL[" << s << "]);\n";
                        else if (derivSp[idx][s] == -1)
                            ss << "        deriv = deriv - poly(spL[" << s << "]);\n";
                        else if (derivSp[idx][s] > 0)
                            ss << "        deriv = deriv + " << derivSp[idx][s] << "*poly(spL[" << s << "]);\n";
                        else
                            ss << "        deriv = deriv - " << -derivSp[idx][s] << "*poly(spL[" << s << "]);\n";
                        hasTerms = true;
                    }
                }
                // Add constant term
                if (derivConst[idx] != "0") {
                    ss << "        deriv = deriv + (" << derivConst[idx] << ");\n";
                    hasTerms = true;
                }
                if (!hasTerms) ss << "        deriv = 0;\n";
                ss << "      }\n";
                ++idx;
            }
        }
    }

    ss << "      deriv_result = deriv_result + string(deriv) + \"|\";\n";
    ss << "    }\n";
    ss << "  }\n";
    ss << "}\n";
    ss << "deriv_result = string(nibp) + \"|\" + deriv_result;\n";

    return ss.str();
}

// Execute the derivative computation via Singular
struct DerivResult {
    bool success = false;
    int nibp;
    std::vector<std::string> derivPolynomials; // one per (i,j,k) in z-basis
};

inline DerivResult runIBPDerivatives(const FamilyDef& fam) {
    auto spb = buildScalarProductBasis(fam);
    std::string script = buildIBPDerivativeScript(spb, fam.propagators, fam.modulus,
                                                   fam.numeric, fam.loopMomenta,
                                                   fam.externalMomenta, fam.kinematicRules);

    DerivResult result;
    result.nibp = fam.nLoop() * (fam.nLoop() + fam.nExt()) * fam.nProp();

    std::map<std::string, std::string> output;
    try {
        output = SingularRunner::runSingular(script, {"deriv_result"}, "/tmp");
    } catch (const std::exception& e) {
        std::cerr << "[IBP-Deriv] Singular error: " << e.what() << std::endl;
        return result;
    }

    auto it = output.find("deriv_result");
    if (it == output.end()) return result;

    std::string data = it->second;
    while (!data.empty() && data[0] == ' ') data.erase(0, 1);
    while (!data.empty() && (data.back() == '\r' || data.back() == '\n' || data.back() == ' '))
        data.pop_back();

    // Parse: nibp|poly1|poly2|...|polyN
    std::vector<std::string> fields;
    {
        std::stringstream lss(data);
        std::string field;
        while (std::getline(lss, field, '|')) {
            while (!field.empty() && field[0] == ' ') field.erase(0, 1);
            while (!field.empty() && (field.back() == '\r' || field.back() == '\n' || field.back() == ' '))
                field.pop_back();
            if (!field.empty()) fields.push_back(field);
        }
    }

    if (fields.empty()) return result;
    int nibpRead = std::stoi(fields[0]);
    if (nibpRead != result.nibp) {
        std::cerr << "[IBP-Deriv] Expected nibp=" << result.nibp << " got " << nibpRead << std::endl;
        return result;
    }

    result.success = true;
    for (size_t i = 1; i < fields.size() && (int)(i-1) < result.nibp; ++i)
        result.derivPolynomials.push_back(fields[i]);

    return result;
}

// ============================================================
// Step 3b — IBP equation assembly from derivatives (in C++)
// ============================================================
//
// MMA formula (genIBP in LIEUtility.wl):
//   eqmon = d * δ_jk - Σ_i n_i * deriv[i,j,k] / z_i
//
// After multiplying by ∏_m z_m (to clear denominators):
//   num = d * δ_jk * ∏_m z_m - Σ_i n_i * deriv[i,j,k] * ∏_{m≠i} z_m
//
// Each monomial  c * ∏_m z_m^{a_m}  maps via mon2F (+ sign flip) to:
//   c * g^{(a_1, ..., a_ne)}    where g^α shifts ν → ν + α
//
// The n_i factor (symbolic exponent) becomes ν_i in the g-operator coefficient.

// z-Monomial: product of z_i^{exp_i}
struct ZMonomial {
    std::vector<int> exps;  // [ne], exponent of z[i+1]
    int64_t coeff = 0;      // coefficient (over finite field)
    int nIdx = 0;           // 0 = no n-factor, i = multiplied by n_i (1-based)
    bool hasD = false;      // true if multiplied by d
};

// Parse a Singular polynomial string into monomials.
// Format: "c1*z(1)^e1*z(2)^e2+..." or "c1*z(1)^e1*z(2)^e2-..."
// The polynomial may contain: numeric coefficients, z(i)^n, d, n(i), s, msq, etc.
// Returns list of ZMonomial structs.
inline std::vector<ZMonomial> parseZPolynomial(const std::string& poly, int ne, int64_t modulus) {
    std::vector<ZMonomial> result;
    if (poly.empty() || poly == "0") return result;

    std::string s = poly;
    // Remove spaces
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());

    // Tokenize: split at '+' or '-' that are NOT inside parentheses
    struct RawTerm { int sign; std::string expr; };
    std::vector<RawTerm> rawTerms;

    int depth = 0;
    size_t start = 0;
    int sign = 1;
    // Check if first char is a sign
    if (!s.empty() && s[0] == '-') { sign = -1; start = 1; }
    else if (!s.empty() && s[0] == '+') { sign = 1; start = 1; }

    std::string current;
    for (size_t i = start; i <= s.size(); ++i) {
        char c = (i < s.size()) ? s[i] : '\0';
        if (c == '(') ++depth;
        else if (c == ')') --depth;
        if (depth == 0 && (c == '+' || c == '-' || c == '\0')) {
            if (!current.empty()) {
                rawTerms.push_back({sign, current});
                current.clear();
            }
            sign = (c == '-') ? -1 : 1;
        } else {
            current += c;
        }
    }
    if (!current.empty()) rawTerms.push_back({sign, current});

    // Parse each raw term
    for (auto& [tsign, expr] : rawTerms) {
        ZMonomial zm;
        zm.exps.assign(ne, 0);
        zm.coeff = tsign;  // start with ±1
        zm.nIdx = 0;
        zm.hasD = false;

        // Parse factors in the term: split by '*'
        std::stringstream fss(expr);
        std::string factor;
        bool hasNumericFactor = false;

        while (std::getline(fss, factor, '*')) {
            if (factor.empty()) continue;

            // Check for numeric coefficient (including negative, which was already handled)
            bool isNumeric = true;
            for (char ch : factor) {
                if (!std::isdigit(ch) && ch != '-') { isNumeric = false; break; }
            }
            if (isNumeric && !factor.empty() && factor != "-") {
                int64_t val = std::stoll(factor);
                zm.coeff *= val;
                hasNumericFactor = true;
                continue;
            }

            // Check for d (dimension)
            if (factor == "d") {
                zm.hasD = true;
                continue;
            }

            // Check for n(i)
            if (factor.size() > 2 && factor[0] == 'n' && factor[1] == '(') {
                size_t rp = factor.find(')');
                if (rp != std::string::npos) {
                    int idx = std::stoi(factor.substr(2, rp - 2));
                    zm.nIdx = idx;
                }
                continue;
            }

            // Check for z(i)^e or z(i)
            if (factor.size() >= 3 && factor[0] == 'z' && factor[1] == '(') {
                size_t rp = factor.find(')');
                if (rp != std::string::npos) {
                    int idx = std::stoi(factor.substr(2, rp - 2));
                    // Check for exponent
                    size_t hat = factor.find('^');
                    int exp = 1;
                    if (hat != std::string::npos) {
                        exp = std::stoi(factor.substr(hat + 1));
                    }
                    if (idx >= 1 && idx <= ne) {
                        zm.exps[idx - 1] += exp;
                    }
                }
                continue;
            }

            // Parameter names (s, msq, etc.) — treat as constant = factor is numeric-ish
            // Skip: parameter names like "s", "msq" are Singular numbers, already reduced mod p
        }

        // Reduce coefficient modulo modulus
        zm.coeff = ((zm.coeff % modulus) + modulus) % modulus;
        if (zm.coeff != 0) {
            result.push_back(zm);
        }
    }

    return result;
}

// mon2F conversion: map z-exponents to g-operator shift.
// In MMA: mon2F maps z_i^{a_i} → f[a_i], then sign flip f[a] → f[-a].
// So f[-a] = g^a (shift ν by +a).
// Therefore z_i^a in the numerator → g^a_i (shift ν_i up by a).
//
// gShift[j] = exponent a_j of z_j (which equals the shift amount for ν_j).
// The shift is always non-negative since we multiplied through by ∏ z_i.
struct GTerm {
    std::vector<int> gShift;  // [ne], g-operator shift vector α
    int64_t coeff = 0;         // polynomial coefficient in ν (includes n_i factors)
    int nIdx = 0;              // which n_i, or 0
    bool hasD = false;
};

inline std::vector<GTerm> mon2F(const std::vector<ZMonomial>& zMons, int ne) {
    std::vector<GTerm> result;
    for (auto& zm : zMons) {
        GTerm gt;
        gt.gShift = zm.exps;  // z_i^{a_i} → g^{a_i} (shift ν_i by +a_i)
        gt.coeff = zm.coeff;
        gt.nIdx = zm.nIdx;
        gt.hasD = zm.hasD;
        result.push_back(gt);
    }
    return result;
}

// Group g-terms by shift vector and nIdx/hasD, summing coefficients
inline std::vector<GTerm> groupGTerms(const std::vector<GTerm>& terms, int ne) {
    // Map key: shift vector + nIdx + hasD
    std::map<std::string, GTerm> grouped;
    for (auto& gt : terms) {
        std::stringstream key;
        for (int e : gt.gShift) key << e << ",";
        key << "|n" << gt.nIdx << "|d" << gt.hasD;
        auto it = grouped.find(key.str());
        if (it != grouped.end()) {
            it->second.coeff = (it->second.coeff + gt.coeff);
        } else {
            grouped[key.str()] = gt;
        }
    }
    std::vector<GTerm> result;
    for (auto& [k, v] : grouped) {
        if (v.coeff != 0) result.push_back(v);
    }
    return result;
}

// Assemble IBP equations from derivative polynomials.
//
// For each (j,k), the IBP numerator is:
//   d * δ_jk * ∏_m z_m + Σ_i n_i * (-deriv[i,j,k] * ∏_{m≠i} z_m)
//
// This multiplies through by ∏ z_i to clear the 1/z_i denominators in the
// standard IBP formula: d*δ_jk - Σ_i n_i * deriv[i,j,k] / z_i
//
// Returns: vector of IBP equations, indexed by (j,k) pair.
// Each equation is a list of GTerms in g-operator form.
struct IBPEquationG {
    int j, k;                     // (j,k) indices (1-based)
    std::vector<GTerm> gTerms;    // g-operator terms
};

inline std::vector<IBPEquationG> assembleIBPFromDerivatives(
    const DerivResult& deriv,
    int ne, int nl, int nE, int64_t modulus,
    int64_t dCoeff = 1)
{
    int nTotal = nl + nE;
    std::vector<IBPEquationG> equations;

    // Product z1 * z2 * ... * zne as exponent vector: all 1's
    std::vector<int> prodAll(ne, 1);

    int idx = 0; // index into deriv.derivPolynomials
    for (int j = 1; j <= nl; ++j) {
        for (int kk = 1; kk <= nTotal; ++kk) {
            IBPEquationG eq;
            eq.j = j;
            eq.k = kk;

            std::vector<GTerm> allTerms;

            // d-term: d * δ_jk * ∏_m z_m
            if (j == kk) {
                GTerm dt;
                dt.gShift = prodAll;
                dt.coeff = dCoeff;
                dt.hasD = true;
                allTerms.push_back(dt);
            }

            // n_i terms: n_i * (-deriv[i,j,k]) * ∏_{m≠i} z_m
            for (int i = 1; i <= ne; ++i) {
                int didx = (i - 1) * nl * nTotal + (j - 1) * nTotal + (kk - 1);
                if (didx >= (int)deriv.derivPolynomials.size()) continue;

                std::string dpoly = deriv.derivPolynomials[didx];
                auto zMons = parseZPolynomial(dpoly, ne, modulus);

                for (auto& zm : zMons) {
                    // Multiply zm by ∏_{m≠i} z_m
                    ZMonomial multiplied = zm;
                    for (int m = 0; m < ne; ++m) {
                        if (m != i - 1) multiplied.exps[m] += 1;
                    }
                    // The n_i factor: negate (MMA formula has -n_i * deriv)
                    multiplied.coeff = (-multiplied.coeff) % modulus;
                    if (multiplied.coeff < 0) multiplied.coeff += modulus;
                    multiplied.nIdx = i;

                    auto gts = mon2F({multiplied}, ne);
                    allTerms.insert(allTerms.end(), gts.begin(), gts.end());
                }
            }

            eq.gTerms = groupGTerms(allTerms, ne);
            equations.push_back(eq);
            ++idx;
        }
    }

    return equations;
}

// ============================================================
// Master script: generate IBP derivative script for a family
// ============================================================

inline std::string buildFullIBPScript(const FamilyDef& fam) {
    auto spb = buildScalarProductBasis(fam);

    std::stringstream ss;
    ss << "// ========================================\n";
    ss << "// Full IBP + LI Equation Generation\n";
    ss << "// Family: " << fam.name << "\n";
    ss << "// nLoop=" << fam.nLoop() << " nExt=" << fam.nExt()
       << " nProp=" << fam.nProp() << "\n";
    ss << "// Modulus: " << fam.modulus << "\n";
    ss << "// ========================================\n\n";

    // Only need the derivative script — it already includes SP2PD setup
    ss << buildIBPDerivativeScript(spb, fam.propagators, fam.modulus, fam.numeric,
                                   fam.loopMomenta, fam.externalMomenta, fam.kinematicRules);

    return ss.str();
}

// ============================================================
// Data structures for IBP equation output
// ============================================================

// One term in an IBP equation: coeff * g^shift, with optional n_i and d factors.
// In full g-operator form: Σ (c + Σ d_k * ν_k) * g^α.
// For now we track nIdx (which ν_i contributes) and hasD separately.
struct IBPTerm {
    std::vector<int> gShift;  // [ne], g-operator shift α
    int64_t coeff = 0;
    int nIdx = 0;    // which ν_i, 0 = no ν factor
    bool hasD = false;
};

struct IBPEquation {
    int j, k;  // (loop momentum index, momentum component index)
    std::vector<IBPTerm> terms;
};

struct IBPEquations {
    int ne, nl, nE, nibp;                     // nibp = nl * (nl+nE)
    std::vector<std::string> Alist, vlist;     // A[i], v[i] symbols
    std::vector<std::vector<int>> sectorlist;  // all subsectors
    std::vector<IBPEquation> equations;        // one per (j,k) pair
};

// ============================================================
// Subsector generation
// ============================================================

// Generate subsectors matching MMA's Subsets[activeIndices, {nl, ne}]
// Only generates subsets with sizes from nl to ne (inclusive), excluding empty set.
inline std::vector<std::vector<int>> generateSubsectors(const std::vector<int>& topSector, int nl) {
    std::vector<int> activeProps;
    for (int i = 0; i < (int)topSector.size(); ++i)
        if (topSector[i] == 1) activeProps.push_back(i);

    int nActive = (int)activeProps.size();
    int ne = (int)topSector.size();
    std::vector<std::vector<int>> result;
    // Generate all combinations of sizes nl..ne
    int minSize = std::max(1, nl); // exclude empty set
    int maxSize = ne;
    std::function<void(int, int, std::vector<int>&)> dfs = [&](int start, int depth, std::vector<int>& chosen) {
        if ((int)chosen.size() >= minSize && (int)chosen.size() <= maxSize) {
            std::vector<int> sector(ne, 0);
            for (int idx : chosen) sector[idx] = 1;
            result.push_back(sector);
        }
        if ((int)chosen.size() >= maxSize) return;
        for (int i = start; i < nActive; ++i) {
            chosen.push_back(activeProps[i]);
            dfs(i + 1, depth + 1, chosen);
            chosen.pop_back();
        }
    };
    std::vector<int> chosen;
    dfs(0, 0, chosen);
    return result;
}

// ============================================================
// Top-level: generate IBP equations for a family
//
// Phase 2a: builds Singular script, but equations loaded from stub.
// Phase 2b: completes Singular integration for real computation.
// ============================================================

inline IBPEquations generateIBPEquations(const FamilyDef& fam) {
    IBPEquations eqs;
    eqs.ne   = fam.nProp();
    eqs.nl   = fam.nLoop();
    eqs.nE   = fam.nExt();
    eqs.nibp = eqs.nl * (eqs.nl + eqs.nE);

    for (int i = 1; i <= eqs.ne; ++i) {
        eqs.Alist.push_back("\"A\"[" + std::to_string(i) + "]");
        eqs.vlist.push_back("\"v\"[" + std::to_string(i) + "]");
    }

    eqs.sectorlist = generateSubsectors(fam.topSector, eqs.nl);

    std::cout << "[IBPEqGenerator] Family " << fam.name
              << ": ne=" << eqs.ne << " nl=" << eqs.nl << " nE=" << eqs.nE
              << " nIBP=" << eqs.nibp
              << " subsectors=" << eqs.sectorlist.size() << std::endl;

    // Step 1: SP2PD — express scalar products in z-basis
    std::cout << "[IBPEqGenerator] Step 1: Running SP2PD..." << std::endl;
    SP2PDResult sp2pd = runSP2PD(fam);
    if (sp2pd.success) {
        std::cout << "[IBPEqGenerator] SP2PD OK: " << sp2pd.nSP
                  << " scalar products → z-basis" << std::endl;
    } else {
        std::cout << "[IBPEqGenerator] SP2PD FAILED" << std::endl;
        return eqs;
    }

    // Step 2: Derivatives — compute deriv[i,j,k] in z-basis
    std::cout << "[IBPEqGenerator] Step 2: Computing derivatives..." << std::endl;
    DerivResult deriv = runIBPDerivatives(fam);
    if (!deriv.success) {
        std::cout << "[IBPEqGenerator] Derivative computation FAILED" << std::endl;
        return eqs;
    }
    std::cout << "[IBPEqGenerator] Derivatives OK: " << deriv.nibp
              << " polynomials computed" << std::endl;

    // Parse d-value from family config (rational string → integer mod p)
    int64_t dCoeff = 1;
    auto it = fam.numeric.find("d");
    if (it != fam.numeric.end()) {
        const std::string& ds = it->second;
        size_t slash = ds.find('/');
        if (slash != std::string::npos) {
            int64_t num = std::stoll(ds.substr(0, slash));
            int64_t den = std::stoll(ds.substr(slash + 1));
            // Modular inverse via Fermat's little theorem (modulus is prime)
            int64_t inv = 1, base = ((den % fam.modulus) + fam.modulus) % fam.modulus;
            int64_t exp = fam.modulus - 2;
            while (exp > 0) {
                if (exp & 1) inv = (inv * base) % fam.modulus;
                base = (base * base) % fam.modulus;
                exp >>= 1;
            }
            dCoeff = ((num % fam.modulus) + fam.modulus) % fam.modulus;
            dCoeff = (dCoeff * inv) % fam.modulus;
        } else {
            dCoeff = std::stoll(ds);
            dCoeff = (dCoeff % fam.modulus + fam.modulus) % fam.modulus;
        }
    }

    // DEBUG: print derivatives for small families
    if (eqs.ne <= 4) {
        std::cout << "\n[DEBUG DERIV] Raw derivatives:" << std::endl;
        for (size_t di=0; di<deriv.derivPolynomials.size(); ++di)
            std::cout << "  deriv[" << di << "] = " << deriv.derivPolynomials[di] << std::endl;
    }
    // Step 3: Assemble IBP equations from derivatives
    std::cout << "[IBPEqGenerator] Step 3: Assembling IBP equations (mon2F)..." << std::endl;
    auto ibpEqsG = assembleIBPFromDerivatives(
        deriv, eqs.ne, eqs.nl, eqs.nE, fam.modulus, dCoeff);

    for (auto& geq : ibpEqsG) {
        IBPEquation eq;
        eq.j = geq.j;
        eq.k = geq.k;
        for (auto& gt : geq.gTerms) {
            IBPTerm term;
            term.gShift = gt.gShift;
            term.coeff = gt.coeff;
            term.nIdx = gt.nIdx;
            term.hasD = gt.hasD;
            eq.terms.push_back(term);
        }
        eqs.equations.push_back(eq);
    }

    // Summary
    int totalTerms = 0;
    for (auto& eq : eqs.equations) totalTerms += eq.terms.size();
    std::cout << "[IBPEqGenerator] Generated " << eqs.equations.size()
              << " IBP equations with " << totalTerms << " total g-operator terms."
              << std::endl;

    return eqs;
}

} // namespace IBPEqGenerator
#endif
