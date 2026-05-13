#ifndef POLY_ARITH_HPP
#define POLY_ARITH_HPP

// Polynomial arithmetic modulo a prime p.
// Used by IBPAnalyzer, RecursionBuilder, and RingBuilder for
// symbolic algebra in the quotient ring F_p[VarIndep] / (MinPoly).

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <cstdint>
#include <cctype>
#include <algorithm>
#include <stdexcept>

namespace PolyArith {

// A monomial: coeff * ∏_i var_i^{exps[i]}, with coefficient mod p.
struct Monomial {
    std::vector<int> exps;
    int64_t coeff = 0;
};

using Polynomial = std::vector<Monomial>;

// ============================================================
// Utilities
// ============================================================

// Lexicographic compare of two exponent vectors (same length).
inline int cmpExps(const std::vector<int>& a, const std::vector<int>& b) {
    for (size_t i = 0; i < a.size() && i < b.size(); ++i) {
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    }
    if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
    return 0;
}

inline bool expsEqual(const std::vector<int>& a, const std::vector<int>& b) {
    return a == b; // vector::operator== does element-wise comparison
}

// Sort a polynomial by exponent vector (lex order) and merge like terms.
inline void canonicalize(Polynomial& p, int64_t mod) {
    if (p.empty()) return;
    // Sort by exponent vector
    std::sort(p.begin(), p.end(), [](const Monomial& a, const Monomial& b) {
        return cmpExps(a.exps, b.exps) < 0;
    });
    // Merge like terms
    Polynomial merged;
    for (size_t i = 0; i < p.size(); ++i) {
        if (p[i].coeff == 0) continue;
        if (!merged.empty() && expsEqual(merged.back().exps, p[i].exps)) {
            merged.back().coeff = (merged.back().coeff + p[i].coeff) % mod;
            if (merged.back().coeff < 0) merged.back().coeff += mod;
            if (merged.back().coeff == 0) merged.pop_back();
        } else {
            int64_t c = p[i].coeff % mod;
            if (c < 0) c += mod;
            if (c != 0) merged.push_back({p[i].exps, c});
        }
    }
    p = std::move(merged);
}

// ============================================================
// Arithmetic operations
// ============================================================

// Add two polynomials
inline Polynomial polyAdd(const Polynomial& a, const Polynomial& b, int64_t mod) {
    Polynomial result;
    result.reserve(a.size() + b.size());
    result.insert(result.end(), a.begin(), a.end());
    result.insert(result.end(), b.begin(), b.end());
    canonicalize(result, mod);
    return result;
}

// Multiply two polynomials
inline Polynomial polyMul(const Polynomial& a, const Polynomial& b, int64_t mod) {
    if (a.empty() || b.empty()) return {};
    int nVars = static_cast<int>(a[0].exps.size());
    Polynomial result;
    result.reserve(a.size() * b.size());
    for (const auto& ma : a) {
        for (const auto& mb : b) {
            Monomial m;
            m.exps.resize(nVars);
            for (int i = 0; i < nVars; ++i)
                m.exps[i] = ma.exps[i] + mb.exps[i];
            m.coeff = (ma.coeff * mb.coeff) % mod;
            if (m.coeff < 0) m.coeff += mod;
            if (m.coeff != 0) result.push_back(m);
        }
    }
    canonicalize(result, mod);
    return result;
}

// Multiply polynomial by a scalar
inline Polynomial polyScale(const Polynomial& a, int64_t scalar, int64_t mod) {
    scalar = scalar % mod;
    if (scalar < 0) scalar += mod;
    if (scalar == 0) return {};
    Polynomial result;
    result.reserve(a.size());
    for (const auto& m : a) {
        int64_t c = (m.coeff * scalar) % mod;
        if (c < 0) c += mod;
        if (c != 0) result.push_back({m.exps, c});
    }
    return result;
}

// Negate a polynomial
inline Polynomial polyNeg(const Polynomial& a, int64_t mod) {
    return polyScale(a, -1, mod);
}

// ============================================================
// monomialRulesPower
// ============================================================

// Extract a map from exponent vector to coefficient.
// This is the C++ equivalent of MMA's monomialRulesPower:
//   expr → Association[exponent_vector → coefficient]
inline std::map<std::vector<int>, int64_t> monomialRulesPower(const Polynomial& p) {
    std::map<std::vector<int>, int64_t> result;
    for (const auto& m : p) {
        if (m.coeff != 0) {
            result[m.exps] = (result[m.exps] + m.coeff);
        }
    }
    // Mod reduction not needed here — caller handles it
    return result;
}

// ============================================================
// Substitution
// ============================================================

// Substitute varIdx → replacement in polynomial p.
// All polynomials use the same variable ordering.
inline Polynomial polySubstitute(const Polynomial& p, int varIdx,
                                  const Polynomial& replacement, int64_t mod) {
    if (p.empty()) return {};
    int nVars = static_cast<int>(p[0].exps.size());

    // Start with constant 1
    Polynomial result;
    result.push_back({std::vector<int>(nVars, 0), int64_t(1)});

    Polynomial total;
    for (const auto& mon : p) {
        // Build the monomial: coeff * ∏_i (var_i)^{mon.exps[i]}
        // For i == varIdx: use replacement^{mon.exps[i]}
        // For i != varIdx: use the monomial with exps only at i
        Polynomial termProd;
        termProd.push_back({std::vector<int>(nVars, 0), mon.coeff});

        for (int i = 0; i < nVars; ++i) {
            if (mon.exps[i] == 0) continue;
            Polynomial factor;
            if (i == varIdx) {
                factor = replacement;
            } else {
                factor.push_back({std::vector<int>(nVars, 0), int64_t(1)});
                factor[0].exps[i] = 1;
            }
            // Raise factor to power mon.exps[i]
            Polynomial powered;
            powered.push_back({std::vector<int>(nVars, 0), int64_t(1)});
            for (int e = 0; e < mon.exps[i]; ++e) {
                powered = polyMul(powered, factor, mod);
            }
            termProd = polyMul(termProd, powered, mod);
        }
        total = polyAdd(total, termProd, mod);
    }
    return total;
}

// Substitute multiple variables at once.
// rules[i] = replacement polynomial for variable i (empty poly = no replacement).
inline Polynomial polySubstituteMulti(const Polynomial& p,
                                       const std::vector<Polynomial>& replacements,
                                       int64_t mod,
                                       int nOutVars = -1) {
    if (p.empty()) return {};
    int nInVars = static_cast<int>(p[0].exps.size());
    if (nOutVars < 0) nOutVars = nInVars;  // default: same variable count
    Polynomial total;

    for (const auto& mon : p) {
        Polynomial term;
        term.push_back({std::vector<int>(nOutVars, 0), mon.coeff});

        for (int i = 0; i < nInVars; ++i) {
            if (mon.exps[i] == 0) continue;
            if (i >= (int)replacements.size()) continue;
            // Raise var_i's replacement to mon.exps[i]
            Polynomial factor = replacements[i];
            if (factor.empty()) {
                if (nOutVars == nInVars) {
                    // Variable not in replacement map → keep as variable
                    factor.push_back({std::vector<int>(nOutVars, 0), int64_t(1)});
                    factor[0].exps[i] = 1;
                } else {
                    // Different variable spaces: missing replacement → term vanishes
                    term.clear();
                    break;
                }
            }
            if (term.empty()) break;
            Polynomial powered;
            powered.push_back({std::vector<int>(nOutVars, 0), int64_t(1)});
            for (int e = 0; e < mon.exps[i]; ++e)
                powered = polyMul(powered, factor, mod);
            term = polyMul(term, powered, mod);
        }
        if (!term.empty())
            total = polyAdd(total, term, mod);
    }
    return total;
}

// ============================================================
// Polynomial reduction modulo a Groebner basis
// ============================================================

// Reduce a single monomial product (exponent vector) modulo a list of
// leading monomials. This is a simple version: we just check if the
// monomial is divisible by any leading monomial.
inline bool isDivisibleByLM(const std::vector<int>& exps,
                             const std::vector<std::vector<int>>& leadingMonomials) {
    for (const auto& lm : leadingMonomials) {
        bool div = true;
        for (size_t i = 0; i < exps.size() && i < lm.size(); ++i) {
            if (exps[i] < lm[i]) { div = false; break; }
        }
        if (div) return true;
    }
    return false;
}

// ============================================================
// Parse / Format utilities
// ============================================================

// Format a polynomial as a Singular-compatible string.
// varNames[i] = variable name for index i.
inline std::string polyToSingularString(const Polynomial& p,
                                         const std::vector<std::string>& varNames) {
    if (p.empty()) return "0";
    std::stringstream ss;
    for (size_t i = 0; i < p.size(); ++i) {
        int64_t c = p[i].coeff;
        bool hasVar = false;
        for (size_t j = 0; j < p[i].exps.size(); ++j)
            if (p[i].exps[j] > 0) hasVar = true;

        // Build the variable part
        std::stringstream varPart;
        bool first = true;
        for (size_t j = 0; j < p[i].exps.size(); ++j) {
            if (p[i].exps[j] == 0) continue;
            if (!first) varPart << "*";
            varPart << varNames[j];
            if (p[i].exps[j] > 1)
                varPart << "^" << p[i].exps[j];
            first = false;
        }
        std::string varStr = varPart.str();

        // Sign handling
        if (i > 0) {
            ss << (c > 0 ? "+" : "");
        }

        // Coefficient + variable
        if (!hasVar) {
            // Constant term
            ss << c;
        } else if (c == 1) {
            ss << varStr;
        } else if (c == -1) {
            ss << "-" << varStr;
        } else {
            ss << c << "*" << varStr;
        }
    }
    return ss.str();
}

// ============================================================
// Parse a Singular polynomial string
//
// Format: terms separated by + or -, each term is:
//   [coeff*]var1[^exp1][*var2[^exp2]...]
// Example: "x1^2+x2-3", "3*A1*A2-2*A1+1", "0"
// ============================================================

inline Polynomial parseSingularPolynomial(const std::string& s,
                                           const std::vector<std::string>& varNames,
                                           int64_t mod) {
    int nVars = static_cast<int>(varNames.size());
if (nVars == 0) {
        // No variable names: the polynomial should be a constant.
        // If it contains variables (shouldn't happen for properly projected data),
        // just return as constant 0 to avoid "unknown variable" crashes.
        std::string trimmed = s;
        trimmed.erase(std::remove_if(trimmed.begin(), trimmed.end(), ::isspace), trimmed.end());
        if (trimmed.empty() || trimmed == "0") return {};
        bool allDigits = true;
        for (char c : trimmed) if (c != '+' && c != '-' && !std::isdigit(c)) { allDigits = false; break; }
        if (allDigits) {
            try {
                int64_t val = std::stoll(trimmed) % mod;
                if (val < 0) val += mod;
                Polynomial p;
                p.push_back({{}, val});
                return p;
            } catch (...) {}
        }
        // Has non-digit characters but no variables defined — treat as zero.
        // This can happen when vargen is empty and the "polynomial" is actually
        // just a constant that wasn't recognized as digit-only.
        return {};
    }
    Polynomial result;

    if (s.empty() || s == "0") return result;

    // Build a lookup from var name to index
    std::map<std::string, int> varIdx;
    for (int i = 0; i < nVars; ++i) varIdx[varNames[i]] = i;

    size_t pos = 0;
    int sign = 1;
    bool firstTerm = true;

    while (pos < s.size()) {
        // Skip whitespace
        while (pos < s.size() && s[pos] == ' ') ++pos;
        if (pos >= s.size()) break;

        // Determine sign
        if (!firstTerm) {
            if (s[pos] == '+') { sign = 1; ++pos; }
            else if (s[pos] == '-') { sign = -1; ++pos; }
        } else {
            if (s[pos] == '-') { sign = -1; ++pos; }
            else if (s[pos] == '+') { sign = 1; ++pos; }
            else sign = 1;
        }
        firstTerm = false;

        // Skip whitespace after sign
        while (pos < s.size() && s[pos] == ' ') ++pos;
        if (pos >= s.size()) break;

        // Parse coefficient and variables
        int64_t coeff = 1;
        std::vector<int> exps(nVars, 0);
        bool hasCoeff = false;
        bool hasVar = false;

        // Check if this term starts with a number (coefficient)
        size_t numStart = pos;
        if (std::isdigit(s[pos])) {
            hasCoeff = true;
            while (pos < s.size() && std::isdigit(s[pos])) ++pos;
            coeff = std::stoll(s.substr(numStart, pos - numStart));
        }

        // Parse variable factors
        while (pos < s.size() && s[pos] != '+' && s[pos] != '-') {
            // Skip whitespace and *
            while (pos < s.size() && (s[pos] == ' ' || s[pos] == '*')) ++pos;
            if (pos >= s.size() || s[pos] == '+' || s[pos] == '-') break;

            // Read variable name (letters followed by optional digits)
            size_t varStart = pos;
            while (pos < s.size() && std::isalpha(s[pos])) ++pos;
            while (pos < s.size() && std::isdigit(s[pos])) ++pos;
            std::string varName = s.substr(varStart, pos - varStart);

            if (varName.empty()) {
                // Shouldn't happen, but skip invalid content
                ++pos;
                continue;
            }

            hasVar = true;
            auto it = varIdx.find(varName);
            if (it == varIdx.end()) {
                std::cerr << "[PARSE ERROR] varName='" << varName << "' not in varNames=[";
                for (size_t vi=0;vi<varNames.size();++vi){if(vi)std::cerr<<",";std::cerr<<varNames[vi];}
                std::cerr << "] (nVars=" << nVars << ") s='" << s << "'" << std::endl;
                throw std::runtime_error("Unknown variable in Singular polynomial: " + varName);
            }

            // Check for exponent: ^N
            int exp = 1;
            if (pos < s.size() && s[pos] == '^') {
                ++pos;
                size_t expStart = pos;
                while (pos < s.size() && std::isdigit(s[pos])) ++pos;
                exp = std::stoi(s.substr(expStart, pos - expStart));
            }
            exps[it->second] += exp;
        }

        if (!hasCoeff && !hasVar) {
            // This shouldn't happen in well-formed input
            continue;
        }

        if (!hasVar) {
            // Constant term
            exps.assign(nVars, 0);
        }

        int64_t c = (coeff * sign) % mod;
        if (c < 0) c += mod;
        if (c != 0) result.push_back({exps, c});
    }

    canonicalize(result, mod);
    return result;
}

// ============================================================
// Evaluate a polynomial at a point (exponent-wise, trivial)
// ============================================================

// Build a monomial in a fresh variable space.
// nVars = total number of variables, activeIndices = which ones get exponent 1,
// coeff = coefficient.
inline Polynomial makeMonomial(const std::vector<int>& exps, int64_t coeff, int64_t mod) {
    int64_t c = coeff % mod;
    if (c < 0) c += mod;
    if (c == 0) return {};
    return {{exps, c}};
}

// Zero polynomial
inline Polynomial zeroPoly() { return {}; }

// Constant polynomial
inline Polynomial constPoly(int64_t c, int nVars, int64_t mod) {
    int64_t cv = c % mod;
    if (cv < 0) cv += mod;
    if (cv == 0) return {};
    return {{std::vector<int>(nVars, 0), cv}};
}

} // namespace PolyArith

#endif
