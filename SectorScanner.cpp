// SectorScanner — scan subsectors of an IBP family and determine which are ZERO
//
// Usage:
//   ./SectorScanner <family.json> <sector_mask> [--verbose]
//
//   <family.json>   Path to a family JSON file (e.g., families/bub11.json)
//   <sector_mask>   Bitstring filter, e.g., "11001" means propagators D0,D1,D4 active
//                   Only subsectors of this mask are tested
//   --verbose       Print G, linear system, etc. for each subsector
//
// Example:
//   ./SectorScanner families/bub11.json 11
//   ./SectorScanner families/Tri.json 111 --verbose
//   ./SectorScanner families/SR212.json 11111
//
// For each subsector with rank >= L (number of loops), calls ZeroSectorQ
// and prints whether the sector is ZERO or NON-ZERO.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "json.hpp"
#include "ZeroSectorQ.hpp"

using json = nlohmann::json;
using Rat = ZeroSectorQ::Rat;

// ============================================================
// Expression evaluator — simple arithmetic over Q
// Supports: +, -, *, /, (), integers, fractions (a/b), variables
// ============================================================
struct ExprEval {
    std::map<std::string, Rat> vars;

    void setVar(const std::string &name, Rat val) { vars[name] = val; }

    // Tokenize
    struct Token {
        enum Type { NUM, VAR, PLUS, MINUS, MUL, DIV, LPAREN, RPAREN, END };
        Type type;
        Rat numVal;
        std::string varName;
    };

    std::vector<Token> tokens;
    size_t tokPos;

    Token nextToken(const std::string &s, size_t &pos) {
        while (pos < s.size() && s[pos] == ' ') ++pos;
        if (pos >= s.size()) return {Token::END};
        char c = s[pos];
        if (c == '+') { ++pos; return {Token::PLUS}; }
        if (c == '-') { ++pos; return {Token::MINUS}; }
        if (c == '*') { ++pos; return {Token::MUL}; }
        if (c == '/') { ++pos; return {Token::DIV}; }
        if (c == '(') { ++pos; return {Token::LPAREN}; }
        if (c == ')') { ++pos; return {Token::RPAREN}; }
        if (std::isdigit(c)) {
            size_t start = pos;
            bool isFrac = false, isDecimal = false;
            while (pos < s.size() && (std::isdigit(s[pos]) || s[pos] == '/' || s[pos] == '.')) {
                if (s[pos] == '/') isFrac = true;
                if (s[pos] == '.') isDecimal = true;
                ++pos;
            }
            std::string numStr = s.substr(start, pos - start);
            Token t;
            t.type = Token::NUM;
            if (isFrac) {
                auto slash = numStr.find('/');
                int64_t n = std::stoll(numStr.substr(0, slash));
                int64_t d = std::stoll(numStr.substr(slash + 1));
                t.numVal = Rat(n, d);
            } else if (isDecimal) {
                double val = std::stod(numStr);
                int64_t num = (int64_t)std::round(val * 1e9);
                int64_t den = (int64_t)1e9;
                int64_t g = std::gcd(std::abs(num), std::abs(den));
                num /= g; den /= g;
                t.numVal = Rat(num, den);
            } else {
                t.numVal = Rat(std::stoll(numStr));
            }
            return t;
        }
        if (std::isalpha(c) || c == '_') {
            size_t start = pos;
            while (pos < s.size() && (std::isalnum(s[pos]) || s[pos] == '_')) ++pos;
            std::string name = s.substr(start, pos - start);
            // Check if it's a known variable or looks like a variable name
            Token t;
            t.type = Token::VAR;
            t.varName = name;
            return t;
        }
        throw std::runtime_error("Unexpected char '" + std::string(1, c) + "' in expression: " + s);
    }

    void tokenize(const std::string &s) {
        tokens.clear();
        size_t pos = 0;
        while (true) {
            Token t = nextToken(s, pos);
            tokens.push_back(t);
            if (t.type == Token::END) break;
        }
        tokPos = 0;
    }

    Token peek() { return tokens[tokPos]; }
    Token consume() { return tokens[tokPos++]; }

    Rat parseExpr(const std::string &s) {
        tokenize(s);
        tokPos = 0;
        Rat r = parseAddSub();
        return r;
    }

    Rat parseAddSub() {
        Rat r = parseMulDiv();
        while (peek().type == Token::PLUS || peek().type == Token::MINUS) {
            if (consume().type == Token::PLUS) r = r + parseMulDiv();
            else r = r - parseMulDiv();
        }
        return r;
    }

    Rat parseMulDiv() {
        Rat r = parseUnary();
        while (true) {
            auto t = peek();
            if (t.type == Token::MUL) { consume(); r = r * parseUnary(); }
            else if (t.type == Token::DIV) { consume(); r = r / parseUnary(); }
            else if (t.type == Token::VAR || t.type == Token::LPAREN || t.type == Token::NUM) {
                // implicit multiplication: 2x, (a)(b), etc.
                r = r * parseUnary();
            } else break;
        }
        return r;
    }

    Rat parseUnary() {
        if (peek().type == Token::MINUS) {
            consume();
            return -parseAtom();
        }
        if (peek().type == Token::PLUS) {
            consume();
            return parseAtom();
        }
        return parseAtom();
    }

    Rat parseAtom() {
        if (peek().type == Token::NUM) {
            return consume().numVal;
        }
        if (peek().type == Token::VAR) {
            std::string name = consume().varName;
            auto it = vars.find(name);
            if (it == vars.end())
                throw std::runtime_error("Unknown variable: " + name);
            return it->second;
        }
        if (peek().type == Token::LPAREN) {
            consume();
            Rat r = parseAddSub();
            if (peek().type != Token::RPAREN)
                throw std::runtime_error("Missing )");
            consume();
            return r;
        }
        throw std::runtime_error("Unexpected token in expression");
    }
};

// ============================================================
// Propagator parser
//   Input:  "-k1^2 + msq",  "-(k1-p1)^2 + msq",  "(l1+p1)^2",  "-l1^2"
//   Output: momentumExpression, massExpression (string)
// ============================================================
struct PropParts {
    std::string momentum; // e.g., "k1-p1"
    std::string massStr;  // e.g., "msq" or ""
    int massSign = 1;      // sign of mass term
};

static PropParts parsePropagator(const std::string &prop) {
    std::string s = prop;

    int massSign = 1;
    std::string massStr;
    int parenLevel = 0;
    size_t splitPos = std::string::npos;
    char splitOp = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '(') ++parenLevel;
        if (s[i] == ')') --parenLevel;
        if (parenLevel == 0 && i + 3 <= s.size()) {
            if ((s[i] == ' ' && s[i+1] == '+' && s[i+2] == ' ') ||
                (s[i] == ' ' && s[i+1] == '-' && s[i+2] == ' ')) {
                if (i >= 2 && s[i-1] == '2' && s[i-2] == '^') continue;
                splitPos = i;
                splitOp = s[i+1];
                break;
            }
        }
    }
    if (splitPos != std::string::npos) {
        massStr = s.substr(splitPos + 3);
        s = s.substr(0, splitPos);
        massSign = (splitOp == '+') ? 1 : -1;
    }

    if (!s.empty() && s[0] == '-')
        s = s.substr(1);

    size_t hatPos = s.rfind("^2");
    if (hatPos != std::string::npos)
        s = s.substr(0, hatPos);

    if (s.size() >= 2 && s[0] == '(' && s.back() == ')')
        s = s.substr(1, s.size() - 2);

    return {s, massStr, massSign};
}

// ============================================================
// Parse momentum expression into coefficients of loop/ext momenta
//   momentum: "k1-p1", "l1+l2+p", "k1-p1-p2", "l1"
//   loopNames: {"k1"}  or  {"l1","l2"}
//   extNames:  {"p1"}  or  {"p1","p2"}
//
//   term: "+- name" where name is a loop or ext momentum
//   Returns: coefficients in L0_row, E0_row (same for L1_row, E1_row)
// ============================================================
struct MomCoeffs {
    std::vector<int> loopCoeff; // per loop momentum
    std::vector<int> extCoeff;  // per external momentum
};

static MomCoeffs parseMomentum(const std::string &momentum,
                                const std::vector<std::string> &loopNames,
                                const std::vector<std::string> &extNames) {
    MomCoeffs mc;
    mc.loopCoeff.assign(loopNames.size(), 0);
    mc.extCoeff.assign(extNames.size(), 0);

    if (momentum.empty()) return mc;

    // Build name→index maps
    std::map<std::string, int> loopIdx, extIdx;
    for (size_t i = 0; i < loopNames.size(); ++i) loopIdx[loopNames[i]] = (int)i;
    for (size_t i = 0; i < extNames.size(); ++i)  extIdx[extNames[i]] = (int)i;

    // Parse as +/- separated terms
    std::string s = momentum;
    // Normalize: replace spaces
    for (auto &c : s) if (c == ' ') c = ' '; // already fine
    
    // Split by +/- operators (keeping the sign)
    // Strategy: iterate character by character
    // A term starts with + or - (or beginning of string), followed by name chars
    
    // First, strip leading spaces
    size_t pos = 0;
    int sign = 1;
    
    while (pos < s.size()) {
        // Skip spaces
        while (pos < s.size() && s[pos] == ' ') ++pos;
        if (pos >= s.size()) break;
        
        // Read sign
        if (s[pos] == '+') { sign = 1; ++pos; }
        else if (s[pos] == '-') { sign = -1; ++pos; }
        // else: first term, sign stays as default (+1)
        
        // Skip spaces after sign
        while (pos < s.size() && s[pos] == ' ') ++pos;
        
        // Read term name (alphanumeric + underscore)
        size_t start = pos;
        while (pos < s.size() && (std::isalnum(s[pos]) || s[pos] == '_')) ++pos;
        std::string name = s.substr(start, pos - start);
        
        if (name.empty()) continue;
        
        // Look up in loop and ext maps
        auto lit = loopIdx.find(name);
        if (lit != loopIdx.end()) {
            mc.loopCoeff[lit->second] += sign;
        } else {
            auto eit = extIdx.find(name);
            if (eit != extIdx.end()) {
                mc.extCoeff[eit->second] += sign;
            } else {
                // Unknown name — could be a parameter, ignore
                // But warn about it? For now just skip.
            }
        }
    }
    
    return mc;
}

// ============================================================
// Kinematic rule parser — build dot product matrix
//   kinRules: {p1^2 -> s, p2^2 -> s2, p1*p2 -> (s3-s1-s2)/2}
//   numeric:  {s -> 3, s2 -> 8, s3 -> 0, ...}
//   extNames: {"p1", "p2"}
//
//   Returns: vector<Rat> of size E*E, kr[e1*E+e2] = pExt[e1]·pExt[e2]
// ============================================================
static std::vector<Rat> buildKinMatrix(
    const std::vector<std::string> &extNames,
    const std::map<std::string, std::string> &kinRules,
    const std::map<std::string, std::string> &numeric)
{
    int E = (int)extNames.size();
    std::vector<Rat> kr(E * E, Rat(0));

    if (E == 0) return kr;

    // Set up expression evaluator with numeric values
    ExprEval eval;
    for (auto &[k, v] : numeric) eval.setVar(k, eval.parseExpr(v));

    // Build name→index map
    std::map<std::string, int> extIdx;
    for (int i = 0; i < E; ++i) extIdx[extNames[i]] = i;

    auto normalize = [](const std::string &s) -> std::string {
        std::string r;
        for (char c : s) if (c != ' ') r += c;
        return r;
    };

    for (auto &[rule, expr] : kinRules) {
        std::string r = normalize(rule);
        std::string e = normalize(expr);

        // Pattern 1: "p1^2" -> p1·p1
        // Extract momentum name before ^2
        size_t hatPos = r.find("^2");
        if (hatPos != std::string::npos && hatPos > 0) {
            std::string pName = r.substr(0, hatPos);
            auto it = extIdx.find(pName);
            if (it != extIdx.end()) {
                Rat val = eval.parseExpr(e);
                kr[it->second * E + it->second] = val;
                continue;
            }
        }

        // Pattern 2: "(p1+p2)^2" -> (p1+p2)·(p1+p2) = p1·p1 + p2·p2 + 2*p1·p2
        if (r.size() > 4 && r[0] == '(') {
            size_t closeParen = r.find(')');
            if (closeParen != std::string::npos && closeParen + 2 < r.size() && r.substr(closeParen, 3) == ")^2") {
                std::string inner = r.substr(1, closeParen - 1);
                // Parse inner as sum of momenta
                // e.g., "p1+p2" means we need the sum's square
                // Compute as: Σ_i p_i^2 + 2*Σ_{i<j} p_i·p_j
                // We know the RHS value, so this gives us dot products
                // But we need individual p_i·p_j. If there are only 2 momenta:
                // (p1+p2)^2 = p1^2 + p2^2 + 2*p1·p2
                // So p1·p2 = ((p1+p2)^2 - p1^2 - p2^2) / 2
                // This is handled by explicit "p1*p2" rules below
                
                // If we already have all the individual p_i^2, we can solve for cross terms
                // But since individual dot products might be given explicitly, skip here
                // unless no explicit dot product rule exists
                
                // For now, store the sum-square value and handle cross terms
                Rat sumSq = eval.parseExpr(e);
                
                // Parse the inner sum to find momenta
                std::vector<std::string> terms;
                std::string cur;
                for (size_t i = 0; i < inner.size(); ++i) {
                    if (inner[i] == '+') {
                        if (!cur.empty()) { terms.push_back(cur); cur.clear(); }
                    } else if (inner[i] != ' ') {
                        cur += inner[i];
                    }
                }
                if (!cur.empty()) terms.push_back(cur);
                
                // If we have all individual squared, and all pairwise dot products,
                // Then (sum)^2 is just a consistency check. Store it only if we need it.
                // For now, compute dot products only if specific p_i·p_j rules are given.
                
                // Special case: known dot products give cross terms directly
                // If only 2 momenta: p1·p2 = (sumSq - p1^2 - p2^2) / 2
                if (terms.size() == 2) {
                    auto it1 = extIdx.find(terms[0]);
                    auto it2 = extIdx.find(terms[1]);
                    if (it1 != extIdx.end() && it2 != extIdx.end()) {
                        // Check if already set via explicit p1*p2 rule
                        int idx1 = it1->second, idx2 = it2->second;
                        if (idx1 != idx2 && kr[idx1 * E + idx2].isZero() && kr[idx2 * E + idx1].isZero()) {
                            Rat p1sq = kr[idx1 * E + idx1];
                            Rat p2sq = kr[idx2 * E + idx2];
                            // If p1^2 and p2^2 are set
                            if (!p1sq.isZero() || !p2sq.isZero()) {
                                Rat cross = (sumSq - p1sq - p2sq) * Rat(1, 2);
                                kr[idx1 * E + idx2] = cross;
                                kr[idx2 * E + idx1] = cross;
                            }
                        }
                    }
                }
                continue;
            }
        }

        // Pattern 3: "p1*p2" -> p1·p2
        size_t starPos = r.find('*');
        if (starPos != std::string::npos) {
            std::string p1 = r.substr(0, starPos);
            std::string p2 = r.substr(starPos + 1);
            auto it1 = extIdx.find(p1);
            auto it2 = extIdx.find(p2);
            if (it1 != extIdx.end() && it2 != extIdx.end()) {
                Rat val = eval.parseExpr(e);
                kr[it1->second * E + it2->second] = val;
                kr[it2->second * E + it1->second] = val;
                continue;
            }
        }
    }

    return kr;
}

// ============================================================
// Evaluate mass string using numeric map
// ============================================================
static Rat evalMass(const std::string &massStr,
                     const std::map<std::string, std::string> &numeric) {
    if (massStr.empty()) return Rat(0);
    
    ExprEval eval;
    for (auto &[k, v] : numeric) eval.setVar(k, eval.parseExpr(v));
    
    return eval.parseExpr(massStr);
}

// ============================================================
// Build coefficient matrices for ALL propagators in the family
// ============================================================
struct FamilyMatrices {
    int L, E, N;
    std::vector<std::vector<int>> L0; // [L][N]
    std::vector<std::vector<int>> L1; // [L][N]
    std::vector<std::vector<int>> E0; // [E][N]
    std::vector<std::vector<int>> E1; // [E][N]
    std::vector<Rat> masses;          // [N]
};

static FamilyMatrices buildFamilyMatrices(
    const std::vector<std::string> &propagators,
    const std::vector<std::string> &loopMomenta,
    const std::vector<std::string> &externalMomenta,
    const std::map<std::string, std::string> &numeric)
{
    int L = (int)loopMomenta.size();
    int E = (int)externalMomenta.size();
    int N = (int)propagators.size();

    FamilyMatrices fm;
    fm.L = L; fm.E = E; fm.N = N;
    fm.L0.assign(L, std::vector<int>(N, 0));
    fm.L1.assign(L, std::vector<int>(N, 0));
    fm.E0.assign(E, std::vector<int>(N, 0));
    fm.E1.assign(E, std::vector<int>(N, 0));
    fm.masses.resize(N, Rat(0));

    for (int p = 0; p < N; ++p) {
        auto parts = parsePropagator(propagators[p]);
        
        // Parse momentum → coefficients (same for mom1 and mom2 since propagator is squared)
        auto mc = parseMomentum(parts.momentum, loopMomenta, externalMomenta);
        
        for (int l = 0; l < L; ++l) {
            fm.L0[l][p] = mc.loopCoeff[l];
            fm.L1[l][p] = mc.loopCoeff[l];
        }
        for (int e = 0; e < E; ++e) {
            fm.E0[e][p] = mc.extCoeff[e];
            fm.E1[e][p] = mc.extCoeff[e];
        }
        
        Rat m = evalMass(parts.massStr, numeric);
        if (parts.massSign < 0) m = -m;
        fm.masses[p] = m;
    }

    return fm;
}

// ============================================================
// Sector utilities
// ============================================================

// Count 1-bits in an integer
static int popcount(uint64_t x) {
    return (int)__builtin_popcountll(x);
}

// Bit → index list
static std::vector<int> bitsToIndices(uint64_t mask, int N) {
    std::vector<int> idx;
    for (int i = 0; i < N; ++i)
        if (mask & (1ULL << i)) idx.push_back(i);
    return idx;
}

// Extract submatrices for a given sector (subset of propagators)
static bool testSector(
    uint64_t sectorMask,
    const FamilyMatrices &fm,
    const std::vector<Rat> &kinRules,
    bool verbose = false)
{
    auto idx = bitsToIndices(sectorMask, fm.N);
    int n = (int)idx.size();

    if (n == 0) {
        if (verbose) std::cout << "  (empty sector) → ZERO\n";
        return true;
    }

    std::vector<std::vector<int>> L0(fm.L, std::vector<int>(n, 0));
    std::vector<std::vector<int>> L1(fm.L, std::vector<int>(n, 0));
    std::vector<std::vector<int>> E0(fm.E, std::vector<int>(n, 0));
    std::vector<std::vector<int>> E1(fm.E, std::vector<int>(n, 0));
    std::vector<Rat> masses(n);

    for (int j = 0; j < n; ++j) {
        int p = idx[j];
        for (int l = 0; l < fm.L; ++l) {
            L0[l][j] = fm.L0[l][p];
            L1[l][j] = fm.L1[l][p];
        }
        for (int e = 0; e < fm.E; ++e) {
            E0[e][j] = fm.E0[e][p];
            E1[e][j] = fm.E1[e][p];
        }
        masses[j] = fm.masses[p];
    }

    return ZeroSectorQ::zeroSectorQ(fm.L, fm.E, n, L0, L1, E0, E1, masses, kinRules, verbose);
}

// ============================================================
// Format sector mask as binary string
// ============================================================
static std::string maskStr(uint64_t mask, int N) {
    std::string s;
    for (int i = 0; i < N; ++i)
        s += (mask & (1ULL << i)) ? '1' : '0';
    return s;
}

// ============================================================
// Main
// ============================================================
int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <family.json> <sector_mask> [--verbose]\n"
                  << "  sector_mask: binary string, e.g., \"11101\"\n"
                  << "  Scans all subsectors of the given mask with rank >= L\n";
        return 1;
    }

    std::string famPath = argv[1];
    std::string maskStrInput = argv[2];
    bool verbose = (argc > 3 && std::string(argv[3]) == "--verbose");

    // Parse the filter mask
    uint64_t filterMask = 0;
    int N_input = (int)maskStrInput.size();
    for (int i = 0; i < N_input; ++i) {
        if (maskStrInput[i] == '1')
            filterMask |= (1ULL << i); // bit i = propagator i = D_{i+1}
    }

    // Read JSON
    std::ifstream f(famPath);
    if (!f) { std::cerr << "Cannot open " << famPath << "\n"; return 1; }
    json j;
    f >> j;

    std::string name = j.value("name", "unnamed");
    std::vector<std::string> propagators = j["propagators"].get<std::vector<std::string>>();
    std::vector<std::string> loopMomenta = j["loopMomenta"].get<std::vector<std::string>>();
    std::vector<std::string> extMomenta = j["externalMomenta"].get<std::vector<std::string>>();

    int N_total = (int)propagators.size();
    int L = (int)loopMomenta.size();
    int E = (int)extMomenta.size();

    // Build numeric map
    std::map<std::string, std::string> numeric;
    if (j.contains("numeric")) {
        for (auto &[k, v] : j["numeric"].items()) {
            if (v.is_number_integer()) numeric[k] = std::to_string(v.get<int64_t>());
            else if (v.is_number_float()) numeric[k] = std::to_string(v.get<double>());
            else numeric[k] = v.get<std::string>();
        }
    }

    // Build kinematic rules map
    std::map<std::string, std::string> kinRules;
    if (j.contains("kinematicRules")) {
        for (auto &[k, v] : j["kinematicRules"].items())
            kinRules[k] = v.get<std::string>();
    }

    // Decode top sector from JSON (if present)
    uint64_t topMask = filterMask; // default = user-provided mask
    if (j.contains("topSector")) {
        auto ts = j["topSector"].get<std::vector<int>>();
        topMask = 0;
        for (int i = 0; i < (int)ts.size() && i < N_total; ++i)
            if (ts[i]) topMask |= (1ULL << i);
    }
    // Intersect with filter
    topMask &= filterMask;

    std::cout << "Family: " << name << "\n";
    std::cout << "  Loops: " << L << ", External: " << E << ", Props: " << N_total << "\n";
    std::cout << "  Filter mask: " << maskStrInput << "\n";
    std::cout << "  Effective top sector: " << maskStr(topMask, N_total) << "\n\n";

    // Build full coefficient matrices
    FamilyMatrices fm = buildFamilyMatrices(propagators, loopMomenta, extMomenta, numeric);

    // Build kinematic matrix
    std::vector<Rat> kinMatrix = buildKinMatrix(extMomenta, kinRules, numeric);

    // Print kinematic matrix if verbose
    if (verbose && E > 0) {
        std::cout << "External momentum dot matrix:\n";
        for (int e1 = 0; e1 < E; ++e1) {
            for (int e2 = 0; e2 < E; ++e2)
                std::cout << "  " << std::setw(6) << kinMatrix[e1*E+e2].str();
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    // Print propagator info
    if (verbose) {
        std::cout << "Propagator coefficients:\n";
        for (int p = 0; p < N_total; ++p) {
            std::cout << "  D" << (p+1) << ": loop=[";
            for (int l = 0; l < L; ++l) std::cout << fm.L0[l][p] << (l+1<L?",":"");
            std::cout << "] ext=[";
            for (int e = 0; e < E; ++e) std::cout << fm.E0[e][p] << (e+1<E?",":"");
            std::cout << "] mass=" << fm.masses[p].str() << "\n";
        }
        std::cout << "\n";
    }

    // Enumerate all subsectors of filterMask
    std::cout << "Sector scan (rank >= " << L << "):\n";
    std::cout << "  mask        rank  result\n";
    std::cout << "  " << std::string(N_total + 8, '-') << "\n";

    // Iterate over all submasks of filterMask
    uint64_t sub = topMask;
    int total = 0, zero = 0, nonzero = 0;
    do {
        int rank = popcount(sub);
        if (rank < L) {
            if (sub == 0) break;
            sub = (sub - 1) & topMask;
            continue;
        }

        bool isZero = testSector(sub, fm, kinMatrix, verbose);

        std::cout << "  " << maskStr(sub, N_total)
                  << "  " << std::setw(3) << rank
                  << "    " << (isZero ? "ZERO" : "NON-ZERO") << "\n";

        total++;
        if (isZero) zero++; else nonzero++;

        if (sub == 0) break;
        sub = (sub - 1) & topMask;
    } while (true);

    std::cout << "\nSummary: " << total << " sectors tested, "
              << zero << " ZERO, " << nonzero << " NON-ZERO\n";

    return 0;
}
