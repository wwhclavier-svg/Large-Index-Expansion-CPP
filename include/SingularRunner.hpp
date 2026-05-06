#ifndef SINGULAR_RUNNER_HPP
#define SINGULAR_RUNNER_HPP

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <unistd.h>

// ==========================================
// SingularRunner — C++ subprocess bridge to Singular CAS
//
// Mirrors MMA's SingularInterface.wl:SingularRun
// ==========================================

namespace SingularRunner {

// ==========================================
// Low-level: execute a Singular script and collect output files
// ==========================================

// Write script to temp file, run `Singular -q -t < script`, read output files.
// Returns map of output_var -> file_content (raw string).
// On error, returns empty map.
inline std::map<std::string, std::string> runSingular(
    const std::string& script,
    const std::vector<std::string>& outVars,
    const std::string& workDir = "/tmp")
{
    // Unique prefix to avoid collisions
    std::string prefix = workDir + "/cpp_sing_" + std::to_string(getpid());

    // Write the Singular script, appending write() commands for each output
    std::string scriptFile = prefix + ".sing";
    {
        std::ofstream f(scriptFile);
        if (!f) throw std::runtime_error("Cannot write Singular script: " + scriptFile);
        f << script << "\n";
        for (size_t i = 0; i < outVars.size(); ++i) {
            f << "write(\"" << prefix << "_out" << i
              << "\", string(" << outVars[i] << "));\n";
        }
    }

    // Execute Singular
    std::string cmd = "Singular -q -t < " + scriptFile + " 2>&1";
    std::string singularOutput;
    {
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            std::remove(scriptFile.c_str());
            throw std::runtime_error("Failed to execute Singular");
        }
        char buf[4096];
        while (fgets(buf, sizeof(buf), pipe))
            singularOutput += buf;
        int rc = pclose(pipe);
        if (rc != 0) {
            std::remove(scriptFile.c_str());
            throw std::runtime_error(
                "Singular exited with code " + std::to_string(rc) +
                ":\n" + singularOutput);
        }
    }

    // Read output files
    std::map<std::string, std::string> results;
    for (size_t i = 0; i < outVars.size(); ++i) {
        std::string outFile = prefix + "_out" + std::to_string(i);
        std::ifstream f(outFile);
        if (f) {
            std::stringstream ss;
            ss << f.rdbuf();
            results[outVars[i]] = ss.str();
        }
        std::remove(outFile.c_str());
    }

    // Cleanup
    std::remove(scriptFile.c_str());

    return results;
}

// ==========================================
// Template-based invocation (mirrors MMA's FileTemplateApply + SingularRun)
//
// `script` is a template with `placeholder` patterns that get replaced.
// ==========================================

inline std::string replaceAll(std::string s, const std::string& from, const std::string& to) {
    size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
    return s;
}

inline std::map<std::string, std::string> runSingularTemplate(
    const std::string& scriptTemplate,
    const std::map<std::string, std::string>& replacements,
    const std::vector<std::string>& outVars,
    const std::string& workDir = "/tmp")
{
    std::string script = scriptTemplate;
    for (const auto& [k, v] : replacements)
        script = replaceAll(script, "<" + k + ">", v);
    return runSingular(script, outVars, workDir);
}

// ==========================================
// Parsing helpers: Singular string format → C++ vectors
// ==========================================

// Split a Singular list-of-lists string like:
//   x1,x2,x3,x1^2-2,x2^2-3
// with nlist[i] = number of elements in sublist i
// Returns vector of sublist strings.
inline std::vector<std::string> parseSingularList(
    const std::string& raw,
    const std::vector<int>& sublistSizes)
{
    // Strip whitespace
    std::string s;
    for (char c : raw) if (c != ' ' && c != '\n' && c != '\r') s += c;

    // Split by comma
    std::vector<std::string> tokens;
    size_t pos = 0, depth = 0;
    std::string token;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '(' || s[i] == '[') ++depth;
        else if (s[i] == ')' || s[i] == ']') --depth;
        else if (s[i] == ',' && depth == 0) {
            tokens.push_back(token);
            token.clear();
            continue;
        }
        token += s[i];
    }
    if (!token.empty()) tokens.push_back(token);

    // Reassemble sublists
    std::vector<std::string> result;
    size_t idx = 0;
    for (int sz : sublistSizes) {
        std::string sublist;
        for (int j = 0; j < sz && idx < tokens.size(); ++j) {
            if (j > 0) sublist += ",";
            sublist += tokens[idx++];
        }
        result.push_back(sublist);
    }
    return result;
}

// Parse a simple flat list like "x3,x2+1,x1-5" → vector of strings
inline std::vector<std::string> parseIdeal(const std::string& raw) {
    // Handle common Singular output: remove outer brackets if present
    std::string s = raw;
    // Trim whitespace
    while (!s.empty() && (s[0] == ' ' || s[0] == '\n')) s.erase(0, 1);
    while (!s.empty() && (s.back() == ' ' || s.back() == '\n')) s.pop_back();

    std::vector<std::string> result;
    size_t pos = 0;
    std::string token;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == ',') {
            if (!token.empty()) { result.push_back(token); token.clear(); }
        } else {
            token += s[i];
        }
    }
    if (!token.empty()) result.push_back(token);
    return result;
}

// Parse dim(GroebnerBasis) output: integer per component
inline std::vector<int> parseDimensionList(const std::string& raw) {
    std::vector<int> result;
    std::string num;
    for (char c : raw) {
        if (c >= '0' && c <= '9' || c == '-') {
            num += c;
        } else if (!num.empty()) {
            result.push_back(std::stoi(num));
            num.clear();
        }
    }
    if (!num.empty()) result.push_back(std::stoi(num));
    return result;
}

// ==========================================
// High-level Singular algebra wrappers
// ==========================================

// Compute Groebner basis of an ideal over Z/p
// `ideal` is a list of polynomial strings in Singular notation
// `varNames` is the variable ordering (first = most significant for dp)
// Returns Groebner basis as list of polynomial strings
inline std::vector<std::string> groebnerBasis(
    const std::vector<std::string>& ideal,
    const std::vector<std::string>& varNames,
    int64_t modulus)
{
    // Build variable string: x1,...,xN
    std::string vars;
    for (size_t i = 0; i < varNames.size(); ++i) {
        if (i > 0) vars += ",";
        vars += varNames[i];
    }

    // Build ideal string
    std::string idealStr;
    for (size_t i = 0; i < ideal.size(); ++i) {
        if (i > 0) idealStr += ",";
        idealStr += ideal[i];
    }

    std::string script = R"(
ring r = <char>, (<vars>), (dp);
ideal I = <ideal>;
ideal gb = groebner(I);
ideal gbs = simplify(gb, 34);
)";

    auto result = runSingularTemplate(script,
        {{"char", std::to_string(modulus)},
         {"vars", vars},
         {"ideal", idealStr}},
        {"gbs"});

    if (result.count("gbs"))
        return parseIdeal(result["gbs"]);
    throw std::runtime_error("Singular groebnerBasis: no output");
}

// Compute Groebner basis + lift matrix
// Returns {gb, liftMatrix} where liftMatrix is flat: dims = {#gens, #gb}
inline std::pair<std::vector<std::string>, std::vector<int64_t>>
groebnerBasisWithLift(
    const std::vector<std::string>& moduleGen,
    const std::vector<std::string>& varNames,
    int64_t modulus)
{
    std::string vars;
    for (size_t i = 0; i < varNames.size(); ++i) {
        if (i > 0) vars += ",";
        vars += varNames[i];
    }

    std::string modStr;
    for (size_t i = 0; i < moduleGen.size(); ++i) {
        if (i > 0) modStr += ",";
        modStr += moduleGen[i];
    }

    std::string script = R"(
ring r = <char>, (<vars>), (dp);
module m = <module>;
module gb = groebner(m);
module gbs = simplify(gb, 34);
matrix lm = lift(m, gb);
)";

    auto result = runSingularTemplate(script,
        {{"char", std::to_string(modulus)},
         {"vars", vars},
         {"module", modStr}},
        {"gbs", "lm"});

    if (!result.count("gbs"))
        throw std::runtime_error("Singular groebnerBasisWithLift: no output");

    std::vector<std::string> gb = parseIdeal(result["gbs"]);
    // Parse lift matrix — for now return empty (complex parsing)
    return {gb, {}};
}

// Minimal associated primes via primdecGTZ
// Returns {primelist, dimlist} where primelist is vector of components,
// each component is a vector of polynomial strings (the associated prime ideal).
//
// NOTE: We use primdecGTZ rather than minAssGTZ because minAssGTZ is buggy
// in some Singular versions (returns positive-dimensional primes for
// zero-dimensional ideals). primdecGTZ returns {primary, prime} pairs,
// and we extract the prime component from each pair.
inline std::pair<std::vector<std::vector<std::string>>, std::vector<int>>
minimalAssPrimes(
    const std::vector<std::string>& ideal,
    const std::vector<std::string>& varNames,
    int64_t modulus,
    const std::string& monomialOrder = "lp")
{
    std::string vars;
    for (size_t i = 0; i < varNames.size(); ++i) {
        if (i > 0) vars += ",";
        vars += varNames[i];
    }

    std::string idealStr;
    for (size_t i = 0; i < ideal.size(); ++i) {
        if (i > 0) idealStr += ",";
        idealStr += ideal[i];
    }

    std::string script = R"(
LIB "primdec.lib";
ring r = <char>, (<vars>), (<mo>);
ideal I = <ideal>;
option(redSB);
ideal gb = std(I);
list primdec = primdecGTZ(gb);
int n = size(primdec);
list plist, rdim, npoly = list(), list(), list();
for (int i = 1; i <= n; i++)
{
  list comp = primdec[i];
  ideal pideal = comp[2];
  plist = plist + list(string(pideal));
  rdim = rdim + list(string(dim(std(pideal))));
  npoly = npoly + list(string(size(pideal)));
}
)";

    auto result = runSingularTemplate(script,
        {{"char", std::to_string(modulus)},
         {"vars", vars},
         {"mo", monomialOrder},
         {"ideal", idealStr}},
        {"plist", "npoly", "rdim"});

    if (!result.count("plist"))
        throw std::runtime_error("Singular minAssPrime: no output");

    // Parse poly counts per component and dimensions
    std::vector<int> npoly = parseDimensionList(result["npoly"]);
    std::vector<int> dims  = parseDimensionList(result["rdim"]);

    // Parse prime ideals: each component's string is a comma-separated list
    std::vector<std::string> tokens;
    {
        std::string s = result["plist"];
        std::string tok;
        size_t depth = 0;
        for (char c : s) {
            if (c == '(' || c == '[') ++depth;
            else if (c == ')' || c == ']') --depth;
            else if (c == ',' && depth == 0) {
                if (!tok.empty()) tokens.push_back(tok);
                tok.clear();
                continue;
            }
            if (c != ' ' && c != '\n' && c != '\r') tok += c;
        }
        if (!tok.empty()) tokens.push_back(tok);
    }

    std::vector<std::vector<std::string>> primelist;
    size_t tidx = 0;
    for (int sz : npoly) {
        std::vector<std::string> comp;
        for (int j = 0; j < sz && tidx < tokens.size(); ++j)
            comp.push_back(tokens[tidx++]);
        primelist.push_back(comp);
    }

    return {primelist, dims};
}

// Primary decomposition with lexicographic output (mirrors SingularPrimaryDecompLex)
// Computes primdecGTZE in dp ring, then imap to lp ring, then groebner each component.
// Returns {primelist (lex GB), dimlist}
inline std::pair<std::vector<std::vector<std::string>>, std::vector<int>>
primaryDecompLex(
    const std::vector<std::string>& ideal,
    const std::vector<std::string>& varNames,
    int64_t modulus)
{
    std::string vars;
    for (size_t i = 0; i < varNames.size(); ++i) {
        if (i > 0) vars += ",";
        vars += varNames[i];
    }

    std::string idealStr;
    for (size_t i = 0; i < ideal.size(); ++i) {
        if (i > 0) idealStr += ",";
        idealStr += ideal[i];
    }

    std::string script = R"(
LIB "primdec.lib";
ring r = <char>, (<vars>), (dp);
ideal I = <ideal>;
list pr0 = primdecGTZE(I);
ring s = <char>, (<vars>), (lp);
list pr = imap(r, pr0);
list prr,rdim,npoly=list(),list(),list();
for(int i=1; i<=size(pr); i=i+1)
{
  pr[i][2] = groebner(pr[i][2]);
  prr = prr+list(string(groebner(pr[i][2])));
  rdim = rdim+list(dim(pr[i][2]));
  npoly = npoly+list(size(pr[i][2]));
};
)";

    auto result = runSingularTemplate(script,
        {{"char", std::to_string(modulus)},
         {"vars", vars},
         {"ideal", idealStr}},
        {"prr", "npoly", "rdim"});

    if (!result.count("prr"))
        throw std::runtime_error("Singular primaryDecompLex: no output");

    std::vector<int> npoly = parseDimensionList(result["npoly"]);
    std::vector<int> dims  = parseDimensionList(result["rdim"]);

    std::vector<std::string> tokens;
    {
        std::string s = result["prr"];
        std::string tok;
        size_t depth = 0;
        for (char c : s) {
            if (c == '(' || c == '[') ++depth;
            else if (c == ')' || c == ']') --depth;
            else if (c == ',' && depth == 0) {
                if (!tok.empty()) tokens.push_back(tok);
                tok.clear();
                continue;
            }
            if (c != ' ' && c != '\n' && c != '\r') tok += c;
        }
        if (!tok.empty()) tokens.push_back(tok);
    }

    std::vector<std::vector<std::string>> primelist;
    size_t tidx = 0;
    for (int sz : npoly) {
        std::vector<std::string> comp;
        for (int j = 0; j < sz && tidx < tokens.size(); ++j)
            comp.push_back(tokens[tidx++]);
        primelist.push_back(comp);
    }

    return {primelist, dims};
}

// Compute dimension of an ideal
inline int idealDim(
    const std::vector<std::string>& ideal,
    const std::vector<std::string>& varNames,
    int64_t modulus)
{
    std::string vars;
    for (size_t i = 0; i < varNames.size(); ++i) {
        if (i > 0) vars += ",";
        vars += varNames[i];
    }

    std::string idealStr;
    for (size_t i = 0; i < ideal.size(); ++i) {
        if (i > 0) idealStr += ",";
        idealStr += ideal[i];
    }

    std::string script = R"(
ring r = <char>, (<vars>), (dp);
ideal I = <ideal>;
ideal gb = groebner(I);
int d = dim(gb);
)";

    auto result = runSingularTemplate(script,
        {{"char", std::to_string(modulus)},
         {"vars", vars},
         {"ideal", idealStr}},
        {"d"});

    if (result.count("d")) {
        std::string dStr = result["d"];
        while (!dStr.empty() && (dStr[0] == ' ' || dStr[0] == '\n')) dStr.erase(0, 1);
        while (!dStr.empty() && (dStr.back() == ' ' || dStr.back() == '\n')) dStr.pop_back();
        return std::stoi(dStr);
    }
    return -1;
}

} // namespace SingularRunner
#endif
