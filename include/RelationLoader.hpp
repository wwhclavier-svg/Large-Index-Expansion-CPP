#ifndef RELATION_LOADER_HPP
#define RELATION_LOADER_HPP

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <iostream>
#include "firefly/FFInt.hpp"

struct RelationEntry {
    std::string family;
    int lev = 0;
    int deg = 0;
    int ne = 0;
    std::string ansatzMode;
    int64_t modulus = 0;
    int order = 0;
    int stableOrder = 0;
    int numRelations = 0;
    int activeVars = 0;
    int totalVars = 0;
    int numVariables = 0;
    int numSolutions = 0;
    bool hasSolution = false;
    std::vector<std::vector<int>> alphas;
    std::vector<std::vector<int>> betas;
    std::vector<std::vector<firefly::FFInt>> coefficients;
    std::vector<std::pair<std::vector<int>, std::vector<int>>> independentPairs;
};

struct RelationData {
    std::string family;
    int ne = 0;
    int64_t modulus = 0;
    int order = 0;
    std::vector<RelationEntry> entries;
};

inline std::string trimString(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

inline std::vector<int> parseMMAIntList(const std::string& s) {
    std::vector<int> result;
    std::string clean;
    for (char c : s) {
        if (c == '{' || c == '}' || c == ' ' || c == '\n' || c == '\r') continue;
        clean += c;
    }
    std::stringstream ss(clean);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (!token.empty()) result.push_back(std::stoi(token));
    }
    return result;
}

inline std::vector<std::vector<int>> parseMMAMatrix(const std::string& s) {
    std::vector<std::vector<int>> result;
    std::regex rowRegex(R"(\{([^}]+)\})");
    auto begin = std::sregex_iterator(s.begin(), s.end(), rowRegex);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        result.push_back(parseMMAIntList(it->str()));
    }
    return result;
}

inline std::vector<firefly::FFInt> parseMMACoeffList(const std::string& s, int64_t modulus) {
    std::vector<firefly::FFInt> result;
    std::string clean;
    for (char c : s) {
        if (c == '{' || c == '}' || c == ' ' || c == '\n' || c == '\r') continue;
        clean += c;
    }
    std::stringstream ss(clean);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (!token.empty()) {
            long long val = std::stoll(token);
            val %= modulus;
            if (val < 0) val += modulus;
            result.emplace_back(static_cast<uint64_t>(val));
        }
    }
    return result;
}

inline std::vector<std::pair<std::vector<int>, std::vector<int>>>
parseMMAIndependentPairs(const std::string& s) {
    std::vector<std::pair<std::vector<int>, std::vector<int>>> result;
    std::regex pairRegex(R"(\{\{([^}]+)\},\{([^}]+)\}\})");
    auto begin = std::sregex_iterator(s.begin(), s.end(), pairRegex);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        auto lhs = parseMMAIntList(it->str(1));
        auto rhs = parseMMAIntList(it->str(2));
        result.emplace_back(lhs, rhs);
    }
    return result;
}

inline std::string extractValue(const std::string& text, const std::string& key) {
    // Find the key position
    std::string searchKey = "\"" + key + "\"";
    size_t pos = text.find(searchKey);
    if (pos == std::string::npos) return "";

    // Find " -> " after the key
    pos = text.find("->", pos + searchKey.size());
    if (pos == std::string::npos) return "";
    pos += 2; // skip "->"

    // Skip whitespace
    while (pos < text.size() && (text[pos] == ' ' || text[pos] == '\t' || text[pos] == '\n' || text[pos] == '\r'))
        ++pos;
    if (pos >= text.size()) return "";

    // Extract value based on first character
    char first = text[pos];
    std::string val;

    if (first == '{') {
        // Brace-balanced extraction
        int depth = 0;
        for (size_t i = pos; i < text.size(); ++i) {
            if (text[i] == '{') ++depth;
            else if (text[i] == '}') { --depth; val += text[i]; if (depth == 0) break; continue; }
            val += text[i];
        }
    } else if (first == '"') {
        // Quoted string
        ++pos;
        for (size_t i = pos; i < text.size(); ++i) {
            if (text[i] == '"') break;
            val += text[i];
        }
    } else {
        // Simple value (number, True/False)
        for (size_t i = pos; i < text.size(); ++i) {
            if (text[i] == ',' || text[i] == '\n' || text[i] == '\r' || text[i] == '}')
                break;
            val += text[i];
        }
        val = trimString(val);
    }

    return val;
}

inline RelationData loadAllRelations(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filepath);

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string text = buffer.str();

    RelationData data;
    std::regex entryRegex(R"(<\|\s*([\s\S]*?)\s*\|>)");
    auto begin = std::sregex_iterator(text.begin(), text.end(), entryRegex);
    auto end = std::sregex_iterator();

    for (auto it = begin; it != end; ++it) {
        std::string block = it->str(1);
        RelationEntry entry;

        auto val = [&](const std::string& k) { return extractValue(block, k); };

        auto readInt = [&](const std::string& k) -> int {
            std::string v = val(k);
            return v.empty() ? 0 : std::stoi(v);
        };

        auto readLong = [&](const std::string& k) -> int64_t {
            std::string v = val(k);
            return v.empty() ? 0 : std::stoll(v);
        };

        entry.family = val("Family");
        entry.lev = readInt("Lev");
        entry.deg = readInt("Deg");
        entry.ne = readInt("NE");
        entry.ansatzMode = val("AnsatzMode");
        entry.modulus = readLong("Modulus");
        entry.order = readInt("Order");
        entry.stableOrder = readInt("StableOrder");
        entry.numRelations = readInt("NumRelations");
        entry.activeVars = readInt("ActiveVars");
        entry.totalVars = readInt("TotalVars");
        entry.numVariables = readInt("NumVariables");
        entry.numSolutions = readInt("NumSolutions");

        std::string hs = val("HasSolution");
        entry.hasSolution = (hs == "True");

        std::string alphasStr = val("Alphas");
        if (!alphasStr.empty() && alphasStr.size() > 2) {
            int depth = 0;
            std::string current;
            for (size_t ci = 0; ci < alphasStr.size(); ++ci) {
                char c = alphasStr[ci];
                if (c == '{') { ++depth; if (depth == 2) current.clear(); }
                else if (c == '}') { if (depth == 2) entry.alphas.push_back(parseMMAIntList(current)); --depth; }
                else if (depth == 2) current += c;
            }
        }

        std::string betasStr = val("Betas");
        if (!betasStr.empty() && betasStr.size() > 2) {
            int depth = 0;
            std::string current;
            for (size_t ci = 0; ci < betasStr.size(); ++ci) {
                char c = betasStr[ci];
                if (c == '{') { ++depth; if (depth == 2) current.clear(); }
                else if (c == '}') { if (depth == 2) entry.betas.push_back(parseMMAIntList(current)); --depth; }
                else if (depth == 2) current += c;
            }
        }

        std::string coeffsStr = val("Coefficients");
        if (!coeffsStr.empty() && entry.modulus > 0 && coeffsStr.size() > 2) {
            // Brace-level parsing: iterate through chars, extract rows at depth 1→0
            std::vector<std::string> rows;
            int depth = 0;
            std::string current;
            for (size_t ci = 0; ci < coeffsStr.size(); ++ci) {
                char c = coeffsStr[ci];
                if (c == '{') {
                    ++depth;
                    if (depth == 2) current.clear(); // start of inner row
                } else if (c == '}') {
                    if (depth == 2) rows.push_back(current); // end of inner row
                    --depth;
                } else if (depth == 2) {
                    current += c;
                }
            }
            for (auto& rowStr : rows) {
                std::vector<firefly::FFInt> ffRow;
                std::stringstream ss(rowStr);
                std::string token;
                while (std::getline(ss, token, ',')) {
                    std::string trimmed;
                    for (char c : token) if (c != ' ' && c != '\n' && c != '\r' && c != '\t') trimmed += c;
                    if (trimmed.empty()) continue;
                    long long v = std::stoll(trimmed);
                    v %= entry.modulus; if (v < 0) v += entry.modulus;
                    ffRow.emplace_back(static_cast<uint64_t>(v));
                }
                entry.coefficients.push_back(ffRow);
            }
        }

        std::string ipStr = val("IndependentPairs");
        if (!ipStr.empty()) entry.independentPairs = parseMMAIndependentPairs(ipStr);

        if (data.family.empty()) data.family = entry.family;
        if (data.ne == 0) data.ne = entry.ne;
        if (data.modulus == 0) data.modulus = entry.modulus;
        if (data.order == 0) data.order = entry.order;

        data.entries.push_back(entry);
    }

    return data;
}

#endif
