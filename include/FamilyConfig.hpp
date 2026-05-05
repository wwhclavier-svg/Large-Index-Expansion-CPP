#ifndef FAMILY_CONFIG_HPP
#define FAMILY_CONFIG_HPP

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <stdexcept>
#include "json.hpp"

// ==========================================
// FamilyDef — IBP integral family definition
// ==========================================
struct FamilyDef {
    std::string name;
    std::string description;
    std::vector<std::string> propagators;
    std::vector<std::string> loopMomenta;
    std::vector<std::string> externalMomenta;
    std::map<std::string, std::string> kinematicRules; // var -> expr
    std::vector<int> topSector;
    std::map<std::string, std::string> numeric;         // param -> value (string preserves fractions)
    int64_t modulus = 0;

    // Derived quantities
    int nProp() const { return static_cast<int>(propagators.size()); }
    int nLoop() const { return static_cast<int>(loopMomenta.size()); }
    int nExt()  const { return static_cast<int>(externalMomenta.size()); }
};

// ==========================================
// Parse family.json → FamilyDef
// ==========================================
inline FamilyDef parseFamilyConfig(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f) throw std::runtime_error("Cannot open family config: " + filepath);

    nlohmann::json j;
    f >> j;

    FamilyDef fam;
    fam.name        = j.value("name", "");
    fam.description = j.value("description", "");
    fam.modulus     = j.value("modulus", 0LL);

    for (const auto& p : j.at("propagators"))
        fam.propagators.push_back(p.get<std::string>());

    for (const auto& lm : j.at("loopMomenta"))
        fam.loopMomenta.push_back(lm.get<std::string>());

    for (const auto& em : j.at("externalMomenta"))
        fam.externalMomenta.push_back(em.get<std::string>());

    for (const auto& [k, v] : j.at("kinematicRules").items())
        fam.kinematicRules[k] = v.get<std::string>();

    for (const auto& s : j.at("topSector"))
        fam.topSector.push_back(s.get<int>());

    for (const auto& [k, v] : j.at("numeric").items()) {
        if (v.is_number_integer())
            fam.numeric[k] = std::to_string(v.get<int64_t>());
        else if (v.is_number_float())
            fam.numeric[k] = std::to_string(v.get<double>());
        else
            fam.numeric[k] = v.get<std::string>();
    }

    return fam;
}
#endif
