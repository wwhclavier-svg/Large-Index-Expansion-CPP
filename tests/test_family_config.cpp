// test_family_config.cpp — Verify FamilyConfig JSON parser
#include <iostream>
#include <exception>
#include "FamilyConfig.hpp"

int main() {
    std::vector<std::string> families = {"SR","bub00","SR3m","bub10","bub11","Tri","Box","SR5m"};
    bool allPass = true;

    for (auto& f : families) {
        try {
            FamilyDef fam = parseFamilyConfig("families/" + f + ".json");
            bool ok = true;

            if (fam.nProp() != (int)fam.topSector.size())  { ok = false; }
            if (fam.modulus != 179424673)                  { ok = false; }
            if (fam.name != f)                             { ok = false; }

            std::cout << f << ": nProp=" << fam.nProp() << " nLoop=" << fam.nLoop()
                      << " nExt=" << fam.nExt() << " mod=" << fam.modulus
                      << " — " << (ok ? "PASS" : "FAIL") << std::endl;

            if (!ok) allPass = false;
        } catch (const std::exception& e) {
            std::cout << f << ": FAIL — " << e.what() << std::endl;
            allPass = false;
        }
    }

    std::cout << "\n=== " << (allPass ? "ALL PASS" : "SOME FAILED") << " ===" << std::endl;
    return allPass ? 0 : 1;
}
