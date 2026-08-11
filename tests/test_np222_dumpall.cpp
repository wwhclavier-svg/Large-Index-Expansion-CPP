// test_np222_dumpall.cpp — dump 全部 792 个候选族的 top-sector AB 系统
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

int main() {
    int64_t mod = 179424673;
    std::vector<std::string> pool = {
        "l1^2", "l2^2", "(l1+l2)^2", "(l1-l2)^2",
        "(l1+p)^2", "(l1-p)^2", "(l2+p)^2", "(l2-p)^2",
        "(l1+l2+p)^2", "(l1+l2-p)^2", "(l1-l2+p)^2", "(l1-l2-p)^2"
    };
    int NP = pool.size();
    std::filesystem::create_directories("tools/hunt_all");
    std::vector<int> idx(5);
    int count = 0, ok = 0;
    for (idx[0] = 0; idx[0] < NP; ++idx[0])
    for (idx[1] = idx[0]+1; idx[1] < NP; ++idx[1])
    for (idx[2] = idx[1]+1; idx[2] < NP; ++idx[2])
    for (idx[3] = idx[2]+1; idx[3] < NP; ++idx[3])
    for (idx[4] = idx[3]+1; idx[4] < NP; ++idx[4]) {
        ++count;
        FamilyDef fam;
        fam.name = "hunt";
        for (int t : idx) fam.propagators.push_back(pool[t]);
        fam.loopMomenta = {"l1", "l2"};
        fam.externalMomenta = {"p"};
        fam.kinematicRules = {{"p^2", "s"}};
        fam.topSector = {1,1,1,1,1};
        fam.numeric = {{"s","1"},{"d","1/13"}};
        fam.modulus = mod;
        try {
            auto ibp = IBPEqGenerator::generateIBPEquations(fam);
            if (ibp.ne != 5 || ibp.nibp != 6) continue;
            auto abEqs = IBPAnalyzer::buildABEquations(ibp, fam.topSector, mod);
            std::string fn = "tools/hunt_all/fam";
            for (int t : idx) fn += "_" + std::to_string(t);
            fn += ".txt";
            std::ofstream ofs(fn);
            for (size_t e = 0; e < abEqs.size(); ++e)
                ofs << (e ? ",\n" : "") << abEqs[e];
            ofs << "\n";
            ++ok;
        } catch (...) { continue; }
    }
    std::cerr << "total=" << count << " dumped=" << ok << std::endl;
    return 0;
}
