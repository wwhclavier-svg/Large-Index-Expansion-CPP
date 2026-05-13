#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>

#include "firefly/FFInt.hpp"
#include "RelationLoader.hpp"
#include "SymbolicRule.hpp"

void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " <allrelations.m> [sector_bits] [level_bound] [--verbose]\n";
    std::cerr << "  allrelations.m : AllRelation_*.m 文件路径\n";
    std::cerr << "  sector_bits    : 二进制 sector (默认全1)\n";
    std::cerr << "  level_bound    : level 上限 (默认 10)\n";
    std::cerr << "  --verbose,-v   : 显示 GB 对角元详情\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) { printUsage(argv[0]); return 1; }

    std::string relPath = argv[1];
    std::vector<int> sector;
    int levelBound = 10;
    bool verbose = false;

    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--verbose" || a == "-v") verbose = true;
        else if (a.find_first_not_of("01") == std::string::npos && sector.empty())
            for (char c : a) sector.push_back(c == '1' ? 1 : 0);
        else if (a.find_first_not_of("0123456789") == std::string::npos)
            levelBound = std::stoi(a);
    }

    try {
        auto t0 = std::chrono::high_resolution_clock::now();

        auto rd = loadAllRelations(relPath);
        if (sector.empty()) sector.resize(rd.ne, 1);

        std::cout << "Loaded " << rd.entries.size() << " entries | Family=" << rd.family
                  << " NE=" << rd.ne << " Modulus=" << rd.modulus << " Order=" << rd.order << "\n";
        std::cout << "Sector: [";
        for (size_t i=0;i<sector.size();++i) { if(i)std::cout<<","; std::cout<<sector[i]; }
        std::cout << "]\n\n";

        auto t1 = std::chrono::high_resolution_clock::now();
        auto res = setupGeneratingCone(rd, sector, levelBound);
        res.printSummary(std::cout);

        if (verbose) {
            auto printConeDetail = [&](const std::string& name, const auto& cones) {
                for (auto& c : cones) {
                    if (c.gbDiagonal.empty() && c.reducedVarIdx.empty()) continue;
                    std::cout << name << " level=" << c.level << ":\n";
                    if (!c.reducedVarIdx.empty()) {
                        std::cout << "  reducedVars:";
                        for (int vi : c.reducedVarIdx) std::cout << " " << vi;
                        std::cout << "\n";
                    }
                    if (!c.gbDiagonal.empty()) {
                        std::cout << "  GB diagonals:\n";
                        for (size_t i=0;i<c.gbDiagonal.size();++i) {
                            std::string d = c.gbDiagonal[i];
                            bool hasNu=false;
                            for (int j=0;j<rd.ne;++j)
                                if(d.find("nu"+std::to_string(j))!=std::string::npos) hasNu=true;
                            std::cout << "    [" << i << "] " << d
                                      << (hasNu ? " ← 含 ISP 变量!" : "")
                                      << "\n";
                        }
                    }
                    if (!c.singularISP.empty()) {
                        std::cout << "  ISP singularities:";
                        for (auto& isp : c.singularISP) {
                            std::cout << " {";
                            for (int v : isp) std::cout << v << ",";
                            std::cout << "}";
                        }
                        std::cout << "\n";
                    }
                    std::cout << "  flags: " << c.flags.toString() << "\n\n";
                }
            };
            printConeDetail("B1", res.B1);
            printConeDetail("B2", res.B2);
            printConeDetail("B3", res.B3);

            // Build symbolic rules with shift info
            std::vector<RelationEntry> pent;
            for (auto& e : rd.entries) if (e.ansatzMode=="Pyramid") pent.push_back(e);
            auto symRules = generateSymbolicRules(res, pent, rd.ne, sector);
            std::cout << "--- Symbolic Rules (" << symRules.size() << "总) 前20条 ---\n";
            for (size_t i = 0; i < symRules.size() && i < 20; ++i) {
                std::cout << "  [" << i << "] " << symRules[i].toString(rd.ne);
                if (symRules[i].hasISPSingularity) {
                    std::cout << "  ⚠ ISP: {";
                    for (int v : symRules[i].ispVars) std::cout << "ν" << v << ",";
                    std::cout << "}";
                }
                std::cout << "\n";
            }
            if (symRules.size() > 10) std::cout << "  ... 共" << symRules.size() << "条\n";
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        auto rpt = runVerification(rd, sector, levelBound);
        rpt.print(std::cout);

        auto t3 = std::chrono::high_resolution_clock::now();
        double t_setup = std::chrono::duration<double>(t1 - t0).count();
        double t_cone  = std::chrono::duration<double>(t2 - t1).count();
        double t_vfy   = std::chrono::duration<double>(t3 - t2).count();
        double t_total = std::chrono::duration<double>(t3 - t0).count();
        std::cout << "\nTiming: load=" << t_setup << "s cone=" << t_cone
                  << "s verify=" << t_vfy << "s total=" << t_total << " s\n";

    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
