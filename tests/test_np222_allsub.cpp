// test_np222_allsub.cpp — 对 7-prop NP222 的所有 C(7,5)=21 个五传播子子扇区
// dump AB 特征系统到 tools/ab_sub/ 目录（每个扇区一个文件）
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdint>

int main() {
    int64_t mod = 179424673;
    auto fam = parseFamilyConfig("families/NP222.json");
    auto ibp = IBPEqGenerator::generateIBPEquations(fam);
    std::cout << "ne=" << ibp.ne << "  nibp=" << ibp.nibp << std::endl;

    std::filesystem::create_directories("tools/ab_sub");
    // 枚举丢掉的传播子对 (i,j)，i<j，0-based
    for (int i = 0; i < 7; ++i) {
        for (int j = i + 1; j < 7; ++j) {
            std::vector<int> sector(7, 1);
            sector[i] = 0;
            sector[j] = 0;
            auto abEqs = IBPAnalyzer::buildABEquations(ibp, sector, mod);
            std::string fname = "tools/ab_sub/drop" + std::to_string(i + 1)
                              + std::to_string(j + 1) + ".txt";
            std::ofstream ofs(fname);
            for (size_t k = 0; k < abEqs.size(); ++k) {
                if (k) ofs << ",\n";
                ofs << abEqs[k];
            }
            ofs << "\n";
            std::cout << fname << "  (" << abEqs.size() << " eqs)" << std::endl;
        }
    }
    return 0;
}
