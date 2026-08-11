// test_np222_hunt.cpp — 在 12 个传播子池中枚举 C(12,5)=792 个 2L1P 族，
// 用 meta 的 4 个 top-sector regime 点（120 种坐标置换）代入 AB 特征系统，
// 找出与 RelationMeta_NP222.m 一致的族。
#include "../include/FamilyConfig.hpp"
#include "../include/IBPEqGenerator.hpp"
#include "../include/IBPAnalyzer.hpp"
#include <array>
#include <fstream>
#include <set>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

static const int64_t MOD = 179424673;

// AB 单项：c + c*A_j + c*B_i*A_j
struct Mon { int64_t c; int b; int a; }; // b=0 表示无 B，a=0 表示无 A

static std::vector<Mon> parseEq(const std::string& s) {
    std::vector<Mon> out;
    size_t i = 0, n = s.size();
    while (i < n) {
        size_t j = s.find('+', i);
        std::string t = s.substr(i, j == std::string::npos ? n : j - i);
        i = (j == std::string::npos) ? n : j + 1;
        Mon m{1, 0, 0};
        // 拆分 * 因子
        std::vector<std::string> facs;
        size_t k = 0;
        while (true) {
            size_t l = t.find('*', k);
            facs.push_back(t.substr(k, l == std::string::npos ? t.size() : l - k));
            if (l == std::string::npos) break;
            k = l + 1;
        }
        for (auto& f : facs) {
            if (f.empty()) continue;
            if (f[0] == 'A') m.a = std::stoi(f.substr(1));
            else if (f[0] == 'B') m.b = std::stoi(f.substr(1));
            else m.c = std::stoll(f) % MOD;
        }
        out.push_back(m);
    }
    return out;
}

int main() {
    // meta top-sector 4 个 NB=1 regime 的 A 点
    std::vector<std::array<int64_t, 5>> regs = {
        {179424671, 2, 179424671, 2, 179424672},
        {179424668, 179424668, 179424668, 179424668, 5},
        {112304868, 157489941, 112304868, 157489941, 134239605},
        {67119801, 21934740, 67119801, 21934740, 45185066}
    };
    // 模逆（B = 1/A）
    auto modinv = [](int64_t x) {
        int64_t a = x % MOD, b = MOD, u = 1, v = 0;
        while (b) { int64_t q = a / b; std::tie(a, b) = std::make_pair(b, a - q * b); std::tie(u, v) = std::make_pair(v, u - q * v); }
        return (u % MOD + MOD) % MOD;
    };

    std::vector<std::string> pool = {
        "l1^2", "l2^2", "(l1+l2)^2", "(l1-l2)^2",
        "(l1+p)^2", "(l1-p)^2", "(l2+p)^2", "(l2-p)^2",
        "(l1+l2+p)^2", "(l1+l2-p)^2", "(l1-l2+p)^2", "(l1-l2-p)^2"
    };
    int NP = pool.size();

    // 120 个置换
    std::vector<std::array<int,5>> perms = {{0,1,2,3,4}};
    {
        std::array<int,5> base = {0,1,2,3,4};
        do { perms.push_back(base); } while (std::next_permutation(base.begin(), base.end()));
        perms.erase(perms.begin());
    }

    int tested = 0, matched = 0;
    std::set<std::vector<int>> dumped;
    std::filesystem::create_directories("tools/hunt");
    std::vector<int> idx(5);
    for (idx[0] = 0; idx[0] < NP; ++idx[0])
    for (idx[1] = idx[0]+1; idx[1] < NP; ++idx[1])
    for (idx[2] = idx[1]+1; idx[2] < NP; ++idx[2])
    for (idx[3] = idx[2]+1; idx[3] < NP; ++idx[3])
    for (idx[4] = idx[3]+1; idx[4] < NP; ++idx[4]) {
        ++tested;
        FamilyDef fam;
        fam.name = "hunt";
        for (int t : idx) fam.propagators.push_back(pool[t]);
        fam.loopMomenta = {"l1", "l2"};
        fam.externalMomenta = {"p"};
        fam.kinematicRules = {{"p^2", "s"}};
        fam.topSector = {1,1,1,1,1};
        fam.numeric = {{"s","1"},{"d","1/13"}};
        fam.modulus = MOD;

        std::vector<std::string> abEqs;
        try {
            auto ibp = IBPEqGenerator::generateIBPEquations(fam);
            if (ibp.ne != 5 || ibp.nibp != 6) continue;
            abEqs = IBPAnalyzer::buildABEquations(ibp, fam.topSector, MOD);
        } catch (...) { continue; }

        std::vector<std::vector<Mon>> eqs;
        for (auto& e : abEqs) eqs.push_back(parseEq(e));

        for (auto& pm : perms) {
            bool ok = true;
            for (auto& r : regs) {
                int64_t A[5], B[5];
                for (int t = 0; t < 5; ++t) { A[t] = r[pm[t]]; B[t] = modinv(A[t]); }
                for (auto& eq : eqs) {
                    int64_t s = 0;
                    for (auto& m : eq) {
                        int64_t term = m.c;
                        if (m.b) term = term * B[m.b-1] % MOD;
                        if (m.a) term = term * A[m.a-1] % MOD;
                        s = (s + term) % MOD;
                    }
                    if (s != 0) { ok = false; break; }
                }
                if (!ok) break;
            }
            if (ok) {
                ++matched;
                std::cout << "MATCH props={";
                for (int t = 0; t < 5; ++t) std::cout << (t?",":"") << pool[idx[t]];
                std::cout << "} perm={";
                for (int t = 0; t < 5; ++t) std::cout << (t?",":"") << pm[t];
                std::cout << "}" << std::endl;
                // dump AB 系统（每个族只 dump 一次）
                if (dumped.insert(idx).second) {
                    std::string fn = "tools/hunt/fam";
                    for (int t : idx) fn += "_" + std::to_string(t);
                    fn += ".txt";
                    std::ofstream ofs(fn);
                    for (size_t e = 0; e < abEqs.size(); ++e)
                        ofs << (e ? ",\n" : "") << abEqs[e];
                    ofs << "\n";
                }
            }
        }
    }
    std::cerr << "tested=" << tested << " matched=" << matched << std::endl;
    return 0;
}
