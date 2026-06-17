#ifndef SYMBOLIC_RULE_HPP
#define SYMBOLIC_RULE_HPP

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <functional>

#include "firefly/FFInt.hpp"
#include "RelationLoader.hpp"
#include "SingularRunner.hpp"
#include "Combinatorics.hpp"

// ═════════════════════════════════════════════════════════════
// 输出数据结构
// ═════════════════════════════════════════════════════════════

// 模块 GB 约化规则: gen(src) + Σ c_j * gen(tgt_j) = 0
// → g(ν-α_src) = Σ (-c_j) * g(ν-α_tgt_j)  (c_j ∈ FF[ν])
struct ModuleReductionRule {
    int srcCol = -1;
    std::string srcCoeffStr;
    std::vector<std::pair<int, std::string>> targets;
};

struct CompletenessFlags {
    bool TopFlag1 = false, TopFlag2 = false;
    bool FullFlag1 = false, FullFlag2 = false, GlobalFlag = false;
    std::string toString() const {
        return std::string(TopFlag1 ? "T1✓ ":"T1✗ ")+(TopFlag2 ? "T2✓ ":"T2✗ ")
             +(FullFlag1 ? "F1✓ ":"F1✗ ")+(FullFlag2 ? "F2✓ ":"F2✗ ")
             +(GlobalFlag ? "G✓":"G✗");
    }
};

struct ConeRule {
    int level = 0;
    std::vector<int> reducedVarIdx;       // 可约化变量在 αShifts 中的索引 (original indices)
    std::vector<std::vector<int>> singularISP;
    std::vector<std::string> gbDiagonal;  // GB 对角元多项式 (POT order)
    CompletenessFlags flags;
    std::string singularLocus;
    // P6: concrete reduction rules from Gaussian elimination
    // Each entry: {pivotVarIdx, [(targetVarIdx, coeff), ...]}
    //   meaning: g(alphaShifts[pivotVarIdx]) = -Σ coeff_j * g(alphaShifts[targetVarIdx_j])
    std::vector<std::vector<std::pair<int,firefly::FFInt>>> reductionRules;
    std::vector<std::vector<int>> alphaShifts;  // copy of BlockMatrix.alphaShifts for mapping
    std::vector<ModuleReductionRule> moduleRules;
};

struct GeneratingConeResult {
    std::vector<ConeRule> B1, B2, B3;
    int levelBound = 0, topRank = 0, ne = 0;
    bool fullCoverage = false;
    void printSummary(std::ostream& os) const {
        os << "=== GeneratingCone Result ===\nNE="<<ne<<" LevelBound="<<levelBound<<" TopRank="<<topRank<<"\n";
        auto pc = [&](const std::string& n, const auto& cv) {
            os<<n<<" ("<<cv.size()<<" cones):\n";
            for (auto& c : cv) os<<"  level="<<c.level<<" redVars="<<c.reducedVarIdx.size()
              <<" singISP="<<c.singularISP.size()<<" gbDiag="<<c.gbDiagonal.size()<<" ["<<c.flags.toString()<<"]\n";
        };
        pc("B1 high-rank", B1); pc("B2 mid-rank", B2); pc("B3 near-corner", B3);
        os << "Full coverage: "<<(fullCoverage?"YES":"NO")<<"\n";
    }
};

struct SymbolicRule {
    std::vector<int> sourceShift;
    std::vector<std::vector<int>> targetShifts;
    std::vector<std::string> coeffStrs;       // ν-多项式系数字符串（已取反: -c_j/c_src）
    bool hasISPSingularity = false;
    std::vector<int> ispVars;
    std::string toString(int ne) const;
};

// ═════════════════════════════════════════════════════════════
// 种子生成
// ═════════════════════════════════════════════════════════════

inline std::vector<std::vector<int>> generateDotSeeds(int ne, int level, bool topRankOnly=false, int minRank=0) {
    std::vector<std::vector<int>> seeds;
    std::vector<int> cur(ne,0);
    std::function<void(int,int)> rec = [&](int p, int r) {
        if (p==ne) { if (!r) { int tr=std::accumulate(cur.begin(),cur.end(),0); if(!topRankOnly||tr>=minRank) seeds.push_back(cur); } return; }
        for (int v=0;v<=r;++v) { cur[p]=v; rec(p+1,r-v); }
    };
    for (int t=0;t<=level;++t) rec(0,t);
    return seeds;
}

// MinusRank seeds (MMA minusranklist, SymbolicBTForm.wl:242-245)
// Flat mode: generate sparse 0/1 vectors with 1's on ISP positions
// from Subsets[ispset, mrnkLevel] — each represents a supersector shift
inline std::vector<std::vector<int>> generateMinusSeeds(
    int ne, int mrnkLevel, const std::vector<int>& sector)
{
    std::vector<std::vector<int>> seeds;
    if (mrnkLevel <= 0) return seeds;

    // Collect ISP positions (sector[i] == 0)
    std::vector<int> ispPos;
    for (int i = 0; i < ne; ++i)
        if (i < (int)sector.size() && sector[i] == 0) ispPos.push_back(i);

    int nIsp = (int)ispPos.size();
    if (mrnkLevel > nIsp || nIsp == 0) return seeds;

    // Generate all subsets of ispPos of size mrnkLevel (MMA: Subsets[ispset, mrnklevel])
    std::vector<int> subset(mrnkLevel);
    std::function<void(int, int)> genSubsets = [&](int start, int depth) {
        if (depth == mrnkLevel) {
            std::vector<int> vec(ne, 0);
            for (int pos : subset) vec[pos] = 1;
            seeds.push_back(vec);
            return;
        }
        for (int i = start; i < nIsp; ++i) {
            subset[depth] = ispPos[i];
            genSubsets(i + 1, depth + 1);
        }
    };
    genSubsets(0, 0);
    return seeds;
}

// ═════════════════════════════════════════════════════════════
// Block 矩阵构建 (高斯消元用简化版)
// 行 = numRelations × seeds，列 = α-shifts 去重
// 系数映射简化但不影响秩
// ν-多项式结构留给模块 GB 独立处理
// ═════════════════════════════════════════════════════════════

struct BlockMatrix {
    std::vector<std::vector<firefly::FFInt>> coeffs;
    std::vector<std::vector<int>> alphaShifts;       // 每列对应的 ν-shift
    std::vector<int> rankPerVar;
    int nRelations = 0, nVars = 0;
};

inline BlockMatrix buildBlockMatrix(
    const std::vector<RelationEntry>& entries,
    const std::vector<std::vector<int>>& seeds,
    const std::vector<int>& sector,
    int ne, int64_t modulus)
{
    BlockMatrix bm;
    std::map<std::vector<int>,int> s2c;
    // Phase 1: column map from α-shifts with sector-aware signs
    // pdset (sector[j]=1): sd = -seed  → shift = sd + α = -seed + α
    // ispset (sector[j]=0): sd = +seed → shift = sd + α = +seed + α
    for (auto& seed : seeds)
        for (auto& entry : entries) {
            if (entry.numRelations==0) continue;
            for (auto& alpha : entry.alphas) {
                std::vector<int> key(ne);
                for (int j=0;j<ne;++j)
                    key[j] = sector[j] ? (-seed[j] + alpha[j]) : (seed[j] + alpha[j]);
                if (!s2c.count(key)) {
                    s2c[key]=s2c.size();
                    bm.alphaShifts.push_back(key);
                    bm.rankPerVar.push_back(std::accumulate(key.begin(),key.end(),0));
                }
            }
        }
    bm.nVars = s2c.size();

    // Phase 2: fill coefficients — one row per (seed, α_i, relation_k)
    // Sum over β_j: coefficient = Σ_j b0[ai*nb+bj][k+1]  (MMA: CoefficientArrays sum over β)
    // This avoids redundant rows since all β for the same (seed,α,k) share a column.
    int rowCap = 50000;
    std::vector<std::vector<firefly::FFInt>> coeffs;
    coeffs.reserve(rowCap);
    for (auto& entry : entries) {
        if (entry.numRelations==0||entry.alphas.empty()) continue;
        int nb = entry.betas.size();
        for (auto& seed : seeds) {
            for (int ai = 0; ai < (int)entry.alphas.size() && (int)coeffs.size() < rowCap; ++ai) {
                // Column shift = sd + α_ai
                std::vector<int> colKey(ne);
                for (int j=0;j<ne;++j)
                    colKey[j] = sector[j] ? (-seed[j] + entry.alphas[ai][j]) : (seed[j] + entry.alphas[ai][j]);
                
                // sectorFilter: skip if pdset component < 0
                bool skip=false;
                for (int j=0;j<ne;++j)
                    if (sector[j] && colKey[j] < 0) { skip=true; break; }
                if (skip) continue;

                auto colIt = s2c.find(colKey);
                if (colIt == s2c.end()) continue;
                int ci = colIt->second;

                for (int k = 0; k < entry.numRelations && (int)coeffs.size() < rowCap; ++k) {
                    std::vector<firefly::FFInt> row(bm.nVars, firefly::FFInt(0));
                    // Sum β-contributions: Σ_j b0[ai*nb+bj][k+1]
                    for (int bj = 0; bj < nb; ++bj) {
                        int ri = ai * nb + bj;
                        if (ri >= (int)entry.coefficients.size()) continue;
                        auto& cr = entry.coefficients[ri];
                        if (k + 1 < (int)cr.size())
                            row[ci] = row[ci] + cr[k + 1];
                    }
                    if (row[ci].n != 0) coeffs.push_back(row);
                }
            }
        }
    }
    bm.nRelations = (int)coeffs.size();
    bm.coeffs = coeffs;
    return bm;
}

// ═════════════════════════════════════════════════════════════
// Helper: extract leading e_i index and its coefficient from a Singular GB poly string
// In POT order, the first e_i in each GB element is the diagonal position.
// Returns {eIndex, coeffIsOne, dependsOnNu, nuVarIndices}
// coeffIsOne=true means the coefficient of e_i is 1 (unit diagonal)
// dependsOnNu=true means the coefficient of e_i contains ν variables
// ═════════════════════════════════════════════════════════════

struct LeadingEInfo {
    int eIndex = -1;
    bool coeffIsOne = false;
    bool dependsOnNu = false;
    std::vector<int> nuVarIndices;
};

inline LeadingEInfo parseLeadingE(const std::string& poly, int ne) {
    LeadingEInfo info;

    // Split polynomial into terms at + and -
    // Keep track of sign for each term
    std::vector<std::pair<std::string, int>> terms; // {term, sign} sign=+1 or -1
    std::string current;
    int sign = 1;
    for (size_t i = 0; i < poly.size(); ++i) {
        char c = poly[i];
        if (c == '+') {
            if (!current.empty()) terms.push_back({current, sign});
            sign = 1;
            current.clear();
        } else if (c == '-') {
            if (!current.empty()) terms.push_back({current, sign});
            sign = -1;
            current.clear();
        } else if (c != ' ') {
            current += c;
        }
    }
    if (!current.empty()) terms.push_back({current, sign});

    // Find the smallest e_i index among all terms
    int minEIndex = 1000000;
    std::string eCoeffStr; // coefficient string for the leading e_i term

    for (auto& [term, sgn] : terms) {
        // Check if term contains e_i
        for (int i = 0; i < 100; ++i) { // scan reasonable range
            std::string eStr = "e" + std::to_string(i);
            if (term.find(eStr) != std::string::npos) {
                if (i < minEIndex) {
                    minEIndex = i;
                    // Extract coefficient: everything before "*e"
                    size_t starPos = term.find("*" + eStr);
                    if (starPos != std::string::npos) {
                        eCoeffStr = term.substr(0, starPos);
                    } else {
                        // Coefficient is the sign (implied 1 or -1)
                        eCoeffStr = (sgn == 1) ? "1" : "-1";
                    }
                }
                break;
            }
        }
    }

    if (minEIndex >= 100) return info; // no e_i found (pure ν relation)

    info.eIndex = minEIndex;

    // Check if coefficient is "1" or "-1" (unit)
    info.coeffIsOne = (eCoeffStr == "1" || eCoeffStr == "-1");

    // Check if coefficient depends on ν variables
    for (int j = 0; j < ne; ++j) {
        std::string ns = "nu" + std::to_string(j);
        if (eCoeffStr.find(ns) != std::string::npos ||
            poly.find(ns) != std::string::npos) {
            // Only mark as ν-dependent if ν appears in the coefficient of e_i
            // or in the entire polynomial (conservative but safer)
        }
    }

    // More precise: check if ν appears in the coefficient expression for e_i
    if (eCoeffStr.find("nu") != std::string::npos) {
        info.dependsOnNu = true;
        for (int j = 0; j < ne; ++j) {
            if (eCoeffStr.find("nu" + std::to_string(j)) != std::string::npos)
                info.nuVarIndices.push_back(j);
        }
    }

    return info;
}

// ═════════════════════════════════════════════════════════════
// 辅助: 展开 ν-多项式系数 c_i(ν-sd) = Σ_j b0_j · (ν-sd)^{β_j}
// 返回 Singular 多项式字符串，如 "3*nu0^2*nu1+5*nu0"
// 其中 sd 的符号由 sector 决定: sd[j] = sector[j] ? -seed[j] : seed[j]
// ═════════════════════════════════════════════════════════════

inline std::string buildNuPolynomial(
    const std::vector<firefly::FFInt>& b0Row,  // b0[β_j] for j=0..nb-1
    const std::vector<std::vector<int>>& betas,
    const std::vector<int>& seed,
    const std::vector<int>& sector,
    int ne, int64_t modulus)
{
    // 累加所有 β-项的单项式系数到统一的 exponent→coeff 映射
    // 展开 b0_j · (ν-sd)^{β_j} = b0_j · Π_m (ν_m - sd_m)^{β_j[m]}
    // 其中 sd_m = sector[m] ? -seed[m] : seed[m]
    // (ν_m - sd_m)^{e} = Σ_{t=0}^{e} C(e, t) · ν_m^t · (-sd_m)^{e-t}

    std::map<std::vector<int>, int64_t> total;

    for (size_t bj = 0; bj < betas.size() && bj < b0Row.size(); ++bj) {
        uint64_t b0val = b0Row[bj].n;
        if (b0val == 0) continue;

        // 从 {ν^0} 开始，逐维度展开 (ν_m - sd_m)^{β_j[m]}
        std::map<std::vector<int>, int64_t> terms;
        terms[std::vector<int>(ne, 0)] = 1;  // start with 1

        for (int m = 0; m < ne; ++m) {
            int e = betas[bj][m];
            if (e == 0) continue;
            int sd_m = sector[m] ? -seed[m] : seed[m];
            int64_t negSd = (-sd_m) % modulus;
            if (negSd < 0) negSd += modulus;

            std::map<std::vector<int>, int64_t> next;
            for (auto& [expV, c] : terms) {
                for (int t = 0; t <= e; ++t) {
                    auto newExp = expV;
                    newExp[m] += t;
                    int64_t binom = comb(e, t);
                    int64_t shiftPow = 1;
                    for (int p = 0; p < e - t; ++p)
                        shiftPow = (shiftPow * negSd) % modulus;
                    int64_t nc = (((c * binom) % modulus) * shiftPow) % modulus;
                    if (nc < 0) nc += modulus;
                    next[newExp] = (next[newExp] + nc) % modulus;
                }
            }
            terms = std::move(next);
        }

        // 乘以 b0val 并累加
        for (auto& [expV, c] : terms) {
            int64_t nc = (c * b0val) % modulus;
            total[expV] = (total[expV] + nc) % modulus;
        }
    }

    // 构建 Singular 多项式字符串
    std::string result;
    for (auto& [expV, coeff] : total) {
        if (coeff == 0) continue;
        if (!result.empty()) result += "+";

        bool hasNu = false;
        for (int v : expV) if (v != 0) { hasNu = true; break; }

        if (coeff != 1 || !hasNu)
            result += std::to_string(coeff);

        if (hasNu) {
            if (coeff != 1 || result.empty() || result.back() == '+')
                result += "*";
        }

        for (int m = 0; m < ne; ++m) {
            if (expV[m] == 0) continue;
            result += "nu" + std::to_string(m);
            if (expV[m] > 1) result += "^" + std::to_string(expV[m]);
            result += "*";
        }
        if (!result.empty() && result.back() == '*') result.pop_back();
    }

    return result.empty() ? "0" : result;
}

// ═════════════════════════════════════════════════════════════
// eliminatedRule (MMA SymbolicBTForm.wl:142-204)
// 使用 Singular POT 模块 GB + firstNonZero 分组提取对角元 + 五级完备性检验
// ═════════════════════════════════════════════════════════════

struct EliminatedRuleResult {
    std::vector<std::string> gb;                  // POT GB polynomials (ideal)
    std::vector<std::vector<int>> singularISP;    // ISP vars per diagonal
    std::vector<int> reducedVarTable;             // 可消去变量索引 (vorder-indexed)
    std::vector<int> vorder;                      // variable ordering
    std::vector<ModuleReductionRule> moduleRules; // gen(src) + Σ c_j * gen(tgt_j) = 0
    std::vector<std::string> gbDiagonal;          // display strings
    int blockStart=0, blockEnd=0, nVarsTotal=0;
    CompletenessFlags flags;
    int elimRank=0;
};

inline EliminatedRuleResult eliminatedRule(
    const BlockMatrix& bm,
    int level,
    int64_t modulus,
    int ne,
    bool topRankOnly=false,
    const std::vector<int>* numericReduce=nullptr,
    // 模块 GB 需要的原始数据（提供则启用 ν-多项式模块 GB 替代理想 GB）
    const std::vector<RelationEntry>* entriesForModule=nullptr,
    const std::vector<std::vector<int>>* seedsForModule=nullptr,
    const std::vector<int>* sectorForModule=nullptr)
{
    EliminatedRuleResult res;
    res.nVarsTotal = bm.nVars;

    // 1. Sort variables by rank descending (MMA IntgSort order-0)
    std::vector<int> vorder(bm.nVars);
    std::iota(vorder.begin(), vorder.end(), 0);
    std::sort(vorder.begin(), vorder.end(), [&](int a, int b) {
        if (bm.rankPerVar[a]!=bm.rankPerVar[b]) return bm.rankPerVar[a]>bm.rankPerVar[b];
        return a<b;
    });
    res.vorder = vorder;

    // 2. Determine block range (MMA L147-155)
    std::vector<int> highRank;
    for (int vi : vorder) if (bm.rankPerVar[vi]>=level) highRank.push_back(vi);
    res.blockStart = topRankOnly ? 0 : 0;
    res.blockEnd = topRankOnly ? (int)highRank.size() : bm.nVars;
    int nCols = res.blockEnd;
    int nRows = bm.nRelations;
    auto& mat = bm.coeffs;

    // 3. Gaussian elimination for TopFlag1 / FullFlag1 (MMA L157)
    if (nCols>0 && nRows>0) {
        std::vector<std::vector<firefly::FFInt>> red(nRows, std::vector<firefly::FFInt>(nCols));
        for (int i=0;i<nRows;++i) for (int j=0;j<nCols;++j) red[i][j]=mat[i][vorder[j]];

        int rnk=0;
        std::vector<int> pivotCols;
        for (int c=0;c<nCols && rnk<nRows;++c) {
            int pv=-1;
            for (int r=rnk;r<nRows;++r) if (red[r][c].n!=0) { pv=r; break; }
            if (pv<0) continue;
            std::swap(red[rnk], red[pv]);
            auto inv = firefly::FFInt(1)/red[rnk][c];
            for (int j=c;j<nCols;++j) red[rnk][j]*=inv;
            for (int r=0;r<nRows;++r) {
                if (r==rnk||red[r][c].n==0) continue;
                auto f=red[r][c];
                for (int j=c;j<nCols;++j) red[r][j]-=f*red[rnk][j];
            }
            pivotCols.push_back(vorder[c]);
            ++rnk;
        }
        res.elimRank = rnk;

        int hrCnt=highRank.size(), hrElim=0;
        for (int vi : highRank)
            if (std::find(pivotCols.begin(),pivotCols.end(),vi)!=pivotCols.end()) ++hrElim;
        res.flags.TopFlag1 = hrElim>=hrCnt;
        res.flags.FullFlag1 = rnk>=nCols;
        for (int c=0;c<(int)pivotCols.size();++c) res.reducedVarTable.push_back(pivotCols[c]);

        // 4. 模块 GB: 用 ν-多项式系数构建 Singular 多分量模块（MMA L160-175）
        // 使用 gen(i) 模块分量（Singular 1-indexed），POT 排序 (lp(nCols), dp(ne))
        // MMA: SingularGroebnerBasis with "PositionOverTerm"->True
        // 只在提供了 entriesForModule 时启用
        bool useModuleGB = (entriesForModule != nullptr && seedsForModule != nullptr && sectorForModule != nullptr && !entriesForModule->empty());
        if (useModuleGB) {
            try {
                // 构建模块生成器: 每个 (entry, seed, relation_k) 一行
                // poly(ν)*gen(1) + poly(ν)*gen(2) + ... + poly(ν)*gen(nCols)
                // Singular gen(i) 是 1-indexed，gen(1) 对应第 0 列
                std::vector<std::string> moduleGens;

                for (auto& entry : *entriesForModule) {
                    if (entry.numRelations == 0 || entry.alphas.empty()) continue;
                    int nb = (int)entry.betas.size();
                    int na = (int)entry.alphas.size();

                    for (auto& seed : *seedsForModule) {
                        for (int k = 0; k < entry.numRelations; ++k) {
                            // 初始化 nCols 个 ν-多项式的字符串数组
                            std::vector<std::string> polyByCol(nCols, "0");

                            // 对每个 α_i，找到对应列，填入 ν-多项式
                            for (int ai = 0; ai < na; ++ai) {
                                // 计算列 key = sd + α_ai
                                std::vector<int> colKey(ne);
                                for (int j = 0; j < ne; ++j)
                                    colKey[j] = (*sectorForModule)[j] ? (-seed[j] + entry.alphas[ai][j]) : (seed[j] + entry.alphas[ai][j]);

                                // sectorFilter: skip if pdset component < 0
                                bool skip = false;
                                for (int j = 0; j < ne; ++j)
                                    if ((*sectorForModule)[j] && colKey[j] < 0) { skip = true; break; }
                                if (skip) continue;

                                // 在 bm.alphaShifts 中找此列
                                auto it = std::find(bm.alphaShifts.begin(), bm.alphaShifts.end(), colKey);
                                if (it == bm.alphaShifts.end()) continue;
                                int col = (int)(it - bm.alphaShifts.begin());

                                // 在 vorder 中找此列位置
                                int vpos = -1;
                                for (int vi = 0; vi < nCols; ++vi)
                                    if (vorder[vi] == col) { vpos = vi; break; }
                                if (vpos < 0) continue;

                                // 构建 ν-多项式系数: Σ_j b0[α_i, β_j] · (ν-sd)^{β_j}
                                std::vector<firefly::FFInt> b0Row;
                                for (int bj = 0; bj < nb; ++bj) {
                                    int ri = ai * nb + bj;
                                    if (ri >= (int)entry.coefficients.size() || k + 1 >= (int)entry.coefficients[ri].size())
                                        b0Row.push_back(firefly::FFInt(0));
                                    else
                                        b0Row.push_back(entry.coefficients[ri][k + 1]);
                                }

                                std::string nuPoly = buildNuPolynomial(b0Row, entry.betas, seed,
                                    *sectorForModule, ne, modulus);

                                // 如果 numericReduce 提供，将 ν 替换为常数（B2: ISP→0, B3: ISP→0+prop→1）
                                if (numericReduce) {
                                    for (size_t ni = 0; ni < (size_t)ne && ni < numericReduce->size(); ++ni) {
                                        long long val = (*numericReduce)[ni];
                                        if (val < 0) continue;  // -1 = no substitution
                                        std::string target = "nu" + std::to_string(ni);
                                        size_t p;
                                        while ((p = nuPoly.find(target)) != std::string::npos) {
                                            size_t end = p + target.size();
                                            int exp = 1;
                                            if (end < nuPoly.size() && nuPoly[end] == '^') {
                                                exp = std::stoi(nuPoly.substr(end + 1));
                                                end = nuPoly.find_first_not_of("0123456789", end + 1);
                                            }
                                            long long powVal = 1;
                                            for (int e = 0; e < exp; ++e) powVal = (powVal * val) % modulus;
                                            nuPoly.replace(p, end - p, std::to_string(powVal));
                                        }
                                    }
                                }

                                polyByCol[vpos] = nuPoly;
                            }

                            // 构建 Singular 模块生成器: (coeff)*gen(1) + (coeff)*gen(2) + ...
                            // Singular gen(i) 是 1-indexed 模块分量
                            std::string modElem;
                            for (int ci = 0; ci < nCols; ++ci) {
                                if (polyByCol[ci] == "0") continue;
                                if (!modElem.empty()) modElem += "+";
                                modElem += "(" + polyByCol[ci] + ")*gen(" + std::to_string(ci + 1) + ")";
                            }
                            if (!modElem.empty()) moduleGens.push_back(modElem);
                        }
                    }
                }

                // 调用 Singular 模块 GB（POT 排序）
                // Ring 仅含 ν 变量，模块分量通过 gen(i) 提供，排序为 (lp(nCols), dp(ne))
                if (!moduleGens.empty()) {
                    std::vector<std::string> nuVars;
                    for (int i = 0; i < ne; ++i) nuVars.push_back("nu" + std::to_string(i));

                    auto gbList = SingularRunner::modulePOTGroebner(moduleGens, nCols, nuVars, modulus);

                    // 解析模块 GB: 每个元素为 coeff_i(ν)*gen(pos_i) + coeff_j(ν)*gen(pos_j) + ...
                    // 按 firstNonZero gen(i) 分组提取对角元（MMA L166-167）
                    // Singular gen() 是 1-indexed；转换为 0-indexed
                    std::map<int, std::string> diagByPos;

                    for (auto& element : gbList) {
                        // 按 gen(N) 拆分，第一个出现的 gen(N) 即 firstNonZero 分量
                        // 解析: "(coeff)*gen(N)" 格式，可能前有符号
                        struct GenTerm { int idx; std::string coeff; };
                        std::vector<GenTerm> gterms;

                        size_t pos = 0;
                        while (pos < element.size()) {
                            // 跳过空格
                            while (pos < element.size() && element[pos] == ' ') ++pos;
                            if (pos >= element.size()) break;

                            // 判断符号
                            int sgn = 1;
                            if (element[pos] == '+') { sgn = 1; ++pos; }
                            else if (element[pos] == '-') { sgn = -1; ++pos; }
                            while (pos < element.size() && element[pos] == ' ') ++pos;

                            // 读取系数（在括号中或直接是数字/1）
                            std::string coeff;
                            if (pos < element.size() && element[pos] == '(') {
                                ++pos;
                                int depth = 1;
                                while (pos < element.size() && depth > 0) {
                                    if (element[pos] == '(') ++depth;
                                    else if (element[pos] == ')') --depth;
                                    if (depth > 0) coeff += element[pos];
                                    ++pos;
                                }
                                // 跳过 "*gen"
                                while (pos < element.size() && element[pos] != 'g') ++pos;
                            } else {
                                // 无括号: 系数是 "1" 或数字
                                while (pos < element.size() && element[pos] != 'g' && element[pos] != '+' && element[pos] != '-') {
                                    coeff += element[pos];
                                    ++pos;
                                }
                                if (coeff.empty() || coeff == "*") coeff = "1";
                            }

                            // 跳过 "gen("
                            while (pos < element.size() && element[pos] != '(') ++pos;
                            if (pos >= element.size()) break;
                            ++pos; // skip '('

                            // 读取 gen 索引
                            std::string idxStr;
                            while (pos < element.size() && element[pos] >= '0' && element[pos] <= '9') {
                                idxStr += element[pos];
                                ++pos;
                            }
                            // skip ')'
                            if (pos < element.size() && element[pos] == ')') ++pos;

                            if (!idxStr.empty()) {
                                int gi = std::stoi(idxStr) - 1;  // 1-indexed → 0-indexed
                                std::string fullCoeff = (sgn == -1) ? "-" + coeff : coeff;
                                if (fullCoeff.empty()) fullCoeff = "1";
                                gterms.push_back({gi, fullCoeff});
                            }
                        }

                        if (gterms.empty()) continue;

                        // firstNonZero = smallest gen-index (POT 排序下首个分量)
                        int firstPos = gterms[0].idx;
                        std::string firstCoeff = gterms[0].coeff;

                        // target: 所有 gen-index > firstPos 的分量
                        std::vector<std::pair<int, std::string>> targets;
                        for (size_t t = 1; t < gterms.size(); ++t) {
                            if (gterms[t].idx > firstPos && gterms[t].idx < res.blockEnd)
                                targets.push_back({gterms[t].idx, gterms[t].coeff});
                        }

                        if (firstPos >= 0 && firstPos < res.blockEnd) {
                            diagByPos[firstPos] = firstCoeff;
                            res.gbDiagonal.push_back("gen("+std::to_string(firstPos)+"): coeff="+firstCoeff);

                            ModuleReductionRule mr;
                            mr.srcCol = firstPos;
                            mr.srcCoeffStr = firstCoeff;
                            mr.targets = targets;
                            res.moduleRules.push_back(mr);

                            // ISP singularity detection: 检查对角系数的分母依赖
                            bool hasNu = false;
                            std::vector<int> nuV;
                            for (int j = 0; j < ne; ++j)
                                if (firstCoeff.find("nu" + std::to_string(j)) != std::string::npos) { hasNu = true; nuV.push_back(j); }
                            if (hasNu) res.singularISP.push_back(nuV);
                        }
                    }

                    // Flags (MMA L188-201)
                    res.flags.TopFlag2 = false;
                    res.flags.FullFlag2 = false;
                    res.flags.GlobalFlag = false;

                    if (!res.gbDiagonal.empty()) {
                        bool anyNu = false, allUnit = true;
                        bool allHighRankCovered = true, allBlockCovered = true;

                        int nHighRank = (int)highRank.size();
                        for (int vi = 0; vi < nHighRank; ++vi)
                            if (diagByPos.find(vi) == diagByPos.end()) { allHighRankCovered = false; break; }
                        for (int vi = 0; vi < res.blockEnd; ++vi)
                            if (diagByPos.find(vi) == diagByPos.end()) { allBlockCovered = false; break; }

                        for (auto& [pos_, coeff] : diagByPos) {
                            bool hn = false;
                            for (int j = 0; j < ne; ++j)
                                if (coeff.find("nu" + std::to_string(j)) != std::string::npos) { hn = true; break; }
                            if (hn) anyNu = true;
                            if (coeff != "1" && coeff != "-1") allUnit = false;
                        }

                        if (allHighRankCovered && !anyNu) res.flags.TopFlag2 = true;
                        if (allBlockCovered && !anyNu) res.flags.FullFlag2 = true;
                        if (allBlockCovered && allUnit) res.flags.GlobalFlag = true;
                    }
                }
            } catch (std::exception& e) {
                std::cerr << "[eliminatedRule] Module POT GB failed: " << e.what() << "\n";
            }
        }

        // Fallback: RREF-based concrete rules (constant coefficient case)
        bool hasRealTargets = false;
        for (auto& mr : res.moduleRules)
            if (!mr.targets.empty()) { hasRealTargets = true; break; }

        if (!hasRealTargets && res.elimRank > 0) {
            for (int r = 0; r < rnk; ++r) {
                int pivotCol = -1;
                for (int c = 0; c < nCols; ++c)
                    if (red[r][c].n != 0) { pivotCol = c; break; }
                if (pivotCol < 0) continue;

                int origSrc = vorder[pivotCol];
                ModuleReductionRule mr;
                mr.srcCol = origSrc; // original alphaShifts index
                mr.srcCoeffStr = "1";
                for (int c = pivotCol + 1; c < nCols; ++c) {
                    if (red[r][c].n == 0) continue;
                    mr.targets.push_back({vorder[c], std::to_string((-red[r][c]).n)});
                }
                res.moduleRules.push_back(mr);
                res.gbDiagonal.push_back("gen("+std::to_string(origSrc)+"): coeff=1");
            }
        }
    }

    return res;
}

inline ConeRule sectorBlockReducible(
    const std::vector<RelationEntry>& entries,
    const std::vector<int>& sector,
    int level, int64_t modulus, int ne,
    bool topRankOnly=false, int sectorUp=0,
    const std::vector<int>* numericReduce=nullptr)
{
    ConeRule cr;
    cr.level = level;

    auto seeds = topRankOnly
        ? generateDotSeeds(ne, level, true, level)
        : generateDotSeeds(ne, level+sectorUp, false, 0);
    auto mSeeds = generateMinusSeeds(ne, sectorUp, sector);
    for (auto& s : mSeeds) seeds.push_back(s);

    auto bm = buildBlockMatrix(entries, seeds, sector, ne, modulus);
    auto elim = eliminatedRule(bm, level, modulus, ne, topRankOnly, numericReduce,
                               &entries, &seeds, &sector);

    cr.flags = elim.flags;
    cr.gbDiagonal = elim.gbDiagonal;
    cr.singularISP = elim.singularISP;
    cr.alphaShifts = bm.alphaShifts;
    for (int vi : elim.reducedVarTable) cr.reducedVarIdx.push_back(vi);
    for (auto& mr : elim.moduleRules) cr.moduleRules.push_back(mr);

    // Build singular locus description
    if (!elim.singularISP.empty()) {
        std::string loc;
        for (auto& isp : elim.singularISP) {
            if (!loc.empty()) loc+=" ∪ ";
            loc += "denom[ν";
            for (size_t i=0;i<isp.size();++i) { if(i) loc+=","; loc+=std::to_string(isp[i]); }
            loc += "]=0";
        }
        cr.singularLocus = loc;
    }

    return cr;
}

// ═════════════════════════════════════════════════════════════
// setupGeneratingCone  (MMA SymbolicBTForm.wl:281-365)
// ═════════════════════════════════════════════════════════════

inline GeneratingConeResult setupGeneratingCone(
    const RelationData& rd,
    const std::vector<int>& sector,
    int levelBound=10)
{
    GeneratingConeResult res;
    res.ne=rd.ne; res.levelBound=levelBound;
    int ne=rd.ne; int64_t mod=rd.modulus;
    int topRank=std::accumulate(sector.begin(),sector.end(),0);
    res.topRank=topRank;

    std::vector<RelationEntry> pent;
    for (auto& e : rd.entries) if (e.ansatzMode=="Pyramid") pent.push_back(e);
    if (pent.empty()) { std::cerr<<"[setupGeneratingCone] No Pyramid entries\n"; return res; }

    // MMA: rgen = Length@rel — number of relation levels
    // B1 starts at rgen-1 (not topRank-1)
    int rgen = 0;
    for (auto& e : pent) rgen = std::max(rgen, (int)e.alphas.size());
    if (rgen == 0) rgen = topRank;

    // ── B1: High-rank cone (MMA L288-300) ──
    {
        bool done=false;
        for (int lvl = std::max(topRank, rgen-1); lvl < levelBound && !done; ++lvl) {
            auto cr = sectorBlockReducible(pent, sector, lvl, mod, ne, true, 0);
            res.B1.push_back(cr);
            done = cr.flags.TopFlag2;
        }
    }

    // ── B2: Mid-rank cone (MMA L312-328) ──
    // B2 uses NumericReduce: ISP vars → 0 (MMA L316)
    // Convergence: Most[B2mistab] === B2mistab0 (L324)
    {
        std::vector<int> nr2(ne, -1); // -1 = no substitution
        for (size_t i = 0; i < (size_t)ne; ++i)
            if (i < sector.size() && sector[i] == 0) nr2[i] = 0; // ISP→0

        int b1lvl = res.B1.empty() ? topRank : res.B1.back().level;

        // Count ISP-only dimensions for mistab generation
        int ne_isp = 0;
        for (size_t i = 0; i < (size_t)ne; ++i)
            if (sector[i] == 0) ++ne_isp;

        std::vector<std::vector<int>> prevMistabB2;
        bool done = false;
        for (int lvl = std::max(1, b1lvl - 1); lvl < levelBound && !done; ++lvl) {
            auto cr = sectorBlockReducible(pent, sector, lvl, mod, ne, false, 0, &nr2);
            res.B2.push_back(cr);

            // Build B2 mistab: per-ISP-rank list of non-reduced integrals (MMA L321-323)
            std::vector<std::vector<int>> curMistab;
            for (int r = 0; r <= lvl; ++r) {
                std::vector<int> nonRed;
                auto rankSeeds = generateDotSeeds(ne_isp, r);
                for (auto& seed : rankSeeds) {
                    std::vector<int> full(ne, 0);
                    int idx = 0;
                    for (int i = 0; i < ne; ++i) {
                        if (sector[i] == 0 && idx < (int)seed.size()) {
                            full[i] = seed[idx];
                            ++idx;
                        }
                    }
                    // Check if integral is in reducedVarIdx
                    auto it = std::find(cr.alphaShifts.begin(), cr.alphaShifts.end(), full);
                    int pos = -1;
                    if (it != cr.alphaShifts.end()) pos = (int)(it - cr.alphaShifts.begin());
                    if (pos >= 0 && std::find(cr.reducedVarIdx.begin(), cr.reducedVarIdx.end(), pos) == cr.reducedVarIdx.end())
                        nonRed.push_back(pos);
                }
                curMistab.push_back(nonRed);
            }

            // MMA convergence: Most[B2mistab] === B2mistab0 (L324)
            // Most drops the last (highest ISP rank) level
            if (!prevMistabB2.empty() && curMistab.size() > 1 && prevMistabB2.size() > 1) {
                bool same = true;
                size_t cmpLen = std::min(curMistab.size() - 1, prevMistabB2.size() - 1);
                for (size_t i = 0; i < cmpLen; ++i)
                    if (curMistab[i] != prevMistabB2[i]) { same = false; break; }
                if (same) done = true;
            }
            if (!done) prevMistabB2 = curMistab;
        }
    }

    // ── B3: Near-corner cone (MMA L336-360) ──
    {
        // B3 NumericReduce: ISP→0, propagator vars→1 (MMA L342)
        std::vector<int> nr3(ne, -1);
        for (size_t i = 0; i < (size_t)ne; ++i) {
            if (i < sector.size()) {
                if (sector[i] == 0) nr3[i] = 0;       // ISP→0
                else nr3[i] = 1;                       // propagator→1
            }
        }

        auto ispRank = [&](const std::vector<int>& shift) -> int {
            int r = 0;
            for (size_t i = 0; i < shift.size() && i < sector.size(); ++i)
                if (sector[i] == 0) r += shift[i];
            return r;
        };

        int ne_tot = ne;
        for (size_t i = 0; i < (size_t)ne; ++i)
            if (sector[i] == 1) ne_tot--; // count ISP-only dimensions

        std::vector<std::vector<int>> prevMistab;
        bool done = false;
        int sup = 1, maxSup = 3;
        int maxIter = 12;

        for (int lvl = std::max(topRank, rgen-1), iter = 0; lvl < levelBound && !done && iter < maxIter; ++lvl, ++iter) {
            auto cr = sectorBlockReducible(pent, sector, lvl, mod, ne, false, sup, &nr3);
            res.B3.push_back(cr);

            // MMA B3 mistab: non-reduced integrals per ISP rank (MMA L344)
            std::vector<std::vector<int>> curMistab;
            for (int r = 0; r <= lvl; ++r) {
                std::vector<int> nonRed;
                // Generate all integrals of ISP rank r
                auto rankSeeds = generateDotSeeds(ne_tot, r);
                for (auto& seed : rankSeeds) {
                    std::vector<int> full(ne, 0);
                    int idx = 0;
                    for (int i = 0; i < ne; ++i) {
                        if (sector[i] == 0) {
                            full[i] = seed[idx];  // ISP index from seed
                            idx++;
                        }
                    }
                    // Check if this integral is reduced
                    // Find it in alphaShifts, then check if it's in reducedVarIdx
                    auto it = std::find(cr.alphaShifts.begin(), cr.alphaShifts.end(), full);
                    int pos = -1;
                    if (it != cr.alphaShifts.end()) pos = (int)(it - cr.alphaShifts.begin());
                    if (pos >= 0 && std::find(cr.reducedVarIdx.begin(), cr.reducedVarIdx.end(), pos) == cr.reducedVarIdx.end())
                        nonRed.push_back(pos);
                }
                curMistab.push_back(nonRed);
            }

            // MMA convergence: Most[B3mistab] === B3mistab0 (MMA L346)
            // Most drops the last (highest ISP rank) level — it naturally changes as level increases
            if (!prevMistab.empty() && curMistab.size() > 1 && prevMistab.size() > 1) {
                bool same = true;
                size_t cmpLen = std::min(curMistab.size() - 1, prevMistab.size() - 1);
                for (size_t i = 0; i < cmpLen; ++i)
                    if (curMistab[i] != prevMistab[i]) { same = false; break; }
                if (same) done = true;
            }
            if (!done) prevMistab = curMistab;

            if (iter >= maxIter - 2 && !done && sup < maxSup) ++sup;
        }
    }

    // Coverage checks
    res.fullCoverage = true;
    for (auto& b1 : res.B1) if (!b1.flags.TopFlag1) res.fullCoverage=false;
    if (res.B2.empty() || !res.B2.back().flags.FullFlag2) res.fullCoverage=false;

    return res;
}

// ═════════════════════════════════════════════════════════════
// P6: 符号规则输出 (from eliminatedRule GB diagonal + lift matrix)
// ═════════════════════════════════════════════════════════════

inline std::string SymbolicRule::toString(int ne) const {
    std::string s = "g(";
    for (int i=0;i<ne;++i) {
        if(i)s+=",";
        if(sourceShift[i]>=0) s+="ν"+std::to_string(i)+"-"+std::to_string(sourceShift[i]);
        else s+="ν"+std::to_string(i)+"+"+std::to_string(-sourceShift[i]);
    }
    if (!targetShifts.empty()) {
        s+=") = ";
        for (size_t ti=0;ti<targetShifts.size();++ti) {
            if(ti) s+=" + ";
            if(ti<coeffStrs.size()) {
                s += coeffStrs[ti];
            } else {
                s += "?";
            }
            s+="*g(";
            for(int i=0;i<ne;++i) {
                if(i)s+=",";
                int v=targetShifts[ti][i];
                if(v>=0) s+="ν"+std::to_string(i)+"-"+std::to_string(v);
                else s+="ν"+std::to_string(i)+"+"+std::to_string(-v);
            }
            s+=")";
        }
    }
    if (hasISPSingularity) {
        if (targetShifts.empty()) s += ")";
        s += "  ⚠ISP：ν";
        for (size_t i=0;i<ispVars.size();++i) { if(i)s+=","; s+=std::to_string(ispVars[i]); }
        s += " 在分母";
    }
    return s;
}

inline std::vector<SymbolicRule> generateSymbolicRules(
    const GeneratingConeResult& gcr,
    const std::vector<RelationEntry>& entries,
    int ne,
    const std::vector<int>& sector)
{
    std::vector<SymbolicRule> rules;
    auto allCones = gcr.B1;
    allCones.insert(allCones.end(), gcr.B2.begin(), gcr.B2.end());
    allCones.insert(allCones.end(), gcr.B3.begin(), gcr.B3.end());

    for (auto& cone : allCones) {
        for (auto& mr : cone.moduleRules) {
            SymbolicRule sr;

            // Source shift
            if (mr.srcCol >= (int)cone.alphaShifts.size()) continue;
            sr.sourceShift = cone.alphaShifts[mr.srcCol];

            // Parse source coefficient for ISP detection
            sr.hasISPSingularity = false;
            for (int j = 0; j < ne; ++j)
                if (mr.srcCoeffStr.find("nu" + std::to_string(j)) != std::string::npos) {
                    sr.hasISPSingularity = true;
                    sr.ispVars.push_back(j);
                }

            // Target shifts: gen(src) + Σ c_j * gen(tgt_j) = 0
            // → g(ν-α_src) = Σ (-c_j(ν)/c_src(ν)) * g(ν-α_tgt_j)
            std::string srcCoeff = mr.srcCoeffStr;
            bool srcIsUnit = (srcCoeff == "1" || srcCoeff == "-1");

            for (auto& [tgtCol, coeffStr] : mr.targets) {
                if (tgtCol >= (int)cone.alphaShifts.size()) continue;
                sr.targetShifts.push_back(cone.alphaShifts[tgtCol]);

                // 构建约化系数: -c_j(ν)/c_src(ν)
                // 当 c_src=1 时就是 -c_j(ν)
                std::string ruleCoeff;
                if (srcIsUnit) {
                    // 检测 coeffStr 是否含 ν（即多项式，非纯常数）
                    bool hasNu = (coeffStr.find("nu") != std::string::npos);
                    if (hasNu) {
                        // ν-多项式：用括号包装
                        ruleCoeff = "-(" + coeffStr + ")";
                    } else if (coeffStr == "1") {
                        ruleCoeff = "-1";
                    } else if (coeffStr == "-1") {
                        ruleCoeff = "1";
                    } else if (!coeffStr.empty() && coeffStr[0] == '-') {
                        // "-N" → "N" (取反)
                        ruleCoeff = coeffStr.substr(1);
                    } else {
                        ruleCoeff = "-" + coeffStr;
                    }
                } else {
                    ruleCoeff = "-(" + coeffStr + ")/(" + srcCoeff + ")";
                }
                sr.coeffStrs.push_back(ruleCoeff);
            }

            // Include rule even if empty targets (reduces to 0)
            rules.push_back(sr);
        }
    }
    return rules;
}

// ═════════════════════════════════════════════════════════════
// reduceTargetsByCone  (MMA SymbolicBTForm.wl:383-446)
// 对目标积分列表逐项尝试 B1→B2→B3 规则进行约化
// nuBase: propagator vars → 1, ISP vars → 0 (MMA L391-392)
// ═════════════════════════════════════════════════════════════

struct ConeReductionResult {
    std::vector<std::vector<int>> sourceShifts;      // 每个目标的原始 shift
    std::vector<std::string> coneLabels;              // "B1"/"B2"/"B3" 或 "NONE"
    std::vector<std::string> reducedExprs;            // 约化结果字符串
    int totalTargets = 0, reducibleTargets = 0;
};

// Evaluate a ν-polynomial string at nuBase values
// This is a simplified numeric evaluator: substitutes nu0, nu1, ... with their values
// and computes the result in the finite field
inline int64_t evaluateNuPoly(const std::string& poly, const std::vector<int>& nuVals, int64_t modulus) {
    // Simple approach: replace nu<N> with its value, then evaluate via Singular
    // For pure constants and simple expressions, do inline evaluation
    // For complex polynomials, we return 0 as a fallback (the target is not checked via this path)

    // Count occurrences of nu variables
    bool hasAnyNu = false;
    for (size_t i = 0; i < nuVals.size(); ++i)
        if (poly.find("nu" + std::to_string(i)) != std::string::npos) { hasAnyNu = true; break; }

    if (!hasAnyNu) {
        // Pure constant — parse directly
        try {
            long long val = std::stoll(poly);
            val = val % modulus;
            if (val < 0) val += modulus;
            return val;
        } catch (...) { return 0; }
    }

    // For ν-dependent polynomials, substitute ν values and evaluate
    std::string evald = poly;
    for (size_t i = 0; i < nuVals.size(); ++i) {
        std::string target = "nu" + std::to_string(i);
        std::string repl = std::to_string(nuVals[i]);
        size_t p = 0;
        while ((p = evald.find(target, p)) != std::string::npos) {
            // Check if followed by exponent: nu<N>^exp
            size_t end = p + target.size();
            int exp = 1;
            if (end < evald.size() && evald[end] == '^') {
                size_t expStart = end + 1;
                size_t expEnd = expStart;
                while (expEnd < evald.size() && evald[expEnd] >= '0' && evald[expEnd] <= '9') ++expEnd;
                exp = std::stoi(evald.substr(expStart, expEnd - expStart));
                end = expEnd;
            }
            // Compute nuVal^exp
            long long powVal = 1;
            for (int e = 0; e < exp; ++e) powVal = (powVal * nuVals[i]) % modulus;
            evald.replace(p, end - p, std::to_string(powVal));
            p += std::to_string(powVal).size();
        }
    }

    // Now evaluate the constant expression: sum of terms like "a*b*c"
    // Simple evaluator for product-of-constants
    try {
        long long total = 0, cur = 0, prod = 1;
        bool inNum = false, inProd = false;
        char lastOp = '+';
        std::string numBuf;

        for (size_t i = 0; i <= evald.size(); ++i) {
            char c = (i < evald.size()) ? evald[i] : '\0';
            if (c == ' ' || c == '*') continue;
            if (c == '+' || c == '-' || c == '\0') {
                if (!numBuf.empty()) {
                    long long v = std::stoll(numBuf) % modulus;
                    if (v < 0) v += modulus;
                    prod = (prod * v) % modulus;
                    numBuf.clear();
                }
                if (lastOp == '+') total = (total + prod) % modulus;
                else if (lastOp == '-') total = (total - prod + modulus) % modulus;
                prod = 1;
                lastOp = (c == '-') ? '-' : '+';
            } else if (c >= '0' && c <= '9') {
                numBuf += c;
            }
        }
        return total;
    } catch (...) { return 0; }
}

inline ConeReductionResult reduceTargetsByCone(
    const GeneratingConeResult& gcr,
    const std::vector<std::vector<int>>& targetShifts,
    const std::vector<int>& sector,
    int ne)
{
    ConeReductionResult result;
    result.totalTargets = (int)targetShifts.size();
    result.sourceShifts = targetShifts;
    result.coneLabels.resize(targetShifts.size(), "NONE");
    result.reducedExprs.resize(targetShifts.size());

    // Build rule lookup: shift → {coneLabel, moduleRules}
    // For each cone in B1, B2, B3, map alphaShifts entries to their module rules
    struct RuleEntry {
        std::string coneName;
        const ModuleReductionRule* rule;
        const std::vector<std::vector<int>>* shifts;
    };
    std::map<std::vector<int>, RuleEntry> ruleMap;

    auto addConeRules = [&](const std::vector<ConeRule>& cones, const std::string& name) {
        for (auto& cone : cones) {
            for (auto& mr : cone.moduleRules) {
                if (mr.srcCol >= (int)cone.alphaShifts.size()) continue;
                auto shift = cone.alphaShifts[mr.srcCol];
                // Don't overwrite B1 with B2/B3 (B1 has highest priority)
                if (ruleMap.find(shift) == ruleMap.end())
                    ruleMap[shift] = {name, &mr, &cone.alphaShifts};
            }
        }
    };
    addConeRules(gcr.B1, "B1");
    addConeRules(gcr.B2, "B2");
    addConeRules(gcr.B3, "B3");

    // nuBase: propagator → 1, ISP → 0 (MMA L391-392)
    std::vector<int> nuBase(ne, 0);
    for (int i = 0; i < ne; ++i)
        if (i < (int)sector.size() && sector[i] == 1) nuBase[i] = 1;

    // For each target, try to find and apply a reduction rule
    for (size_t ti = 0; ti < targetShifts.size(); ++ti) {
        auto& ts = targetShifts[ti];

        // The g-key for rule lookup: g(ν_base - target) → rule maps g(shift)
        // In MMA: sh = g@@(nuBase - List@@jx) where jx is the target integral
        // j[x] = g(ν_base - x), so g_key = ν_base - target
        std::vector<int> gKey(ne);
        for (int i = 0; i < ne; ++i)
            gKey[i] = nuBase[i] - ts[i];

        auto it = ruleMap.find(gKey);
        if (it != ruleMap.end()) {
            result.coneLabels[ti] = it->second.coneName;
            ++result.reducibleTargets;

            // Build reduction expression: source = Σ coeff * target
            auto& mr = *it->second.rule;
            std::ostringstream expr;
            expr << "g(v-[" << gKey[0];
            for (int i = 1; i < ne; ++i) expr << "," << gKey[i];
            expr << "]) =";
            for (size_t j = 0; j < mr.targets.size(); ++j) {
                if (j > 0) expr << " + ";
                auto& tgt = mr.targets[j];
                expr << "(" << tgt.second << ")*g(v-[";
                for (int i = 0; i < ne; ++i) {
                    if (i) expr << ",";
                    expr << (*it->second.shifts)[tgt.first][i];
                }
                expr << "])";
            }
            result.reducedExprs[ti] = expr.str();
        }
    }

    return result;
}

// ═════════════════════════════════════════════════════════════
// 验证管线
// ═════════════════════════════════════════════════════════════

struct VerificationReport {
    std::string family;
    int totalTargets=0, reducibleTargets=0;
    int symbolicRules=0, concreteRules=0;
    bool allConesPass=false;
    std::vector<std::string> warnings;
    void print(std::ostream& os) const {
        os << "=== Verification Report ===\nFamily: "<<family<<"\nTotal targets: "<<totalTargets
           <<"\nReducible targets: "<<reducibleTargets<<"/"<<totalTargets
           <<"\nSymbolic rules: "<<symbolicRules<<"\nConcrete rules: "<<concreteRules
           <<"\nAll cones pass: "<<(allConesPass?"YES":"NO")<<"\n";
        for (auto& w : warnings) os<<"Warning: "<<w<<"\n";
    }
};

inline VerificationReport runVerification(
    const RelationData& rd,
    const std::vector<int>& sector,
    int levelBound=10)
{
    VerificationReport rpt;
    rpt.family=rd.family;
    auto gcr = setupGeneratingCone(rd, sector, levelBound);
    std::vector<RelationEntry> pent;
    for (auto& e : rd.entries) if (e.ansatzMode=="Pyramid") pent.push_back(e);

    auto srs = generateSymbolicRules(gcr, pent, rd.ne, sector);
    rpt.symbolicRules = srs.size();

    // Build target list and check reducibility via reduceTargetsByCone
    std::vector<std::vector<int>> allTargets;
    for (int l = 0; l <= levelBound; ++l) {
        auto ss = generateDotSeeds(rd.ne, l);
        for (auto& s : ss) allTargets.push_back(s);
    }
    rpt.totalTargets = (int)allTargets.size();
    auto cr = reduceTargetsByCone(gcr, allTargets, sector, rd.ne);
    rpt.reducibleTargets = cr.reducibleTargets;

    if (!gcr.fullCoverage) rpt.warnings.push_back("Not all targets covered");
    if (rpt.reducibleTargets < rpt.totalTargets)
        rpt.warnings.push_back(std::to_string(rpt.totalTargets - rpt.reducibleTargets) + " targets not reducible");

    rpt.allConesPass = gcr.fullCoverage;
    for (auto& b1 : gcr.B1) if (!b1.flags.TopFlag1) rpt.allConesPass=false;
    for (auto& b2 : gcr.B2) if (!b2.flags.FullFlag2) rpt.allConesPass=false;

    return rpt;
}

#endif
