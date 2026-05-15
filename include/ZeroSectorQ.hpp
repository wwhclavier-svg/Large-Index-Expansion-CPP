#ifndef ZERO_SECTOR_Q_HPP
#define ZERO_SECTOR_Q_HPP

// ZeroSectorQ — C++ implementation of Mathematica ZeroSectorQ algorithm
// Determines whether a Feynman integral sector is zero by the criterion
//   Σ k_i · x_i · ∂G/∂x_i = G   where G = U + F (Symanzik polynomials)
//
// Translated from ZeroSectorQ_v3.wl
//
// Convention: returns true  (= ZERO sector)  when the criterion equation has
//                        a solution (system is consistent).
//                        returns false (= NON-ZERO sector) when inconsistent.
//             This matches Mathematica's ZeroSectorQ[] return convention.
//
// Usage:
//   #include "ZeroSectorQ.hpp"
//   ...
//   bool isZero = zeroSectorQ(L, E, N, L0, L1, E0, E1, masses, kinRules);

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace ZeroSectorQ {

// ------------------------------------------------------------
// Rational number
// ------------------------------------------------------------
struct Rat {
  int64_t num, den;
  Rat(int64_t n = 0, int64_t d = 1) : num(n), den(d) {
    if (den < 0) { num = -num; den = -den; }
    int64_t g = std::gcd(std::abs(num), std::abs(den));
    if (g) { num /= g; den /= g; }
  }
  Rat operator+(const Rat &o) const { return {num * o.den + o.num * den, den * o.den}; }
  Rat operator-(const Rat &o) const { return {num * o.den - o.num * den, den * o.den}; }
  Rat operator*(const Rat &o) const { return {num * o.num, den * o.den}; }
  Rat operator/(const Rat &o) const { return {num * o.den, den * o.num}; }
  Rat operator-() const { return {-num, den}; }
  bool operator==(const Rat &o) const { return num * o.den == o.num * den; }
  bool operator!=(const Rat &o) const { return !(*this == o); }
  bool isZero() const { return num == 0; }
  std::string str() const {
    if (den == 1) return std::to_string(num);
    return std::to_string(num) + "/" + std::to_string(den);
  }
};

// ------------------------------------------------------------
// Multivariate polynomial in x[0]..x[n-1]
// ------------------------------------------------------------
using Exp = std::vector<int>;

struct ExpCmp {
  bool operator()(const Exp &a, const Exp &b) const {
    if (a.size() != b.size()) return a.size() < b.size();
    for (int i = (int)a.size() - 1; i >= 0; --i)
      if (a[i] != b[i]) return a[i] < b[i];
    return false;
  }
};

struct Poly {
  int nv;
  std::map<Exp, Rat, ExpCmp> terms;

  Poly(int n = 0) : nv(n) {}
  Poly(int n, const Exp &e, Rat c) : nv(n) { if (!c.isZero()) terms[e] = c; }

  bool isZero() const { return terms.empty(); }

  Poly operator+(const Poly &o) const {
    assert(nv == o.nv); Poly r(nv);
    for (auto &kv : terms)   r.terms[kv.first] = r.terms[kv.first] + kv.second;
    for (auto &kv : o.terms) r.terms[kv.first] = r.terms[kv.first] + kv.second;
    for (auto it = r.terms.begin(); it != r.terms.end(); )
      if (it->second.isZero()) it = r.terms.erase(it); else ++it;
    return r;
  }
  Poly operator-(const Poly &o) const {
    assert(nv == o.nv); Poly r(nv);
    for (auto &kv : terms)   r.terms[kv.first] = r.terms[kv.first] + kv.second;
    for (auto &kv : o.terms) r.terms[kv.first] = r.terms[kv.first] - kv.second;
    for (auto it = r.terms.begin(); it != r.terms.end(); )
      if (it->second.isZero()) it = r.terms.erase(it); else ++it;
    return r;
  }
  Poly operator-() const { Poly r(nv); for (auto &kv : terms) r.terms[kv.first] = -kv.second; return r; }
  Poly operator*(const Poly &o) const {
    assert(nv == o.nv); Poly r(nv);
    for (auto &ka : terms)
      for (auto &kb : o.terms) {
        Exp e = ka.first;
        if (kb.first.size() > e.size()) e.resize(kb.first.size(), 0);
        for (size_t i = 0; i < kb.first.size(); ++i) e[i] += kb.first[i];
        r.terms[e] = r.terms[e] + ka.second * kb.second;
      }
    for (auto it = r.terms.begin(); it != r.terms.end(); )
      if (it->second.isZero()) it = r.terms.erase(it); else ++it;
    return r;
  }
  Poly operator*(Rat c) const { Poly r(nv); for (auto &kv : terms) r.terms[kv.first] = kv.second * c; return r; }
  friend Poly operator*(Rat c, const Poly &p) { return p * c; }

  Poly deriv(int var) const {
    Poly r(nv);
    for (auto &kv : terms) {
      if ((int)kv.first.size() <= var || kv.first[var] == 0) continue;
      Exp e = kv.first;
      Rat c = kv.second * Rat(e[var]);
      e[var] -= 1;
      r.terms[e] = r.terms[e] + c;
    }
    for (auto it = r.terms.begin(); it != r.terms.end(); )
      if (it->second.isZero()) it = r.terms.erase(it); else ++it;
    return r;
  }
};

// ------------------------------------------------------------
// Polynomial matrix
// ------------------------------------------------------------
struct PolyMat {
  int rows, cols, nv;
  std::vector<Poly> data;
  PolyMat(int r, int c, int nv_) : rows(r), cols(c), nv(nv_), data(r * c, Poly(nv_)) {}
  Poly &at(int i, int j) { return data[i * cols + j]; }
  const Poly &at(int i, int j) const { return data[i * cols + j]; }
};

// determinant via Laplace expansion
static Poly detPoly(const PolyMat &m) {
  assert(m.rows == m.cols);
  int n = m.rows;
  if (n == 0) return Poly(0, Exp{}, Rat(1));
  if (n == 1) return m.at(0, 0);
  Poly r(m.nv);
  for (int j = 0; j < n; ++j) {
    PolyMat sub(n - 1, n - 1, m.nv);
    for (int si = 1; si < n; ++si)
      for (int sj = 0; sj < n; ++sj) {
        if (sj == j) continue;
        sub.at(si - 1, sj < j ? sj : sj - 1) = m.at(si, sj);
      }
    Poly term = m.at(0, j) * detPoly(sub);
    if (j % 2 == 0) r = r + term; else r = r - term;
  }
  return r;
}

// adjugate: adj(A)_{ji} = (-1)^{i+j} * det(minor removing row i, col j)
static PolyMat adjPoly(const PolyMat &m) {
  assert(m.rows == m.cols);
  int n = m.rows;
  PolyMat adj(n, n, m.nv);
  if (n == 1) { adj.at(0, 0) = Poly(m.nv, Exp{}, Rat(1)); return adj; }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      PolyMat sub(n - 1, n - 1, m.nv);
      for (int si = 0; si < n; ++si) {
        if (si == i) continue;
        for (int sj = 0; sj < n; ++sj) {
          if (sj == j) continue;
          sub.at(si < i ? si : si - 1, sj < j ? sj : sj - 1) = m.at(si, sj);
        }
      }
      Poly d = detPoly(sub);
      if ((i + j) % 2 == 1) d = Poly(m.nv) - d;
      adj.at(j, i) = d;
    }
  return adj;
}

// ------------------------------------------------------------
// Gaussian elimination
// ------------------------------------------------------------
static bool solveLinearSystem(const std::vector<std::vector<Rat>> &aug,
                               std::vector<Rat> &solution) {
  if (aug.empty()) { solution.clear(); return true; }
  int nEq = (int)aug.size();
  int nV = (int)aug[0].size() - 1;
  auto m = aug;
  int row = 0;
  for (int col = 0; col < nV && row < nEq; ++col) {
    int pivot = row;
    while (pivot < nEq && m[pivot][col].isZero()) ++pivot;
    if (pivot == nEq) continue;
    std::swap(m[row], m[pivot]);
    Rat pivVal = m[row][col];
    for (int j = col; j <= nV; ++j) m[row][j] = m[row][j] / pivVal;
    for (int i = 0; i < nEq; ++i) {
      if (i == row || m[i][col].isZero()) continue;
      Rat factor = m[i][col];
      for (int j = col; j <= nV; ++j)
        m[i][j] = m[i][j] - factor * m[row][j];
    }
    ++row;
  }
  for (int i = 0; i < nEq; ++i) {
    bool allZero = true;
    for (int j = 0; j < nV; ++j)
      if (!m[i][j].isZero()) { allZero = false; break; }
    if (allZero && !m[i][nV].isZero()) return false;
  }
  solution.assign(nV, Rat(0));
  for (int i = 0; i < nEq; ++i) {
    int pivotCol = -1;
    for (int j = 0; j < nV; ++j)
      if (!m[i][j].isZero()) { pivotCol = j; break; }
    if (pivotCol >= 0) solution[pivotCol] = m[i][nV];
  }
  return true;
}

// ------------------------------------------------------------
// ZeroSectorQ — main algorithm
//
// Input:
//   L  = number of loop momenta
//   E  = number of external momenta
//   N  = number of propagators in this sector
//   L0[L][N] = coeff of each loop momentum in mom1 of each propagator
//   L1[L][N] = coeff of each loop momentum in mom2 of each propagator
//   E0[E][N] = coeff of each external momentum in mom1
//   E1[E][N] = coeff of each external momentum in mom2
//   masses[N] = mass-squared values (Rat for exact rational)
//   kinRules[E*E] = dot products pExt[e1]·pExt[e2]
//   verbose = print debug info
//
// Output: true  = ZERO sector (equation system is consistent)
//         false = NON-ZERO sector (system inconsistent)
// ------------------------------------------------------------
inline bool zeroSectorQ(
    int L, int E, int N,
    const std::vector<std::vector<int>> &L0,
    const std::vector<std::vector<int>> &L1,
    const std::vector<std::vector<int>> &E0,
    const std::vector<std::vector<int>> &E1,
    const std::vector<Rat> &masses,
    const std::vector<Rat> &kinRules,
    bool verbose = false)
{
    if (N == 0) return true;

    // ---- Step 1: Build A = L0·diag(x)·L1^T (L×L) ----
    PolyMat A(L, L, N);
    for (int a = 0; a < L; ++a)
        for (int b = 0; b < L; ++b) {
            Poly sum(N);
            for (int p = 0; p < N; ++p) {
                int c0 = L0[a][p], c1 = L1[b][p];
                if (c0 == 0 || c1 == 0) continue;
                Exp e(N, 0); e[p] = 1;
                sum = sum + Poly(N, e, Rat(c0 * c1));
            }
            A.at(a, b) = sum;
        }

    // ---- Step 2: U = det(A) ----
    Poly U = detPoly(A);

    // ---- Step 3: Build F ----
    Poly F(N);
    Poly massSum(N);
    bool hasMass = false;
    for (int p = 0; p < N; ++p) {
        if (masses[p].isZero()) continue;
        hasMass = true;
        Exp ep(N, 0); ep[p] = 1;
        massSum = massSum + Poly(N, ep, masses[p]);
    }

    if (L > 0) {
        Poly detA = U;
        PolyMat adjA = adjPoly(A);

        if (hasMass) F = F - massSum * detA;

        if (E > 0) {
            // Braw[a][e] = Σ_p (E0[e][p]*L1[a][p] + E1[e][p]*L0[a][p]) * x[p]
            std::vector<std::vector<Poly>> Braw(L, std::vector<Poly>(E, Poly(N)));
            for (int a = 0; a < L; ++a)
                for (int e = 0; e < E; ++e)
                    for (int p = 0; p < N; ++p) {
                        int ce0 = E0[e][p], ce1 = E1[e][p];
                        int cl0 = L0[a][p], cl1 = L1[a][p];
                        Exp ep(N, 0); ep[p] = 1;
                        if (ce0 != 0 && cl1 != 0) Braw[a][e] = Braw[a][e] + Poly(N, ep, Rat(ce0 * cl1));
                        if (ce1 != 0 && cl0 != 0) Braw[a][e] = Braw[a][e] + Poly(N, ep, Rat(ce1 * cl0));
                    }

            // Craw[e1][e2] = Σ_p E0[e1][p] * x[p] * E1[e2][p]
            std::vector<std::vector<Poly>> Craw(E, std::vector<Poly>(E, Poly(N)));
            for (int e1 = 0; e1 < E; ++e1)
                for (int e2 = 0; e2 < E; ++e2)
                    for (int p = 0; p < N; ++p) {
                        int c0 = E0[e1][p], c1 = E1[e2][p];
                        if (c0 == 0 || c1 == 0) continue;
                        Exp ep(N, 0); ep[p] = 1;
                        Craw[e1][e2] = Craw[e1][e2] + Poly(N, ep, Rat(c0 * c1));
                    }

            // F += C_ext * detA
            for (int e1 = 0; e1 < E; ++e1)
                for (int e2 = 0; e2 < E; ++e2) {
                    if (Craw[e1][e2].isZero()) continue;
                    Rat kr = kinRules[e1 * E + e2];
                    if (kr.isZero() && e1 != e2) continue;
                    F = F + Craw[e1][e2] * detA * kr;
                }

            // F -= 0.25 * B · adj(A) · B
            for (int a1 = 0; a1 < L; ++a1)
                for (int a2 = 0; a2 < L; ++a2) {
                    Poly adjEl = adjA.at(a1, a2);
                    if (adjEl.isZero()) continue;
                    for (int e1 = 0; e1 < E; ++e1) {
                        if (Braw[a1][e1].isZero()) continue;
                        for (int e2 = 0; e2 < E; ++e2) {
                            if (Braw[a2][e2].isZero()) continue;
                            Rat kr = kinRules[e1 * E + e2];
                            if (kr.isZero()) continue;
                            F = F - Braw[a1][e1] * Braw[a2][e2] * adjEl * kr * Rat(1, 4);
                        }
                    }
                }
        }
    }

    Poly G = U + F;

    // ---- Step 4: Criterion Σ k[i]·x[i]·∂G/∂x[i] - G = 0 ----
    std::vector<Poly> xdG(N);
    for (int i = 0; i < N; ++i) {
        Poly d = G.deriv(i);
        for (auto &kv : d.terms) {
            Exp e = kv.first;
            if ((int)e.size() <= i) e.resize(i + 1, 0);
            e[i] += 1;
            xdG[i].terms[e] = xdG[i].terms[e] + kv.second;
        }
    }

    std::map<Exp, std::vector<Rat>, ExpCmp> monEqs;
    for (int i = 0; i < N; ++i)
        for (auto &kv : xdG[i].terms) {
            auto &vec = monEqs[kv.first];
            if (vec.empty()) vec.resize(N + 1, Rat(0));
            vec[i] = vec[i] + kv.second;
        }
    for (auto &kv : G.terms) {
        auto &vec = monEqs[kv.first];
        if (vec.empty()) vec.resize(N + 1, Rat(0));
        vec[N] = vec[N] - kv.second;
    }

    std::vector<std::vector<Rat>> aug;
    for (auto &kv : monEqs) {
        bool allZero = true;
        for (auto &v : kv.second) if (!v.isZero()) { allZero = false; break; }
        if (allZero) continue;
        aug.push_back(kv.second);
    }

    if (aug.empty()) return true;

    std::vector<Rat> solution;
    return solveLinearSystem(aug, solution);
}

} // namespace ZeroSectorQ

#endif // ZERO_SECTOR_Q_HPP
