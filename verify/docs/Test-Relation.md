# Relation Verification Manual

Unified verification for C++ `AllRelations_*.m` output. Single entry point, four methods.

---

## Entry Point

```bash
cd /home/ykm/Large-Index-Expansion-CPP
wolframscript -file VerifyRelation.wl <famname>                 # all modules
wolframscript -file VerifyRelation.wl <famname> skip M1,M2       # skip modules

# Individual module invocation (can also be loaded via Get[] from VerifyRelation):
cd verify/VerifyUtility
wolframscript -file Verify-Blade.wl    <famname> [nuSample]
wolframscript -file Verify-Series.wl   <famname>
wolframscript -file Verify-Kira.wl     <famname>
wolframscript -file Verify-MMACompare.wl <famname>
```

Each module auto-skips if its data files are missing. VerifyRelation runs all 4 in-process, sharing a single AllRelations load.

---

## Method 1: BladeVerify (Blade IBP)

**Purpose:** Verify C++ nullspace vectors vanish under Blade IBP reduction at a test nu point.
```
Sum_{alpha,beta} b_{alpha,beta} * nu^beta * BLReduce[ BL[family, nu-alpha] ] == 0 (mod p)
```

**Run:** Always runs when you call `VerifyRelation.wl <famname>`. No extra flags.

**Prerequisites:**
| File | Purpose |
|------|---------|
| `/home/ykm/blade-workspace/Blade.wl` | Blade IBP library |
| `verify/FamilyDatabase/FamilyDatabase.wl` | Family definitions |
| `AllRelations_<fam>_k<N>.m` (project root) | C++ relations (auto-picks highest k) |

**Output:**
```
=== MODULE 1: BladeVerify =============================
  nu: {2, 5}  eps: 27/14
  [1/3] Top sector reduction...
  Masters: 1  (2.3s)
  [2/3] Reducing integrals...
  Unique integrals: 5
  Reduced in 3.1s
  [3/3] Verifying...
  (lev=0  deg=0 ) PASS                pass=1/1
  BladeVerify: 10/10 PASSED
```

**Pitfalls:**
| Symptom | Cause | Fix |
|---------|-------|-----|
| `BLReduce[{...}]` returns nothing | Batch overload | Already handled: per-integral calls with cache |
| Single-master family shows `BL[family,{...}]` | Valid master, not error | `Coefficient` check handles this |
| Results differ between runs | Blade disk cache | `rm -rf build/cache/<family>*` |
| `results/usints` missing | maxdots too low | Delete cache, re-run; BlackBoxIBP auto-computes |

---

## Method 2: SeriesVerify (Expansion Substitution)

**Purpose:** Verify C++ relations cancel order-by-order when substituted into asymptotic expansion h_k coefficients.
```
Sum b_{alpha,beta} * nu^beta * A^{-alpha} * h_k(nu-alpha) == 0   (each k, mod p)
```

**Run:** Auto-runs if `ExpansionMMA_<famname>.m` and `RelationMeta_<famname>.m` exist.

**Prerequisites:**
| File | Source |
|------|--------|
| `ExpansionMMA_<famname>.m` (project root) | C++ test_relationFF auto-export |
| `RelationMeta_<famname>.m` (project root) | C++ test_relationFF auto-export |
| `verify/<fam>/Compare-CPPResult-<fam>.m` | C++ test_relationFF auto-export |

**Output:**
```
=== MODULE 2: SeriesVerify ============================
  (lev=0, deg=1) Series: PASS
  (lev=1, deg=2) Series: PASS
  SeriesVerify: 4/4 PASS
```
Or `SKIP: ExpansionMMA_*.m not found` if file missing.

**Pitfalls:**
| Symptom | Cause | Fix |
|---------|-------|-----|
| All PASS at m=0, FAIL at m>=1 | C++ solver `stable_order` truncation | Only m <= stable_order guaranteed; use higher k_max |
| `SKIP: no expansion data` for some (lev,deg) | No matching config in expansion file | Normal — only configs in both sources are tested |

**Specialist variant:** `verify/VerifyUtility/Verify-Series.wl` — detailed per-(lev,deg) verification with region-aware pAlpha/theta.

---

## Method 3: KiraVerify (Kira Cross-Validation)

**Purpose:** Substitute C++ relations at a test nu point into Kira IBP rules and verify they reduce to 0.
```
jExpr = Sum coeff * nu^beta * j[nu-alpha]
kiraReduce[jExpr] == 0 (mod p)
```

**Run:** Auto-runs if Kira rules and loader exist in expected locations.
```bash
cd /home/ykm/Large-Index-Expansion-CPP
wolframscript -file VerifyRelation.wl <famname> skip BladeVerify,SeriesVerify,MMACompare
```

**Prerequisites:**
| File | Source |
|------|--------|
| `verify/<fam>/kira_integrals.m` | Kira IBP output |
| `verify/VerifyUtility/KiraRuleLoader.wl` | Rule loader library |

### Generating Kira Rules

**Step 1: Kira input via M2Kira.wl**
```mathematica
<< "verify/VerifyUtility/M2Kira.wl"
KiraInputGenerate[famname, props, loopMom, kin,
  "KiraWorkingDir" -> "/tmp/kira_tests",
  "RMax" -> 6, "SMax" -> 2, "Sector" -> {1, 1}];
```
**Step 2: Run Kira**
```bash
cd /tmp/kira_tests/<fam>/<fam>
/usr/local/bin/kira jobs.yaml --parallel=physical
```
**Step 3: Copy rules**
```bash
cp /tmp/kira_tests/<fam>/results/<fam>/kira_integrals.m verify/<fam>/
```

### RMax/SMax Selection

| Parameter | Meaning | How to Determine |
|-----------|---------|------------------|
| `RMax` | Max dot product (sum of indices) | Max index in relations (max ν + max α-offset) |
| `SMax` | Max rank (number of non-zero indices) | Must cover top sector; do NOT increase for zero-component indices |

**Key insight:** Integrals with zero components (e.g. `j[4,0]`, `j[0,5]`) belong to **trivial sectors** and are set to zero by `KiraRuleLoader.wl` using sector mapping rules — they do NOT need Kira reduction. Do not increase SMax to cover them.

### Trivial Sector & Zero-Index Rules

`KiraRuleLoader.wl` automatically loads Kira explicit rules and adds two categories:
| Rule Type | Example | Sector | Result |
|-----------|---------|--------|--------|
| Zero-index | `j[a]/; all a ≤ 0` | `{0,0}` | → 0 |
| Trivial sector | `j[a,b]/; a>0 && b≤0` | `{1,0}` | → 0 |
| Trivial sector | `j[a,b]/; a≤0 && b>0` | `{0,1}` | → 0 |
| Kira explicit | `j[2,3]` etc. | Non-trivial | → master integrals |

### Index Coverage Verification

Select test ν points such that after substitution `ν → ν-α`, all resulting `j[ν-α]` indices are covered by Kira rules. The minimum ν requirement per relation depends on the α-offsets:
- `α={-1,0}` → needs `ν1 ≥ 2`
- `α={0,-2}` → needs `ν2 ≥ 3`

**Output:**
```
=== MODULE 3: KiraVerify ==============================
  Rules loaded: 83 symbolic rules
  (lev=2, deg=0) Kira: PASS  pass=1/1
  (lev=2, deg=1) Kira: PASS  pass=3/3
  (lev=2, deg=2) Kira: PASS  pass=7/7
  KiraVerify: 11/11 PASSED
```
Or `SKIP: Kira rules or KiraRuleLoader.wl not found`.

**Pitfalls:**
| Symptom | Cause | Fix |
|---------|-------|-----|
| All FAIL with non-zero residual | `j[alpha]` wrong format | Kira uses `j[a,b]` not `j[{a,b}]` — handled by loader |
| `INCONCLUSIVE` for some integrals | Missing trivial sector rules | `KiraRuleLoader` adds zero-index + trivial sector rules |
| Only high-lev configs tested | Low-lev NumSolutions=0 (stable_order not reached) | Run C++ with higher k_max |
| Uncovered integrals like `j[4,0]` at some ν | ν point creates indices in trivial sectors | Use ν where ν_i > max α-offset; loader auto-zeros trivial sectors |

---

## Method 4: MMACompare (MMA Relation Comparison)

**Purpose:** Compare C++ AllRelations coefficients against MMA `LIEReconstruct` output. Structural match, not numerical.

**Run:** Auto-runs if `MMARelations_<famname>.m` exists.

**Prerequisites:**
| File | Source |
|------|--------|
| `MMARelations_<famname>.m` (project root) | MMA `LIEReconstruct` export |

**Status:** Placeholder — prints config listings, coefficient comparison logic TBD.

---

## Standalone: EquationVerify (C++ Internal)

**Purpose:** M(nu) × coeffs = 0. Confirms C++ nullspace computation is internally consistent.

**Run:**
```bash
cd /home/ykm/Large-Index-Expansion-CPP
./build/test_relationFF bub00 4 2 2
```
The C++ binary also runs EquationVerify internally and prints per-config residual checks.

**Prerequisites:** C++ project built with `cmake --build .`

**Output:** Per-config residual check. All should be 0.

---

## Method Comparison

| Method | Backend | Verifies | Needs External Data | Status |
|--------|---------|----------|:---:|:------:|
| BladeVerify | Blade IBP | Relation ≡ 0 at nu point | No | Active |
| SeriesVerify | Expansion h_k | Relation vanishes order-by-order | Yes (ExpansionMMA) | Active |
| KiraVerify | Kira IBP | Relation ≡ 0 under Kira rules | Yes (kira rules) | Active |
| MMACompare | MMA LIEReconstruct | Coefficient match | Yes (MMA relations) | Placeholder |
| EquationVerify | C++ internal | M(nu)·coeffs = 0 | No | In C++ test |

---

## Verified Families

| Family | Blade | Series | Kira | Notes |
|--------|:-----:|:------:|:----:|-------|
| bub00 | PASS | PASS | 11/11 | 1L bubble |
| bub10 | PASS | — | — | |
| bub11 | PASS | — | — | |
| Tri | PASS | — | — | 1L triangle |
| SR212 | PASS | — | — | 2L sunrise |
| SR212-3m | PASS | — | — | Single master |
| SR212-5m | PASS | — | — | |
| NP222 | PASS | — | — | |
| TB123 | PASS* | — | — | All NumSolutions=0 |
| DB313 | PASS | — | — | |
| NP322 | PASS | — | — | |
| Box | FAIL | — | — | BLMaximalCutMasters fails; use Kira |

---

## Key Files

| File | Role |
|------|------|
| `VerifyRelation.wl` | Unified entry point (project root) |
| `verify/VerifyUtility/Verify-Blade.wl` | BladeVerify module |
| `verify/VerifyUtility/Verify-Series.wl` | SeriesVerify module |
| `verify/VerifyUtility/Verify-Kira.wl` | KiraVerify module |
| `verify/VerifyUtility/Verify-MMACompare.wl` | MMACompare module (placeholder) |
| `verify/FamilyDatabase/FamilyDatabase.wl` | Family configs (19 families) |
| `verify/VerifyUtility/KiraRuleLoader.wl` | Kira rule loader |
| `verify/docs/Test-Expand.md` | Expansion verification docs |
| `/home/ykm/blade-workspace/Blade.wl` | Blade IBP library |
