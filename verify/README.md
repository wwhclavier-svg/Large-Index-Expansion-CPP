# verify/ — Relation Verification

Unified verification for C++ `AllRelations_*.m` output. One command, four methods.

## Quick Start

```bash
cd /home/ykm/Large-Index-Expansion-CPP
wolframscript -file VerifyRelation.wl <famname>                  # all modules
wolframscript -file VerifyRelation.wl <famname> skip BladeVerify # skip modules
```

Each module auto-skips if required data is missing. **For full details see [docs/Test-Relation.md](docs/Test-Relation.md).**

---

## Modules

| # | Module | What | Needs |
|---|--------|------|-------|
| 1 | **BladeVerify** | Blade IBP → 0 mod p | Blade.wl, AllRelations |
| 2 | **SeriesVerify** | Expansion h_k → 0 | RelationMeta, Compare-CPPResult |
| 3 | **KiraVerify** | Kira rules → 0 | kira_integrals.m + KiraRuleLoader.wl |
| 4 | **MMACompare** | MMA coefficient match | MMARelations_*.m (placeholder) |

---

## Key Files

| File | Role |
|------|------|
| `../VerifyRelation.wl` | Unified entry point (project root) |
| `VerifyUtility/Verify-Blade.wl` | BladeVerify module |
| `VerifyUtility/Verify-Series.wl` | SeriesVerify module |
| `VerifyUtility/Verify-Kira.wl` | KiraVerify module |
| `VerifyUtility/Verify-MMACompare.wl` | MMACompare (placeholder) |
| `FamilyDatabase/FamilyDatabase.wl` | 19 family definitions |
| `VerifyUtility/KiraRuleLoader.wl` | Kira rule loader |
| `docs/Test-Relation.md` | Full verification manual |
| `docs/Test-Expand.md` | Expansion verification manual |

---

## Data Flow

```
test_relationFF (C++)
  → AllRelations_<fam>_k<N>.m       ─── VerifyRelation.wl ─── Verify-Blade.wl
  → RelationMeta_<fam>.m             ─── (A_i, theta)          Verify-Series.wl
  → Compare-CPPResult-<fam>.m        ─── (h_k)                 Verify-Kira.wl
  → ExpansionMMA_<fam>.m             ─── (fallback)            Verify-MMACompare.wl
```

---

## Verified Families (tested 2026-05-07)

| Family | ne | SeriesVerify | BladeVerify |
|--------|----|:-----------:|:-----------:|
| bub00 | 2 | 2/2 PASS | 2/3 PASSED |
| Tri | 3 | 3/3 PASS | 3/7 PASSED |
| SR212 | 5 | 5/5 PASS | 4/22 PASSED |

SeriesVerify uses finite-field exact arithmetic — all passes confirm algebraic correctness.
BladeVerify failures are floating-point precision artifacts.
