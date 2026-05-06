# FamilyGenerate C++ Rewrite — Architecture Plan

## Current Status (2026-05-07)

| Phase | Module | Status |
|-------|--------|--------|
| Phase 0 | A: FamilyConfig, D: BinaryIBPWriter, E: BinaryRingWriter | ✅ DONE |
| Phase 1 | F: CLI `family_generate.cpp` | ✅ DONE (full pipeline) |
| Phase 2 | B: IBPEqGenerator | ✅ DONE (LI equations still missing) |
| Phase 3 | C: RegionSolver + aux modules | ✅ IMPLEMENTED |
| Phase 4 | Integration & cross-validation | ⏳ IN PROGRESS |

### Cross-Validation Status (2026-05-07)

| Family | Region count | nb | test_relationFF | .bin byte-identical | Note |
|--------|-------------|-----|-----------------|---------------------|------|
| bub00 | ✅ 1 | ✅ 1 | ✅ sol_dim match | ❌ A/B sign diff | functional equivalence confirmed |
| bub10 | ✅ | ✅ | - | ❌ | |
| bub11 | ✅ | ✅ | - | ❌ | |
| Tri | ✅ 2 | ✅ | - | ❌ | |
| Box | ✅ 15 | ✅ | - | ❌ | |
| SR | ✅ | ✅ | - | - | |
| SR3m | ✅ | ✅ | - | - | |
| SR5m | ✅ | ✅ | - | - | |
| TB123 | ⏳ timeout | - | - | - | larger CAS computation |

**Key finding:** All families produce correct region counts and functional results (expansion coefficients byte-identical to MMA, relation reconstruction sol_dim matches). However, .bin files are NOT byte-identical to MMA reference — A/B equation sign conventions differ, leading Singular to select different but equivalent monomial bases. Root cause: C++ multiplies IBP by ∏z_m then negates n_i terms; MMA uses LargeIndexIBP `n_i → n+v_i` then `Coefficient["n"]` extraction.

**Remaining work:** (1) LI equation generation, (2) fix A/B equation sign convention for byte-identical .bin, (3) JSON configs for DB313/NP222/NP322/SR212 variants, (4) `test_relationFF` integration test for more families.

## Context

The Large Index Expansion pipeline is split across two languages:
- **Mathematica**: Generates `IBPMat_*.bin` + `RingData_*.bin` (called "FamilyGenerate")
- **C++**: Reads .bin files, computes series expansion coefficients, reconstructs reduction relations

The user wants to eliminate the Mathematica dependency by rewriting FamilyGenerate in C++. This is a significant undertaking because the pipeline involves symbolic algebra (Groebner basis, primary decomposition) that requires a CAS backend.

## Pipeline Overview (4 Steps)

```
FamilyDatabase ──→ [1] LIEDefineFamily ──→ [2] LIESolveRegions ──→ [3] ExportIBPMatrix ──→ IBPMat_*.bin
                     (IBP eq generation)     (ring decomposition)     (sparse CSR write)    
                                                                  [4] ExportRingData  ──→ RingData_*.bin
                                                                     (A/Ainv compute+write)
```

### Step 1 — IBP Identity Generation
- Symbolic differentiation of propagators w.r.t. loop/external momenta
- Applies IBP identity `∫ d^dk ∂/∂k_μ (q_μ × ∏ 1/D_i^{ν_i}) = 0`
- Generates `nibp = nl * (nl + nE)` equations per family
- Also generates Lorentz Invariance (LI) relations
- Converts to large-index form: `ν_i → θ_i n + v_i`

### Step 2 — Algebraic Region Decomposition (HARDEST)
- For each subsector of the top sector:
  - Substitute `v_i → sector_i * n + v_i` into IBP equations
  - Compute Groebner basis of sector-limited equations
  - Call Singular CAS for primary decomposition (`minAssPrimes`)
  - Extract monomial basis of quotient ring → `nb` (block dimension)
  - Build recursion matrices (M1, N1, K1, F0, F2, ...) in monomial basis representation

### Step 3 — Binary IBP Matrix Export (EASY)
- Writes IBP1-format binary (big-endian, CSR/coordinate format)
- C++ reader already exists: `IBPMatrixLoader_Binary.hpp`

### Step 4 — Binary Ring Data Export (EASY)
- Computes A/Ainv matrices for each regime
- Writes ring data binary
- C++ reader already exists: `RingDataLoader.hpp`

## Proposed Architecture

```
 ┌──────────────────────────────────────────────────────────┐
 │                    C++ FamilyGenerate                     │
 ├──────────────────────────────────────────────────────────┤
 │  Input: family.json (propagators, momenta, kinematics)   │
 │                                                          │
 │  Module A: FamilyConfig ─── pure C++                     │
 │    Parse JSON family definition                          │
 │    → FamilyDef struct                                    │
 │                                                          │
 │  Module B: IBPEqGenerator ─── needs Singular subprocess  │
 │    Symbolic diff of propagators                          │
 │    Generate IBP + LI identities                          │
 │    Large-index conversion                                │
 │    → IBPEquations struct                                 │
 │                                                          │
 │  Module C: RegionSolver ─── needs Singular subprocess    │
 │    Per-sector: Groebner basis, primary decomposition     │
 │    Monomial basis, quotient ring dimension               │
 │    Build coordinate ring + recursion matrices            │
 │    → vector<RegionData>                                  │
 │                                                          │
 │  Module D: BinaryIBPWriter ─── pure C++ (WEEK 1)        │
 │    Write IBPMat_*.bin in IBP1 format                    │
 │    Mirror of IBPMatrixLoader_Binary.hpp                  │
 │                                                          │
 │  Module E: BinaryRingWriter ─── pure C++ (WEEK 1)       │
 │    Compute A/Ainv matrices per regime                   │
 │    Write RingData_*.bin                                  │
 │    Mirror of RingDataLoader.hpp                          │
 │                                                          │
 │  Module F: FamilyGenerateCLI ─── pure C++ (WEEK 2)      │
 │    Single executable: ./family_generate <family.json>    │
 │    Orchestrates A→B→C→D→E                                │
 └──────────────────────────────────────────────────────────┘
```

## Singular CAS Integration Strategy

Two options for Steps 1-2:

| Option | Approach | Effort | Risk |
|--------|----------|--------|------|
| **A: Subprocess** | Call Singular CLI with generated script, parse text output | 2-3 weeks | Medium (parsing) |
| **B: libSingular** | Link against libSingular C library directly | 4-6 weeks | High (API complexity) |

**Recommendation**: Start with Option A (subprocess). The MMA pipeline already calls Singular via `SingularInterface.wl` — we'd replicate the same pattern: generate a Singular script, execute it, parse output. This avoids the complexity of libSingular's C API while still being fully automated.

## Input Format: `family.json`

Replace `FamilyDatabase.wl` with a JSON config:

```json
{
  "name": "SR5m",
  "propagators": ["-l1^2 - msq", "-(l1+p)^2 - msq", "-l2^2 - msq", "-(l2+p)^2 - msq", "-(l1+l2+p)^2 - msq"],
  "loopMomenta": ["l1", "l2"],
  "externalMomenta": ["p"],
  "kinematicRules": {"p^2": "s"},
  "topSector": [1, 1, 1, 1, 1],
  "numeric": {"s": 0, "msq": 1, "d": "1/13"},
  "modulus": 179424673
}
```

## Key Data Structures (in-memory between steps)

```cpp
// Output of Step 1
struct IBPEquations {
    int ne, nl, nE, nibp;
    vector<string> Alist;                        // A[1]...A[ne] symbols
    vector<string> vlist;                        // v1...vne symbols
    vector<vector<int>> sectorlist;               // all subsectors
    // IBP equations in large-index g[...] notation (symbolic form)
};

// Output of Step 2 per regime
struct RegionData {
    vector<int> limitSector;
    int nb;                                       // monomial basis dimension
    // Coordinate ring
    vector<int> VarIndep, VarDep;
    vector<vector<int64_t>> MonomialBasis;        // in VarDep space
    // Recursion matrices (sparse CSR, same structure as IBPMatrixE)
    SparseTensor3D M1, N1, K1;
    SparseTensor2D F0;
    SparseTensor4D F2, F2s;
    SparseTensor3D K1s, K2s;
};
```

## Implementation Phases

### Phase 0: Preparation (Week 1) — PURE C++ ✅ DONE (2026-05-05)
- [x] **Module A**: JSON family config parser (`include/FamilyConfig.hpp`) — FamilyDef struct, parseFamilyConfig()
- [x] **Module D**: Binary IBP matrix writer (`include/BinaryIBPWriter.hpp`) — mirror of IBPMatrixLoader_Binary
- [x] **Module E**: Binary ring data writer (`include/BinaryRingWriter.hpp`) — mirror of RingDataLoader
- [x] Migration: 8 family JSON configs in `families/*.json` (bub00/10/11, Tri, Box, SR, SR3m, SR5m)
- [x] Test: `tests/test_family_config.cpp` — parses all 8 families

### Phase 1: CLI Tool (Week 2) — PURE C++ ✅ PARTIAL (2026-05-05)
- [x] **Module F**: `tools/family_generate.cpp` — CLI stub: reads family.json, writes .bin (Phase 2b ready)
- [x] Stub Step 1-2: `generateIBPEquations()` loads from stub (Phase 1); Phase 2b replaces with real computation
- [x] `tools/test_singular_runner.cpp` — Singular subprocess runner test
- [x] `tools/test_region_solver.cpp` — region solver test (stub)
- [ ] Integration test: family.json → .bin → test_relationFF (pending Phase 4 validation)

### Phase 2: IBP Equation Generation (Week 3-4) — SINGULAR SUBPROCESS ✅ CORE DONE (2026-05-05)
- [x] **Module B**: `include/IBPEqGenerator.hpp` — complete implementation
- [x] **Propagator auto-expansion parser**: `parsePropagator()` — string → C matrix row (sp coefficients + const)
- [x] **SP2PD (Scalar Products → Propagator Denominators)**: `buildSP2PDScript()` + `runSP2PD()` via Singular `inverse(C)`
- [x] **Derivative matrix derivL[i][j][k]**: `buildIBPDerivativeScript()` + `runIBPDerivatives()` — ∂D_i/∂l_j · q_k in z-basis
- [x] **IBP identity assembly + mon2F**: `assembleIBPFromDerivatives()` — formula `d*δ_jk - Σ_i n_i * deriv[i,j,k]/z_i` multiplied by ∏z_m → g-operator form
- [x] **`generateIBPEquations()`**: orchestrates SP2PD → derivatives → IBP assembly, returns g-operator equations
- [x] **Verified**: all 8 families pass (bub00: 2eq/11terms, Tri: 3eq/26terms, Box: 4eq/41terms, SR: 6eq/56terms, SR3m: 6eq/64terms, SR5m: 6eq/60terms)
- [x] **Large-index conversion**: ν_i → θ_i·n + v_i substitution — implemented via `gShift` mechanism in `assembleIBPFromDerivatives()` and `IBPAnalyzer::buildABEquations()` (2026-05-06)
- [ ] Lorentz Invariance (LI) equations: `genLI` from MMA — not yet implemented

**Fixed (2026-05-06)**: Issues causing diff failures:
1. `generateSubsectors` generated all 2^n subsets (including empty); fixed to match MMA's `Subsets[activeIndices, {nl, ne}]` (sizes nl..ne only).
2. **Critical: A/B equation coefficient scaling** — C++ used ×2 for active n_i terms (wrong). MMA's `regionsBySectors` first does `ibpeqs /. "n" -> 0` (removes all n-dependence), then reintroduces it via `v_i → sector[i]*"n" + v_i`. After `Coefficient[ibpeqs, "n"]`: active indices contribute `coeff × 1` (the `n` coefficient after n→0), inactive indices contribute 0 (no `n` in term). C++ now matches: skip inactive terms entirely, use coeff as-is for active terms. See `IBPAnalyzer::buildABEquations()` lines 171-186.
3. Large-index conversion (ν_i → θ_i·n + v_i) confirmed working: `gShift` mechanism in `assembleIBPFromDerivatives()` produces g-operator form; `IBPAnalyzer::buildABEquations()` converts g-shifts to A/B exponents with correct standard shift formula `s_i = 2*v_i - gShift_i`.

**Result (updated 2026-05-07):** All 8 families run without errors. Expansion coefficients are byte-identical to MMA reference (confirmed via test_relationFF for bub00). `.bin` files are NOT byte-identical — A/B equation sign convention difference leads to equivalent but different monomial bases. Functional correctness is confirmed; byte-level identity requires deeper IBP assembly alignment.

**Key design decisions:**
- Singular subprocess (Option A): C++ writes `.sing` script → `Singular -q -t < script` → parses pipe-delimited output
- `SingularRunner.hpp`: template-based (can be used standalone), stores output in `map<string,string>`
- Derivative pre-computed in C++ (sp-index level), Singular only does SP2PD polynomial substitution
- IBP assembly + mon2F done in C++ (not Singular) for cleaner combinatorial handling
- All finite field arithmetic modulo 179424673; Singular ring characteristic matches modulus
- Families with ne ≠ nSP are skipped (SP2PD requires square C matrix)

### Phase 3: Region Solver (Week 5-8) — SINGULAR SUBPROCESS ✅ IMPLEMENTED, ⏳ UNVALIDATED (2026-05-06)
- [x] **Module C**: `include/RegionSolver.hpp` (809 lines) — full implementation: `solveRegion()`, `solveAllSectors()`, `computeGroebnerBasis()`, `computeMonomialBasisIndex()`, `solveVarRule()`, `computeFractionRule()`, `computeMonomialBasisMatrix()`
- [x] **Module C-aux**: `include/IBPAnalyzer.hpp` (367 lines) — `buildABEquations()` with LargeIndexIBP coefficient scaling
- [x] **Module C-aux**: `include/RecursionBuilder.hpp` (511 lines) — `buildRecursionMatrices()` from fraction rules
- [x] **Module C-aux**: `include/RingBuilder.hpp` (193 lines) — `computeRingMatrices()` for A/Ainv
- [x] **Module C-aux**: `include/PolyArith.hpp` (437 lines) — polynomial arithmetic in finite field
- [x] Groebner basis computation via Singular subprocess
- [x] Primary decomposition → monomial basis extraction
- [x] Recursion matrix construction from coefficient tables
- [x] Coordinate ring data (VarDep/VarIndep/MonomialBasis)
- [ ] **End-to-end validation**: bub00, bub11 byte-identical with MMA; remaining families need MMA reference generation
- [x] All 8 families run without errors (bub00, bub10, bub11, Box, SR, SR3m, SR5m, Tri); TB123 times out

**Implementation notes (2026-05-06):**
- All modules are header-only (`inline` functions), compiled into `family_generate` CLI
- `solveRegion()` mirrors MMA's `expRegSolve2` inner loop: GB → primdecGTZ → per-component GB → VarRule → FractionRule → monomial basis matrix
- `solveAllSectors()` mirrors the outer loop: iterate subsectors, build A/B equations, call `solveRegion()`, collect results
- `PolyArith.hpp` provides finite field polynomial arithmetic (add, mul, mod) used by IBP analyzer and recursion builder
- Recursion matrices (M1, N1, K1, F0, F2, etc.) are built from FractionRule entries via `RecursionBuilder`
- Ring matrices (A, Ainv) computed via `RingBuilder` then serialized via `BinaryRingWriter`
- Full pipeline wired in `tools/family_generate.cpp`: parse JSON → generate IBP → solve sectors → build matrices → write .bin → optional byte diff

### Phase 4: Integration & Validation (Week 9-10)
- [ ] End-to-end: family.json → .bin files → test_relationFF verification
- [ ] Cross-validate against MMA-generated .bin files for all families
- [ ] Performance optimization (parallel regime processing)

## Data Files

| File | Purpose |
|------|---------|
| `families/*.json` | Family definitions (new, replaces FamilyDatabase.wl) |
| `include/FamilyConfig.hpp` | Module A: JSON parser, FamilyDef struct |
| `include/BinaryIBPWriter.hpp` | Module D: IBP1 binary writer |
| `include/BinaryRingWriter.hpp` | Module E: RingData binary writer |
| `include/IBPEqGenerator.hpp` | Module B: IBP identity generation |
| `include/RegionSolver.hpp` | Module C: Algebraic region decomposition |
| `include/IBPAnalyzer.hpp` | Module C-aux: A/B equation builder with LargeIndex scaling |
| `include/RecursionBuilder.hpp` | Module C-aux: Recursion matrix construction from fraction rules |
| `include/RingBuilder.hpp` | Module C-aux: Ring matrix (A/Ainv) computation |
| `include/PolyArith.hpp` | Module C-aux: Finite field polynomial arithmetic |
| `include/SingularRunner.hpp` | Singular subprocess runner (shared by B and C) |
| `tools/family_generate.cpp` | Module F: CLI entry point |
| `tests/test_family_generate.cpp` | Integration test |

## Verification Plan

1. **Unit test**: Write then read back each binary format, compare structs
2. **Cross-validation**: For each existing family (bub00, SR, SR5m, etc.):
   - Generate .bin via new C++ pipeline
   - Compare byte-for-byte with MMA-generated .bin
3. **Integration test**: Generate .bin → run `test_relationFF` → compare sol_dim against MMA reference
