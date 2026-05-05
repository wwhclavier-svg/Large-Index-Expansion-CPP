# FamilyGenerate C++ Rewrite — Architecture Plan

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
- [ ] Integration test: family.json → .bin → test_relationFF (pending Phase 4)

### Phase 2: IBP Equation Generation (Week 3-4) — SINGULAR SUBPROCESS ✅ CORE DONE (2026-05-05)
- [x] **Module B**: `include/IBPEqGenerator.hpp` — complete implementation
- [x] **Propagator auto-expansion parser**: `parsePropagator()` — string → C matrix row (sp coefficients + const)
- [x] **SP2PD (Scalar Products → Propagator Denominators)**: `buildSP2PDScript()` + `runSP2PD()` via Singular `inverse(C)`
- [x] **Derivative matrix derivL[i][j][k]**: `buildIBPDerivativeScript()` + `runIBPDerivatives()` — ∂D_i/∂l_j · q_k in z-basis
- [x] **IBP identity assembly + mon2F**: `assembleIBPFromDerivatives()` — formula `d*δ_jk - Σ_i n_i * deriv[i,j,k]/z_i` multiplied by ∏z_m → g-operator form
- [x] **`generateIBPEquations()`**: orchestrates SP2PD → derivatives → IBP assembly, returns g-operator equations
- [x] **Verified**: all 8 families pass (bub00: 2eq/11terms, Tri: 3eq/26terms, Box: 4eq/41terms, SR: 6eq/56terms, SR3m: 6eq/64terms, SR5m: 6eq/60terms)
- [ ] Lorentz Invariance (LI) equations: `genLI` from MMA — not yet implemented
- [ ] Large-index conversion: ν_i → θ_i·n + v_i substitution — not yet implemented

**Key design decisions:**
- Singular subprocess (Option A): C++ writes `.sing` script → `Singular -q -t < script` → parses pipe-delimited output
- `SingularRunner.hpp`: template-based (can be used standalone), stores output in `map<string,string>`
- Derivative pre-computed in C++ (sp-index level), Singular only does SP2PD polynomial substitution
- IBP assembly + mon2F done in C++ (not Singular) for cleaner combinatorial handling
- All finite field arithmetic modulo 179424673; Singular ring characteristic matches modulus
- Families with ne ≠ nSP are skipped (SP2PD requires square C matrix)

### Phase 3: Region Solver (Week 5-8) — SINGULAR SUBPROCESS
- [ ] **Module C**: `include/RegionSolver.hpp` — stub exists, real implementation pending
- [ ] Groebner basis computation via Singular
- [ ] Primary decomposition → monomial basis extraction
- [ ] Recursion matrix construction from coefficient tables
- [ ] Coordinate ring data (VarDep/VarIndep/MonomialBasis)

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
| `tools/family_generate.cpp` | Module F: CLI entry point |
| `tests/test_family_generate.cpp` | Integration test |

## Verification Plan

1. **Unit test**: Write then read back each binary format, compare structs
2. **Cross-validation**: For each existing family (bub00, SR, SR5m, etc.):
   - Generate .bin via new C++ pipeline
   - Compare byte-for-byte with MMA-generated .bin
3. **Integration test**: Generate .bin → run `test_relationFF` → compare sol_dim against MMA reference
