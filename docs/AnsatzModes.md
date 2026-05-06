# Ansatz Modes in RelationSolver

## Overview

The C++ RelationSolver supports four ansatz modes for generating the unknown coefficient variables `b[α, β]` in the LIE (Large Index Expansion) reduction relation reconstruction. Each mode defines how the α (seed/level) multi-indices are enumerated.

## Modified Files

| File | Changes |
|------|---------|
| [include/RelationSolver.hpp](../include/RelationSolver.hpp) | `AnsatzMode` enum, `alphaLevel()`, `generateAlphaSeeds()`, `negatedKernel_` propagation through `RegimeEvaluator`/`GlobalEquationAssembler`/`AdaptiveEquationBuilder`, updated `reconstructAllRelations` signature and logic |
| [tests/test_relationFF.cpp](../tests/test_relationFF.cpp) | `--mode` and `--sector` CLI flags, export metadata, mode-specific config output |

## CLI Usage

```
./test_relationFF <family> [order] [lev_min] [lev_max] [deg_max] [--topsector] [--mode <n>] [--sector <s>]
```

| Flag | Description |
|------|-------------|
| `--mode 0` | Pyramid (default) |
| `--mode 1` | DotPyramid |
| `--mode 2` | Star |
| `--mode 3` | ExtendedPyramid |
| `--sector 101` | ExtendedPyramid sector override (string of 0/1, default = auto top sector) |

Examples:
```bash
./test_relationFF bub 4 1 2 2                  # Pyramid (default)
./test_relationFF bub 4 1 2 2 --mode 1          # DotPyramid
./test_relationFF bub 4 1 2 2 --mode 2          # Star
./test_relationFF bub 4 0 2 2 --mode 3          # ExtendedPyramid
./test_relationFF bub 4 0 2 2 --mode 3 --sector 110  # ExtendedPyramid, only dims 0,1
```

## Mode Specifications

### Common Properties

- **β (degree) indices**: always non-negative integer vectors with `|β| ≤ deg`, generated via `generateAllIndices` in MMA `GenerateSeeds` order (graded by L1 then colexicographic descending).
- **MMA ordering**: all α lists are sorted to match MMA `GenerateSeeds` / `SectorMO` order, enabling direct column-by-column comparison of exported `.m` files.
- **RemoveSolvedVariables**: `filterVariablePairs` applies componentwise dominance filtering identically across all modes.
- **Exported `.m` files** include an `"AnsatzMode"` metadata field (`"Pyramid"` / `"DotPyramid"` / `"Star"` / `"ExtendedPyramid"`).

---

### Mode 0: Pyramid

**Kernel**: `g(ν − α)` — `negatedKernel = false`

**α generation**: all non-negative integer vectors with `|α| ≤ lev`. Same as MMA `GenerateAnsatz` Mode 0 base ansatz.

**Level filter**: `Total(α)` (sum of components).

**Level range**: `lev_min .. lev_max`, can start from 0.

**Example** (ne=2, lev=2): 6 alphas
```
|α|=0:  {0,0}
|α|=1:  {0,1}, {1,0}
|α|=2:  {0,2}, {1,1}, {2,0}
```

**Boundary terms** (MMA Pyramid boundary for ISP positions) are **not yet implemented** in C++. This is a deferred feature.

---

### Mode 1: DotPyramid

**Kernel**: `g(ν + α)` — `negatedKernel = true`

**α generation**: all non-negative integer vectors with `|α| ≤ lev+1` (one extra level vs Pyramid). Same as MMA `GenerateAnsatz` Mode 1.

**Level filter**: `Total(α)`.

**Level range**: `max(1, lev_min) .. lev_max` (no level 0 — force-starts from 1).

**Example** (ne=2, lev=2): 10 alphas (includes |α|=3 shell)
```
|α|=0:  {0,0}
|α|=1:  {0,1}, {1,0}
|α|=2:  {0,2}, {1,1}, {2,0}
|α|=3:  {0,3}, {1,2}, {2,1}, {3,0}
```

---

### Mode 2: Star

**Kernel**: `g(ν + α)` — `negatedKernel = true`

**α generation**: mixed-sign seeds via `DotRankSeeds` per-sector logic. Iterates all `2^ne` sector subsets. For each subset:
- **pdSet** (propagator-like, sector=1): seed[i] = 1 + dotSeed[i] ≥ 1
- **ispSet** (ISP-like, sector=0): seed[i] = -rankSeed[i] ≤ 0

Seeds are sorted by MMA `SectorMO` ordering, then by L1 norm. Duplicates across sectors are removed.

**Level filter**: `Total[Abs(α)]` (L1 norm of signed components).

**Level range**: `max(1, lev_min) .. lev_max`.

**Sector parameter**: uses `sector[0]` from the ring data to determine pdSet/ispSet split. Falls back to all-1s if no sector data is available.

**Example** (ne=2, lev=2): 30 alphas (all subsets, mixed signs)

---

### Mode 3: ExtendedPyramid

**Kernel**: `g(ν − α)` — `negatedKernel = false`

**α generation**: take `Pyramid(lev+1)` (all non-negative with `|rk| ≤ lev+1`), and for each active direction `i` where `sector[i] == 1`, create `α = rk − e_i`. Take the union over all such shifts.

This produces alphas where at most one component is `−1` (only in active sector directions), and all other components are ≥ 0.

**Level filter**: `Total(α)` (sum of components, may be as low as −1).

**Level range**: `lev_min .. lev_max`, can start from 0.

**Sector parameter**: `--sector` CLI flag sets which directions are active for the `−e_i` shift.
- Default (no flag): auto-detects the top sector (maximum sum) from the ring data — typically all-1s.
- `--sector 10`: only shift in dimension 0.
- `--sector 01`: only shift in dimension 1.
- `--sector 11`: explicit all directions (equivalent to default for ne=2).

**Example** (ne=2, lev=2, sector=[1,1]): 14 alphas
```
Total=-1:  {-1,0}, {0,-1}
Total= 0:  {-1,1}, {0,0}, {1,-1}
Total= 1:  {-1,2}, {0,1}, {1,0}, {2,-1}
Total= 2:  {-1,3}, {0,2}, {1,1}, {2,0}, {3,-1}
```

**Example** (ne=2, lev=2, sector=[1,0]): 10 alphas (only dim 0 gets −1)
```
Total=-1:  {-1,0}
Total= 0:  {-1,1}, {0,0}
Total= 1:  {-1,2}, {0,1}, {1,0}
Total= 2:  {-1,3}, {0,2}, {1,1}, {2,0}
```

## Implementation Architecture

### Data Flow

```
generateAlphaSeeds(mode, ne, lev, sector) → alphas
         ↓
reconstructAllRelations(CTable, ..., mode, ext_sector)
         ↓
    for lev in [eff_lev_min .. lev_max]:
      filter alphas by alphaLevel(mode, α) ≤ lev
      for deg in [0 .. deg_max]:
        generate betas
        filterVariablePairs (RemoveSolvedVariables)
        builder.build(regimes, nimax_lists, ne, kernel_negated, &alphas, &betas)
```

### Key Changes

1. **`AnsatzMode` enum** (line 37): `Pyramid=0, DotPyramid=1, Star=2, ExtendedPyramid=3`
2. **`alphaLevel()`**: mode-specific level function — `Total` for Pyramid/DotPyramid/ExtendedPyramid, `Total[Abs]` for Star
3. **`generateAlphaSeeds()`**: per-mode alpha list generation with MMA-compatible ordering
4. **`negatedKernel_`**: boolean propagated through `RegimeEvaluator → GlobalEquationAssembler → AdaptiveEquationBuilder`. Controls `ν−α` vs `ν+α` in `step2_computeG`
5. **`AdaptiveEquationBuilder::build()`**: accepts optional `ext_alphas`/`ext_betas` pointers — when provided (non-Pyramid modes), uses caller-generated ansatz instead of generating Pyramid from config hints
6. **`reconstructAllRelations()`**: new `mode` and `ext_sector` parameters; per-mode lev_min adjustment, kernel sign determination, alpha seed generation
7. **`test_relationFF.cpp`**: `--mode` and `--sector` CLI parsing; mode metadata in `.m` exports

## Verification

All four modes compile and run successfully on the `bub` family. The exported `.m` files include `"AnsatzMode"` metadata for MMA-side verification scripts.

### Relation counts (bub, k=4, lev=2, deg=2)

| Mode | Active vars (max) | Relations at lev=2,deg=1 |
|------|------------------|--------------------------|
| Pyramid | 26/36 | 4 |
| DotPyramid | 30/60 | 4 |
| Star | 65/130 | — |
| ExtendedPyramid | 24/84 | 4 |
