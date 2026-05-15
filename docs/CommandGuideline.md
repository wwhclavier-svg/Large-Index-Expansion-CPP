# Command Guideline

## region-solve-cpp (Stage 1 — IBPMat + RingData generation)

### Build

| Step | (a) Read from | (b) Command | (c) Output to |
|------|-------------|------------|--------------|
| Configure | `CMakeLists.txt`, `families/*.json` | `mkdir -p build && cd build && cmake ..` | `build/` |
| Compile | `tools/family_generate.cpp`, `include/*.hpp` | `cmake --build . --target family_generate && cp family_generate ..` | `build/family_generate` → project root[^1] |

[^1]: `family_scheduler.py` expects the binary at project root.

### Run — Option A: Tiered Scheduler (recommended for large families)

| Step | (a) Read from | (b) Command | (c) Output to |
|------|-------------|------------|--------------|
| Launch | `families/<family>.json` | `python3 tools/family_scheduler.py <family> --workers 4` | Per-sector: `cache/<family>_scheduled/SectorCache/IBPMat_<fam>_sector_*.bin` + `RingData_*` |
| Assembly (auto) | per-sector `.bin` files | auto-runs `family_generate --merge-sectors` | `cache/<family>_scheduled/output/IBPMat_<family>.bin` + `RingData_<family>.bin` |
| Resume after crash | `cache/<family>_scheduled/status.json` | `python3 tools/family_scheduler.py <family> --resume --workers 2` | same as above |
| Force re-assembly | per-sector `.bin` files | `python3 tools/family_scheduler.py <family> --assemble` | `cache/.../output/IBPMat_<family>.bin` + `RingData_<family>.bin` |
| Custom tiers | families JSON | `--tiers "Quick=120 Fast=1200 Retry=86400"` | same as above |

### Run — Option B: Manual Per-Sector

| Step | (a) Read from | (b) Command | (c) Output to |
|------|-------------|------------|--------------|
| All sectors | `families/<family>.json` | `./build/family_generate families/<family>.json` | `data/IBPMat_<family>.bin` + `data/RingData_<family>.bin` |
| Single sector | `families/<family>.json` | `./build/family_generate families/<family>.json --sector 11 --cache-dir cache/sectors` | `cache/sectors/IBPMat_<fam>_sector_11.bin` |
| Merge manually | per-sector `.bin` files | `./build/family_generate families/<family>.json --merge-sectors cache/sectors --output data/` | `data/IBPMat_<family>.bin` + `data/RingData_<family>.bin` |
| Diff against ref | families JSON + `data/` | `./build/family_generate families/<family>.json --diff --source-dir .` | stdout: IDENTICAL / DIFFER |

> Must run from project root. `--sector` key = binary string (e.g. `"110"`, `"1011"`).

---

## relation-test (Stage 3 — Relation solving)

### Build

| Step | (a) Read from | (b) Command | (c) Output to |
|------|-------------|------------|--------------|
| Configure | same CMakeLists | `cd build && cmake ..` | `build/` |
| Compile | `tests/test_relationFF.cpp`, `include/RelationSolver.hpp` | `cmake --build . --target test_relationFF` | `build/test_relationFF` |

### Run

| Step | (a) Read from | (b) Command | (c) Output to |
|------|-------------|------------|--------------|
| Solve relations | `data/IBPMat_<family>.bin` + `data/RingData_<family>.bin` | `./build/test_relationFF <family> [order] [lev_min] [lev_max] [deg_max]` | `relations/AllRelations_<family>_k<order>.m` + `relations/RelationMeta_<family>.m` |
| Custom output dir | same `.bin` files | add `--output <dir>` | `<dir>/AllRelations_<family>_k<order>.m` |
| Top sector only | same `.bin` files | add `--topsector` | same as above |
| Change ansatz | same `.bin` files | add `--mode <0|1|2|3>` | same as above (0=Pyramid, 1=DotPyramid, 2=Star, 3=ExtendedPyramid) |

> Must run from project root. The binary also exports expansion coefficients automatically (see below).

### Verification

| Method | (a) Read from | (b) Command | (c) Output to |
|--------|-------------|------------|--------------|
| All-in-one | `AllRelations_*.m`, `RelationMeta_*.m` | `wolframscript -file VerifyRelation.wl <fam>` | stdout pass/fail per config |
| BladeVerify | `AllRelations_*.m` + Blade.wl | `wolframscript -file workspace/shared/VerifyUtility/Verify-Blade.wl <fam>` | stdout |
| SeriesVerify | `AllRelations_*.m` + `RelationMeta_*.m` + `ExpansionMMA_*.m` | `wolframscript -file workspace/shared/VerifyUtility/Verify-Series.wl <fam>` | stdout |
| KiraVerify | `AllRelations_*.m` + `verify/<fam>/kira_integrals.m` | `wolframscript -file workspace/shared/VerifyUtility/Verify-Kira.wl <fam>` | stdout |

---

### Expansion Series Coefficient Generation

Both `test_relationFF` and `test_expandFF` generate the intermediate expansion series coefficients (Stage 2 — Layer Recursion). Neither one is standalone-only — they share the same `batchProcessRecursion<FFInt>()` backend:

| Binary | Is the coefficient file generated? | File path |
|--------|------------------------------------|-----------|
| `test_relationFF` | **Yes** (always, as part of its pipeline before relation solving) | Cached to `data/ExpansionCache_<family>_k<order>.bin`; exported to `ExpansionMMA_<family>.m` and `verify/<family>/Compare-CPPResult-<family>.m` |
| `test_expandFF` | **Yes** (standalone Stage 2 only, no relation solving) | Cached to `data/ExpansionCache_<family>_k<order>.bin`; exported to `verify/<family>/Compare-CPPResult-<family>.m` and `verify/<family>/Compare-CPPMeta-<family>.m` |

**When to use `test_expandFF` instead of `test_relationFF`:**
- You only need the expansion coefficients and do NOT want relation solving.
- You want the `Compare-CPPMeta-<family>.m` metadata export (which includes region metadata + k=0,1 coefficients in a compact summary format).

**When `test_relationFF` is sufficient:**
- You need relations anyway (it computes the expansion as a prerequisite).
- The expansion is cached to `ExpansionCache_<family>_k<order>.bin` on first run, and re-used on subsequent runs (unless you delete the cache file).

### Combined Example

```bash
# Build everything once
cd build && cmake .. && cmake --build . --target family_generate && cmake --build . --target test_relationFF && cp family_generate .. && cd ..

# Full pipeline: regions → relations
python3 tools/family_scheduler.py bub00 --workers 4
./build/test_relationFF bub00 4 1 2 2

# Verify
wolframscript -file VerifyRelation.wl bub00
```

## See Also

| Document | Content |
|----------|---------|
| [`docs/README.md`](./README.md) | Documentation index, all executables mapped to pipeline stages |
| [`CLAUDE.md`](../CLAUDE.md) | Project overview, build dependencies, architecture (4-stage pipeline) |
| [`AGENTS.md`](../AGENTS.md) | Bilingual project guide, type system, integration |
| [`DEV_STATUS.md`](../DEV_STATUS.md) | Engineering status per family |
| [`docs/algorithms/RegionSolverAlgorithm.md`](./algorithms/RegionSolverAlgorithm.md) | Stage 1 algorithm: MMA vs C++ comparison |
| [`docs/algorithms/ReconstructAlgorithm.md`](./algorithms/ReconstructAlgorithm.md) | Stage 3 algorithm: relation reconstruction theory |
| [`docs/algorithms/LayerRecursion_Algorithm.md`](./algorithms/LayerRecursion_Algorithm.md) | Stage 2 algorithm: layer recursion expansion |
| [`docs/components/RelationSolver_ComponentGuide.md`](./components/RelationSolver_ComponentGuide.md) | RelationSolver full API reference |
| [`docs/components/RelationSolver_QuickReference.md`](./components/RelationSolver_QuickReference.md) | RelationSolver quick lookup |
| [`docs/components/RelationSolver_Documentation_Hub.md`](./components/RelationSolver_Documentation_Hub.md) | Central documentation hub for RelationSolver |
| [`docs/components/AnsatzModes.md`](./components/AnsatzModes.md) | Multi-mode ansatz guide with `--mode` CLI reference |
| [`docs/plans/RegionSolver-Debug-Plan.md`](./plans/RegionSolver-Debug-Plan.md) | C++ RegionSolver debug plan |
| [`workspace/shared/verify-docs/Test-Relation.md`](../workspace/shared/verify-docs/Test-Relation.md) | Relation verification (4 methods) |
| [`workspace/shared/verify-docs/Test-Expand.md`](../workspace/shared/verify-docs/Test-Expand.md) | Expansion consistency verification |

### Data flow summary

```
families/<fam>.json
  → [region-solve-cpp] family_generate / family_scheduler.py
  → data/IBPMat_<fam>.bin + data/RingData_<fam>.bin
  → [relation-test] test_relationFF
     ├── (internal) batchProcessRecursion → ExpansionCache + ExpansionMMA_*.m + Compare-CPPResult-*.m
     └── (internal) reconstructAllRelations → AllRelations_<fam>_k<N>.m + RelationMeta_<fam>.m
  → [verify] VerifyRelation.wl  →  pass/fail per (lev,deg) config
```
