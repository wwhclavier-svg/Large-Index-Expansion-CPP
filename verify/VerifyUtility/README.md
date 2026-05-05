# VerifyUtility

Verification toolchain for C++ vs Mathematica cross-validation of IBP expansion coefficients and LIE relation reconstruction.

**Size**: ~292K total (20 `.wl` scripts).

## Directory Layout

```
verify/VerifyUtility/
├── README.md                        # This file
│
├── M2Kira.wl                        # Mathematica↔Kira interface (external dependency)
│
├── lib/                             # (flat — core libraries load each other via Get[])
│   ├── LIEWorkflow.wl               # LIE workflow orchestrator
│   ├── LIEReconstruct.wl            # Relation reconstruction algorithm
│   ├── LIEExpand.wl                 # Series expansion engine
│   ├── LIEFamilyDefine.wl           # IBP family definition
│   ├── LIERegions.wl                # Region/chracteristic equation solving
│   ├── LIECoreAlgebra.wl            # Coordinate ring & algebra core
│   ├── LIEUtility.wl                # Shared utilities
│   ├── SingularInterface.wl         # SINGULAR CAS interface
│   ├── KiraRuleLoader.wl            # Kira reduction rule loader
│   └── ExportBinary_IBPMatrix.wl    # Binary .bin IBP matrix export
│
├── VerifyExpand-Prepare.wl          # [Expand] Step 1: Generate .bin + checkpoint
├── VerifyExpand-MMAExpand.wl        # [Expand] Step 2: MMA series expansion
├── VerifyExpand-Compare.wl          # [Expand] Step 3: C++ vs MMA diff + VerifyLog
├── VerifyExpand-SeriesVerify.wl     # [Expand] Optional: IBP series substitution
│
├── Compare-FamilyGenerate.wl        # [Relation] Step 1: Define family + export .bin
├── Compare-Expand.wl                # [Relation] Step 2: LIE expansion
├── Compare-Results.wl               # [Relation] Step 3: C++ vs MMA comparison
├── Compare-VerifyLog.wl             # [Relation] Step 4: Generate verification log
└── Compare-Reconstruct-bub00.wl     # [Relation] Full pipeline (bub00): reconstruct + Kira verify
```

## Quick Start

### Test-Expand: C++ vs MMA coefficient comparison

```bash
cd verify/VerifyUtility

# Step 1: Prepare — generate .bin + checkpoint
wolframscript -file VerifyExpand-Prepare.wl bub00

# Step 2: C++ expansion
cd ../..
./build/test_expandFF bub00

# Step 3: MMA expansion
cd verify/VerifyUtility
wolframscript -file VerifyExpand-MMAExpand.wl bub00

# Step 4: Compare + generate VerifyLog
wolframscript -file VerifyExpand-Compare.wl bub00
```

### Test-Relation: LIE relation reconstruction + Kira verification

```bash
cd verify/VerifyUtility

# Full pipeline (define → expand → reconstruct → Kira verify)
wolframscript -file Compare-Reconstruct-bub00.wl

# C++ relation reconstruction
cd ../..
./build/test_relationFF bub00 4 2 2
```

### Other families

All scripts accept `<famname>` as first argument. Example:

```bash
wolframscript -file VerifyExpand-Prepare.wl Box
wolframscript -file VerifyExpand-Compare.wl Tri
wolframscript -file Compare-Expand.wl SR 5
```

List available families:

```bash
wolframscript -file Compare-FamilyGenerate.wl
```

## Paths

Scripts auto-detect paths relative to their own location (`$InputFileName`). Key path variables:

| Variable | Where | Default |
|----------|-------|---------|
| `$CPPVerifyRoot` | `VerifyExpand-Prepare.wl` | `ParentDirectory[$InputFileName]` = `verify/` |
| `$KiraBaseDir` | `Compare-Reconstruct-bub00.wl` | `/root/kira_tests/bub00_new/bub00` |

`$KiraBaseDir` can be overridden in the calling script before `Get[]`:

```mathematica
$KiraBaseDir = "/path/to/my/kira/project";
Get["VerifyUtility/Compare-Reconstruct-bub00.wl"];
```

Result files are written to `verify/<fam>/` (which is `$CPPVerifyRoot <> famname <> "/"`).

## External Dependencies

These are NOT included in this directory and must be installed separately:

| Dependency | Purpose | Typical location |
|------------|---------|-----------------|
| Mathematica / wolframscript | Run `.wl` scripts | System path |
| SINGULAR | Grobner basis computation | System path |
| Kira | IBP reduction engine | `/root/kira/` |
| Kira reduction rules | Pre-computed IBP rules per family | `/root/kira_tests/<fam>/` |
| FireFly | Finite field arithmetic (C++ side) | `/root/firefly-2.0.3/` |

## Script Roles

### Core Libraries (loaded via `Get[]`, not run directly)

| Script | Provides |
|--------|----------|
| `LIEWorkflow.wl` | `LIEDefineFamily`, `LIESolveRegions`, `LIEExpandSeries`, `LIEGetRelations` |
| `LIEReconstruct.wl` | Low-level relation reconstruction functions |
| `KiraRuleLoader.wl` | `LoadKiraRules[]` — loads Kira explicit + trivial sector + zero-index rules |
| `ExportBinary_IBPMatrix.wl` | `ExportBinaryIBPMatrix`, `ExportBinaryRingData` |
| `M2Kira.wl` | `KiraInputGenerate`, `KiraRun`, `KiraRelationImport` |

### Pipeline Scripts (run via `wolframscript -file`)

| Script | Input | Output |
|--------|-------|--------|
| `VerifyExpand-Prepare.wl <fam>` | FamilyDatabase config | `IBPMat_*.bin`, `RingData_*.bin`, `PrepareCheckpoint-*.wdx`, timing `.m` |
| `VerifyExpand-MMAExpand.wl <fam> [order]` | Checkpoint from Prepare | `VerifyExpansion-MMAExpansion.m` |
| `VerifyExpand-Compare.wl <fam>` | C++ `.m` + MMA `.m` | `VerifyLog-*.md` |
| `Compare-FamilyGenerate.wl <fam>` | FamilyDatabase config | `.bin` files (legacy location) |
| `Compare-Expand.wl <fam> [order]` | `.bin` files | `Compare-MMAResult-*.m` |
| `Compare-Results.wl <fam>` | C++ + MMA result `.m` | Terminal diff output |
| `Compare-VerifyLog.wl <fam>` | All intermediate files | `VerifyLog-*.md` |
| `Compare-Reconstruct-bub00.wl` | FamilyDatabase + Kira rules | Relations + Kira verification |

### Key Difference: `VerifyExpand-*` vs `Compare-*`

- **VerifyExpand-*** — Unified 3-step pipeline (Prepare → MMAExpand → Compare). Uses checkpoint to avoid re-computation. Recommended for routine verification.
- **Compare-*** — Legacy per-step scripts with individual file search. Still functional for one-off comparisons.

## Adding a New Family to the Workflow

1. Register in `verify/FamilyDatabase/FamilyDatabase.wl` (see [FamilyDatabase/README.md](../FamilyDatabase/README.md))
2. Run Prepare: `wolframscript -file VerifyExpand-Prepare.wl <newfam>`
3. C++ test: `./build/test_expandFF <newfam>`
4. MMA expand: `wolframscript -file VerifyExpand-MMAExpand.wl <newfam>`
5. Compare: `wolframscript -file VerifyExpand-Compare.wl <newfam>`

No new scripts needed — the existing ones read config from FamilyDatabase.

## Troubleshooting

### "No matrices loaded" (C++ side)
Check that `.bin` files exist in the working directory or `verify/<fam>/`.

### "FamilyDatabase not found"
The scripts use `ParentDirectory[$InputFileName]` to find `verify/FamilyDatabase/FamilyDatabase.wl`. Ensure the directory structure is intact.

### "Kira rules file not found"
Check `$KiraBaseDir` points to a directory containing `results/<fam>/kira_integrals.m`. For bub00, default is `/root/kira_tests/bub00_new/bub00`.

### Path issues when running from a different directory
Scripts always `SetDirectory[DirectoryName[$InputFileName]]` at startup, so they can be invoked from anywhere:

```bash
wolframscript -file /absolute/path/to/verify/VerifyUtility/VerifyExpand-Prepare.wl bub00
```
