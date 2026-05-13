# VerifyUtility

C++ vs Mathematica cross-validation scripts for IBP expansion coefficients and LIE relation reconstruction.

**Size**: ~22 `.wl` scripts.

## Directory Layout

```
verify/VerifyUtility/

  Core Libraries (8):
    LIECoreAlgebra.wl       LIEExpand.wl      LIEFamilyDefine.wl
    LIEReconstruct.wl       LIERegions.wl     LIEUtility.wl
    LIEWorkflow.wl          SingularInterface.wl

  Export/Bridge (2):
    ExportBinary_IBPMatrix.wl         M2Kira.wl
    KiraRuleLoader.wl

  Expand Verification Pipeline (4):
    VerifyExpand-Prepare.wl           Step 1: Generate .bin + RegionInfo
    VerifyExpand-MMAExpand.wl         Step 2: MMA expansion
    VerifyExpand-Compare.wl           Step 3: C++ vs MMA diff
    VerifyExpand-SeriesVerify.wl      Independent series substitution

  Relation Verification Modules (4):
    Verify-Blade.wl                   Blade IBP reduction
    Verify-Series.wl                  Expansion substitution
    Verify-Kira.wl                    Kira cross-validation
    Verify-MMACompare.wl              MMA comparison (placeholder)

  Auxiliary (3):
    Compare-VerifyLog.wl              Markdown log generator
    Debug-DumpFamilyGenerate-v2.wl    Pipeline debugging
```

## Quick Start: Expand Verification

```bash
cd /home/ykm/Large-Index-Expansion-CPP/verify/VerifyUtility

# Step 1: Generate binary data
wolframscript -file VerifyExpand-Prepare.wl <famname>

# Step 2: C++ expansion (from project root)
cd ../.. && ./build/test_expandFF <famname>

# Step 3: MMA expansion
cd verify/VerifyUtility
wolframscript -file VerifyExpand-MMAExpand.wl <famname>

# Step 4: Compare + log (from any directory, paths relative to script)
wolframscript -file VerifyExpand-Compare.wl <famname>
```

## Quick Start: Relation Verification

```bash
# From project root — unified entry point
cd /home/ykm/Large-Index-Expansion-CPP
./build/test_relationFF <famname> <order> <lev_min> <lev_max> <deg_max>
wolframscript -file VerifyRelation.wl <famname>
wolframscript -file VerifyRelation.wl <famname> skip BladeVerify

# Individual modules
cd verify/VerifyUtility
wolframscript -file Verify-Blade.wl    <famname>
wolframscript -file Verify-Series.wl   <famname>
wolframscript -file Verify-Kira.wl     <famname>
```

## Prerequisites

| Dependency | Needed By |
|------------|-----------|
| FamilyDatabase.wl | All scripts |
| Blade.wl + FiniteFlow | Verify-Blade.wl |
| Kira binary | Verify-Kira.wl |
| Singular 4.3.2 | VerifyExpand-Prepare.wl |
| test_relationFF (C++) | Generates AllRelations_*.m + RelationMeta_*.m |

## Data Flow

```
C++ tests/                                        MMA scripts/
  test_expandFF                                     VerifyExpand-Prepare.wl
    → Compare-CPPResult-*.m        ──compare──        → RegionInfo, .bin, checkpoint
  test_relationFF                                   VerifyExpand-MMAExpand.wl
    → AllRelations_*.m             ──verify──         → MMA expansions
    → RelationMeta_*.m                               VerifyExpand-Compare.wl
    → ExpansionMMA_*.m               Verify-Kira.wl (needs Kira rules)
      Verify-Series.wl               Verify-MMACompare.wl (placeholder)
      Verify-Blade.wl
```

## Key Files

| File | Location | Purpose |
|------|----------|---------|
| VerifyRelation.wl | project root | Relation verification unified entry |
| Verify-Blade.wl | VerifyUtility/ | Blade IBP reduction |
| Verify-Series.wl | VerifyUtility/ | Expansion substitution |
| Verify-Kira.wl | VerifyUtility/ | Kira cross-validation |
| Verify-MMACompare.wl | VerifyUtility/ | MMA comparison (placeholder) |
| VerifyExpand-*.wl | VerifyUtility/ | Expand pipeline (4 scripts) |
| KiraRuleLoader.wl | VerifyUtility/ | Kira rule loader library |
| Test-Expand.md | verify/docs/ | Expand verification docs |
| Test-Relation.md | verify/docs/ | Relation verification docs |
