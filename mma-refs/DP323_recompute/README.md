# DP323 Cascading-Timeout Sector Recompute

> **Note**: This directory has been refactored to use `ScheduledRegionSolver` (see `verify/VerifyUtility/ScheduledRegionSolver.wl`). The original `run.wl` has been simplified to a thin wrapper around the generic scheduler.

## Directory Structure

```
verify/DP323_recompute/
├── run.wl                 ← Thin wrapper calling ScheduledRegionSolver
├── README.md              ← This document
├── cache/
│   └── status.wdx         ← Per-sector state tracking (WDX format)
├── output/
│   ├── IBPMat_DP323.bin         (159MB)  ✓ Generated
│   ├── RingData_DP323.bin       (2.1MB)  ✓ Generated
│   └── PrepareCheckpoint-DP323.wdx  (43MB)  ✓
├── SingularTempFile/      ← Singular runtime temp files
└── archive/               ← Historical scripts (cross-family tests, debug scripts, logs)
    ├── compare/           ← Cross-family comparison scripts
    ├── debug/             ← One-off debug/inspect scripts
    ├── export/            ← Export scripts for other families
    └── tests/             ← Test scripts
```

## Usage

### Run

```bash
cd verify/DP323_recompute
wolframscript -file run.wl
```

### Resume after interruption

`ScheduledRegionSolver` automatically saves state to `cache/status.wdx` atomically. Safe to `kill` and restart — the scheduler will resume from the last checkpoint.

### Instance lock

A file lock (`cache/instance.lock`) prevents multiple concurrent instances. If a stale lock is detected (dead PID), it is automatically removed.

## Architecture

```
run.wl (thin wrapper)
  ↓
ScheduledRegionSolver`ScheduledRegionSolve["DP323", config, scheduleConfig]
  ↓
┌─────────────────────────────────────────┐
│  InstanceLock      — prevents concurrent runs    │
│  SectorStateManager — atomic state persistence   │
│  TieredScheduler    — 1200s → 12000s → failed    │
│  RegionComputeEngine — Singular with timeout      │
└─────────────────────────────────────────┘
```

## Timeout Architecture

Singular subprocess timeout is now **properly propagated** through the full call chain:

```
run.wl: ComputeSector[..., timeout]
  → regionsBySectors[..., "Timeout" -> timeout]
    → expRegSolve2[..., "Timeout" -> timeout]
      → SingularGroebnerBasisDefault[..., "SingularTimeout" -> timeout]
        → SingularRun[..., "SingularTimeout" -> timeout]
          → timeout --signal=KILL {timeout}s Singular
```

Previously, timeout was set via a non-existent variable (`LIECoreAlgebra`Private`$SingularTimeout`) and had no effect.

## Output Files

| File | Description | Size |
|------|-------------|------|
| `output/IBPMat_DP323.bin` | IBP matrix binary | ~159MB |
| `output/RingData_DP323.bin` | Ring data binary | ~2.1MB |
| `output/PrepareCheckpoint-DP323.wdx` | Full checkpoint (with region data) | ~43MB |

Byte-level comparison with `verify/DP323/` original files verified as IDENTICAL.

## Dependencies

- Mathematica (WolframScript)
- `LIEWorkflow.wl` (loads `ScheduledRegionSolver.wl`)
- Singular (for Groebner basis computation)
