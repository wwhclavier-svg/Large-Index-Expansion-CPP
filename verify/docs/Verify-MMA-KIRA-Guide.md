# Verify-LIEMMA-to-KiraReduce-Guide

Guide to verify LIE (Large Index Expansion) MMA relations using Kira reduction rules.

---

## Parameters

| Parameter | Value |
|-----------|-------|
| Family | `bub00` |
| s | 3 |
| msq | 0 |
| d | 1/3 |
| Modulus | 179424673 (Prime[10000000]) |
| NE | 2 |

---

## Workflow

```
Step 1: Define integral family
        ↓
Step 2: LIE series expansion
        ↓
Step 3: Relation reconstruction
        ↓
Step 4: Generate Kira reduction rules
        ↓
Step 5: Verify MMA relations via Kira
```

---

## Step 1: Define Integral Family

**File:** `Compare-FamilyGenerate-bub00.wl`

```mathematica
bub00Config = <|
    "Propagators" -> ({-k1^2 + msq, -(k1 - p1)^2 + msq} /. numericRules),
    "LoopMomenta" -> {k1},
    "ExternalMomenta" -> {p1},
    "KinematicRules" -> ({p1^2 -> s} /. numericRules),
    "TopSector" -> {1, 1},
    "Numeric" -> numericRules,
    "Modulus" -> modulus
|>;

data = LIEDefineFamily[
    bub00Config["Propagators"],
    bub00Config["LoopMomenta"],
    bub00Config["ExternalMomenta"],
    bub00Config["KinematicRules"],
    bub00Config["TopSector"],
    "Numeric" -> bub00Config["Numeric"],
    Modulus -> modulus
];
```

---

## Step 2: LIE Expansion

**File:** `Compare-Expand-bub00.wl`

```mathematica
data = LIESolveRegions[data, Verbose -> False];
data = LIEExpandSeries[data,
    "Order" -> 4,
    Modulus -> modulus,
    "Increment" -> 2,
    "LayerByLayer" -> True,
    Verbose -> False
];
```

---

## Step 3: Relation Reconstruction

**File:** `Compare-Reconstruct-bub00.wl`

```mathematica
rank = 2;  (* max level |alpha| *)
maxDeg = 2;  (* max degree |beta| *)
data = LIEGetRelations[data, Verbose -> False, "MaxCoefDeg" -> maxDeg];
relations = data["Relations", "Relations"];
```

**Output:** `relations[[lev+1, deg+1]]` — LIE g-form relations

**Example relation (lev=0, deg=1):**
```
-119616448*g[-1+v1, v2] + v1*g[-1+v1, v2] + 59808226*g[v1, -1+v2]
+ 179424670*v1*g[v1, -1+v2] + 179424671*v2*g[v1, -1+v2]
+ 3*g[v1, v2] + 3*v1*g[v1, v2] + 179424667*v2*g[v1, v2] = 0
```

**MMA Relations:**

| Level | Degree | Expression |
|-------|--------|------------|
| lev=0 | deg=1 | `-119616448*g[-1+v1,v2] + v1*g[-1+v1,v2] + 59808226*g[v1,-1+v2] + 179424670*v1*g[v1,-1+v2] + 179424671*v2*g[v1,-1+v2] + 3*g[v1,v2] + 3*v1*g[v1,v2] + 179424667*v2*g[v1,v2] = 0` |
| lev=0 | deg=2 | `169456636*g[v1,-1+v2] + 89712335*v1*g[v1,-1+v2] + v1^2*g[v1,-1+v2] + 89712335*v2*g[v1,-1+v2] - 179424671*v1*v2*g[v1,-1+v2] + v2^2*g[v1,-1+v2] + 89712337*g[v1,v2] + 89712333*v2*g[v1,v2] + 3*v2^2*g[v1,v2] = 0` |
| lev=1 | deg=2 | `g[-2+v1,v2] - 57270444*g[-1+v1,v2] - 72594664*g[v1,-2+v2] - 122472592*v1*g[v1,-2+v2] - 170087892*v2*g[v1,-2+v2] + 116595564*g[v1,-1+v2] + 11118049*v1*g[v1,-1+v2] + 175140458*v2*g[v1,-1+v2] + 16984649*g[v1,v2] + 146678964*v1*g[v1,v2] + 67867852*v2*g[v1,v2] = 0` |
| lev=2 | deg=2 | `-168270948*g[-1+v1,v2] + v1*g[v1,-3+v2] - 88181728*g[v1,-2+v2] - 35120916*v1*g[v1,-2+v2] - 134422015*v2*g[v1,-2+v2] + 43849420*g[v1,-1+v2] + 151461548*v1*g[v1,-1+v2] + 149171382*v2*g[v1,-1+v2] + 70877199*g[v1,v2] + 160095842*v1*g[v1,v2] + 141122181*v2*g[v1,v2] = 0` |

---

## Step 4: Generate Kira Reduction Rules

### 4.1 Design Kira Input

**Step 1: Analyze indices appearing in MMA relations after ν substitution**

After substituting ν values into MMA relations and converting `g → j`, the following indices appear:

| Index Type | Examples | Count |
|-----------|----------|-------|
| Negative | `j[-1,4]`, `j[-1,6]`, `j[4,-1]`, `j[6,-1]` | 4 |
| Zero components | `j[0,4]`, `j[0,5]`, `j[0,6]`, `j[4,0]`, `j[5,0]`, `j[6,0]` | 6 |
| Positive (covered) | `j[1,2]`, `j[2,1]`, ..., `j[5,1]`, `j[4,2]` | 12 |
| **Missing positive** | `j[0,4]`, `j[0,5]`, `j[0,6]`, `j[4,0]`, `j[5,0]`, `j[6,0]` | 6 |

**Issue:** RMax=6, SMax=2 only covers `j[a,b]` with `a≥1 AND b≥1`. Indices with `a=0` or `b=0` belong to **trivial sectors** (`{1,0}` or `{0,1}`) and are not reduced by Kira — they must be set to zero via trivial sector rules.

**Step 2: Determine required RMax/SMax**

For the non-trivial sector `{1,1}`, RMax=6 covers all indices up to `j[6,6]`:

```mathematica
<< "/root/M2Kira.wl"

fam = "bub00";
props = {-k1^2 + msq, -(k1 - p1)^2 + msq};
loopMom = {k1};
kin = <|
  "Incoming" -> {p1},
  "Outgoing" -> {},
  "MomentumConservation" -> {},
  "Invariants" -> {{s, 2}},
  "ScalarRules" -> {{p1, p1} -> s},
  "SymbolToReplaceByOne" -> s
|>;

(* RMax=6 suffices for non-trivial sector {1,1} *)
(* SMax does NOT need to cover a=0 or b=0 — those are trivial sectors *)
workDir = KiraInputGenerate[
  fam, props, loopMom, kin,
  "KiraWorkingDir" -> "/root/kira_tests/bub00",
  "RMax" -> 6,
  "SMax" -> 2,
  "Sector" -> {1, 1},
  Verbose -> True
];
```

**Key insight:** Do NOT increase `SMax` to cover `a=0` or `b=0` indices. Integrals like `j[4,0]` belong to **trivial sector `{1,0}`** and should be set to zero by trivial sector rules, not by Kira reduction. The `KiraRuleLoader.wl` module automatically extracts trivial sector information from Kira's `sectormappings` output and generates the corresponding zero-rules.

### 4.2 Run Kira

```bash
cd /root/kira_tests/bub00_new/bub00
/root/kira/builddir/src/kira/kira jobs.yaml --parallel=physical
```

**Runtime:** 6.4s
**Master integrals:** 8 — `bub00[1,1], [1,7], [2,6], [3,5], [4,4], [5,3], [6,2], [7,1]`

### 4.3 Kira Output — Reduction Rules

```
j[2,1] -> (d-3) j[1,1]
j[1,2] -> (d-3) j[1,1]
j[3,1] -> ((d^2-7d+12)/2) j[1,1]
j[1,3] -> ((d^2-7d+12)/2) j[1,1]
j[2,2] -> (d^2-9d+18) j[1,1]
j[4,1] -> ((d^3-12d^2+47d-60)/6) j[1,1]
j[1,4] -> ((d^3-12d^2+47d-60)/6) j[1,1]
j[3,2] -> ((d^3-16d^2+79d-120)/2) j[1,1]
j[2,3] -> ((d^3-16d^2+79d-120)/2) j[1,1]
j[5,1] -> ((d^4-18d^3+119d^2-342d+360)/24) j[1,1]
j[1,5] -> ((d^4-18d^3+119d^2-342d+360)/24) j[1,1]
j[4,2] -> ((d^4-24d^3+203d^2-720d+900)/6) j[1,1]
j[2,4] -> ((d^4-24d^3+203d^2-720d+900)/6) j[1,1]
j[3,3] -> ((d^4-26d^3+239d^2-910d+1200)/4) j[1,1]
j[6,1] -> ((d^5-25d^4+245d^3-1175d^2+2754d-2520)/120) j[1,1]
j[1,6] -> ((d^5-25d^4+245d^3-1175d^2+2754d-2520)/120) j[1,1]
j[5,2] -> ((d^5-33d^4+413d^3-2463d^2+7002d-7560)/24) j[1,1]
j[2,5] -> ((d^5-33d^4+413d^3-2463d^2+7002d-7560)/24) j[1,1]
j[4,3] -> ((d^5-37d^4+521d^3-3467d^2+10830d-12600)/12) j[1,1]
j[3,4] -> ((d^5-37d^4+521d^3-3467d^2+10830d-12600)/12) j[1,1]
```

**Master Integral:** `j[1,1]` — all other integrals reduce to it

### 4.4 Key Configuration

| Parameter | Description | How to Determine |
|-----------|-------------|------------------|
| `RMax` | Max dot product | Max index in MMA relations (max ν + max offset). For ν∈[1,6] and offsets up to 3, RMax=6 suffices |
| `SMax` | Max rank | Must cover indices with `a=0` or `b=0` (e.g., `j[0,4]`, `j[4,0]`). Set SMax ≥ RMax to ensure full coverage |
| `Sector` | Top-level sector | From LIE config TopSector |

### 4.5 Index Coverage Summary

**Zero-Index Rule (Exact):**

All integrals with **all** indices non-positive are zero:

```mathematica
j[a__]/;!Or@@({a}/.{b_/;b>0->True,b_/;b<=0->False})->0
```

**Trivial Sector Rules:**

Integrals in trivial sectors (where only a subset of propagators are present) also evaluate to zero. For bub00 with 2 propagators, sectors `{1,0}` and `{0,1}` are trivial:

```mathematica
(* Sector {1,0}: first index > 0, second <= 0 *)
j[a_, b_] /; a > 0 && b <= 0 -> 0

(* Sector {0,1}: first index <= 0, second > 0 *)
j[a_, b_] /; a <= 0 && b > 0 -> 0
```

These rules are automatically generated by `KiraRuleLoader.wl` from Kira's `sectormappings/bub00/trivialsector` file.

| Index Pattern | Example | Sector | Reduction Rule | Result |
|-------------|---------|--------|----------------|--------|
| `a ≤ 0 AND b ≤ 0` | `j[-1,0]`, `j[0,0]`, `j[-2,-1]` | `{0,0}` | Zero-index rule | → 0 |
| `a > 0 AND b ≤ 0` | `j[4,0]`, `j[5,0]`, `j[4,-1]` | `{1,0}` | Trivial sector rule | → 0 |
| `a ≤ 0 AND b > 0` | `j[0,4]`, `j[0,5]`, `j[-1,4]` | `{0,1}` | Trivial sector rule | → 0 |
| `a ≥ 1 AND b ≥ 1` | `j[2,3]`, `j[3,2]` | `{1,1}` | Kira explicit rules | → `j[1,1]` |

---

## Step 5: Verify MMA Relations via Kira

### 5.1 Verification Function

```mathematica
(* Load Kira rules from file *)
kiraFile = "/root/kira_tests/bub00_new/bub00/results/bub00/kira_integrals.m";
kiraRaw = Get[kiraFile];
kiraRulesSymbolic = kiraRaw /. {bub00[a_, b_] :> j[a, b]};

(* Exact zero-index rule: all non-positive indices -> 0 *)
zeroIndexRule = j[a__]/;!Or@@({a}/.{b_/;b>0->True,b_/;b<=0->False})->0;

(* Combine and finite-field-ize *)
allRules = Join[kiraRulesSymbolic, {zeroIndexRule}] /. d -> 1/3;
allRules = PolynomialMod[allRules, modulus];

(* Kira reduction *)
kiraReduce[expr_] := Module[
  {result},
  result = FixedPoint[(# /. allRules) &, expr, 100];
  PolynomialMod[result, modulus]
];

(* ν sampling verification *)
verifyAtNu[rel_, nu1_, nu2_] := Module[
  {substituted, jForm, reduced},
  (* Substitute ν values *)
  substituted = Evaluate[rel /. {v1 -> nu1, v2 -> nu2}];
  (* Convert g -> j *)
  jForm = substituted /. {g[a_, b_] :> j[a, b]};
  (* Kira reduction *)
  reduced = kiraReduce[jForm];
  reduced === 0
];
```

### 5.2 Verification Script

```mathematica
(* Load M2Kira and generate rules *)
<< "/root/M2Kira.wl"
workDir = KiraInputGenerate["bub00", props, loopMom, kin,
  "KiraWorkingDir" -> "/root/kira_tests/bub00", "RMax" -> 6, "SMax" -> 2];
KiraRun["bub00", "KiraExecutable" -> "/root/kira/builddir/src/kira/kira"];
{kiraRules, kiraMasters} = KiraRelationImport["bub00", "KiraPath" -> workDir];

(* Verify MMA relations *)
seeds = {{3, 3}, {4, 1}, {1, 4}, {5, 2}, {2, 5}, {6, 1}, {1, 6}};
rel1 = -119616448*g[-1+v1,v2] + v1*g[-1+v1,v2] + ...;

Do[
  result = verifyAtNu[rel1, seeds[[i,1]], seeds[[i,2]]];
  Print["nu={", seeds[[i,1]], ",", seeds[[i,2]], "}: ",
        If[result, "PASS", "FAIL"]],
  {i, Length[seeds]}
]
```

### 5.3 Notes

1. **Kira rules must be complete**: If Kira rules do not cover some index combinations, verification will fail. Increase `RMax`/`SMax` and re-run Kira if needed.
2. **Master Integral**: In bub00 case, master is `bub00[1,1]`
3. **d parameter**: Kira rules contain symbol d; substitute d=1/3 during verification

### 5.4 Verification Status (2026-04-30)

#### The definitive verification strategy

**Key insight:** To avoid `INCONCLUSIVE` results from uncovered integrals (e.g. `j[4,0]`), **select test points ν such that every term in the relation produces j-indices fully covered by Kira rules.**

For a relation containing `g[v1 + a, v2 + b]`, substituting ν yields `j[ν1+a, ν2+b]`. To ensure Kira coverage:
1. Both indices must be > 0 (otherwise zero-index rule applies, trivializing the test)
2. The resulting index pair must be in Kira's target integral set

From the relations:
| Relation | g-offsets | Minimum ν requirement |
|----------|-----------|----------------------|
| Rel1 | `{-1,0}, {0,-1}` | `ν1≥2, ν2≥2` |
| Rel2 | `{0,-1}` | `ν1≥1, ν2≥2` |
| Rel3 | `{-2,0}, {-1,0}, {0,-2}` | `ν1≥3, ν2≥3` |
| Rel4 | `{-1,0}, {0,-3}, {0,-2}` | `ν1≥2, ν2≥4` |

**Global requirement: `ν1≥3, ν2≥4`** (and upper-bounded by Kira's RMax/SMax coverage).

This yields **6 valid test points**: `{3,4}, {3,5}, {3,6}, {4,4}, {4,5}, {5,4}`.

#### Results at Kira-covered points

| Relation | Result at 6 valid points | Conclusion |
|----------|-------------------------|------------|
| Rel1 | 0/6 PASS | **FAIL** (e.g. `ν={3,4}` → `13045393*j[1,1]`) |
| Rel2 | 0/6 PASS | **FAIL** (e.g. `ν={3,4}` → `7630983*j[1,1]`) |
| Rel3 | 0/6 PASS | **FAIL** (e.g. `ν={3,4}` → `137269204*j[1,1]`) |
| Rel4 | 0/6 PASS | **FAIL** (e.g. `ν={3,4}` → `167117323*j[1,1]`) |

**All failures are definitive `FAIL_RELATION`** — every term reduces to Kira masters with nonzero coefficients. There are no uncovered integrals to create ambiguity.

#### Summary of other methods

| Test Method | Result | Assessment |
|-------------|--------|------------|
| Method A (arbitrary ν) | Mixed FAIL/INCONCLUSIVE | Some points have uncovered integrals |
| Method B (symbolic offset) | 4/4 PASS | **Trivial** — equivalent to ν={0,0} where zero-index rule zeros everything |
| sd={1,1}−seeds | Mostly INCONCLUSIVE | Uncovered integrals like `j[1,0]` prevent judgment |

**Conclusion:** The MMA relations for d=1/3, s=3 are **not valid IBP identities** for general ν.

**Zero-index rule (exact):**
```mathematica
j[a__]/;!Or@@({a}/.{b_/;b>0->True,b_/;b<=0->False})->0
```

---

## Run Commands

```bash
cd /root/Large-Index-Expansion-MMA-Mini

# Steps 1-3: Full workflow (family -> expand -> reconstruct)
wolframscript -file Compare-Reconstruct-bub00.wl

# Step 4: Kira verification (standalone)
wolframscript -file Compare-KiraVerify-bub00.wl
```

---

## Key Files

| File | Purpose |
|------|---------|
| `Compare-FamilyGenerate-bub00.wl` | Define integral family |
| `Compare-Expand-bub00.wl` | LIE expansion |
| `Compare-Reconstruct-bub00.wl` | Relation reconstruction + Kira verification |
| `Compare-KiraVerify-bub00.wl` | Kira verification (standalone) |
| `/root/M2Kira.wl` | Kira interface |
| `/root/kira_tests/bub00_new/bub00/results/bub00/kira_integrals.m` | Kira reduction rules |
