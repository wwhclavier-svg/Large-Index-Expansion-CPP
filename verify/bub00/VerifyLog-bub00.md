# Verify Log: bub00
**Date:** 2026-05-02 19:17:19
**Modulus:** 179424673
**Description:** 1-loop bubble, msq=0

---

## 1. Characteristic Equation: Minimal Associated Primes

### Summary

| Sector | #Regions | Total #Solutions |
|--------|:--------:|:----------------:|
| {1, 1} | 1 | 1 |

### Sector {1, 1}


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{59808223, 59808223}` |
| **#Solutions** | **1** |
| Solution | `A[2] -> 59808223, A[1] -> 59808223` |
| Min.Poly. | `(none -- rational)` |

## 2. Determined Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `incre` | 2 | growth per order |
| `nimax` | 9 | solution space dim |
| `ne` | 2 | #external variables |
| `nb` | 1 | monomial basis dim |

## 3. Computation Time

| Step | MMA (s) | C++ (s) |
|------|:------:|:------:|
| Family Definition | 0.556 | -- |
| Region Solving | 1.128 | -- |
| Expansion | 0.419 | 0.000204 |
| **Total** | **2.103** | **0.000204** |

## 4. Series Expansion: Orders 0 and 1

### Order k=0

| | Expression |
|---|-----------|
| MMA | `1` |
| C++ | `1` |
| Diff | `0` **[MATCH]** |

### Order k=1

| | Expression |
|---|-----------|
| MMA | `105689379*v1 + 89712336*v1^2 + 105689379*v2 + v1*v2 + 89712336*v2^2` |
| C++ | `105689379*v1 + 89712336*v1^2 + 105689379*v2 + v1*v2 + 89712336*v2^2` |
| Diff | `0` **[MATCH]** |

## 5. Full Comparison Summary

| k | MMA | C++ | Diff | Result |
|:--|------|------|------|:------:|
| 0 | `1` | `1` | `0` | MATCH |
| 1 | `105689379*v1 + 89712336*v1^2 + 105689379*v2 + ... (5 terms)` | `105689379*v1 + 89712336*v1^2 + 105689379*v2 + ... (5 terms)` | `0` | MATCH |
| 2 | `81536841*v1 + 112499694*v1^2 + 21573929*v1^3 + ... (11 terms)` | `81536841*v1 + 112499694*v1^2 + 21573929*v1^3 + ... (11 terms)` | `0` | MATCH |
| 3 | `132236424*v1 + 26494910*v1^2 + 37148041*v1^3 + ... (17 terms)` | `132236424*v1 + 26494910*v1^2 + 37148041*v1^3 + ... (17 terms)` | `0` | MATCH |
| 4 | `129267832*v1 + 25376144*v1^2 + 157513150*v1^3 + ... (23 terms)` | `129267832*v1 + 25376144*v1^2 + 157513150*v1^3 + ... (23 terms)` | `0` | MATCH |

**Verdict: [PASS] -- All 5 orders match.**

