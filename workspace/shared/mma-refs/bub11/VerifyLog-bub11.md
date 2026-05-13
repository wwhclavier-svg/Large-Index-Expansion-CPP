# Verify Log: bub11
**Date:** 2026-05-10 10:21:32
**Modulus:** 179424673
**Description:** 1-loop bubble, msq=1

---

## 1. Characteristic Equation: Minimal Associated Primes

### Summary

| Sector | #Regions | Total #Solutions |
|--------|:--------:|:----------------:|
| {0, 1} | 1 | 1 |
| {1, 0} | 1 | 1 |
| {1, 1} | 3 | 3 |

### Sector {0, 1}


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{89712336, 1}` |
| **#Solutions** | **1** |
| Solution | `A[2] -> 1, A[1] -> 89712336` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 0}


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{1, 89712336}` |
| **#Solutions** | **1** |
| Solution | `A[2] -> 89712336, A[1] -> 1` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 1}


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{4513913, 174910761}` |
| **#Solutions** | **1** |
| Solution | `A[2] -> 174910761, A[1] -> 4513913` |
| Min.Poly. | `(none -- rational)` |


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{174910761, 4513913}` |
| **#Solutions** | **1** |
| Solution | `A[2] -> 4513913, A[1] -> 174910761` |
| Min.Poly. | `(none -- rational)` |


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{4, 4}` |
| **#Solutions** | **1** |
| Solution | `A[2] -> 4, A[1] -> 4` |
| Min.Poly. | `(none -- rational)` |

## 2. Determined Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `incre` | ? | growth per order |
| `nimax` | ? | solution space dim |
| `ne` | ? | #external variables |
| `nb` | ? | monomial basis dim |

## 3. Computation Time

| Step | MMA (s) | C++ (s) |
|------|:------:|:------:|
| Family Definition | 0.750 | -- |
| Region Solving | 2.272 | -- |
| Expansion | 1.194 | ? |
| **Total** | **4.216** | **?** |

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
| MMA | `164472617*v1 + 38448144*v1^2 + 164472617*v2 + 102528385*v1*v2 + 38448144*v2^2` |
| C++ | `164472617*v1 + 38448144*v1^2 + 164472617*v2 + 102528385*v1*v2 + 38448144*v2^2` |
| Diff | `0` **[MATCH]** |

## 5. Full Comparison Summary

| k | MMA | C++ | Diff | Result |
|:--|------|------|------|:------:|
| 0 | `1` | `1` | `0` | MATCH |
| 1 | `164472617*v1 + 38448144*v1^2 + 164472617*v2 + ... (5 terms)` | `164472617*v1 + 38448144*v1^2 + 164472617*v2 + ... (5 terms)` | `0` | MATCH |
| 2 | `120163165*v1 + 80643384*v1^2 + 163404613*v1^3 + ... (14 terms)` | `120163165*v1 + 80643384*v1^2 + 163404613*v1^3 + ... (14 terms)` | `0` | MATCH |

**Verdict: [PASS] -- All 3 orders match.**

