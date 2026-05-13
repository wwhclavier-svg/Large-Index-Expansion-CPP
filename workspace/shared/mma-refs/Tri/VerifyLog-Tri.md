# Verify Log: Tri
**Date:** 2026-05-10 10:21:48
**Modulus:** 179424673
**Description:** 1-loop triangle, 3 props, s1=3,s2=8,s3=0

---

## 1. Characteristic Equation: Minimal Associated Primes

### Summary

| Sector | #Regions | Total #Solutions |
|--------|:--------:|:----------------:|
| {1, 0, 1} | 1 | 1 |
| {1, 1, 0} | 1 | 1 |
| {1, 1, 1} | 2 | 2 |

### Sector {1, 0, 1}


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{89712336, 2, 89712336}` |
| **#Solutions** | **1** |
| Solution | `A[3] -> 89712336, A[2] -> 2, A[1] -> 89712336` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 1, 0}


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{59808223, 59808223, 124217081}` |
| **#Solutions** | **1** |
| Solution | `A[3] -> 124217081, A[2] -> 59808223, A[1] -> 59808223` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 1, 1}


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{44856166, 53827401, 148025355}` |
| **#Solutions** | **1** |
| Solution | `A[3] -> 148025355, A[2] -> 53827401, A[1] -> 44856166` |
| Min.Poly. | `(none -- rational)` |


| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{89712336, 89712337, 134568504}` |
| **#Solutions** | **1** |
| Solution | `A[3] -> 134568504, A[2] -> 89712337, A[1] -> 89712336` |
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
| Family Definition | 0.749 | -- |
| Region Solving | 2.575 | -- |
| Expansion | 3.038 | ? |
| **Total** | **6.362** | **?** |

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
| MMA | `12205760*v1 + 113208424*v1^2 + 7323456*v2 + 34176129*v1*v2 + 49128184*v2^2 + 10374896*v3 + 98256369*v1*v3 + 46992176*v2*v3 + 17088064*v3^2` |
| C++ | `12205760*v1 + 113208424*v1^2 + 7323456*v2 + 34176129*v1*v2 + 49128184*v2^2 + 10374896*v3 + 98256369*v1*v3 + 46992176*v2*v3 + 17088064*v3^2` |
| Diff | `0` **[MATCH]** |

## 5. Full Comparison Summary

| k | MMA | C++ | Diff | Result |
|:--|------|------|------|:------:|
| 0 | `1` | `1` | `0` | MATCH |
| 1 | `12205760*v1 + 113208424*v1^2 + 7323456*v2 + ... (9 terms)` | `12205760*v1 + 113208424*v1^2 + 7323456*v2 + ... (9 terms)` | `0` | MATCH |
| 2 | `133159372*v1 + 162149786*v1^2 + 19730224*v1^3 + ... (34 terms)` | `133159372*v1 + 162149786*v1^2 + 19730224*v1^3 + ... (34 terms)` | `0` | MATCH |

**Verdict: [PASS] -- All 3 orders match.**

