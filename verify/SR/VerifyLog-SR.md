# Verify Log: SR
**Date:** 2026-05-02 00:41:35
**Modulus:** 179424673
**Description:** Sunrise 2L1P, massless, s=1

---

## 1. Characteristic Equation: Minimal Associated Primes

### Summary

| Sector | #Regions | Total #Solutions |
|--------|:--------:|:----------------:|
| {0, 1, 0, 1, 1} | 1 | 1 |
| {0, 1, 1, 1, 1} | 1 | 1 |
| {1, 0, 1, 0, 1} | 1 | 1 |
| {1, 0, 1, 1, 1} | 1 | 1 |
| {1, 1, 0, 1, 1} | 1 | 1 |
| {1, 1, 1, 0, 1} | 1 | 1 |

### Sector {0, 1, 0, 1, 1}

#### Primary Component (1)

| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`, `A[4], 1`, `A[5], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{44856166, 179424664, 44856166, 179424664, 179424664}` |
| **#Solutions** | **1** (= 1x1x1x1x1) |
| Solution | `A[5] -> 179424664, A[4] -> 179424664, A[3] -> 44856166, A[2] -> 179424664, A[1] -> 44856166` |
| Min.Poly. | `(none -- rational)` |

### Sector {0, 1, 1, 1, 1}

#### Primary Component (1)

| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`, `A[4], 1`, `A[5], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{129185762, 79744292, 19936073, 179424657, 79744292}` |
| **#Solutions** | **1** (= 1x1x1x1x1) |
| Solution | `A[5] -> 79744292, A[4] -> 179424657, A[3] -> 19936073, A[2] -> 79744292, A[1] -> 129185762` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 0, 1, 0, 1}

#### Primary Component (1)

| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`, `A[4], 1`, `A[5], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{179424664, 44856166, 179424664, 44856166, 179424664}` |
| **#Solutions** | **1** (= 1x1x1x1x1) |
| Solution | `A[5] -> 179424664, A[4] -> 44856166, A[3] -> 179424664, A[2] -> 44856166, A[1] -> 179424664` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 0, 1, 1, 1}

#### Primary Component (1)

| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`, `A[4], 1`, `A[5], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{79744292, 129185762, 179424657, 19936073, 79744292}` |
| **#Solutions** | **1** (= 1x1x1x1x1) |
| Solution | `A[5] -> 79744292, A[4] -> 19936073, A[3] -> 179424657, A[2] -> 129185762, A[1] -> 79744292` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 1, 0, 1, 1}

#### Primary Component (1)

| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`, `A[4], 1`, `A[5], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{19936073, 179424657, 129185762, 79744292, 79744292}` |
| **#Solutions** | **1** (= 1x1x1x1x1) |
| Solution | `A[5] -> 79744292, A[4] -> 79744292, A[3] -> 129185762, A[2] -> 179424657, A[1] -> 19936073` |
| Min.Poly. | `(none -- rational)` |

### Sector {1, 1, 1, 0, 1}

#### Primary Component (1)

| Field | Value |
|-------|-------|
| Leading Degree | `{{`A[1], 1`, `A[2], 1`, `A[3], 1`, `A[4], 1`, `A[5], 1`}}` |
| Generator (indep) | `{}` |
| Parametrized (dep) | `{179424657, 19936073, 79744292, 129185762, 79744292}` |
| **#Solutions** | **1** (= 1x1x1x1x1) |
| Solution | `A[5] -> 79744292, A[4] -> 129185762, A[3] -> 79744292, A[2] -> 19936073, A[1] -> 179424657` |
| Min.Poly. | `(none -- rational)` |

## 2. Determined Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `incre` | 2 | nu-polynomial expansion growth per order |
| `nimax` | 495 | Solution space dimension (independent 1/n series) |
| `ne` | 5 | #external variables (A_i) |
| `nb` | 1 | Monomial basis dimension (coordinate ring) |

### Per-Region Details

| Region | NE | NB | NIBP | Incre | Nimax |
|--------|:--:|:--:|:----:|:-----:|:-----:|
| 1 | 5 | 1 | 6 | 2 | 495 |
| 2 | 5 | 1 | 6 | 2 | 495 |
| 3 | 5 | 1 | 6 | 2 | 495 |
| 4 | 5 | 1 | 6 | 2 | 495 |
| 5 | 5 | 1 | 6 | 2 | 495 |
| 6 | 5 | 1 | 6 | 2 | 495 |

## 3. Computation Time

| Step | MMA (s) | C++ (s) |
|------|:------:|:------:|
| Family Definition | 0.502 | -- |
| Region Solving (GB + minAssGTZ) | 4.366 | -- |
| Expansion | 26.130 | 0.381388 |
| **Total** | **31.000** | **0.381388** |

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
| MMA | `170304580*v1 + 119616448*v1^2 + 81992745*v2 + 59808225*v1*v2 + 29904112*v2^2 + 170304580*v3 + 59808225*v1*v3 + 59808224*v2*v3 + 119616448*v3^2 + 81992745*v4 + 59808224*v1*v4 + 149520561*v2*v4 + 59808225*v3*v4 + 29904112*v4^2 + 170304580*v5 + 59808225*v1*v5 + 59808224*v2*v5 + 59808225*v3*v5 + 59808224*v4*v5 + 119616448*v5^2` |
| C++ | `170304580*v1 + 119616448*v1^2 + 81992745*v2 + 59808225*v1*v2 + 29904112*v2^2 + 170304580*v3 + 59808225*v1*v3 + 59808224*v2*v3 + 119616448*v3^2 + 81992745*v4 + 59808224*v1*v4 + 149520561*v2*v4 + 59808225*v3*v4 + 29904112*v4^2 + 170304580*v5 + 59808225*v1*v5 + 59808224*v2*v5 + 59808225*v3*v5 + 59808224*v4*v5 + 119616448*v5^2` |
| Diff | `0` **[MATCH]** |

## 5. Full Comparison Summary

| k | MMA | C++ | Diff | Result |
|:--|------|------|------|:------:|
| 0 | `1` | `1` | `0` | MATCH |
| 1 | `170304580*v1 + 119616448*v1^2 + 81992745*v2 + ... (20 terms)` | `170304580*v1 + 119616448*v1^2 + 81992745*v2 + ... (20 terms)` | `0` | MATCH |
| 2 | `66745271*v1 + 19936075*v1^4 + 136079496*v2 + ... (59 terms)` | `66745271*v1 + 19936075*v1^4 + 136079496*v2 + ... (59 terms)` | `0` | MATCH |
| 3 | `69572032*v1 + 95250135*v1^6 + 53877456*v2 + ... (215 terms)` | `69572032*v1 + 95250135*v1^6 + 53877456*v2 + ... (215 terms)` | `0` | MATCH |
| 4 | `35113864*v1 + 73837314*v1^8 + 86915377*v2 + ... (500 terms)` | `35113864*v1 + 73837314*v1^8 + 86915377*v2 + ... (500 terms)` | `0` | MATCH |

**Verdict: [PASS] -- All 5 orders match.**

