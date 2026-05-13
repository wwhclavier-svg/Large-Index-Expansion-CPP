# Per-Region Expansion Orders

**Date**: 2026-05-11
**Status**: Design (approved)

## Problem

When `incre > 2 && nb > 1`, the `seriesCoefficient<T>` cache in layer recursion grows extremely large because:

- `size ∝ Σ_k Σ_{l=0}^{incre·k} BINOM[l+ne-1][ne-1] × nb × (nimax+1)`
- All regions use the same global `order`, even though some regions (nb=1) can easily handle higher orders while others (nb>1) produce excessive cache

## Solution

Per-region expansion orders: each IBP matrix region gets its own effective `order` based on IBP matrix properties (nb, incre), computed automatically inside `layerRecursion`.

---

## 1. Order Inference (inside `layerRecursion`)

The per-region effective order is computed inside `layerRecursion`, not in the caller:

```cpp
template<typename T>
auto layerRecursion(const IBPMatrixE<T>& ibpmat, int ne, int nb, int nibp,
                    int base_order, int order_padding = 1, int incre_hint = 2)
{
    int detected_incre = detectIncrement(ibpmat);
    int incre = (detected_incre == 1 || detected_incre == 3) ? detected_incre : incre_hint;
    int effective_order = (detected_incre <= 2 && nb <= 1)
                          ? base_order + order_padding
                          : base_order;
    // Proceed with effective_order and incre (rest unchanged)
}
```

**Rule**: `incre_detected <= 2 && nb <= 1` → `effective_order = base_order + order_padding`, else `effective_order = base_order`.

---

## 2. `batchProcessRecursion` Interface

```cpp
// Before
auto batchProcessRecursion(const vector<IBPMatrixE<T>>& mats, int order, int incre_hint = 2);
// After
auto batchProcessRecursion(const vector<IBPMatrixE<T>>& mats, int base_order,
                           int incre_hint = 2, int order_padding = 1);
```

---

## 3. `test_expandFF` Changes

```cpp
// New CLI: ./test_expandFF <family> [base_order] [order_padding]
int base_order = 4, order_padding = 1;
auto allResults = batchProcessRecursion<FFInt>(ibpmatlist, base_order, incre, order_padding);
```

---

## 4. `test_relationFF` / RelationSolver Changes

### 4a. RegimeEvaluator: Per-Regime k_max
```cpp
void init(...) {
    k_max_ = min(k_max, C_->getKmax());  // per-regime cap
}
```

### 4b. GlobalEquationAssembler: Variable Row Counts
```cpp
size_t totalRowsPerNu() const {
    size_t rows = 0;
    for (auto& eval : evaluators_) rows += eval.rowsPerNu();
    return rows;
}
```

### 4c. stable_order Verification
After solving with `stable_order = N`:
- Regimes with `C_->getKmax() > N` contribute extra rows at `k = N + 1`
- Verify solved relation coefficients satisfy these extra rows

---

## 5. Data Flow Summary

```
test_expandFF
  └─ batchProcessRecursion(mats, base_order, incre, order_padding)
       └─ for each matrix i:
            layerRecursion(mat[i], ..., base_order, order_padding, incre)
              ├─ detectIncrement → incre_i
              └─ effective_order = (incre_i<=2 && nb<=1) ? base+padding : base
              └─ seriesCoefficient<T>(effective_order, ...)

test_relationFF
  └─ load allResults (per-region C, varied getKmax())
  └─ for each (lev, deg):
       └─ RegimeEvaluator::init: k_max_ = min(k_max, C_->getKmax())
       └─ GlobalEquationAssembler: per-regime row counts
       └─ solve → stable_order = N
       └─ verify: extra rows from C_->getKmax() > N regimes
```

## Files Changed

| File | Change |
|------|--------|
| `include/LayerRecursion.hpp` | Signatures; internal order inference |
| `include/RelationSolver.hpp` | Per-regime k_max; variable rows; stable_order verify |
| `tests/test_expandFF.cpp` | New CLI args; updated API call |
| `tests/test_relationFF.cpp` | Per-regime k_max propagation; stable_order verification |
