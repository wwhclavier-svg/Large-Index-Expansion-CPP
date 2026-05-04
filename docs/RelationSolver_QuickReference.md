# RelationSolver Quick Reference

## At a Glance

**What**: Linear relation reconstruction for IBP expansion coefficients  
**Where**: `include/RelationSolver.hpp` (namespace `RelationSolver`)  
**Types**: `double`, `firefly::FFInt`  
**Main Entry**: `reconstructAllRelations<T>()` — multi-(lev,deg) with RemoveSolvedVariables + order-stability detection  
**Legacy Entry**: `reconstructReductionRelation<T>()` — single (lev,deg)

---

## Basic Usage (3 Steps)

### 1️⃣ Prepare Data
```cpp
// Collect from layer recursion
std::vector<std::vector<seriesCoefficient<FFInt>>> CTable;
std::vector<std::vector<int>> sector;
std::vector<std::vector<std::vector<FFInt>>> A_list, Ainv_list;
int ne = 4;  // dimension
```

### 2️⃣ Configure Sampling
```cpp
RelationSolver::AdaptiveSamplingConfig config;
config.min_nu = 3;
config.max_nu = 100;
config.nullity_stable_threshold = 3;
config.plateau_size = 1;      // MMA "PlateauSize": extra orders to confirm stability
```

### 3️⃣ Solve for (lev, deg)
```cpp
auto [linear_result, relation_coeff] = 
    RelationSolver::reconstructReductionRelation<FFInt>(
        CTable, sector, A_list, Ainv_list, 
        ne, lev, deg, config);

if (linear_result.hasSolution) {
    // Access: relation_coeff({0,1}, {2,0})
    // Get: std::vector<FFInt> with solutions
}
```

---

## Key Data Types

| Type | Purpose |
|------|---------|
| `RegimeData<T>` | One sector with C, theta, A-operators |
| `RelationCoefficient<T>` | Index-based coefficient accessor |
| `LinearSystemResult<T>` | Solution: hasSolution, Mext, S |
| `AdaptiveSamplingConfig` | Sampling strategy & convergence |

---

## Important Methods

```cpp
// Generate all alpha with |α|_1 ≤ lev
std::vector<std::vector<int>> alphas;
std::vector<int> temp;
RelationSolver::generateAllIndices(ne, lev, temp, alphas, false);

// Access relation coefficient for (α, β)
const auto& coeff_vec = relation_coeff(alpha, beta);
// coeff_vec[0] = particular solution
// coeff_vec[i] (i>0) = basis elements of solution space
```

---

## Common Patterns

### Pattern 1: Complete Loop (From test_relationFF.cpp)
```cpp
for (int lev = 0; lev <= 3; ++lev) {
    for (int deg = 0; deg <= 3; ++deg) {
        auto config = adjustSamplingConfig(base, lev, deg, ...);
        auto [res, coeff] = RelationSolver::reconstructReductionRelation<FFInt>(
            CTable, sector, A_list, Ainv_list, ne, lev, deg, config);
        
        cout << "Nullity at (" << lev << "," << deg << "): " 
             << res.S.size() << endl;
    }
}
```

### Pattern 2: Single Regime
```cpp
std::vector<std::vector<seriesCoefficient<double>>> CTable_one = {
    {C_single}
};
auto [res, coeff] = RelationSolver::reconstructReductionRelation<double>(
    CTable_one, {theta}, {A}, {Ainv}, ne, lev, deg);
```

### Pattern 3: Multi-Index Enumeration
```cpp
std::vector<int> temp;
std::vector<std::vector<int>> all_alphas;
RelationSolver::generateAllIndices(ne, max_lev, temp, all_alphas, false);
// all_alphas = {(0,0,...,0), (1,0,...,0), ..., (0,0,...,max_lev)}
```

---

## Algorithm Outline

```
Input: CTable, sector, A/Ainv, ne, lev_max, deg_max
  ↓
Build RegimeData list (attach C, theta, A, nb)
  ↓
For lev = 0..lev_max, deg = 0..deg_max:
  Generate α with |α|₁ ≤ lev, β with |β|₁ ≤ deg
  ↓
  RemoveSolvedVariables: filter dominated (α,β) pairs
  ↓
  Phase 1: AdaptiveEquationBuilder — nu-sampling loop
    Sample random nu vector
    Build equations: for each regime, α, β, k
    Split rows by expansion order r (zero-cost grouping)
    Solve homogeneous system, track nullity
    Converge when nullity stable across nu-points
  ↓
  Phase 2: Order-stability analysis (zero redundant evaluation)
    Incrementally add rows by order r = 0..k_max
    Track nullity at each order
    Detect plateau: nullity unchanged for plateau_size+1 orders
    OR nullity=0 → definitive (all variables determined)
    Record stable_order (or -2 if never stable)
  ↓
  Convert result to RelationCoefficient
  ↓
Output: stability matrix + relation counts per (lev, deg)
```

---

## Config Tuning

| Goal | Setting |
|------|---------|
| Fast (small system) | `max_nu=10`, nullity_stable_threshold=2 |
| Accurate (large system) | `max_nu=200`, nullity_stable_threshold=5 |
| Memory-constrained | `max_nu=50`, process one (lev,deg) at a time |
| Finite field (recommended) | Use `FFInt`, `max_nu=100` |
| Order stability (MMA-like) | `plateau_size=1` (default), checks -2 flag for insufficient order |

**Adaptive Adjustment Function** (From test_relationFF.cpp):
```cpp
auto config = adjustSamplingConfig(
    base_config, lev, deg, ne, nb, k_max, num_regimes);
// Automatically scales min_nu/max_nu based on variable count
```

---

## Diagnostic Output

```
  --- (lev=1, deg=1) --- vars=9 active=9
    [ν-sampling] vars=9 eq/sample=50 min_req=1 eff_max=10
    ν points: {1,0}, {0,1}, {1,1}, {1,2}, ...
    sol_dim=2 independent=2 stable_order=3

=== Stability bounds (lev x deg) ===
           0       1       2
  0        0       1       2
  1        1       3      -2
  2        2       4      -2

=== Relation counts (lev x deg) ===
           0       1       2
  0        0       0       0
  1        0       2       1
  2        0       0       4
```

**stable_order** values: non-negative = stable at that expansion order, `-2` = not stable within available orders (need higher k_max). Only stable solutions should be considered reliable.

If not converging:
- Increase `max_nu` for nu-space convergence
- Increase expansion order (`k_max`) for order-space convergence (fixes -2 entries)
- Use FFInt instead of double
- Check data validity (A matrices invertible?)

---

## Error Handling

```cpp
try {
    auto [res, coeff] = RelationSolver::reconstructReductionRelation<FFInt>(...);
    
    if (!res.hasSolution) {
        std::cout << "No relations found (full rank system)" << std::endl;
    } else {
        std::cout << "Relations: " << res.S.size() << std::endl;
        // Use relation_coeff to access details
    }
} catch (const std::exception& e) {
    std::cerr << "Solver error: " << e.what() << std::endl;
}
```

---

## Integration Points

| Component | Role | File |
|-----------|------|------|
| **Input** | Layer recursion coefficients | `seriesCoefficient<T>` |
| **L-Solver** | Linear system solver | `LinearSolver_Eigen.hpp`, `_FF.hpp` |
| **Output** | Relation coefficients | `RelationCoefficient<T>` (to next processing stage) |

---

## Full Documentation

See [RelationSolver Component Guide](RelationSolver_ComponentGuide.md) for:
- Complete API reference
- Algorithm details
- Performance considerations
- Troubleshooting guide

---

## Example Files

- **Usage**: [tests/test_relationFF.cpp](../tests/test_relationFF.cpp)
- **Config helper**: `adjustSamplingConfig()` in test_relationFF.cpp
- **Advanced**: [include/IncrementalRelationSolver.hpp](../include/IncrementalRelationSolver.hpp)
- **Algorithm theory**: [Reconstruct_Algorithm.md](Reconstruct_Algorithm.md) — MMA vs C++ comparison, per-point evaluation details
