# RelationSolver Quick Reference

## At a Glance

**What**: Linear relation reconstruction for IBP expansion coefficients  
**Where**: `include/RelationSolver.hpp` (namespace `RelationSolver`)  
**Types**: `double`, `firefly::FFInt`  
**Main Entry**: `reconstructReductionRelation<T>()`

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
config.convergence_threshold = 3;
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

### Pattern 1: Complete Loop (From test_RelationFF.cpp)
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
Input: CTable, sector, A/Ainv, ne, (lev, deg)
  ↓
Build RegimeData list (attach C, theta, A, nb)
  ↓
Generate all α vectors with |α|_1 ≤ lev
Generate all β vectors with |β|_1 ≤ deg
  ↓
Prepare each regime (pre-compute P(α))
  ↓
AdaptiveEquationBuilder:
  Loop until nullity converges:
    Sample random nu vector
    Build equations: for each regime, α, β, k
    Solve homogeneous system
    Track nullity
  ↓
Convert result to RelationCoefficient
  ↓
Output: (LinearSystemResult, RelationCoefficient)
```

---

## Config Tuning

| Goal | Setting |
|------|---------|
| Fast (small system) | `max_nu=10`, convergence_threshold=2 |
| Accurate (large system) | `max_nu=200`, convergence_threshold=5 |
| Memory-constrained | `max_nu=50`, process one (lev,deg) at a time |
| Finite field (recommended) | Use `FFInt`, `max_nu=100` |

**Adaptive Adjustment Function** (From test_RelationFF.cpp):
```cpp
auto config = adjustSamplingConfig(
    base_config, lev, deg, ne, nb, k_max, num_regimes);
// Automatically scales min_nu/max_nu based on variable count
```

---

## Diagnostic Output

```
Adaptive sampling config:
  Variables: 120
  Equations per sample: 240
  Actual sampling range: [5, 80]

Equation building iteration 0: nullity=3
Equation building iteration 1: nullity=3     ← Stabilized
Equation building iteration 2: nullity=3
Converged! Final nullity: 3
```

If not converging:
- Increase `max_nu`
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

See [RelationSolver Component Guide](./docs/RelationSolver_ComponentGuide.md) for:
- Complete API reference
- Algorithm details
- Performance considerations
- Troubleshooting guide

---

## Example Files

- **Usage**: [tests/test_RelationFF.cpp](../tests/test_RelationFF.cpp)
- **Config helper**: `adjustSamplingConfig()` in test_RelationFF.cpp
- **Advanced**: [include/IncrementalRelationSolver.hpp](../include/IncrementalRelationSolver.hpp)
