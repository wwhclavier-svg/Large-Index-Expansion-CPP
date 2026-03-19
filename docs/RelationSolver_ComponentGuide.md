# RelationSolver Component Guide

## Overview

The **RelationSolver** module reconstructs linear reduction relations between IBP matrix expansion coefficients. It works in conjunction with the layer recursion engine to solve homogeneous linear systems over both floating-point (`double`) and finite fields (`firefly::FFInt`).

**Key Responsibility**: Given expansion coefficient data at a specified level (`lev`) and degree (`deg`), find all linear dependencies among the basis elements.

---

## Architecture

### Type System

| Type | Use Case | Module |
|------|----------|--------|
| `double` | Floating-point validation, numerical testing | `LinearSolver_Eigen.hpp` |
| `firefly::FFInt` | Finite field computations (modular arithmetic) | `LinearSolver_FF.hpp` |

The entire `RelationSolver` namespace is templated on `typename T` to support both types transparently.

### Core Data Structures

#### `RegimeData<T>` (Lines 24–130)
Encapsulates a single **regime** (sector with specific algebraic constraints).

```cpp
template<typename T>
struct RegimeData {
    const seriesCoefficient<T>* C;              // Coefficient container reference
    std::vector<int> theta;                     // Sector constraint vector
    std::vector<std::vector<T>> A_ops;          // A_i matrices
    std::vector<std::vector<T>> A_inv_ops;      // A_i inverse matrices
    int nb;                                     // Basis block dimension
};
```

**Methods:**
- `prepare(lev, alphas)` — Pre-compute P(alpha) matrices for all alpha at given lev
- `computePRecursive(alpha, ...)` — Compute P(alpha) using recursion
- `getP(alpha)` — Retrieve cached P(alpha)

#### `RelationCoefficient<T>` (Lines 1480–1530)
Stores the **solution** to the relation problem: all (alpha, beta) pairs and their corresponding coefficient vectors.

```cpp
template<typename T>
class RelationCoefficient {
    std::vector<std::vector<int>> alphas_;      // All alpha multi-indices
    std::vector<std::vector<int>> betas_;       // All beta multi-indices
    std::vector<std::vector<T>> coeffs_;        // 1 + nullity coefficients per pair
    
    // Fast lookup by (alpha, beta)
    const std::vector<T>& operator()(
        const std::vector<int>& alpha,
        const std::vector<int>& beta) const;
};
```

#### `AdaptiveSamplingConfig` (Lines 180–220)
Controls how sampling points are generated and when to stop.

```cpp
struct AdaptiveSamplingConfig {
    int min_nu = 3;                             // Minimum sampling points
    int max_nu = 50;                            // Maximum sampling points
    int convergence_threshold = 3;              // Nullity stability plateau
    int lev_hint = 0;                           // Current (lev, deg) hints
    int deg_hint = 0;
    // ... additional parameters
};
```

#### `AdaptiveEquationBuilder<T>` (Lines 230–450)
**Core engine**: Incrementally builds and solves the linear system with convergence monitoring.

```cpp
template<typename T>
class AdaptiveEquationBuilder {
    // Builds equations at current sampling points
    BuildResult build(
        const std::vector<RegimeData<T>>& regimes,
        const std::vector<std::vector<int>>& nimaxLists,
        int ne);
    
    // Checks if nullity has stabilized
    bool hasConverged() const;
};
```

---

## Main API: `reconstructReductionRelation<T>()`

### Function Signature (Lines 1560–1660)

```cpp
template<typename T>
static std::pair<LinearSystemResult<T>, RelationCoefficient<T>> 
reconstructReductionRelation(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne, int lev, int deg,
    const AdaptiveSamplingConfig& config = {});
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `CTable` | `vector<vector<seriesCoefficient<T>>>` | Expansion coefficients: `CTable[r][s]` = regime (r, s) |
| `sector` | `vector<vector<int>>` | Sector constraints: `sector[r]` = theta vector for regime r |
| `A_list` | `vector<vector<vector<T>>>` | A matrices: `A_list[r][i]` = matrix A_i for regime r |
| `Ainv_list` | `vector<vector<vector<T>>>` | Inverse A matrices |
| `ne` | `int` | Dimension of multi-index vectors (alpha, beta) |
| `lev` | `int` | Max &#124;alpha&#124; to explore |
| `deg` | `int` | Max &#124;beta&#124; to explore |
| `config` | `AdaptiveSamplingConfig` | Sampling and convergence settings (optional) |

### Return Value

```cpp
std::pair<LinearSystemResult<T>, RelationCoefficient<T>> result;

// result.first: Solution structure
LinearSystemResult<T> lsr = result.first;
bool hasSolution = lsr.hasSolution;
vector<vector<T>> Mext = lsr.Mext;      // [cols] x [1+nullity]
vector<int> S = lsr.S;                  // Free variable indices

// result.second: Indexed coefficient access
RelationCoefficient<T> coeff = result.second;
const vector<T>& c = coeff({0,1}, {2,0});  // Get coeff for alpha=(0,1), beta=(2,0)
```

---

## Algorithm Flow

### Step-by-Step Execution

1. **Regime Construction** (Lines 1593–1610)
   ```cpp
   // For each regime (r,s), create RegimeData<T>
   // Attach: C pointer, theta vector, A/Ainv matrices, nb
   std::vector<RegimeData<T>> regimes;
   for (size_t r = 0; r < CTable.size(); ++r) {
       for (size_t s = 0; s < CTable[r].size(); ++s) {
           RegimeData<T> reg;
           reg.C = &CTable[r][s];
           reg.theta = sector[r];
           reg.A_ops = A_list[r];
           reg.A_inv_ops = Ainv_list[r];
           reg.nb = reg.C->basis_size();
           regimes.push_back(std::move(reg));
       }
   }
   ```

2. **Multi-Index Generation** (Lines 1610–1620)
   ```cpp
   // Generate all alpha with |alpha| <= lev
   std::vector<std::vector<int>> alphas, betas;
   RelationSolver::generateAllIndices(ne, lev, /* */ alphas, false);
   RelationSolver::generateAllIndices(ne, deg, /* */ betas, false);
   // alphas.size() * betas.size() = total variables
   ```

3. **Prepare Regime Data** (Lines 1620–1630)
   ```cpp
   // Pre-compute P(alpha) matrices for each regime
   for (auto& reg : regimes) {
       reg.prepare(lev, alphas);
   }
   ```

4. **Adaptive Equation Building** (Lines 1630–1650)
   ```cpp
   // Incrementally sample and build until convergence
   AdaptiveEquationBuilder<T> builder(config);
   auto build_result = builder.build(regimes, nimaxLists, ne);
   
   // build_result.converged: Has nullity stabilized?
   // build_result.nullspace: Solution matrix and free variables
   ```

5. **Result Conversion** (Lines 1650–1665)
   ```cpp
   // Convert BuildResult -> LinearSystemResult<T>
   LinearSystemResult<T> lsr_result;
   lsr_result.hasSolution = build_result.converged || 
                            build_result.nullspace.is_valid;
   lsr_result.Mext = /* transpose nullspace to [cols] x [1+nullity] */;
   lsr_result.S = build_result.nullspace.free_vars;
   ```

6. **Wrap in RelationCoefficient** (Lines 1665–1670)
   ```cpp
   // Create indexed accessor
   RelationCoefficient<T> coeff(alphas, betas, lsr_result);
   return {lsr_result, coeff};
   ```

---

## Usage Example (From test_RelationFF.cpp)

### Setup Phase

```cpp
// 1. Load all data structures
vector<vector<seriesCoefficient<FFInt>>> allResults;
vector<vector<int>> sector_list;
vector<vector<vector<FFInt>>> A_list, Ainv_list;
// ... populate from IBP/ring data loaders

// 2. Configure adaptive sampling
RelationSolver::AdaptiveSamplingConfig base_config;
base_config.min_nu = 3;
base_config.max_nu = 100;
base_config.convergence_threshold = 3;
```

### Solving Phase (For Each (lev, deg) Pair)

```cpp
for (int lev = lev_min; lev <= lev_max; ++lev) {
    for (int deg = deg_min; deg <= deg_max; ++deg) {
        
        // Dynamic sampling adjustment
        auto config = adjustSamplingConfig(
            base_config, lev, deg, ne, nb, k_max, num_regimes);
        
        // Solve for this level/degree
        auto [lsr_result, rel_coeff] = 
            RelationSolver::reconstructReductionRelation<FFInt>(
                allResults,
                sector_list,
                A_list,
                Ainv_list,
                ne, lev, deg,
                config);
        
        // Check success
        if (lsr_result.hasSolution) {
            cout << "Relations found at (lev=" << lev 
                 << ", deg=" << deg << ")" << endl;
            cout << "Nullity: " << lsr_result.S.size() << endl;
            
            // Access coefficients
            for (const auto& alpha : alphas) {
                for (const auto& beta : betas) {
                    const auto& coeff_vec = rel_coeff(alpha, beta);
                    // coeff_vec[0] = particular solution
                    // coeff_vec[1..] = basis of solution space
                }
            }
        } else {
            cout << "No relations found (system rank = full)" << endl;
        }
    }
}
```

---

## Key Algorithms

### 1. Multi-Index Generation: `generateAllIndices()`

Enumerates all integer vectors of given dimension with bounded ℓ₁ norm.

```cpp
// Generate all alpha with |alpha|_1 <= lev and alpha_i >= 0
std::vector<std::vector<int>> alphas;
std::vector<int> temp;
generateAllIndices(ne, lev, temp, alphas, false);
// Result size = C(lev + ne, ne) = binomial(lev+ne, ne)
```

### 2. P Matrix Computation: `RegimeData<T>::prepare()`

For each multi-index alpha, precompute:
$$P(\alpha) = \prod_{i=1}^{ne} A_i^{\alpha_i}$$

**Optimization**: Uses memoization to avoid redundant matrix multiplications.

### 3. Adaptive Sampling & Convergence

The `AdaptiveEquationBuilder` samples random `nu` vectors and tracks **nullity** (dimension of solution space):

```
Iteration 0: nu = nu_min → nullity = n0
Iteration 1: nu = nu_min+1 → nullity = n1
...
Iteration k: nullity = n_k

Converged when: n0 == n1 == ... == n_{threshold}
```

This ensures the found relations are **stable across random samplings**.

---

## Integration with Other Modules

### Dependencies

```
RelationSolver.hpp
    ├── LinearSolver.hpp (solveLinearSystem<T>)
    ├── LinearSolver_Eigen.hpp (double solver)
    ├── LinearSolver_FF.hpp (FFInt solver)
    ├── SeriesCoefficient.hpp (coefficient storage)
    ├── Combinatorics.hpp (index generation)
    └── UnifiedStorage.hpp (memory management)
```

### Typical Workflow

```
IBPMatrixLoader → LayerRecursion → seriesCoefficient
                                         ↓
                              RelationSolver
                                         ↓
                              RelationCoefficient (output)
```

---

## Performance Considerations

### Computational Complexity

For a given (lev, deg):
- **Variables**: $O(\binom{lev+ne}{ne} \times \binom{deg+ne}{ne})$
- **Equations per sample**: $O(\text{numRegimes} \times nb \times \text{kmax})$
- **Sampling iterations**: $O(\log(\text{variables}))$ in optimal case

### Memory Usage

- `RegimeData` stores P-matrix cache: $O(nb^2 \times \binom{lev+ne}{ne})$
- Global system matrix: $O(\text{equations} \times \text{variables})$

### Optimization Tips

1. **Reuse A matrices** across multiple (lev, deg) pairs (already done in test_RelationFF.cpp)
2. **Tune base_config** dynamically based on (lev, deg) size
3. **Use FFInt for large systems** (more numerically stable than double)

---

## Diagnostic Output

RelationSolver produces console output for debugging:

```
Adaptive sampling config:
  Variables: 120
  Equations per sample: 240
  Requested sampling range: [3, 50]
  Actual sampling range: [5, 80]
  
Equation building iteration 0: nullity=3
Equation building iteration 1: nullity=3
Equation building iteration 2: nullity=3
Converged! Final nullity: 3
System solved successfully.
```

---

## Common Pitfalls & Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| Convergence fails | Too few sampling points | Increase `max_nu` in config |
| Wrong nullity | Numerical instability | Use `FFInt` instead of `double` |
| Slow execution | Too many regimes/high (lev,deg) | Reduce `max_nu` or split into stages |
| Memory overflow | Large P-matrix cache | Process one regime at a time |

---

## Testing

### Unit Tests

```bash
cd build
./test_RelationFF    # Finite field relation solving
./test_RelationNew   # Extended test with multi-regimes
```

### Expected Output

```
Extracted 4 sectors from ring data.
=== Iterative Relation Solving ===
lev range: [0, 3]
deg range: [0, 3]

[lev=0, deg=0]: Relations found. Nullity=2
[lev=0, deg=1]: Relations found. Nullity=1
[lev=1, deg=0]: No relations (full rank)
...
Total time: 5.23 seconds
```

---

## See Also

- [LayerRecursion](./LayerRecursion_ComponentGuide.md) — Coefficient generation
- [LinearSolver](./LinearSolver_ComponentGuide.md) — Linear system solvers
- [AGENTS.md](../AGENTS.md#algorithm-flow) — High-level architecture
