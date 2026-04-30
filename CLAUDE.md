# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **C++17 scientific computing project** focused on **IBP (Integration By Parts) matrix series expansion** and **finite field arithmetic**. It implements algorithms for computing IBP matrix expansion coefficients and reconstructing linear reduction relations.

Core functionality:
- **Layer Recursion**: Computes IBP matrix series expansion coefficients order-by-order
- **Finite Field Arithmetic**: Uses the FireFly library for computations over finite fields
- **Linear Relation Reconstruction**: Solves linear relations among expansion coefficients

The project uses heavy template metaprogramming to support both `double` (floating-point) and `firefly::FFInt` (finite field integers) simultaneously.

## Build Instructions

```bash
# Create and enter build directory
mkdir -p build && cd build

# Configure (requires Eigen3, FireFly, GMP installed)
cmake ..

# Build all targets
cmake --build .

# Build and run a single test target
cmake --build . --target test_expandFF
cd .. && ./build/test_expandFF
```

### Dependencies

**Required:**
- C++17 compatible compiler
- CMake 3.10+
- Eigen3 (linear algebra)
- GMP (multi-precision arithmetic)
- FireFly (finite field, pre-built at `/root/firefly-2.0.3`)

**Optional:**
- Mathematica (for `.wl` data generation scripts)

### Key Build Notes

- FireFly is pre-built at `/root/firefly-2.0.3` — not installed via package manager
- Tests are **standalone executables**, not managed by `ctest` or `add_test()`
- Binary data files (`IBPMat_*.bin`, `RingData_*.bin`) must exist in the working directory for tests

### Alternate Build Config

`CMakeLists_test.txt` uses a local `include/firefly/FFInt.hpp` stub instead of the external FireFly library. Use this for quick compile checks when FireFly is unavailable:

```bash
cp CMakeLists_test.txt CMakeLists.txt
cd build && cmake .. && cmake --build .
```

## Architecture

### Type System

The project uses **template polymorphism** with C++17 `if constexpr` for type dispatch:

```cpp
if constexpr (std::is_same_v<T, firefly::FFInt>) {
    // Finite field path
} else {
    // Floating-point path
}
```

| Type | Use Case | Header |
|------|----------|--------|
| `double` | Floating-point testing, numerical validation | Native |
| `firefly::FFInt` | Finite field computation (prime modulus) | `<firefly/FFInt.hpp>` |

### Module Hierarchy

```
LayerRecursion.hpp                    RelationSolver.hpp
       ↓                                      ↓
LayerRecursionCore.hpp              RegimeData<T>
       ↓                               RegimeEvaluator<T>
LayerRecursionCore.tpp              AdaptiveEquationBuilder<T>
       ↓                                      ↓
LinearSolver.hpp ←────────────────→ LinearSolver_FF.hpp
       ↓                               LinearSolver_Eigen.hpp
SeriesCoefficient.hpp              UnifiedStorage.hpp
IBPMatrixLoader_Binary.hpp         IncrementalRelationSolver.hpp
IBPMatrixLoader.hpp (JSON)         Combinatorics.hpp
RingDataLoader.hpp                 Utilities.hpp
                                   binomial.hpp / binomial2.hpp
```

### Core Data Structures

**`IBPMatrixE<T>`** (`IBPMatrixLoader_Binary.hpp` / `IBPMatrixLoader.hpp`):
- IBP matrix operator storage with 3D operators `N1, K1, M1, K1s, K2s`, 2D `F0`, and 4D `F2, F2s`
- Dimensions: `nibp` (IBP equations), `ne` (external variables), `nb` (block size)
- Two loader implementations: binary (`IBPMatrixLoader_Binary.hpp`) with magic `"IBP1"`, and JSON (`IBPMatrixLoader.hpp`)

**`seriesCoefficient<T>`** (`SeriesCoefficient.hpp`):
- 5D coefficient storage indexed as `(k, l, cid, j, i)` where:
  - `k`: expansion order
  - `l`: layer level
  - `cid`: seed/combinatorial index
  - `j`: basis index (0 to nb-1)
  - `i`: solution index (0 to nimax), i=0 is particular solution, i>0 are homogeneous solutions

**`LinearSystemResult<T>`** (`LinearSolver.hpp`):
```cpp
template<typename T>
struct LinearSystemResult {
    bool hasSolution;
    vector<vector<T>> Mext;  // Solution matrix (particular + nullspace)
    vector<int> S;           // Free variable indices
    vector<int> pivot_cols;  // Pivot column indices
};
```

### Layer Recursion Algorithm

The layer recursion has four layers:
1. **Entry**: `LayerRecursion.hpp` - `layerRecursion<T>()` template wrapper
2. **Core Logic**: `LayerRecursionCore.hpp` - function declarations
3. **Implementation**: `LayerRecursionCore.cpp` - non-template implementation
4. **Template Implementation**: `LayerRecursionCore.tpp` - template implementation

The core computation in `LayerRecursionCore::inhomogTerms<T>` builds:
```
inhomog = NMinus + NZero + NPluMi + NPlus + M1_contrib + MPlus
```

For each `(k, l, seed)`, assembles linear equations `M1 * x = -inhomog` and solves via `solveLinearSystem()`.

**Algorithm flow:**
1. Load IBP matrices from binary files via `loadAllIBPMatricesBinary<T>()`
2. Execute layer recursion via `layerRecursion<T>()` to compute expansion coefficients
3. Serialize results via `SeriesIO::saveAllResults()` to cache computed coefficients
4. Reconstruct relations via `RelationSolver::reconstructReductionRelation<T>()`

### RelationSolver Module

Key types in `include/RelationSolver.hpp`:
- **`RegimeData<T>`**: Encapsulates a sector with its coefficients and A operators
- **`RelationCoefficient<T>`**: Provides indexed access to relation coefficients by `(alpha, beta)` multi-indices
- **`AdaptiveSamplingConfig`**: Controls adaptive sampling and convergence detection
```cpp
struct AdaptiveSamplingConfig {
    int min_nu = 3;                // Minimum sample points
    int max_nu = 50;               // Maximum sample points
    int convergence_threshold = 3; // Nullity stability threshold
    int lev_hint = 0;              // Current (lev, deg) hint
    int deg_hint = 0;
};
```

Main API:
```cpp
template<typename T>
std::pair<LinearSystemResult<T>, RelationCoefficient<T>>
reconstructReductionRelation(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne, int lev, int deg,
    const AdaptiveSamplingConfig& config = {});
```

## Test Executables

| Test | Purpose | Command |
|------|---------|---------|
| `test_expandFF` | Finite field expansion coefficients | `./build/test_expandFF` |
| `test_relationFF` | Finite field linear relation reconstruction | `./build/test_relationFF` |
| `test_load_bub` | Load bub-format IBP matrices | `./build/test_load_bub` |
| `test_expand_family` | Expansion family testing | `./build/test_expand_family` |
| `test_IBPVerification` | IBP matrix verification | `./build/test_IBPVerification` |
| `test_ff_verify` | FireFly library verification | `./build/test_ff_verify` |

### Root-Level Standalone Tools

These files reside in the project root and are **not** part of the main `CMakeLists.txt`. Compile them manually when needed:

| File | Purpose |
|------|---------|
| `test_ff_verify.cpp` | Verify FireFly library basic functionality |
| `test_check_matrix.cpp` | Check matrix consistency |
| `test_matrix_dump.cpp` | Dump matrix contents for debugging |
| `test_solver.cpp` | Linear solver standalone test |
| `test_firefly_simple.cpp` | Minimal FireFly functionality test |
| `test_ffint.cpp` | Local FFInt stub test |
| `test_ff_div_debug.cpp` | Finite field division debugging |

**Note**: No `ctest` or `add_test()` — each test is a standalone executable with its own `main()`. Run from the project root directory where binary data files reside.

## Coding Conventions

### Naming

| Element | Convention | Example |
|---------|------------|---------|
| Files | PascalCase headers | `LayerRecursion.hpp` |
| Classes | PascalCase | `class seriesCoefficient` |
| Functions | camelCase | `layerRecursion()`, `getIndex()` |
| Constants | UPPER_CASE | `MAX_VAL`, `BINOM` |
| Namespaces | PascalCase | `namespace LayerRecursionCore` |

### Code Patterns

1. **Type dispatch with `if constexpr`**:
```cpp
if constexpr (std::is_same_v<T, firefly::FFInt>) {
    // FF path
} else {
    // double path
}
```

2. **Namespace organization**:
```cpp
namespace LayerRecursionCore { /* core functions */ }
namespace AlgebraData { /* data loading */ }
namespace RelationSolver { /* relation solving */ }
namespace SeriesIO { /* serialization */ }
```

### File Organization

- Templates: header files (`.hpp`) or `.tpp` files for large template implementations
- Implementations: `.cpp` files in `src/`
- Template implementations in `.tpp` files are included by corresponding `.hpp` files
- Files prefixed with `#` and ending in `.txt` (e.g., `#RelationRecon.cpp.txt`) are **archived legacy code**, not part of the current build

### Comment Language

- **Chinese (Simplified)**: Used for mathematical and algorithmic explanations
- **English**: Used for API documentation and brief inline comments

### Adding New Tests

1. Create `tests/test_MyFeature.cpp` with a `main()` function
2. Add to `CMakeLists.txt`:
```cmake
add_executable(test_MyFeature tests/test_MyFeature.cpp ${COMMON_SOURCES})
target_link_libraries(test_MyFeature ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)
```
3. Rebuild: `cd build && cmake --build .`

## Global State

- **`BINOM[MAX_VAL][MAX_VAL]`**: Precomputed binomial coefficients (initialized via `initBinomial()`)
- **`FFInt::p`**: Finite field modulus (set via `FFInt::set_new_prime()`)

### Memory Notes

`seriesCoefficient` preallocates large contiguous memory blocks. Be cautious with high expansion orders — they can cause significant memory pressure.

## Binary Data Format

IBP matrix files use a custom binary format with magic number `"IBP1"`. See `IBPMatrixLoader_Binary.hpp` for details.

Required data files in working directory:
- `IBPMat_*.bin` — IBP matrix binary data
- `RingData_*.bin` — Ring/algebra data
- `resCache_*.bin` — Expansion result cache

**Important**: Tests must be run from the project root directory (where `.bin` files reside), not from the `build/` directory.

## Documentation

- `docs/LayerRecursion_Algorithm.md` — Layer recursion algorithm details
- `docs/RelationSolver_ComponentGuide.md` — RelationSolver component guide
- `ReconstructReductionRelation_Documentation.md` — Mathematica package docs

## External References

- FireFly Library: https://github.com/firefly-library/firefly
- Eigen Documentation: https://eigen.tuxfamily.org/
