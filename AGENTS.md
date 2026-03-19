# AGENTS.md - Project Guide for AI Coding Agents

## Project Overview

This is a **C++17 scientific computing project** focused on **IBP (Integration By Parts) matrix expansion** and **finite field operations**. The codebase implements algorithms for:

- **Layer Recursion**: Computing series expansion coefficients for IBP matrices
- **Finite Field Arithmetic**: Computing over finite fields using the FireFly library
- **Linear Relation Reconstruction**: Solving for linear relations between expansion coefficients

The project uses heavy template metaprogramming to support both `double` (floating-point) and `firefly::FFInt` (finite field integer) data types.

---

## Technology Stack

| Component | Purpose | Required |
|-----------|---------|----------|
| C++17 | Core language | Yes |
| CMake 3.10+ | Build system | Yes |
| Eigen3 | Linear algebra (floating-point) | Yes |
| FireFly | Finite field arithmetic | Yes |
| GMP | Multiple precision arithmetic (FireFly dependency) | Yes |
| nlohmann/json | JSON parsing (bundled in `include/json.hpp`) | Bundled |
| MATLAB | Data generation scripts (`.m` files) | Optional |

---

## Project Structure

```
.
├── CMakeLists.txt           # Main build configuration
├── #CMakeLists_full.txt     # Extended build config (backup/reference)
├── include/                 # Header files (templates + interfaces)
│   ├── LayerRecursion.hpp       # Main recursion algorithm wrapper
│   ├── LayerRecursionCore.hpp   # Core recursion functions
│   ├── LayerRecursionCore.tpp   # Template implementations
│   ├── IBPMatrixLoader_Binary.hpp  # Binary IBP matrix loader
│   ├── SeriesCoefficient.hpp     # Coefficient storage class
│   ├── SeriesCoefficientIO.hpp   # Binary serialization for coefficients
│   ├── LinearSolver.hpp          # Unified solver dispatcher
│   ├── LinearSolver_Eigen.hpp    # Eigen-based floating-point solver
│   ├── LinearSolver_FF.hpp       # FireFly finite field solver
│   ├── RelationSolver.hpp        # Linear relation reconstruction
│   ├── RingDataLoader.hpp        # Ring/algebraic data loader
│   ├── Combinatorics.hpp         # Index utilities & seed generation
│   ├── binomial.hpp              # Precomputed binomial coefficients
│   ├── UnifiedStorage.hpp        # Memory management utilities
│   ├── Utilities.hpp             # General utilities
│   └── json.hpp                  # Third-party JSON library
├── src/                     # Source files (non-template implementations)
│   ├── LayerRecursionCore.cpp   # Core function implementations
│   ├── layerRecursion.cpp       # Legacy implementations
│   └── main.cpp                 # Legacy driver (not actively used)
├── tests/                   # Test executables (each with main())
│   ├── test_expandFF.cpp        # Finite field expansion test
│   ├── test_RelationFF.cpp      # Finite field relation solving test
│   ├── test_expand.cpp          # Double-precision expansion test
│   ├── test_recons.cpp          # Reconstruction test
│   └── *_Test.cpp               # Component-specific tests
├── build/                   # Build output directory (generated)
├── .vscode/                 # VS Code configuration
│   ├── tasks.json           # Build tasks
│   └── c_cpp_properties.json # IntelliSense config
└── .github/
    └── copilot-instructions.md  # Additional coding guidelines

# Data files (binary and JSON formats)
├── IBPMat_*.bin             # IBP matrix binary data files
├── RingData_*.bin           # Ring data binary files
├── resCache_*.bin           # Expansion result cache files
└── *.json                   # Matrix data in JSON format

# MATLAB scripts for data generation
├── ibpmatSR5m.m
├── ibpmatDB.m
├── ibpmatNP.m
└── ibpmatE.m
```

---

## Build Instructions

### Standard CMake Build

```bash
# 1. Create and enter build directory
mkdir -p build && cd build

# 2. Configure (requires Eigen3, FireFly, GMP installed)
cmake ..

# 3. Build
cmake --build .

# 4. Run tests
./test_expandFF
./test_RelationFF
```

### Build Options

The `CMakeLists.txt` provides:
- `BUILD_FF_TESTS` (ON by default): Build finite field tests requiring FireFly

### Dependencies Installation

**Ubuntu/Debian:**
```bash
sudo apt-get install libeigen3-dev libgmp-dev
# FireFly must be built from source: https://github.com/firefly-library/firefly
```

**macOS:**
```bash
brew install eigen gmp
# FireFly must be built from source
```

---

## Code Organization & Architecture

### Type System

The project uses **template polymorphism** to support multiple numeric types:

| Type | Use Case | Header |
|------|----------|--------|
| `double` | Floating-point testing, numerical validation | Native |
| `firefly::FFInt` | Finite field computations (prime modulus) | `<firefly/FFInt.hpp>` |

### Key Data Structures

#### `IBPMatrixE<T>` (IBPMatrixLoader_Binary.hpp)
Storage for IBP matrix operators:
```cpp
template<typename T>
struct IBPMatrixE {
    vector<vector<vector<T>>> N1, K1, M1, K1s, K2s;  // 3D operators
    vector<vector<T>> F0;                              // 2D operator
    vector<vector<vector<vector<T>>>> F2, F2s;        // 4D operators
    int nibp, ne, nb;  // dimensions
};
```

#### `seriesCoefficient<T>` (SeriesCoefficient.hpp)
5-dimensional coefficient storage indexed by `(k, l, cid, j, i)`:
- `k`: expansion order
- `l`: layer level
- `cid`: seed/composition index
- `j`: basis index (0 to nb-1)
- `i`: solution index (0 to nimax)

#### `LinearSystemResult<T>` (LinearSolver_*.hpp)
Unified result structure for linear solvers:
```cpp
template<typename T>
struct LinearSystemResult {
    bool hasSolution;
    vector<vector<T>> Mext;  // Solution matrix (particular + nullspace)
    vector<int> S;           // Free variable indices
};
```

---

## RelationSolver Specification

### 📍 Location & Purpose
- **File**: `include/RelationSolver.hpp`
- **Namespace**: `RelationSolver`
- **Role**: Reconstructs linear reduction relations between IBP matrix expansion coefficients
- **Supported Types**: Template-based on `T` (`double`, `firefly::FFInt`)

### Key Data Structures

#### `RegimeData<T>`
Encapsulates a single sector with its coefficients and A-operators.
```cpp
template<typename T>
struct RegimeData {
    const seriesCoefficient<T>* C;           // Coefficient reference
    std::vector<int> theta;                  // Sector identifier
    std::vector<std::vector<T>> A_ops;       // A_i matrices
    std::vector<std::vector<T>> A_inv_ops;   // A_i^{-1} matrices
    int nb;                                  // Basis dimension
};
```

#### `RelationCoefficient<T>`
Result structure providing indexed access to relation coefficients.
```cpp
template<typename T>
class RelationCoefficient {
    // Access coefficient for given multi-indices
    const std::vector<T>& operator()(
        const std::vector<int>& alpha,  // Multi-index constraint
        const std::vector<int>& beta)   // Multi-index power
    const;
};
```

#### `AdaptiveSamplingConfig`
Configuration for adaptive sampling and convergence checking.
```cpp
struct AdaptiveSamplingConfig {
    int min_nu = 3;                    // Minimum sampling points
    int max_nu = 50;                   // Maximum sampling points
    int convergence_threshold = 3;     // Nullity stability detection
    int lev_hint = 0;                  // Current (lev, deg) level hints
    int deg_hint = 0;
};
```

### Main API: `reconstructReductionRelation<T>()`

**Purpose**: Solves for all linear relations at a specific (lev, deg) level.

```cpp
template<typename T>
std::pair<LinearSystemResult<T>, RelationCoefficient<T>> 
reconstructReductionRelation(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne,                            // Dimension of multi-indices
    int lev, int deg,                  // Current level/degree
    const AdaptiveSamplingConfig& config = {});
```

### Usage Pattern (From test_RelationFF.cpp)

```cpp
// For each (lev, deg) pair:
auto [linear_result, relation_coeff] = 
    RelationSolver::reconstructReductionRelation<FFInt>(
        allResults,        // Coefficients from layer recursion
        sector_list,       // Sector identifiers
        A_list, Ainv_list, // A-operators
        ne, lev, deg,      // Dimensions
        config);           // Adaptive sampling config

// Check if relations were found
if (linear_result.hasSolution) {
    int num_relations = linear_result.S.size();
    
    // Access specific relation coefficients
    for (const auto& alpha : alphas) {
        for (const auto& beta : betas) {
            const auto& coeff_vec = relation_coeff(alpha, beta);
            // coeff_vec[0] = particular solution
            // coeff_vec[1..] = solution space basis
        }
    }
}
```

### Algorithm: Adaptive Multipoint Sampling

1. **Regime Preparation** — Compute P(α) matrices for all multi-indices α at given lev
2. **Iterative Sampling** — Generate random nu vectors; each produces |regimes| × nb × (k_max+1) equations
3. **System Solving** — Solve the homogeneous linear system; track nullity dimension
4. **Convergence Check** — Stop when nullity stabilizes over convergence_threshold iterations
5. **Solution Extraction** — Package result as RelationCoefficient indexed by (α, β)

### Related Files

For detailed implementation and usage examples, see:
- [RelationSolver Component Guide](./docs/RelationSolver_ComponentGuide.md)
- [test_RelationFF.cpp](./tests/test_RelationFF.cpp) — Practical usage example
- [Incremental Solver](./include/IncrementalRelationSolver.hpp) — Advanced multi-level solver

---

### Algorithm Flow

1. **Data Loading**: `loadAllIBPMatricesBinary<T>()` loads IBP matrices from binary files
2. **Layer Recursion**: `layerRecursion<T>()` computes expansion coefficients order-by-order
3. **Serialization**: `SeriesIO::saveAllResults()` caches computed coefficients
4. **Relation Solving** ⭐:
   - Configure adaptive sampling with (lev, deg) hints
   - Call `reconstructReductionRelation<T>()` for each (lev, deg)
   - Get back solutions as `RelationCoefficient<T>` for downstream processing

---

## Testing Strategy

### Test Executables

| Test | Purpose | Data Type | Command |
|------|---------|-----------|---------|
| `test_expandFF` | Expansion coefficient calculation | `FFInt` | `./build/test_expandFF` |
| `test_RelationFF` | Linear relation reconstruction | `FFInt` | `./build/test_RelationFF` |
| `test_expand` | Double-precision expansion | `double` | `./build/test_expand` |
| `test_recons` | Reconstruction algorithms | `double` | `./build/test_recons` |

**Note**: This project does not use `ctest` or `add_test()`. Tests are standalone executables.

### Running Tests

```bash
cd build

# Finite field expansion test (requires IBPMat_DPpart_QuadriScale.bin)
./test_expandFF

# Relation solving test (requires RingData_DBtop.bin, IBPMat_DBtop.bin)
./test_RelationFF
```

### Test Data Files

Tests expect specific binary data files in the working directory:
- `IBPMat_DBtop.bin`, `IBPMat_DPpart_QuadriScale.bin` - IBP matrices
- `RingData_DBtop.bin`, `RingData_DPpart_QuadriScale.bin` - Ring data
- `ExpansionResults_cache.bin` - Coefficient cache (auto-generated)

---

## Coding Conventions

### Naming Style

| Element | Convention | Example |
|---------|------------|---------|
| Files | PascalCase for headers | `LayerRecursion.hpp` |
| Classes | PascalCase | `class seriesCoefficient` |
| Functions | camelCase | `layerRecursion()`, `getIndex()` |
| Variables | snake_case (local), camelCase (members) | `int num_regs;`, `int numRegs;` |
| Constants | UPPER_CASE | `MAX_VAL`, `BINOM` |
| Templates | PascalCase | `typename T`, `typename Field` |

### Comment Language

- **Primary**: Chinese (简体中文) - for mathematical/algorithmic explanations
- **Secondary**: English - for API documentation and brief notes

### File Organization

- Templates: Header-only (`.hpp`) or `.tpp` files
- Implementations: `.cpp` files in `src/`
- Each header should have include guards: `#ifndef FILENAME_HPP`

### Code Patterns

1. **Template Specialization** for type dispatch:
```cpp
if constexpr (std::is_same_v<T, firefly::FFInt>) {
    // Finite field path
} else {
    // Floating point path
}
```

2. **Namespace organization**:
```cpp
namespace LayerRecursionCore { /* core functions */ }
namespace AlgebraData { /* data loading */ }
namespace RelationSolver { /* relation solving */ }
namespace SeriesIO { /* serialization */ }
```

---

## Common Tasks

### Adding a New Test

1. Create `tests/test_MyFeature.cpp` with a `main()` function
2. Add to `CMakeLists.txt`:
```cmake
add_executable(test_MyFeature tests/test_MyFeature.cpp ${COMMON_SOURCES})
target_link_libraries(test_MyFeature ${FIREFLY_LIBRARY} ${GMP_LIBRARY} Eigen3::Eigen)
```
3. Rebuild: `cd build && cmake --build .`

### Adding a New Header

1. Place in `include/MyHeader.hpp`
2. Use include guards and namespace
3. Include from `src/` or other headers as needed

### Modifying the Core Algorithm

The layer recursion algorithm has three layers:
1. **Entry**: `LayerRecursion.hpp` - `layerRecursion<T>()` template
2. **Core Logic**: `LayerRecursionCore.hpp` - function declarations
3. **Implementation**: `LayerRecursionCore.cpp` - non-template implementations
4. **Template Impl**: `LayerRecursionCore.tpp` - template implementations

Edit the appropriate layer based on your changes.

---

## Important Notes

1. **Binary Data Format**: IBP matrix files use a custom binary format (magic number `"IBP1"`). See `IBPMatrixLoader_Binary.hpp` for specification.

2. **Global State**: 
   - `BINOM[MAX_VAL][MAX_VAL]` - precomputed binomial coefficients (initialized via `initBinomial()`)
   - `FFInt::p` - finite field modulus (set via `FFInt::set_new_prime()`)

3. **Memory Management**: `seriesCoefficient` preallocates large contiguous memory blocks. Be cautious with large `order` values.

4. **Build Directory**: Do not edit files in `build/` directly. They are generated by CMake.

5. **MATLAB Scripts**: The `.m` files are for data generation in MATLAB/Mathematica. They are not required for C++ builds.

---

## References

- **FireFly Library**: https://github.com/firefly-library/firefly
- **Eigen Documentation**: https://eigen.tuxfamily.org/
- **IBP Method**: Integration-by-parts identities for Feynman integral reduction (physics)
