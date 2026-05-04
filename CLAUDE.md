# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Behavioral Guidelines

**1. Think Before Coding** — State assumptions explicitly. If uncertain, ask. Present tradeoffs when multiple approaches exist.
**2. Simplicity First** — Minimum code that solves the problem. No abstractions for single-use code. No speculative features.
**3. Surgical Changes** — Touch only what the task requires. Don't refactor adjacent code. Match existing style.
**4. Goal-Driven Execution** — Define verifiable success criteria. Loop until verified. Write tests to confirm correctness.

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
| `firefly::FFInt` | Finite field computation (prime modulus) | `<firefly/FFInt.hpp>` (external FireFly library) |

A local minimal FFInt stub exists at `include/firefly/FFInt.hpp` for quick compile checks when the external FireFly library is unavailable (used by `CMakeLists_test.txt`).

### Module Hierarchy

```
LayerRecursion.hpp                    RelationSolver.hpp
       ↓                              (RegimeData<T>, RegimeEvaluator<T>,
LayerRecursionCore.hpp                 AdaptiveEquationBuilder<T>,
       ↓                               RelationCoefficient<T>)
LayerRecursionCore.tpp                        ↓
       ↓                              IncrementalRelationSolver.hpp
LinearSolver.hpp                              ↓
  (dispatches to ↓)                  UnifiedStorage.hpp
LinearSolver_FF.hpp                  Combinatorics.hpp
LinearSolver_Eigen.hpp               Utilities.hpp
       ↓                              binomial.hpp / binomial2.hpp
SeriesCoefficient.hpp
SeriesCoefficientIO.hpp
IBPMatrixLoader_Binary.hpp
IBPMatrixLoader.hpp (JSON)
RingDataLoader.hpp
```

**Note**: `RegimeData<T>`, `RegimeEvaluator<T>`, `AdaptiveEquationBuilder<T>`, and `RelationCoefficient<T>` are all defined within `RelationSolver.hpp` — they are not separate files.

### Core Data Structures

**`IBPMatrixE<T>`** (`IBPMatrixLoader_Binary.hpp` / `IBPMatrixLoader.hpp`):
- IBP matrix operator storage with 3D operators `N1, K1, M1, K1s, K2s`, 2D `F0`, and 4D `F2, F2s`
- Dimensions: `nibp` (IBP equations), `ne` (external variables), `nb` (block size)
- Two loader implementations: binary (`IBPMatrixLoader_Binary.hpp`) with magic `"IBP1"`, and JSON (`IBPMatrixLoader.hpp`)

**`seriesCoefficient<T>`** (`SeriesCoefficient.hpp`):
- 5D coefficient storage indexed as `(k, l, cid, j, i)` where:
  - `k`: expansion order
  - `l`: layer level, range `[0, incre*k]` — **NOT** `[0, k]`
  - `cid`: seed/combinatorial index
  - `j`: basis index (0 to nb-1)
  - `i`: solution index (0 to nimax), i=0 is particular solution, i>0 are homogeneous solutions
- `incre`: level increment per order, typically 2 (lmax = incre * k at each order k)

**`SeriesIO`** (`SeriesCoefficientIO.hpp`):
- Binary serialization of `seriesCoefficient<T>` with magic `"SERCOEF"` and versioning
- Provides `writeCoefficient()`, `readCoefficient()`, `saveAllResults()`, `loadAllResults()`

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
    int min_nu = 0;                    // Minimum sample points (0=auto)
    int max_nu = 200;                  // Maximum sample points (hard cap)
    double safety_factor = 1.2;        // Oversampling factor
    int nullity_stable_threshold = 3;  // Stability count for convergence
    int verification_points = 3;       // Extra verification points
    int lev_hint = 2;                  // |alpha| bound hint
    int deg_hint = 2;                  // |beta| bound hint
};
```

Main APIs:

**Single-configuration** (low-level):
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

**Multi-configuration with RemoveSolvedVariables** (recommended):
```cpp
template<typename T>
std::vector<LevDegResult<T>> reconstructAllRelations(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne, int lev_max, int deg_max,
    const AdaptiveSamplingConfig& config = {});
```

### RemoveSolvedVariables Strategy

The high-level `reconstructAllRelations()` implements the MMA-inspired **RemoveSolvedVariables** strategy to reduce equation redundancy across `(lev, deg)` configurations:

1. **Nested loop**: outer loop over seed levels `lev = 0..lev_max`, inner over coefficient degrees `deg = 0..deg_max`
2. **Variable filtering** before each solve: eliminates `(alpha, beta)` pairs that are "dominated" by previously-solved independent variables:
   - **Same-level, lower-degree**: if `b[α_s, β_s]` solved at `deg_s < deg`, eliminate `b[α, β]` where `α == α_s` and `β >= β_s` componentwise
   - **Cross-level**: if `b[α_s, β_s]` solved at `lev_s < lev`, eliminate `b[α, β]` where `α >= α_s` and `β >= β_s` componentwise
3. **Independent variable tracking**: after each solve, the free (non-pivot) columns map to `(alpha, beta)` pairs that parametrize the nullspace. These feed into subsequent filtering.
4. **Column compression**: the equation matrix is compressed to only active columns before Gaussian elimination, then solutions are expanded back to the full variable space.

This reduces variable counts by up to ~50% (e.g., for bub00 at lev=2,deg=2: 36→19 active variables).

## Test Executables

| Test | Purpose | Command |
|------|---------|---------|
| `test_expandFF` | Finite field expansion coefficients | `./build/test_expandFF` |
| `test_relationFF` | Finite field linear relation reconstruction | `./build/test_relationFF` |
| `test_load_bub` | Load bub-format IBP matrices | `./build/test_load_bub` |
| `test_expand_family` | Expansion family testing | `./build/test_expand_family` |
| `test_IBPVerification` | IBP matrix verification | `./build/test_IBPVerification` |
| `test_ff_verify` | FireFly library verification | `./build/test_ff_verify` |

**Note**: No `ctest` or `add_test()` — each test is a standalone executable with its own `main()`. Run from the project root directory where binary data files reside (tests load data files via **relative paths**).

Historical test files are archived in `tests/archive/` (`test_expand.cpp`, `test_recons.cpp`, `test_RelationNew.cpp`, etc.) — not part of the current build but useful for reference.

Additional standalone tools live in `tools/` (e.g., `test_ff_verify.cpp` is in the CMake build) — compile manually for ad-hoc debugging.

## Common Pitfalls

### FFInt Type Safety: Never Cast Negative Integers

`firefly::FFInt` only has a `FFInt(uint64_t)` constructor — no signed integer constructor. This means:

```cpp
// BROKEN: int(-1) silently converts to uint64_t(2^64-1), then mod p = garbage (~4M for p=179424673)
FFInt x = static_cast<FFInt>(-1);

// CORRECT: Use operator-() on a positive FFInt
FFInt x = -FFInt(1);
```

Any `static_cast<FFInt>(negative_int)` is a bug — C++ silently converts the negative int to a huge uint64_t before the FFInt constructor sees it. No compiler warning is generated.

This affected `sgn()` in `Utilities.hpp` and any code that constructs FFInt from negative literals. If adding a signed constructor, use `int64_t` to avoid ambiguity with the existing `uint64_t` overload for `long long` arguments.

### l-Loop Bound: Always `incre * k`, Never Just `k`

`seriesCoefficient` stores data for `l ∈ [0, incre * k]`, not `[0, k]`. Any loop summing over `l` for a given `k` must use `l <= incre * k` as the upper bound. Using `l <= k` silently drops all coefficients at `l > k`, producing wrong but self-consistent results that pass EquationVerify. With the default `incre = 2`, this means roughly half the coefficients are missed.

Affected functions: `step2_computeG`, `step4_computeF2` (convolution), and any loop iterating `C(k, l, cid, j, i)` over `l`.

### Self-Consistency Checks Are Not Independent Verification

Tests like EquationVerify (`M1*C + Total == 0`) reuse the same equation-building code that generates the coefficients. If the bug is in the equation construction (e.g., a wrong sign in `sgn()`), both generator AND verifier produce wrong but mutually consistent results. Always pair self-consistency checks with an independent validation method — reference implementation comparison (MMA) or substituting results back into the original (pre-decomposition) IBP equations.

## Coding Conventions

### Naming

| Element | Convention | Example |
|---------|------------|---------|
| Files | PascalCase headers | `LayerRecursion.hpp` |
| Classes | PascalCase | `class seriesCoefficient` |
| Functions | camelCase | `layerRecursion()`, `getIndex()` |
| Variables (local) | snake_case | `int num_regs;` |
| Variables (member) | camelCase | `int numRegs;` |
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
- `tools/` — standalone C++ test/utility source files not in main CMake build
- `mma/` — Mathematica `.wl` scripts (IBP expansion, relation reconstruction, bubble generation)
- `archive/` — archived CMake configs and retired build files
- Include guards use `#ifndef FILENAME_HPP` / `#define FILENAME_HPP` pattern
- Third-party `nlohmann/json` is bundled as `include/json.hpp`

### Comment Language

- **Chinese (Simplified)**: Used for mathematical and algorithmic explanations
- **English**: Used for API documentation and brief inline comments

### Extended Documentation

- `AGENTS.md` — Extended bilingual (Chinese/English) guide with exhaustive module specs, testing strategy, and dependency details. Consult for deeper context.
- `docs/ComprehensiveReport.tex` — Canonical theory write-up of the Large Index Expansion method: asymptotic-solution completeness theorem, block-recursive structure, geometric classification of solution spaces, and C++ finite-field implementation. Consult for theoretical principles.

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

- `docs/Core_Working.md` — Module index: theory, workflow, and test cross-reference
- `docs/LayerRecursion_Algorithm.md` — Layer recursion algorithm details
- `docs/RelationSolver_ComponentGuide.md` — RelationSolver component guide
- `docs/RelationSolver_QuickReference.md` — Quick reference for RelationSolver API
- `docs/RelationSolver_Documentation_Hub.md` — Documentation hub/index
- `docs/ReconstructReductionRelation_Documentation.md` — Mathematica package docs
- `docs/Reconstruct_Algorithm.md` — Reconstruction algorithm: MMA vs C++ implementation comparison, verification methods

### Verification (`verify/`)

The `verify/` directory is a structured cross-validation framework comparing C++ output against Mathematica (MMA) and Kira IBP reduction:

```
verify/
├── README.md                    # Full verification workflow with step-by-step commands
├── FamilyDatabase/
│   ├── FamilyDatabase.wl        # Unified integral family definitions (1L + 2L + 3L, 19 families, sorted by L then E)
│   └── README.md
├── docs/
│   ├── Test-Expand.md           # Expansion verification results
│   ├── Test-Relation.md         # Relation reconstruction results
│   ├── IBPVerification.md       # Three verification methods: CompareVerify, EquationVerify, SeriesVerify
│   ├── Verify-MMA-KIRA-Guide.md # Kira cross-validation guide
│   └── CPP-KiraVerify-Debugger.md
├── scripts/
│   └── NuVerify-Relations.wl    # ν-sampling Kira verification of C++ relations
└── results/bub00/               # Verified result snapshots (.m files for MMA consumption)
```

## External References

- FireFly Library: https://github.com/firefly-library/firefly
- Eigen Documentation: https://eigen.tuxfamily.org/
