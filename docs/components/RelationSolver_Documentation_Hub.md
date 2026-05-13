# RelationSolver Documentation Hub

## 📚 Overview

This is the central documentation hub for the **RelationSolver** module - the linear relation reconstruction engine of the IBP matrix expansion project.

---

## 📖 Documentation Hierarchy

### 1. **Quick Reference** ([RelationSolver_QuickReference.md](./RelationSolver_QuickReference.md))
**Start here** for:
- 3-step basic usage pattern
- Key data types table
- Common code patterns
- Configuration tuning guide
- Diagnostic output examples

**Best for**: Quick lookups, code examples, troubleshooting

### 2. **Component Guide** ([RelationSolver_ComponentGuide.md](./RelationSolver_ComponentGuide.md))
Comprehensive documentation covering:
- Architecture and design decisions
- Complete API reference
- Step-by-step algorithm flow
- Full usage examples from test suite
- Performance considerations
- Common pitfalls and solutions

**Best for**: In-depth understanding, implementation details, algorithm explanation

### 3. **Algorithm Theory** ([algorithms/ReconstructAlgorithm.md](algorithms/ReconstructAlgorithm.md))
Mathematical framework and implementation comparison:
- Core expansion formula and relation equations
- MMA symbolic vs C++ numerical approaches
- Per-point evaluation details (computeG, computeF1, computeF2)
- ν-sampling verification theory

**Best for**: Understanding the math behind the solver, algorithm correctness

> **Note**: `../.github/copilot-instructions.md` (previously referenced) does not exist.

### 4. **Architecture Guide** ([../AGENTS.md](../AGENTS.md))
High-level project overview:
- Type system and data structures
- RelationSolver specification section
- Full algorithm flow
- Testing strategy
- Integration with other modules

**Best for**: System design, module relationships, overall project context

### 5. **Ansatz Modes** ([AnsatzModes.md](./AnsatzModes.md))
Multi-mode α-seed generation for LIE reduction:
- Four ansatz modes: Pyramid, DotPyramid, Star, ExtendedPyramid
- Per-mode kernel sign (ν−α vs ν+α), α generation, level filter
- CLI usage with `--mode` and `--sector` flags
- Test results and verification checklist

**Best for**: Understanding ansatz types, choosing the right mode, CLI reference

### 6. **Test Suite** ([../tests/test_relationFF.cpp](../tests/test_relationFF.cpp))
**Living documentation** - practical, working examples:
- Multi-regime setup
- Adaptive config tuning (`adjustSamplingConfig`)
- Complete solve loop
- Result processing and analysis

**Best for**: Learning by example, real-world patterns

---

## 🔄 Recommended Reading Order

### For First-Time Users
1. **Quick Reference** (5 min) → Get basic usage
2. **Component Guide - Usage Example** (10 min) → See full context
3. **algorithms/ReconstructAlgorithm.md** (15 min) → Understand the math and algorithm
4. **test_relationFF.cpp** (15 min) → Study real code

### For Implementers
1. **Component Guide - Algorithm Flow** (15 min) → Understand internals
2. **RelationSolver.hpp** source code (30 min) → Review implementation
3. **Architecture Guide - Integration** (10 min) → See dependencies
4. **Common Tasks** section (varies) → Choose your task

### For Maintainers
1. **Architecture Guide** (full) → System design
2. **Component Guide** (full) → Component details
3. **Performance Considerations** section (10 min)
4. **Test strategy** and diagnostics

---

## 💡 Key Concepts (At a Glance)

### Main Function
```cpp
reconstructReductionRelation<T>()
```
Solves for linear relations at a single (lev, deg) level.

### Input
- `CTable`: Expansion coefficients from layer recursion
- `sector`: Sector/regime identifiers
- `A_list`, `Ainv_list`: Operator matrices
- `ne`, `lev`, `deg`: Dimensions and scope
- `config`: Adaptive sampling configuration

### Output
- `LinearSystemResult<T>`: Solution structure (hasSolution, Mext, S)
- `RelationCoefficient<T>`: Indexed coefficient accessor

### Algorithm
1. Prepare P(α) matrices for each regime
2. Iteratively sample and build equations
3. Solve homogeneous linear system
4. Check nullity convergence
5. Package result as RelationCoefficient

---

## 🚀 Common Tasks

### Task 1: Use RelationSolver in Your Code
📖 **Reference**: Quick Reference → Basic Usage (3 Steps)  
📝 **Example**: test_relationFF.cpp — `main()` solver loop at ~line 438, `adjustSamplingConfig()` at ~line 22  
⏱️ **Time**: 15 minutes

### Task 2: Debug Convergence Issues
📖 **Reference**: Quick Reference → Diagnostic Output  
📖 **Reference**: Component Guide → Common Pitfalls  
💻 **Code**: Look for `nullity_stable_threshold` tuning  
⏱️ **Time**: 10 minutes

### Task 3: Optimize Performance
📖 **Reference**: Component Guide → Performance Considerations  
💻 **Code**: Tune `min_nu`, `max_nu`, `safety_factor` in config  
🔍 **Profile**: Run with diagnostic output enabled  
⏱️ **Time**: 30 minutes

### Task 4: Add New Feature to RelationSolver
📖 **Reference**: Component Guide (full)  
📖 **Reference**: Architecture Guide → Code Organization  
💻 **Code**: Modify `include/RelationSolver.hpp`  
🧪 **Test**: Update `test_relationFF.cpp` or create new test  
⏱️ **Time**: 1-2 hours

### Task 5: Understand the Algorithm
📖 **Reference**: Component Guide → Algorithm Flow (Step-by-Step)  
📖 **Reference**: Component Guide → Key Algorithms  
📖 **Reference**: AGENTS.md → RelationSolver Specification → Algorithm  
⏱️ **Time**: 45 minutes

---

## 🔧 File Structure

```
RelationSolver Module Files:
├── include/
│   ├── RelationSolver.hpp                (Source code - 1700+ lines)
│   └── IncrementalRelationSolver.hpp     (Advanced multi-level solver)
├── tests/
│   ├── test_relationFF.cpp               (Main test - uses FFInt)
│   └── tests/archive/test_RelationNew.cpp (Extended test, archived)
└── docs/
    ├── RelationSolver_ComponentGuide.md  (Comprehensive API + algorithm)
    ├── RelationSolver_QuickReference.md  (Quick lookup)
    ├── AnsatzModes.md                    (Multi-mode ansatz guide + test results)
    └── algorithms/ReconstructAlgorithm.md          (MMA vs C++ algorithm comparison)
```

---

## ⚙️ Configuration Reference

### `AdaptiveSamplingConfig` Structure

```cpp
struct AdaptiveSamplingConfig {
    int min_nu = 3;                    // Minimum sampling points
    int max_nu = 50;                   // Maximum sampling points
    int check_interval = 1;            // Rank check interval (every N points)
    int nullity_stable_threshold = 3;  // Nullity stability count (ν-space)
    int verification_points = 3;       // Extra verification points
    int plateau_size = 1;              // Extra orders for stability confirmation
    double tolerance = 1e-10;          // Numerical tolerance
    bool use_special_points = true;    // Use unit vectors, all-1s, etc.
    double random_min = 3.0;           // Random point min (avoid 0)
    double random_max = 100.0;         // Random point max
    int lev_hint = 2;                  // (lev, deg) hints for
    int deg_hint = 2;                  // dynamic adjustment
};
```

### Common Settings

| Use Case | min_nu | max_nu | nullity_stable_threshold |
|----------|--------|--------|--------------------------|
| Small systems (lev ≤ 2) | 3 | 20 | 2 |
| Medium systems (lev ≤ 4) | 5 | 100 | 3 |
| Large systems (lev ≥ 5) | 10 | 200+ | 4 |
| Memory-constrained | 3 | 50 | 3 |
| Finite field (recommended) | 5 | 150 | 3 |

---

## 🧪 Testing

### Run the Tests
```bash
cd build
./test_relationFF bub00 4 1 2 2          # Pyramid (default)
./test_relationFF bub00 4 1 2 2 --mode 1 # DotPyramid
./test_relationFF bub00 4 1 2 2 --mode 2 # Star
./test_relationFF bub00 4 0 2 2 --mode 3 # ExtendedPyramid
# test_RelationNew                         -- archived in tests/archive/, not in current build
```

### Expected Output
```
Extracted 4 sectors from ring data.
=== Iterative Relation Solving ===
lev range: [0, 3]
deg range: [0, 3]

[lev=0, deg=0]: Relations found. Nullity=2
[lev=0, deg=1]: No relations (full rank)
...
Total time: 5.23 seconds
```

### Troubleshooting Tests
- **Crash on startup**: Check data file paths (IBPMat_*.bin, RingData_*.bin)
- **Slow convergence**: Increase `max_nu` in config
- **Incorrect results**: Verify A matrices are invertible
- **Memory issues**: Process (lev, deg) one at a time

---

## 📞 Integration Guide

### With Layer Recursion
```
LayerRecursion.hpp
  ↓ produces
SeriesCoefficient<T> (CTable)
  ↓ fed into
RelationSolver::reconstructReductionRelation<T>()
```

### With Linear Solvers
```
RelationSolver builds matrix → LinearSolver::solve() → Solution
Config determines size of matrix
```

### With Incremental Solver
```
RelationSolver::reconstructReductionRelation()
  ↓ called by
IncrementalRelationSolver (for multi-level problems)
```

---

## 🔗 Cross-References

- **Linear System Solving**: See LinearSolver_Eigen.hpp, LinearSolver_FF.hpp
- **Coefficient Storage**: See SeriesCoefficient.hpp
- **Multi-Index Tools**: See Combinatorics.hpp
- **Data Loading**: See IBPMatrixLoader_Binary.hpp, RingDataLoader.hpp

---

## 📋 Quick Checklist: Before Using RelationSolver

- [ ] Have expansion coefficients from layer recursion ready
- [ ] Know your `ne` (multi-index dimension)
- [ ] Have A-operators and their inverses
- [ ] Decide on `lev` and `deg` scope
- [ ] Choose `AdaptiveSamplingConfig` settings
- [ ] For FFInt: Set prime via `FFInt::set_new_prime()`
- [ ] For testing: Verify binary data files exist

---

## 📚 Additional Resources

- **Algorithm Theory**: [algorithms/ReconstructAlgorithm.md](algorithms/ReconstructAlgorithm.md) — MMA vs C++ comparison, per-point evaluation
- **Ansatz Modes**: [AnsatzModes.md](AnsatzModes.md) — Pyramid/DotPyramid/Star/ExtendedPyramid guide
- **Project Architecture**: [AGENTS.md](../AGENTS.md)
- **Build System**: [CMakeLists.txt](../CMakeLists.txt)
- **Example Code**: [tests/test_relationFF.cpp](../tests/test_relationFF.cpp)

---

## 🎯 Summary

This documentation hub provides multiple perspectives on RelationSolver:
- **Quick Reference**: For fast lookups and code examples
- **Component Guide**: For deep dives into API and algorithm
- **Test Suite**: For real-world working examples
- **Architecture Guide**: For system-wide context

**Choose your starting point based on your needs** and follow the recommended reading order for your use case.

---

*Last updated: 2026-05-07 (multi-mode ansatz + test results added)*  
*Documentation version: 1.3*  
*Covers: RelationSolver.hpp, test_relationFF.cpp, AnsatzModes.md*

## 补充阅读

- `SymbolicRuleAlgorithm.md` — 符号规则生成与完整性标志
- `Benchmark_Results.md` — 全族 benchmark 时序/内存数据
- `Strategy-Comparison.md` — Strategy 0 vs 4 性能对比
