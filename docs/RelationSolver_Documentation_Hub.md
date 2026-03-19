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

### 3. **Workspace Instructions** ([../.github/copilot-instructions.md](../.github/copilot-instructions.md))
Project-level guidance including:
- Build and run commands
- File organization
- Coding conventions
- RelationSolver-specific prompts
- Where to start looking

**Best for**: Project context, conventions, build system

### 4. **Architecture Guide** ([../AGENTS.md](../AGENTS.md))
High-level project overview:
- Type system and data structures
- RelationSolver specification section
- Full algorithm flow
- Testing strategy
- Integration with other modules

**Best for**: System design, module relationships, overall project context

### 5. **Test Suite** ([../tests/test_RelationFF.cpp](../tests/test_RelationFF.cpp))
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
3. **test_RelationFF.cpp** (15 min) → Study real code
4. **Architecture Guide** (10 min) → Understand overall design

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
📝 **Example**: test_RelationFF.cpp lines 215-240  
⏱️ **Time**: 15 minutes

### Task 2: Debug Convergence Issues
📖 **Reference**: Quick Reference → Diagnostic Output  
📖 **Reference**: Component Guide → Common Pitfalls  
💻 **Code**: Look for `convergence_threshold` tuning  
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
🧪 **Test**: Update `test_RelationFF.cpp` or create new test  
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
│   └── RelationSolver.hpp               (Source code - 1700+ lines)
├── include/
│   └── IncrementalRelationSolver.hpp    (Advanced multi-level solver)
├── tests/
│   ├── test_RelationFF.cpp              (Main test - uses FFInt)
│   └── test_RelationNew.cpp             (Extended test)
└── docs/
    ├── RelationSolver_ComponentGuide.md (THIS DOCUMENTATION)
    ├── RelationSolver_QuickReference.md (Quick lookup)
    └── RelationSolver_Documentation.md  (Linked from Component Guide)
```

---

## ⚙️ Configuration Reference

### `AdaptiveSamplingConfig` Structure

```cpp
struct AdaptiveSamplingConfig {
    int min_nu = 3;                    // Minimum sampling points
    int max_nu = 50;                   // Maximum sampling points  
    int convergence_threshold = 3;     // Nullity stability count
    int lev_hint = 0;                  // (lev, deg) hints for
    int deg_hint = 0;                  // dynamic adjustment
    // ... additional convergence params
};
```

### Common Settings

| Use Case | min_nu | max_nu | convergence_threshold |
|----------|--------|--------|----------------------|
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
./test_RelationFF          # Main finite field test
./test_RelationNew         # Extended multi-regime test
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

- **Workspace Instructions**: [.github/copilot-instructions.md](../.github/copilot-instructions.md)
- **Project Architecture**: [AGENTS.md](../AGENTS.md)
- **Build System**: [CMakeLists.txt](../CMakeLists.txt)
- **Example Code**: [tests/test_RelationFF.cpp](../tests/test_RelationFF.cpp)

---

## 🎯 Summary

This documentation hub provides multiple perspectives on RelationSolver:
- **Quick Reference**: For fast lookups and code examples
- **Component Guide**: For deep dives into API and algorithm
- **Test Suite**: For real-world working examples
- **Architecture Guide**: For system-wide context

**Choose your starting point based on your needs** and follow the recommended reading order for your use case.

---

*Last updated: 2026-03-19*  
*Documentation version: 1.0*  
*Covers: RelationSolver.hpp, test_RelationFF.cpp*
