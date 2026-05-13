# 会话存档：large-index-cpp-20260510

**生成时间**: 2026-05-10 12:32  
**工作目录**: /home/ykm/Large-Index-Expansion-CPP

---

## 前提

IBP矩阵(bub00等)已生成

## 目标

完成有限域展开和关系重建验证

## 计划

layerRecursion→reconstructAllRelations→验证

## 目录设计

```
/:
  [.agents/]
  [.claude/]
  [archive/]
  [docs/]
  [families/]
  [include/]
  [mma/]
  [scripts/]
  [src/]
  [tests/]
  [tools/]
  [verify/]

tools:
  diagnose_ibp.cpp
  diagnose_ibp_reference.cpp
  family_generate.cpp
  test_ff_verify.cpp
  test_firefly_simple.cpp
  test_region_solver.cpp
  test_singular_runner.cpp

src:
  LayerRecursionCore.cpp
  layerRecursion.cpp
  main.cpp

.claude:
  [skills/]

.claude/skills:
  [blade/]
  [large-index-expansion/]

include:
  [firefly/]
  BinaryIBPWriter.hpp
  BinaryRingWriter.hpp
  Combinatorics.hpp
  FamilyConfig.hpp
  IBPAnalyzer.hpp
  IBPEqGenerator.hpp
  IBPMatrixLoader.hpp
  IBPMatrixLoader_Binary.hpp
  IncrementDetector.hpp
  IncrementalRelationSolver.hpp
  LayerRecursion.hpp
  LayerRecursionCore.hpp
  LayerRecursionCore.tpp
  LinearSolver.hpp
  LinearSolver_Eigen.hpp
  LinearSolver_FF.hpp
  ParallelSolver.hpp
  PolyArith.hpp
  RecursionBuilder.hpp
  RegionSolver.hpp
  RelationSolver.hpp
  RingBuilder.hpp
  RingDataLoader.hpp
  SeriesCoefficient.hpp

include/firefly:
  FFInt.hpp

verify:
  [Box/]
  [DB313/]
  [DP323/]
  [DP323_recompute/]
  [FamilyDatabase/]
  [NP222/]
  [NP322/]
  [Penta1L/]
  [SR212/]
  [SR212-3m/]
  [SR212-5m/]
  [TB123/]
  [TB123m/]
  [Tri/]
  [VerifyUtility/]
  [bub00/]
  [bub10/]
  [bub11/]
  [docs/]
  [logs/]
  [timing/]
  benchmark_prepare.sh

verify/DP323:
  [SingularTempFile/]
  [backup_before_recompute/]

verify/VerifyUtility:
  [SingularTempFile/]

verify/DP323_recompute:
  [SingularTempFile/]
  [output/]

tests:
  IBPVerification.hpp
  compare_ringdata_detailed.cpp
  debug_region_solver.cpp
  test_IBPVerification.cpp
  test_compare_MMA.cpp
  test_expandFF.cpp
  test_expand_family.cpp
  test_family_config.cpp
  test_family_generate.cpp
  test_ibp_analyzer.cpp
  test_load_bub.cpp
  test_poly_arith.cpp
  test_recursion_builder.cpp
  test_recursion_builder_tri.cpp
  test_region_solver_full.cpp
  test_relationFF.cpp
  test_ring_builder.cpp
  verify_MMA_lev1_deg1.cpp
  verify_ring_builder.cpp

.agents:
  [skills/]

.agents/skills:
  [expand-verify/]
  [large-index-expansion/]
  [relation-verify/]

```

## 完成状态

核心完成，验证框架搭建中

## 测试结果

bub00 k2 lev1/deg1 PASS

## 下一步

IncrementDetector联调+NP322大矩阵

## 预期测试

EquationVerify+CompareVerify PASS

---

## 自动提取信息

**Git分支**: master  
**最近提交**:
```
e90e9a2 Add IncrementDetector (rank-based incre auto-detection) + verification pipeline restructure
81e124c Update DEV_STATUS and verification docs: reflect May 5-7 project state
99ec194 Add DEV_STATUS.md: comprehensive project status for May 5-7 work
8c98ba2 Add multi-mode ansatz support to RelationSolver
9d5330e Update Benchmark_Results: add ne=7 families (TB123/NP222) timing, expand-vs-solve crossover analysis, and lev=3 full-sector results
8391228 Refactor IncrementalNullspaceSolver to incremental RREF maintenance and simplify nimax_lists
baa0534 Fix .gitignore and add Compare-* utility scripts
8f0df6b Fix VerifyExpand-Prepare: binary export shadowed by package context
```

**变更文件**:
```
M CMakeLists.txt
 M include/RegionSolver.hpp
 M include/SingularRunner.hpp
 M verify/VerifyUtility/LIECoreAlgebra.wl
 M verify/VerifyUtility/LIERegions.wl
 M verify/VerifyUtility/VerifyExpand-Compare.wl
?? .aether/
?? .hermes/
?? AllRelations_Box_k2.m
?? AllRelations_NP222_k2.m
```

**构建目标**: test_expandFF, test_IBPVerification, test_relationFF, test_load_bub, test_expand_family, test_ff_verify, test_family_generate, test_family_config, family_generate, test_singular_runner, test_region_solver, test_ibp_analyzer, test_recursion_builder, test_recursion_builder_tri, test_ring_builder  
**测试文件**: test_IBPVerification.cpp, test_compare_MMA.cpp, test_expandFF.cpp, test_expand_family.cpp, test_family_config.cpp, test_family_generate.cpp, test_ibp_analyzer.cpp, test_load_bub.cpp, test_poly_arith.cpp, test_recursion_builder.cpp, test_recursion_builder_tri.cpp, test_region_solver_full.cpp, test_relationFF.cpp, test_ring_builder.cpp  
**文档**: Plan-FamilyGenerate-CPP-Rewrite.md, RelationSolver_Documentation_Hub.md, LayerRecursion_Algorithm.md, AnsatzModes.md, RelationSolver_ComponentGuide.md, RelationSolver_QuickReference.md, Strategy-Comparison.md, RegionSolver-Debug-Plan.md, Benchmark_Results.md, ReconstructAlgorithm.md, ReconstructReductionRelation_Documentation.md

