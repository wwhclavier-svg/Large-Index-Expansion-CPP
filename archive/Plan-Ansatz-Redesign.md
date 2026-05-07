Multi-Mode Ansatz Support in C++ RelationSolver
Context
The C++ RelationSolver currently hardcodes a single ansatz mode (Pyramid). The MMA reference implementation LIEReconstruct.wl supports three modes via the "AnsatzMode" option. Adding these to C++ enables cross-validation against MMA for all modes and supports families where DotPyramid or Star may produce better results.

Files to Modify
include/RelationSolver.hpp — all core changes (~2449 lines)
tests/test_relationFF.cpp — CLI and export changes
Implementation Steps
Step 1: Add AnsatzMode enum
Insert after namespace RelationSolver { (~line 33):

enum class AnsatzMode {
    Pyramid    = 0,
    DotPyramid = 1,
    Star       = 2
};
Step 2: Add alphaLevel() helper
New inline function after generateAllIndices (~line 88): returns the mode-specific "level" of an alpha vector — Total for Pyramid/DotPyramid, Total[Abs] for Star.

Step 3: Add negatedKernel_ to RegimeEvaluator
In step2_computeG (line 784), change:

T diff = nu[m] - static_cast<T>(alpha[m]);
to:

T diff = negatedKernel_ ? nu[m] + static_cast<T>(alpha[m])
                         : nu[m] - static_cast<T>(alpha[m]);
Add bool negatedKernel_ = false member and setNegatedKernel(bool) to RegimeEvaluator<T>.

Step 4: Propagate negatedKernel through GlobalEquationAssembler
Add kernel_negated_ member to GlobalEquationAssembler<T>
Add bool kernel_negated = false parameter to init() (line 1362)
In addRegime() (line 1374), call evaluator.setNegatedKernel(kernel_negated_)
Step 5: Propagate negatedKernel through AdaptiveEquationBuilder::build
Add bool kernel_negated = false parameter to build() (line 1527). Pass it to assembler.init() call at line 1576.

Step 6: Add generateAlphaSeeds() helper
New function after generateAllIndices: encapsulates per-mode alpha list generation.

Pyramid: generateAllIndices(ne, lev, ..., false) — non-negative, L1 ≤ lev
DotPyramid: generateAllIndices(ne, lev+1, ..., false) — non-negative, L1 ≤ lev+1
Star: Iterate all 2^ne sector subsets. For each subset, split positions into pdSet (1s) and ispSet (0s). Generate dotSeeds (L1 ≤ lev, length |pdSet|) and rankSeeds (L1 ≤ lev, length |ispSet|). Combine: seed[pd_i] = 1 + dotSeed, seed[isp_j] = -rankSeed. Deduplicate via set. Sort by L1 norm then lexicographic.
generateAllIndices's existing allow_neg parameter stays unused (Star uses two separate non-negative generator calls, not a single signed generator).

Step 7: Modify reconstructAllRelations signature and logic
Add AnsatzMode mode = AnsatzMode::Pyramid parameter (line 2175).

Inside the function:

Alpha pre-generation: Replace lines 2202-2205 with call to generateAlphaSeeds(mode, ne, lev_max, first_sector, alphas_max).
Alpha per-level filtering (line 2224-2228): Replace sum(a) <= lev with alphaLevel(mode, a) <= lev.
Level range: For DotPyramid, lev_min = max(1, lev_min). For Star, lev_min = max(1, lev_min) (no zero seed).
Kernel sign: Determine bool kernel_negated = (mode != AnsatzMode::Pyramid). Pass to builder.build().
Star-specific: The sector parameter (already available as sector_list) is needed. For Star, use sector[0] to determine pdSet/ispSet split for the top-level regime (C++ currently uses reconstructAllRelations with the same sector_list as input).
Step 8: Modify test_relationFF.cpp
Add --mode <0|1|2> CLI parsing in main() (line 470-495)
Pass mode to reconstructAllRelations() (line 641)
Add AnsatzMode field to the .m export metadata in exportAllResultsToMMA_SingleFile() (line 376-382)
Print mode in configuration output (line 504-509)
Step 9: Boundary Terms for Pyramid (DEFERRED)
The MMA Pyramid mode generates boundary terms for ISP positions (where limitSector[i] == 0) at |alpha| = lev+1. These modify coefficient indices: b[rk, beta] ← b[rk+e_i, beta-e_i]. This is complex and is ONLY needed when ISPs exist (not in topsector-only mode). Defer this to a follow-up PR.

Verification
Backward compatibility: ./test_relationFF bub 4 1 2 2 (default Pyramid) produces identical .m output as before
Mode 1: ./test_relationFF bub 4 1 2 2 --mode 1 — verify nullity/relation counts are sensible, compare against MMA DotPyramid run
Mode 2: ./test_relationFF bub 4 1 2 2 --mode 2 — same verification
Topsector interaction: --topsector --mode 1 / --topsector --mode 2 with bub family
Build: cd build && cmake --build . --target test_relationFF
Risk Notes
Star mode seed count: For ne=5, lev=2, Star generates thousands of alphas (all subsets × dot/rank combinations). The matrix can be large. Verify RemoveSolvedVariables still reduces variables effectively.
RemoveSolvedVariables with mixed-sign: Componentwise dominance (alpha >= alpha_s) works mathematically for mixed-sign vectors but the Star mode may filter differently than MMA expects. Verification will catch discrepancies.
p-factor sign: For DotPyramid/Star, the kernel uses nu + alpha instead of nu - alpha, but the p-factor p(alpha) = A^{-alpha} computation is unchanged. Since equations are homogeneous, any overall sign difference is absorbed into unknown coefficients.