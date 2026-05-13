# Strategy Optimization Implementation Plan

## Goal
Add `--strategy <n>` CLI option supporting 5 equation-solving strategies for efficiency optimization.

## Strategies

| ID | Name | Mechanic | Cost |
|----|------|----------|------|
| 0 | Default | All regimes at once, all rows (current) | Baseline |
| 1 | BlockAware | Per-basis-block rank check: try block 0 first, skip regime if rank unchanged | ~30 lines |
| 2 | PerOrderOnline | Per-ν, add orders 0→k_max incrementally; stop after 2 consecutive nullity-stable orders | ~20 lines |
| 3 | Combined | Strategy 1 + Strategy 2 | ~10 lines |
| 4 | PerRegimeRankCheck | Sort regimes by |θ| desc, per-regime rank check (original Strategy 4) | ~40 lines |

## Files Changed

### `include/RelationSolver.hpp`
1. Add `Strategy` enum + `strategyToString()`
2. Add `evaluateRegimeAtNu()` to GlobalEquationAssembler
3. Add `strategy_` field + `setStrategy()` to AdaptiveEquationBuilder
4. Add `θWeight()` helper (counts 1s in sector vector)
5. Modify `build()` with strategy branching

### `tests/test_relationFF.cpp`
1. Parse `--strategy <n>` argument
2. Thread strategy enum through `reconstructAllRelations()`

## Strategy Details

### Strategy 1: BlockAware
```
for each ν:
    for each regime r (any order):
        rows = assembler.evaluateRegimeAtNu(r, nu)
        block0 = rows[0:k_max+1]  // first basis only
        null_before = solver.getNullity()
        solver.addRows(block0)
        null_after = solver.getNullity()
        if null_after == null_before: continue  // skip regime
        // block0 contributed → add remaining blocks
        for b = 1..nb-1:
            solver.addRows(rows[b*(k_max+1):(b+1)*(k_max+1)])
```

### Strategy 2: PerOrderOnline
```
for each ν:
    rows = assembler.evaluateAtNu(nu)
    rows_by_order = assembler.splitRowsByOrder(rows)
    stable_count = 0
    for order = 0..k_max:
        null_before = solver.getNullity()
        solver.addRows(rows_by_order[order])
        null_after = solver.getNullity()
        if null_after == null_before:
            stable_count++
            if stable_count >= 2: break  // skip remaining orders
        else:
            stable_count = 0
```

### Strategy 4: PerRegimeRankCheck
```
pre-sort regimes by |θ| descending
for each ν:
    for each regime r (|θ| desc):
        rows = assembler.evaluateRegimeAtNu(r, nu)
        null_before = solver.getNullity()
        solver.addRows(rows)
        null_after = solver.getNullity()
        // rows always added; rank check is informational
        // but the |θ| sort ensures high-value regimes come first
```
Note: Unlike the original Strategy 4 design, we can't easily "discard" rows after addRows(). The benefit comes from |θ|-sorted regime ordering, which maximizes early rank determination, reducing the number of ν needed for convergence.
