#!/bin/bash
# Kira benchmark runner for LIE families.
# Generates kira input from family JSON, runs kira, captures timing.
#
# Usage:
#   ./tools/run_kira_benchmark.sh <family>              # single family
#   ./tools/run_kira_benchmark.sh --all-1l2p             # batch 1L + 2L2P
#
# Environment:
#   KIRA_WORKDIR   - working directory (default: kira_tests)
#   KIRA_PARALLEL  - kira --parallel N (default: 4)

set -euo pipefail

KIRA_WORKDIR="${KIRA_WORKDIR:-$(dirname "$0")/../kira_tests}"
KIRA_PARALLEL="${KIRA_PARALLEL:-4}"
FAILED=0

# 1L + 2L2P families
ONELOOP=(bub00 bub10 bub11 Tri Box Penta1L Penta1L-5m)
TWOLOOP_2POINT=(SR212 SR212-3m SR212-5m)
ALL_1L2P=("${ONELOOP[@]}" "${TWOLOOP_2POINT[@]}")

usage() {
    echo "Usage: $0 [family|--all-1l2p]"
    echo "  KIRA_WORKDIR=$KIRA_WORKDIR  KIRA_PARALLEL=$KIRA_PARALLEL"
    exit 1
}

run_kira() {
    local family="$1"
    local dir="$KIRA_WORKDIR/$family"

    if [[ ! -f "$dir/jobs.yaml" ]]; then
        echo "  [SKIP] $family  (no kira input, run kira_input_gen.py first)"
        return 1
    fi

    echo "  [KIRA] $family  (parallel=$KIRA_PARALLEL, cwd=$dir)"

    mkdir -p "$dir/results"

    local start=$(date +%s%N)
    cd "$dir"
    export FERMATPATH=/usr/local/bin/Ferl7/fer64
    if kira --parallel="$KIRA_PARALLEL" jobs.yaml 2>&1 | tee kira.log; then
        local end=$(date +%s%N)
        local elapsed_ms=$(( (end - start) / 1000000 ))
        echo "  [DONE] $family  (${elapsed_ms}ms)"

        local masters_file="results/$family/masters.final"
        local output_file="results/$family/kira_integrals.m"
        local num_masters=0
        [[ -f "$masters_file" ]] && num_masters=$(wc -l < "$masters_file")

        # Write timing JSON
        cat > "results/timing.json" <<- TIMINGEOF
        {
          "family": "$family",
          "kira_wallclock_ms": $elapsed_ms,
          "num_masters": $num_masters,
          "output_file": "$output_file",
          "timestamp": "$(date -Iseconds)"
        }
TIMINGEOF
        cd - > /dev/null
        return 0
    else
        local end=$(date +%s%N)
        local elapsed_ms=$(( (end - start) / 1000000 ))
        echo "  [FAIL] $family  (exit=$?, ${elapsed_ms}ms)"
        cd - > /dev/null
        return 1
    fi
}

# Parse args
TARGETS=()
case "${1:-}" in
    --all-1l2p)
        TARGETS=("${ALL_1L2P[@]}")
        ;;
    "")
        usage
        ;;
    *)
        TARGETS=("$1")
        ;;
esac

echo "=== Kira Benchmark Batch Run ==="
echo "Workdir: $KIRA_WORKDIR"
echo "Parallel: $KIRA_PARALLEL"
echo "Families: ${TARGETS[*]}"
echo ""

for family in "${TARGETS[@]}"; do
    run_kira "$family" || ((FAILED++))
    echo ""
done

echo "=== Batch Complete ==="
echo "Failed: $FAILED / ${#TARGETS[@]}"
exit $FAILED
