#!/bin/bash
# LIE benchmark runner — 1L + 2L2P families only
# Usage: ./bench/run_lie_bench.sh [order=4] [lev_max=2] [deg_max=2]

ORDER=${1:-4}
LEV_MAX=${2:-2}
DEG_MAX=${3:-2}
OUTDIR="results/lie"
FAILED=0

# 1L (L=1) + 2L2P (L=2, 2-point, #extmom=1)
FAMILIES=(
    "bub00" "bub10" "bub11" "Tri" "Box"
    "Penta1L" "Penta1L-5m"
    "SR212" "SR212-3m" "SR212-5m"
)

mkdir -p "$OUTDIR"

echo "=== LIE Benchmark Batch Run ==="
echo "Order=$ORDER  lev_max=$LEV_MAX  deg_max=$DEG_MAX"
echo "Output: $OUTDIR/"
echo ""

for family in "${FAMILIES[@]}"; do
    ibp="data/IBPMat_${family}.bin"
    ring="data/RingData_${family}.bin"
    if [[ ! -f "$ibp" || ! -f "$ring" ]]; then
        echo "  [SKIP] $family  (missing $ibp or $ring)"
        continue
    fi
    echo "  [RUN]  $family  (order=$ORDER, lev=1..$LEV_MAX, deg=0..$DEG_MAX)"
    start=$(date +%s%N)
    ./build/bench_lie "$family" "$ORDER" 1 "$LEV_MAX" "$DEG_MAX" --output "$OUTDIR" 2>&1
    rc=$?
    end=$(date +%s%N)
    elapsed_ms=$(( (end - start) / 1000000 ))
    if [[ $rc -eq 0 ]]; then
        echo "  [DONE] $family  (${elapsed_ms}ms)"
    else
        echo "  [FAIL] $family  (exit=$rc, ${elapsed_ms}ms)"
        ((FAILED++))
    fi
    echo ""
done

echo "=== Batch Complete ==="
echo "Failed: $FAILED / ${#FAMILIES[@]}"
exit $FAILED
