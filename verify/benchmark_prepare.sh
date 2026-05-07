#!/bin/bash
# benchmark_prepare.sh — Run VerifyExpand-Prepare.wl for all families
# Captures timing breakdown (FamilyDefinition, RegionSolving) and memory (Max RSS)
#
# Usage: cd verify/VerifyUtility && bash ../benchmark_prepare.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
UTIL_DIR="$SCRIPT_DIR/VerifyUtility"
RESULTS_FILE="$SCRIPT_DIR/benchmark_results.md"

# All families from FamilyDatabase.wl, ordered by expected complexity (ne, L)
FAMILIES=(
    # 1-Loop, small ne
    "bub00"
    "bub10"
    "bub11"
    "Tri"
    "Box"
    # 2-Loop, medium ne
    "SR212"
    "SR212-3m"
    "SR212-5m"
    "TB123"
    "TB123m"
    # 2-Loop, large ne
    "NP222"
    "DB313"
    "NP322"
    "NP322m"
    "DP323"
    # 3-Loop
    "BN3L"
    "BN3L1P"
    "BN3L1Pm"
)

TIMEOUT_SEC=3600  # 1 hour per family

mkdir -p "$SCRIPT_DIR"/{logs,timing}

echo "============================================="
echo "Benchmark: VerifyExpand-Prepare for all families"
echo "Started: $(date)"
echo "Families: ${#FAMILIES[@]}"
echo "Timeout per family: ${TIMEOUT_SEC}s"
echo "============================================="
echo

# Write results header
cat > "$RESULTS_FILE" << 'HEADER'
# .bin Generation Benchmark — VerifyExpand-Prepare.wl

**Date**: $(date)
**Environment**: wolframscript + Singular 4.3.2, modulus=179424673

## Results Summary

| Family | L | ne | #Props | TopSector | Step1: FamilyDef (s) | Step2: RegionSolve (s) | Total MMA (s) | Wall Clock (s) | Max RSS (MB) | #Sectors | #Regions | Exit |
|--------|---|----|--------|-----------|---------------------|----------------------|---------------|----------------|-------------|----------|----------|------|
HEADER

run_family() {
    local fam="$1"
    local logfile="$SCRIPT_DIR/logs/${fam}.log"
    local timefile="$SCRIPT_DIR/timing/${fam}.time"

    echo "--- $fam ---"
    echo "  Started at $(date)"

    cd "$UTIL_DIR"

    # Run with timeout, capture both stdout/stderr and /usr/bin/time
    local start_ts=$(date +%s)
    timeout $TIMEOUT_SEC /usr/bin/time -v -o "$timefile" \
        wolframscript -file VerifyExpand-Prepare.wl "$fam" \
        > "$logfile" 2>&1
    local exit_code=$?
    local end_ts=$(date +%s)
    local wall_clock=$((end_ts - start_ts))

    # Extract timing from log
    local time_define=""
    local time_regions=""
    local ne_val=""
    local nibp_val=""
    local sectors_val=""
    local regions_val=""

    if [ -f "$logfile" ]; then
        time_define=$(grep -oP 'Time \(Family Definition\):\s*\K[\d.]+' "$logfile" || echo "N/A")
        time_regions=$(grep -oP 'Time \(Region Solving\):\s*\K[\d.]+' "$logfile" || echo "N/A")
        ne_val=$(grep -oP 'NE\s*=\s*\K\d+' "$logfile" || echo "?")
        nibp_val=$(grep -oP 'NIBP\s*=\s*\K\d+' "$logfile" || echo "?")
        sectors_val=$(grep -oP 'Sectors with regions:\s*\K\d+' "$logfile" || echo "?")
        regions_val=$(grep -oP 'Exported\s+\K\d+' "$logfile" || echo "?")
    fi

    # Extract Max RSS from /usr/bin/time output (in KB, convert to MB)
    local max_rss=""
    if [ -f "$timefile" ]; then
        max_rss=$(grep -oP 'Maximum resident set size \(kbytes\):\s*\K\d+' "$timefile" || echo "0")
        if [ "$max_rss" != "0" ] && [ -n "$max_rss" ]; then
            max_rss=$(echo "scale=1; $max_rss / 1024" | bc)
        else
            max_rss="N/A"
        fi
    fi

    # Determine exit status string
    local exit_str="OK"
    if [ "$exit_code" -eq 124 ]; then
        exit_str="TIMEOUT"
    elif [ "$exit_code" -ne 0 ]; then
        exit_str="ERR($exit_code)"
    fi

    # Total MMA time (sum of two steps if both available)
    local total_mma="N/A"
    if [ "$time_define" != "N/A" ] && [ "$time_regions" != "N/A" ]; then
        total_mma=$(echo "scale=2; $time_define + $time_regions" | bc)
    fi

    echo "  NE=$ne_val NIBP=$nibp_val sectors=$sectors_val regions=$regions_val"
    echo "  FamilyDef: ${time_define}s  RegionSolve: ${time_regions}s  Total: ${total_mma}s"
    echo "  Wall: ${wall_clock}s  MaxRSS: ${max_rss}MB  Exit: $exit_str"
    echo

    # Append to results table
    # Get L from family name heuristics (could parse JSON but this is simpler)
    local L="?"
    case "$fam" in
        bub*|Tri|Box) L=1 ;;
        SR*|TB*|NP*|DB*|DP*) L=2 ;;
        BN*) L=3 ;;
    esac

    # ne and props from known data
    local props="?"
    local ne_known="$ne_val"
    local ts="?"
    case "$fam" in
        bub00) props=2; ts="[1,1]" ;;
        bub10) props=2; ts="[1,1]" ;;
        bub11) props=2; ts="[1,1]" ;;
        Tri)   props=3; ts="[1,1,1]" ;;
        Box)   props=4; ts="[1,1,1,1]" ;;
        SR212) props=5; ts="[1,1,1,1,1]" ;;
        SR212-3m) props=5; ts="[1,0,1,0,1]" ;;
        SR212-5m) props=5; ts="[1,1,1,1,1]" ;;
        TB123) props=7; ts="[1,1,1,1,1,1,0]" ;;
        TB123m) props=7; ts="[1,1,1,1,1,1,0]" ;;
        NP222) props=7; ts="[1,1,1,1,1,1,1]" ;;
        DB313) props=9; ts="[1,1,1,1,1,1,1,0,0]" ;;
        NP322) props=9; ts="[1,1,1,1,1,1,1,1,1]" ;;
        NP322m) props=9; ts="[1,1,1,1,1,1,1,1,1]" ;;
        DP323) props=11; ts="[1,1,1,1,1,1,1,1,0,0,0]" ;;
        BN3L) props=6; ts="[1,1,1,1,1,1]" ;;
        BN3L1P) props=9; ts="[1,1,1,1,1,1,1,1,0]" ;;
        BN3L1Pm) props=9; ts="[1,1,1,1,1,1,1,1,0]" ;;
    esac

    printf "| %-10s | %d | %s | %d | %-19s | %-19s | %-20s | %-13s | %-11s | %-13s | %-8s | %-8s | %s |\n" \
        "$fam" "$L" "$ne_known" "$props" "$ts" \
        "$time_define" "$time_regions" "$total_mma" \
        "${wall_clock}s" "$max_rss" "$sectors_val" "$regions_val" "$exit_str" \
        >> "$RESULTS_FILE"
}

# Run each family
for fam in "${FAMILIES[@]}"; do
    run_family "$fam"
done

# Close results file
echo "" >> "$RESULTS_FILE"
echo "## Notes" >> "$RESULTS_FILE"
echo "" >> "$RESULTS_FILE"
echo "- **Step 1 (Family Definition)**: IBP equation generation in MMA + Singular SP2PD/derivative computation" >> "$RESULTS_FILE"
echo "- **Step 2 (Region Solving)**: Singular Groebner basis + primary decomposition (primdecGTZE) per sector" >> "$RESULTS_FILE"
echo "- **Wall Clock**: includes MMA startup (~5-6s overhead), file I/O" >> "$RESULTS_FILE"
echo "- **Max RSS**: peak physical memory usage (includes MMA + Singular subprocesses)" >> "$RESULTS_FILE"
echo "- Modulus: 179424673 (Prime[10000000])" >> "$RESULTS_FILE"

echo ""
echo "============================================="
echo "Benchmark complete: $(date)"
echo "Results: $RESULTS_FILE"
echo "Per-family logs: $SCRIPT_DIR/logs/"
echo "============================================="
