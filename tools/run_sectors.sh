#!/bin/bash
# run_sectors.sh <family> <sector_list_file> <skip_key> [per_sector_timeout_s]
# 逐 sector 跑 family_generate --sector, 结果写入 cache/<fam>_sectors/
# 跳过 skip_key(由整体运行在攻坚的难 sector), 超时自动放弃该 sector 继续下一个
set -u
FAM="$1"; LIST="$2"; SKIP="$3"; TMO="${4:-1200}"
cd /home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP
CACHE="cache/${FAM}_sectors"
mkdir -p "$CACHE"
LOG="logs/${FAM,,}_sector_loop.log"
echo "=== $(date '+%F %T') start $FAM, skip=$SKIP, timeout=${TMO}s ===" >> "$LOG"
while read -r KEY; do
    [ -z "$KEY" ] && continue
    if [ "$KEY" = "$SKIP" ]; then
        echo "$(date '+%T') SKIP $KEY (hard, handled elsewhere)" >> "$LOG"
        continue
    fi
    if [ -s "$CACHE/IBPMat_${FAM}_sector_${KEY}.bin" ]; then
        echo "$(date '+%T') HAVE $KEY, skip" >> "$LOG"
        continue
    fi
    echo "$(date '+%T') RUN $KEY" >> "$LOG"
    timeout "$TMO" build/family_generate "families/${FAM}.json" --sector "$KEY" --cache-dir "$CACHE" >> "$LOG" 2>&1
    RC=$?
    if [ $RC -eq 124 ]; then
        echo "$(date '+%T') TIMEOUT $KEY (${TMO}s)" >> "$LOG"
    elif [ $RC -ne 0 ]; then
        echo "$(date '+%T') FAIL $KEY rc=$RC" >> "$LOG"
    else
        echo "$(date '+%T') OK $KEY" >> "$LOG"
    fi
done < "$LIST"
echo "=== $(date '+%F %T') loop done $FAM ===" >> "$LOG"
