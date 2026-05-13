#!/bin/bash
set -e

echo "=== Testing Box relation reconstruction (order=6, lev=2, deg=2) ==="

# Save current state
cp data/IBPMat_Box.bin data/IBPMat_Box-current.bin 2>/dev/null || true
cp data/RingData_Box.bin data/RingData_Box-current.bin 2>/dev/null || true

# MMA version
cp verify/Box/IBPMat_Box-MMA.bin data/IBPMat_Box.bin
cp verify/Box/RingData_Box-MMA.bin data/RingData_Box.bin

echo "--- MMA version ---"
time ./build/test_relationFF Box 6 2 2 2 > /tmp/test_relation_box_k6_MMA.txt 2>&1
cp relations/AllRelations_Box_k6.m relations/AllRelations_Box_k6_MMA.m 2>/dev/null || true

# Use Singular version
cp verify/DP323_recompute/IBPMat_Box-Singular.bin data/IBPMat_Box.bin
cp verify/DP323_recompute/RingData_Box-Singular.bin data/RingData_Box.bin

echo ""
echo "--- Singular version ---"
time ./build/test_relationFF Box 6 2 2 2 > /tmp/test_relation_box_k6_Singular.txt 2>&1
cp relations/AllRelations_Box_k6.m relations/AllRelations_Box_k6_Singular.m 2>/dev/null || true

# Restore
cp data/IBPMat_Box-current.bin data/IBPMat_Box.bin 2>/dev/null || true
cp data/RingData_Box-current.bin data/RingData_Box.bin 2>/dev/null || true

echo ""
echo "=== Results ==="
echo "MMA:"
grep "sol_dim=" /tmp/test_relation_box_k6_MMA.txt
echo ""
echo "Singular:"
grep "sol_dim=" /tmp/test_relation_box_k6_Singular.txt
