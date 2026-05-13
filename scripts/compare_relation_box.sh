#!/bin/bash
set -e

echo "=== Testing Box relation reconstruction ==="
echo ""

# Save current state
cp data/IBPMat_Box.bin data/IBPMat_Box-current.bin 2>/dev/null || true
cp data/RingData_Box.bin data/RingData_Box-current.bin 2>/dev/null || true

# Use MMA version from verify/Box
cp verify/Box/IBPMat_Box-MMA.bin data/IBPMat_Box.bin
cp verify/Box/RingData_Box-MMA.bin data/RingData_Box.bin

echo "--- MMA version ---"
time ./build/test_relationFF Box 4 1 2 2 > /tmp/test_relation_box_MMA.txt 2>&1
cp relations/AllRelations_Box_k4.m relations/AllRelations_Box_k4_MMA.m 2>/dev/null || true

# Use Singular version
cp verify/DP323_recompute/IBPMat_Box-Singular.bin data/IBPMat_Box.bin
cp verify/DP323_recompute/RingData_Box-Singular.bin data/RingData_Box.bin

echo ""
echo "--- Singular version ---"
time ./build/test_relationFF Box 4 1 2 2 > /tmp/test_relation_box_Singular.txt 2>&1
cp relations/AllRelations_Box_k4.m relations/AllRelations_Box_k4_Singular.m 2>/dev/null || true

# Restore
cp data/IBPMat_Box-current.bin data/IBPMat_Box.bin 2>/dev/null || true
cp data/RingData_Box-current.bin data/RingData_Box.bin 2>/dev/null || true

echo ""
echo "=== Comparing outputs ==="
echo "MMA:"
grep "sol_dim=" /tmp/test_relation_box_MMA.txt | head -10
echo ""
echo "Singular:"
grep "sol_dim=" /tmp/test_relation_box_Singular.txt | head -10

echo ""
echo "=== Relation file comparison ==="
if [ -f relations/AllRelations_Box_k4_MMA.m ] && [ -f relations/AllRelations_Box_k4_Singular.m ]; then
  diff relations/AllRelations_Box_k4_MMA.m relations/AllRelations_Box_k4_Singular.m > /dev/null && echo "FILES IDENTICAL" || echo "FILES DIFFERENT"
  lines_MMA=$(wc -l < relations/AllRelations_Box_k4_MMA.m)
  lines_Sing=$(wc -l < relations/AllRelations_Box_k4_Singular.m)
  echo "Lines: MMA=$lines_MMA Sing=$lines_Sing"
else
  echo "Relation files not found (may have failed to generate)"
fi
