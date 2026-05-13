#!/bin/bash
set -e

echo "=== Testing Box relation reconstruction (order=2) ==="
echo ""

# Save current state
cp data/IBPMat_Box.bin data/IBPMat_Box-current.bin 2>/dev/null || true
cp data/RingData_Box.bin data/RingData_Box-current.bin 2>/dev/null || true

# Use MMA version
cp verify/Box/IBPMat_Box-MMA.bin data/IBPMat_Box.bin
cp verify/Box/RingData_Box-MMA.bin data/RingData_Box.bin

echo "--- MMA version ---"
time ./build/test_relationFF Box 2 1 2 2 > /tmp/test_relation_box_k2_MMA.txt 2>&1
cp relations/AllRelations_Box_k2.m relations/AllRelations_Box_k2_MMA.m 2>/dev/null || true

# Use Singular version
cp verify/DP323_recompute/IBPMat_Box-Singular.bin data/IBPMat_Box.bin
cp verify/DP323_recompute/RingData_Box-Singular.bin data/RingData_Box.bin

echo ""
echo "--- Singular version ---"
time ./build/test_relationFF Box 2 1 2 2 > /tmp/test_relation_box_k2_Singular.txt 2>&1
cp relations/AllRelations_Box_k2.m relations/AllRelations_Box_k2_Singular.m 2>/dev/null || true

# Restore
cp data/IBPMat_Box-current.bin data/IBPMat_Box.bin 2>/dev/null || true
cp data/RingData_Box-current.bin data/RingData_Box.bin 2>/dev/null || true

echo ""
echo "=== Relation counts ==="
echo "MMA:"
grep "sol_dim=" /tmp/test_relation_box_k2_MMA.txt
echo ""
echo "Singular:"
grep "sol_dim=" /tmp/test_relation_box_k2_Singular.txt

echo ""
echo "=== File comparison ==="
if [ -f relations/AllRelations_Box_k2_MMA.m ] && [ -f relations/AllRelations_Box_k2_Singular.m ]; then
  diff relations/AllRelations_Box_k2_MMA.m relations/AllRelations_Box_k2_Singular.m > /dev/null && echo "FILES IDENTICAL" || echo "FILES DIFFERENT"
  lines_MMA=$(wc -l < relations/AllRelations_Box_k2_MMA.m)
  lines_Sing=$(wc -l < relations/AllRelations_Box_k2_Singular.m)
  echo "Lines: MMA=$lines_MMA Sing=$lines_Sing"
  
  # Compare NumRelations per config
  echo ""
  echo "MMA NumRelations:"
  grep "NumRelations" relations/AllRelations_Box_k2_MMA.m
  echo ""
  echo "Singular NumRelations:"
  grep "NumRelations" relations/AllRelations_Box_k2_Singular.m
else
  echo "Relation files not found"
fi
