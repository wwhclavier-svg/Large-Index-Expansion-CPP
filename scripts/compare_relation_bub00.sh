#!/bin/bash
set -e

echo "=== Testing bub00 relation reconstruction ==="
echo ""

# Backup MMA binaries
mv data/IBPMat_bub00.bin data/IBPMat_bub00-MMA.bin
mv data/RingData_bub00.bin data/RingData_bub00-MMA.bin

# Run with MMA version
cp data/IBPMat_bub00-MMA.bin data/IBPMat_bub00.bin
cp data/RingData_bub00-MMA.bin data/RingData_bub00.bin

echo "--- MMA version ---"
time ./build/test_relationFF bub00 4 1 2 2 > /tmp/test_relation_bub00_MMA.txt 2>&1

# Run with Singular version
cp verify/DP323_recompute/IBPMat_bub00-Singular.bin data/IBPMat_bub00.bin
cp verify/DP323_recompute/RingData_bub00-Singular.bin data/RingData_bub00.bin

echo ""
echo "--- Singular version ---"
time ./build/test_relationFF bub00 4 1 2 2 > /tmp/test_relation_bub00_Singular.txt 2>&1

# Restore default binaries
cp data/IBPMat_bub00-MMA.bin data/IBPMat_bub00.bin
cp data/RingData_bub00-MMA.bin data/RingData_bub00.bin

echo ""
echo "=== Comparing outputs ==="
diff /tmp/test_relation_bub00_MMA.txt /tmp/test_relation_bub00_Singular.txt && echo "OUTPUT IDENTICAL" || echo "OUTPUT DIFFERENT"

echo ""
echo "=== Relation counts ==="
echo "MMA:"
grep "sol_dim=" /tmp/test_relation_bub00_MMA.txt
echo ""
echo "Singular:"
grep "sol_dim=" /tmp/test_relation_bub00_Singular.txt
