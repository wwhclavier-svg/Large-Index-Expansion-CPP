#!/bin/bash
set -e

echo "=== Comparing bub00 relation coefficients ==="

# Backup
mv data/IBPMat_bub00.bin data/IBPMat_bub00-MMA.bin
mv data/RingData_bub00.bin data/RingData_bub00-MMA.bin

# MMA version
cp data/IBPMat_bub00-MMA.bin data/IBPMat_bub00.bin
cp data/RingData_bub00-MMA.bin data/RingData_bub00.bin
./build/test_relationFF bub00 4 1 2 2 > /dev/null 2>&1
cp relations/AllRelations_bub00_k4.m relations/AllRelations_bub00_k4_MMA.m

# Singular version
cp verify/DP323_recompute/IBPMat_bub00-Singular.bin data/IBPMat_bub00.bin
cp verify/DP323_recompute/RingData_bub00-Singular.bin data/RingData_bub00.bin
./build/test_relationFF bub00 4 1 2 2 > /dev/null 2>&1
cp relations/AllRelations_bub00_k4.m relations/AllRelations_bub00_k4_Singular.m

# Restore
cp data/IBPMat_bub00-MMA.bin data/IBPMat_bub00.bin
cp data/RingData_bub00-MMA.bin data/RingData_bub00.bin

echo "MMA relations file: relations/AllRelations_bub00_k4_MMA.m"
echo "Singular relations file: relations/AllRelations_bub00_k4_Singular.m"
echo ""

# Compare line counts
lines_MMA=$(wc -l < relations/AllRelations_bub00_k4_MMA.m)
lines_Sing=$(wc -l < relations/AllRelations_bub00_k4_Singular.m)
echo "Lines: MMA=$lines_MMA Sing=$lines_Sing"

# Compare sizes
size_MMA=$(stat -c%s relations/AllRelations_bub00_k4_MMA.m)
size_Sing=$(stat -c%s relations/AllRelations_bub00_k4_Singular.m)
echo "Bytes: MMA=$size_MMA Sing=$size_Sing"

# Full diff
diff relations/AllRelations_bub00_k4_MMA.m relations/AllRelations_bub00_k4_Singular.m > /tmp/relations_diff.txt && echo "FILES IDENTICAL" || echo "FILES DIFFERENT (see /tmp/relations_diff.txt)"
