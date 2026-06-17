#!/bin/bash
# test_relation_solver_runner.sh
# 运行 RelationSolver 验证全套诊断测试
# 用法: bash test_relation_solver_runner.sh [test_name|all]

set -e

MMA="/usr/bin/wolframscript"
MMA_OPTS="-parallel Kernels -> 4"

CPP_BUILD="/root/Large-Index-Expansion-CPP/build"
VERIFY_DIR="/root/Large-Index-Expansion-CPP/verify"
MMA_DIR="/root/Large-Index-Expansion-MMA-Mini"
TESTS_DIR="$VERIFY_DIR/tests"

CYAN='\033[0;36m'
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN} RelationSolver 验证测试套件${NC}"
echo -e "${CYAN}========================================${NC}"

# 默认运行全部测试
TEST="${1:-all}"

# Pre-check: 确认 C++ 关系已生成
if [ ! -f "$CPP_BUILD/test_relationFF" ]; then
    echo -e "${RED}[ERROR] test_relationFF 未构建${NC}"
    echo "  运行: cd /root/Large-Index-Expansion-CPP && cmake --build build --target test_relationFF"
    exit 1
fi

# Pre-check: 确认 MMA 关系已生成
if [ ! -d "$VERIFY_DIR/bub00" ]; then
    echo -e "${YELLOW}[WARN] verify/bub00 不存在，先运行 MMA 工作流${NC}"
fi

run_mma_test() {
    local name="$1"
    local file="$2"
    echo ""
    echo -e "${CYAN}--- $name ---${NC}"
    if [ ! -f "$file" ]; then
        echo -e "${RED}[SKIP] 文件不存在: $file${NC}"
        return 1
    fi
    $MMA -code "AppendTo[\$Path, \"$MMA_DIR\"]; Get[\"$file\"]" 2>&1 | \
        grep -v -i 'warn\|shadow\|auf\|wiedersehen\|singular\|traceback\|//' | \
        tail -50
}

echo ""
echo -e "${CYAN}[0] Pre-check: C++ RelationSolver MatrixBuild${NC}"
echo -e "${CYAN}========================================${NC}"
cd "$CPP_BUILD"
./test_relationFF bub00 4 2 2 2>&1 | tail -20

if [ "$TEST" = "all" ] || [ "$TEST" = "diagnostic" ]; then
    run_mma_test "TEST 1: KiraVerify 根因诊断" \
        "$TESTS_DIR/test_kira_verify_diagnostic.m"
fi

if [ "$TEST" = "all" ] || [ "$TEST" = "nullspace" ]; then
    run_mma_test "TEST 2: 零空间自洽性" \
        "$TESTS_DIR/test_cpp_nullspace_consistency.m"
fi

if [ "$TEST" = "all" ] || [ "$TEST" = "gtoj" ]; then
    run_mma_test "TEST 3: g→j 转换诊断" \
        "$TESTS_DIR/test_g_to_j_conversion.m"
fi

if [ "$TEST" = "all" ] || [ "$TEST" = "expand" ]; then
    echo ""
    echo -e "${CYAN}--- TEST 4: C++ vs MMA 展开系数对比 ---${NC}"
    cd "$CPP_BUILD"
    ./test_expandFF bub00 2>&1 | tail -20
    
    echo ""
    echo -e "${CYAN}--- MMA SeriesVerify (IBP 方程代入) ---${NC}"
    cd "$MMA_DIR"
    $MMA -file VerifyExpand-SeriesVerify.wl bub00 2>&1 | \
        grep -v -i 'warn\|shadow\|auf\|wiedersehen\|singular\|traceback\|//' | \
        grep -E 'PASS|FAIL|Series|Residual|equation|verif|summary' | \
        tail -30
fi

echo ""
echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN} 测试完成${NC}"
echo -e "${CYAN}========================================${NC}"
