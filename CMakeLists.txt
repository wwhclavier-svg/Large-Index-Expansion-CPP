cmake_minimum_required(VERSION 3.10)
project(IBP_FF_Test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 默认启用 Release 优化，除非显式指定
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ============================================================
# 依赖项查找（可移植设计）
# ============================================================

# --- Eigen3 ---
find_package(Eigen3 REQUIRED)

# --- GMP ---
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)
if(NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP not found. Install libgmp-dev (Debian/Ubuntu) or gmp-devel (RHEL/CentOS).")
endif()

# --- FireFly ---
# 优先级: 1) -DFIREFLY_ROOT=...  2) 环境变量 FIREFLY_ROOT  3) 默认路径
if(NOT DEFINED FIREFLY_ROOT)
    if(DEFINED ENV{FIREFLY_ROOT})
        set(FIREFLY_ROOT "$ENV{FIREFLY_ROOT}")
    else()
        set(FIREFLY_ROOT "/home/ykm/firefly")
    endif()
endif()
set(FIREFLY_ROOT "${FIREFLY_ROOT}" CACHE PATH "FireFly library root directory")

# Try source-tree layout first, then system-installed layout
set(_FF_LIB "${FIREFLY_ROOT}/build/libfirefly.a")
if(NOT EXISTS "${_FF_LIB}")
    find_library(_FF_LIB_FOUND NAMES firefly)
    if(_FF_LIB_FOUND)
        set(_FF_LIB "${_FF_LIB_FOUND}")
    endif()
endif()
set(FIREFLY_LIBRARY "${_FF_LIB}")
unset(_FF_LIB)

if(NOT DEFINED FIREFLY_INCLUDE_DIRS)
    # Check for system-installed headers
    find_path(_FF_INC_FOUND NAMES firefly/FFInt.hpp)
    if(_FF_INC_FOUND)
        set(FIREFLY_INCLUDE_DIRS "${_FF_INC_FOUND}")
        message(STATUS "Using system-installed FireFly: ${_FF_LIB_FOUND}")
    else()
        set(FIREFLY_INCLUDE_DIRS
            "${FIREFLY_ROOT}/source/include"
            "${FIREFLY_ROOT}/build/include"
        )
    endif()
endif()

if(NOT EXISTS "${FIREFLY_LIBRARY}")
    message(FATAL_ERROR
        "FireFly library not found at: ${FIREFLY_LIBRARY}\n"
        "Set FIREFLY_ROOT to the FireFly source directory:\n"
        "  cmake -DFIREFLY_ROOT=/path/to/firefly ..\n"
        "  or: export FIREFLY_ROOT=/path/to/firefly && cmake ..")
endif()

# ============================================================
# 构建配置
# ============================================================

# 包含目录
include_directories(
    include
    ${EIGEN3_INCLUDE_DIR}
    ${FIREFLY_INCLUDE_DIRS}
)

# 源文件
set(COMMON_SOURCES src/LayerRecursionCore.cpp)

# --- 必选测试 ---
add_executable(test_expandFF tests/test_expandFF.cpp ${COMMON_SOURCES})
target_link_libraries(test_expandFF ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)

add_executable(test_IBPVerification tests/test_IBPVerification.cpp ${COMMON_SOURCES})
target_link_libraries(test_IBPVerification ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)

add_executable(test_relationFF tests/test_relationFF.cpp ${COMMON_SOURCES})
target_link_libraries(test_relationFF ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)

add_executable(test_load_bub tests/test_load_bub.cpp ${COMMON_SOURCES})
target_link_libraries(test_load_bub ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)

add_executable(test_expand_family tests/test_expand_family.cpp ${COMMON_SOURCES})
target_link_libraries(test_expand_family ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)

# FireFly 验证测试
add_executable(test_ff_verify tools/test_ff_verify.cpp)
target_link_libraries(test_ff_verify ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} pthread)

# 二进制文件读写往返验证测试
add_executable(test_family_generate tests/test_family_generate.cpp ${COMMON_SOURCES})
target_link_libraries(test_family_generate ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)

# JSON 家族配置解析器测试
add_executable(test_family_config tests/test_family_config.cpp)
target_link_libraries(test_family_config pthread)

# FamilyGenerate CLI 工具 (Phase 1: stub intermediate data)
add_executable(family_generate tools/family_generate.cpp)
target_link_libraries(family_generate ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} pthread)

# Singular CAS 集成测试
add_executable(test_singular_runner tools/test_singular_runner.cpp)
target_link_libraries(test_singular_runner pthread)

add_executable(test_region_solver tools/test_region_solver.cpp)
target_link_libraries(test_region_solver pthread)

add_executable(test_ibp_analyzer tests/test_ibp_analyzer.cpp)
target_link_libraries(test_ibp_analyzer pthread)

add_executable(test_recursion_builder tests/test_recursion_builder.cpp)
target_link_libraries(test_recursion_builder pthread)

add_executable(test_recursion_builder_tri tests/test_recursion_builder_tri.cpp)
target_link_libraries(test_recursion_builder_tri pthread)

add_executable(test_ring_builder tests/test_ring_builder.cpp)
target_link_libraries(test_ring_builder pthread)

add_executable(verify_ring_builder tests/verify_ring_builder.cpp)
target_link_libraries(verify_ring_builder pthread)

add_executable(test_region_solver_full tests/test_region_solver_full.cpp)
target_link_libraries(test_region_solver_full pthread)

add_executable(debug_region_solver tests/debug_region_solver.cpp)
target_link_libraries(debug_region_solver pthread)

add_executable(compare_ringdata_detailed tests/compare_ringdata_detailed.cpp)
target_link_libraries(compare_ringdata_detailed pthread)

add_executable(diagnose_ibp tools/diagnose_ibp.cpp)
target_link_libraries(diagnose_ibp pthread)

# add_executable(diagnostic_dump tools/diagnostic_dump.cpp)
# target_link_libraries(diagnostic_dump pthread)

# --- 可选测试（仅在源文件存在时构建）---
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/test_compare_MMA.cpp")
    add_executable(test_compare_MMA tests/test_compare_MMA.cpp ${COMMON_SOURCES})
    target_link_libraries(test_compare_MMA ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/tests/verify_MMA_lev1_deg1.cpp")
    add_executable(verify_MMA_lev1_deg1 tests/verify_MMA_lev1_deg1.cpp ${COMMON_SOURCES})
    target_link_libraries(verify_MMA_lev1_deg1 ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)
endif()

# T005: GeneratingCone 约化规则生成
add_executable(generating_cone tools/generating_cone.cpp)
target_link_libraries(generating_cone pthread)

add_executable(test_generating_cone tests/test_generating_cone.cpp)
target_link_libraries(test_generating_cone pthread)
