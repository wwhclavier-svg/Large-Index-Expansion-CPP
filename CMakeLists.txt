cmake_minimum_required(VERSION 3.10)
project(IBP_FF_Test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
        set(FIREFLY_ROOT "/root/firefly-2.0.3")
    endif()
endif()
set(FIREFLY_ROOT "${FIREFLY_ROOT}" CACHE PATH "FireFly library root directory")

set(FIREFLY_INCLUDE_DIRS
    "${FIREFLY_ROOT}/source/include"
    "${FIREFLY_ROOT}/build/include"
)
set(FIREFLY_LIBRARY "${FIREFLY_ROOT}/build/libfirefly.a")

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

# --- 可选测试（仅在源文件存在时构建）---
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/test_compare_MMA.cpp")
    add_executable(test_compare_MMA tests/test_compare_MMA.cpp ${COMMON_SOURCES})
    target_link_libraries(test_compare_MMA ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/tests/verify_MMA_lev1_deg1.cpp")
    add_executable(verify_MMA_lev1_deg1 tests/verify_MMA_lev1_deg1.cpp ${COMMON_SOURCES})
    target_link_libraries(verify_MMA_lev1_deg1 ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)
endif()
