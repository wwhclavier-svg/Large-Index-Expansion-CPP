cmake_minimum_required(VERSION 3.10)
project(IBP_FF_Test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找 Eigen3（必需）
find_package(Eigen3 REQUIRED)

# 查找 GMP（必需）
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)

# FireFly 路径（已构建但未安装）
set(FIREFLY_ROOT "/root/firefly-2.0.3")
set(FIREFLY_INCLUDE_DIRS 
    "${FIREFLY_ROOT}/source/include"
    "${FIREFLY_ROOT}/build/include"
)
set(FIREFLY_LIBRARY "${FIREFLY_ROOT}/build/libfirefly.a")

# 包含目录
include_directories(
    include          # 项目自身的头文件
    ${EIGEN3_INCLUDE_DIR}
    ${FIREFLY_INCLUDE_DIRS}
)

# 源文件
set(COMMON_SOURCES src/LayerRecursionCore.cpp)

# 生成可执行文件
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
add_executable(test_ff_verify test_ff_verify.cpp)
target_link_libraries(test_ff_verify ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} pthread)
