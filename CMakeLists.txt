cmake_minimum_required(VERSION 3.10)
project(IBP_FF_Test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找 Eigen3（必需）
find_package(Eigen3 REQUIRED)

# 查找 FireFly 和 GMP（必需）
find_path(FIREFLY_INCLUDE_DIR NAMES firefly/FFInt.hpp REQUIRED)
find_library(FIREFLY_LIBRARY NAMES firefly REQUIRED)
find_library(GMP_LIBRARY NAMES gmp gmpxx REQUIRED)

# 包含目录
include_directories(
    include          # 项目自身的头文件
    ${EIGEN3_INCLUDE_DIR}
    ${FIREFLY_INCLUDE_DIR}
)

# 源文件：只包含非模板的实现文件
set(COMMON_SOURCES src/LayerRecursionCore.cpp)
set(TEST_SOURCE tests/test_expandFF.cpp)

# 生成可执行文件
add_executable(test_expandFF ${TEST_SOURCE} ${COMMON_SOURCES})
target_link_libraries(test_expandFF ${FIREFLY_LIBRARY} ${GMP_LIBRARY} Eigen3::Eigen)

add_executable(test_RelationFF tests/test_RelationFF.cpp ${COMMON_SOURCES})
target_link_libraries(test_RelationFF ${FIREFLY_LIBRARY} ${GMP_LIBRARY} Eigen3::Eigen)

# 以上测试：有限域的读取，展开系数的计算.