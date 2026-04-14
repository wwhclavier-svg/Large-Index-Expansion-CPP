# AGENTS.md - Project Guide for AI Coding Agents

## Project Overview

本项目是一个**C++17科学计算项目**，专注于**IBP（Integration By Parts）矩阵的级数展开**和**有限域运算**。该代码库实现了用于计算IBP矩阵展开系数和重建线性约化关系的算法。

### 核心功能

- **层递归（Layer Recursion）**：逐阶计算IBP矩阵的级数展开系数
- **有限域运算（Finite Field Arithmetic）**：使用FireFly库在有限域上进行计算
- **线性关系重建（Linear Relation Reconstruction）**：求解展开系数之间的线性关系

项目使用重度模板元编程来同时支持 `double`（浮点数）和 `firefly::FFInt`（有限域整数）两种数据类型。

---

## Technology Stack

| 组件 | 用途 | 必需性 |
|------|------|--------|
| C++17 | 核心语言 | 必需 |
| CMake 3.10+ | 构建系统 | 必需 |
| Eigen3 | 线性代数（浮点运算） | 必需 |
| FireFly | 有限域运算 | 必需 |
| GMP | 多精度算术（FireFly依赖） | 必需 |
| nlohmann/json | JSON解析（内置于`include/json.hpp`） | 内置 |
| Mathematica | 数据生成脚本（`.wl`文件） | 可选 |

---

## Project Structure

```
.
├── CMakeLists.txt           # 主构建配置
├── CMakeLists_full.txt      # 扩展构建配置（参考/备份）
├── include/                 # 头文件（模板+接口）
│   ├── LayerRecursion.hpp           # 层递归算法入口封装
│   ├── LayerRecursionCore.hpp       # 核心递归函数声明
│   ├── LayerRecursionCore.tpp       # 模板实现
│   ├── IBPMatrixLoader_Binary.hpp   # 二进制IBP矩阵加载器
│   ├── IBPMatrixLoader.hpp          # JSON格式IBP矩阵加载器
│   ├── SeriesCoefficient.hpp        # 系数存储类
│   ├── SeriesCoefficientIO.hpp      # 系数二进制序列化
│   ├── LinearSolver.hpp             # 统一求解器调度器
│   ├── LinearSolver_Eigen.hpp       # 基于Eigen的浮点求解器
│   ├── LinearSolver_FF.hpp          # 基于FireFly的有限域求解器
│   ├── RelationSolver.hpp           # 线性关系重建求解器
│   ├── IncrementalRelationSolver.hpp # 增量式多层级求解器
│   ├── RingDataLoader.hpp           # 环/代数数据加载器
│   ├── Combinatorics.hpp            # 索引工具与种子生成
│   ├── binomial.hpp                 # 预计算组合数表
│   ├── UnifiedStorage.hpp           # 内存管理工具
│   ├── Utilities.hpp                # 通用工具函数
│   └── json.hpp                     # 第三方JSON库
├── src/                     # 源文件（非模板实现）
│   ├── LayerRecursionCore.cpp   # 核心函数实现
│   ├── layerRecursion.cpp       # 遗留实现
│   └── main.cpp                 # 遗留驱动程序（未活跃使用）
├── tests/                   # 测试可执行文件（每个含main()）
│   ├── test_expandFF.cpp        # 有限域展开测试
│   ├── test_RelationFF.cpp      # 有限域关系求解测试
│   ├── test_RelationNew.cpp     # 扩展关系求解测试
│   ├── test_expand.cpp          # 双精度展开测试
│   ├── test_recons.cpp          # 重建算法测试
│   ├── IBPMatrixLoader_Test.cpp # 矩阵加载器测试
│   └── RingDataLoader_test.cpp  # 环数据加载器测试
├── build/                   # 构建输出目录（自动生成）
├── docs/                    # 文档目录
│   ├── RelationSolver_ComponentGuide.md    # RelationSolver组件指南
│   ├── RelationSolver_QuickReference.md    # 快速参考
│   └── RelationSolver_Documentation_Hub.md # 文档中心
├── ReconstructReductionRelation.wl          # Mathematica关系重建包
└── ReconstructReductionRelation_Documentation.md # Mathematica包文档

# 数据文件（二进制和JSON格式，工作目录中）
├── IBPMat_*.bin             # IBP矩阵二进制数据文件
├── RingData_*.bin           # 环数据二进制文件
├── resCache_*.bin           # 展开结果缓存文件
└── *.json                   # JSON格式的矩阵数据
```

---

## Build Instructions

### 标准CMake构建

```bash
# 1. 创建并进入构建目录
mkdir -p build && cd build

# 2. 配置（需要Eigen3、FireFly、GMP已安装）
cmake ..

# 3. 构建
cmake --build .

# 4. 运行测试
./test_expandFF
./test_RelationFF
```

### 依赖安装

**Ubuntu/Debian:**
```bash
sudo apt-get install libeigen3-dev libgmp-dev
# FireFly需从源码构建：https://github.com/firefly-library/firefly
```

**macOS:**
```bash
brew install eigen gmp
# FireFly需从源码构建
```

### 构建选项

`CMakeLists.txt` 提供：
- `BUILD_FF_TESTS`（默认开启）：构建需要FireFly的有限域测试

---

## Code Organization & Architecture

### 类型系统

项目使用**模板多态性**支持多种数值类型：

| 类型 | 用例 | 头文件 |
|------|------|--------|
| `double` | 浮点测试、数值验证 | 原生 |
| `firefly::FFInt` | 有限域计算（素数模） | `<firefly/FFInt.hpp>` |

类型分发使用C++17的 `if constexpr` 实现：
```cpp
if constexpr (std::is_same_v<T, firefly::FFInt>) {
    // 有限域分支
} else {
    // 浮点分支
}
```

### 核心数据结构

#### `IBPMatrixE<T>` (IBPMatrixLoader_Binary.hpp)
IBP矩阵算子存储结构：
```cpp
template<typename T>
struct IBPMatrixE {
    vector<vector<vector<T>>> N1, K1, M1, K1s, K2s;  // 3D算子
    vector<vector<T>> F0;                              // 2D算子
    vector<vector<vector<vector<T>>>> F2, F2s;        // 4D算子
    int nibp, ne, nb;  // 维度
};
```

#### `seriesCoefficient<T>` (SeriesCoefficient.hpp)
5维系数存储，索引为 `(k, l, cid, j, i)`：
- `k`: 展开阶数
- `l`: 层级别
- `cid`: 种子/组合索引
- `j`: 基索引（0到nb-1）
- `i`: 解索引（0到nimax）

#### `LinearSystemResult<T>` (LinearSolver_*.hpp)
统一的线性求解器结果结构：
```cpp
template<typename T>
struct LinearSystemResult {
    bool hasSolution;
    vector<vector<T>> Mext;  // 解矩阵（特解+零空间）
    vector<int> S;           // 自由变量索引
};
```

### 模块架构

```
LayerRecursion.hpp                    RelationSolver.hpp
       ↓                                      ↓
LayerRecursionCore.hpp              RegimeData<T>
       ↓                               RegimeEvaluator<T>
LayerRecursionCore.tpp              AdaptiveEquationBuilder<T>
       ↓                                      ↓
LinearSolver.hpp ←────────────────→ LinearSolver_FF.hpp
       ↓                               LinearSolver_Eigen.hpp
       ↓                                      ↓
SeriesCoefficient.hpp              UnifiedStorage.hpp
IBPMatrixLoader_Binary.hpp         Combinatorics.hpp
RingDataLoader.hpp
```

---

## RelationSolver 模块规范

### 位置与用途
- **文件**: `include/RelationSolver.hpp`
- **命名空间**: `RelationSolver`
- **作用**: 重建IBP矩阵展开系数之间的线性约化关系
- **支持类型**: 模板参数 `T`（`double`、`firefly::FFInt`）

### 核心数据结构

#### `RegimeData<T>`
封装单个sector及其系数和A算子。
```cpp
template<typename T>
struct RegimeData {
    const seriesCoefficient<T>* C;           // 系数引用
    std::vector<int> theta;                  // Sector标识符
    std::vector<std::vector<T>> A_ops;       // A_i矩阵组
    std::vector<std::vector<T>> A_inv_ops;   // A_i^{-1}矩阵组
    int nb;                                  // 基维度
};
```

#### `RelationCoefficient<T>`
结果结构，提供关系系数的索引访问。
```cpp
template<typename T>
class RelationCoefficient {
    // 访问给定多重指标的系数
    const std::vector<T>& operator()(
        const std::vector<int>& alpha,  // 多重指标约束
        const std::vector<int>& beta)   // 多重指标幂次
    const;
};
```

#### `AdaptiveSamplingConfig`
自适应采样和收敛检测配置。
```cpp
struct AdaptiveSamplingConfig {
    int min_nu = 3;                    // 最小采样点数
    int max_nu = 50;                   // 最大采样点数
    int convergence_threshold = 3;     // 零化度稳定性检测阈值
    int lev_hint = 0;                  // 当前(lev, deg)级别提示
    int deg_hint = 0;
};
```

### 主API：`reconstructReductionRelation<T>()`

**用途**: 求解特定(lev, deg)级别的所有线性关系。

```cpp
template<typename T>
std::pair<LinearSystemResult<T>, RelationCoefficient<T>> 
reconstructReductionRelation(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne,                            // 多重指标维度
    int lev, int deg,                  // 当前级别/次数
    const AdaptiveSamplingConfig& config = {});
```

### 使用模式（来自test_RelationFF.cpp）

```cpp
// 对每个(lev, deg)对：
auto [linear_result, relation_coeff] = 
    RelationSolver::reconstructReductionRelation<FFInt>(
        allResults,        // 层递归的系数
        sector_list,       // Sector标识符
        A_list, Ainv_list, // A算子
        ne, lev, deg,      // 维度
        config);           // 自适应采样配置

// 检查是否找到关系
if (linear_result.hasSolution) {
    int num_relations = linear_result.S.size();
    
    // 访问特定关系系数
    for (const auto& alpha : alphas) {
        for (const auto& beta : betas) {
            const auto& coeff_vec = relation_coeff(alpha, beta);
            // coeff_vec[0] = 特解
            // coeff_vec[1..] = 解空间基
        }
    }
}
```

### 算法：自适应多点采样

1. **Regime准备** — 计算给定lev下所有多重指标alpha的P(alpha)矩阵
2. **迭代采样** — 生成随机nu向量；每个产生|regimes| × nb × (k_max+1)个方程
3. **系统求解** — 求解齐次线性系统；追踪零化度维度
4. **收敛检测** — 当零化度在convergence_threshold次迭代中稳定时停止
5. **解提取** — 将结果打包为按(alpha, beta)索引的RelationCoefficient

---

## Algorithm Flow

```
1. 数据加载：loadAllIBPMatricesBinary<T>() 从二进制文件加载IBP矩阵
       ↓
2. 层递归：layerRecursion<T>() 逐阶计算展开系数
       ↓
3. 序列化：SeriesIO::saveAllResults() 缓存计算的系数
       ↓
4. 关系求解：
   - 用(lev, deg)提示配置自适应采样
   - 对每个(lev, deg)调用reconstructReductionRelation<T>()
   - 以RelationCoefficient<T>形式返回解供下游处理
```

---

## Testing Strategy

### 测试可执行文件

| 测试 | 用途 | 数据类型 | 命令 |
|------|------|----------|------|
| `test_expandFF` | 展开系数计算 | `FFInt` | `./build/test_expandFF` |
| `test_RelationFF` | 线性关系重建 | `FFInt` | `./build/test_RelationFF` |
| `test_expand` | 双精度展开 | `double` | `./build/test_expand` |
| `test_recons` | 重建算法 | `double` | `./build/test_recons` |
| `test_RelationNew` | 扩展多regime测试 | `FFInt` | `./build/test_RelationNew` |

**注意**: 本项目不使用`ctest`或`add_test()`。测试是独立的可执行文件。

### 运行测试

```bash
cd build

# 有限域展开测试（需要IBPMat_DPpart_QuadriScale.bin）
./test_expandFF

# 关系求解测试（需要RingData_*.bin, IBPMat_*.bin）
./test_RelationFF
```

### 测试数据文件

测试需要工作目录中的特定二进制数据文件：
- `IBPMat_DBtop.bin`, `IBPMat_DPpart_QuadriScale.bin` - IBP矩阵
- `RingData_DBtop.bin`, `RingData_DPpart_QuadriScale.bin` - 环数据
- `ExpansionResults_cache.bin` - 系数缓存（自动生成）

---

## Coding Conventions

### 命名风格

| 元素 | 约定 | 示例 |
|------|------|------|
| 文件 | PascalCase（头文件） | `LayerRecursion.hpp` |
| 类 | PascalCase | `class seriesCoefficient` |
| 函数 | camelCase | `layerRecursion()`, `getIndex()` |
| 变量 | snake_case（局部）, camelCase（成员） | `int num_regs;`, `int numRegs;` |
| 常量 | UPPER_CASE | `MAX_VAL`, `BINOM` |
| 模板 | PascalCase | `typename T`, `typename Field` |
| 命名空间 | PascalCase | `namespace LayerRecursionCore` |

### 注释语言

- **主要**: 中文（简体中文）- 用于数学/算法解释
- **次要**: 英文 - 用于API文档和简短注释

### 文件组织

- 模板：头文件（`.hpp`）或`.tpp`文件
- 实现：`.cpp`文件放在`src/`中
- 每个头文件应包含包含守卫：`#ifndef FILENAME_HPP`

### 代码模式

1. **模板特化**用于类型分发：
```cpp
if constexpr (std::is_same_v<T, firefly::FFInt>) {
    // 有限域路径
} else {
    // 浮点路径
}
```

2. **命名空间组织**：
```cpp
namespace LayerRecursionCore { /* 核心函数 */ }
namespace AlgebraData { /* 数据加载 */ }
namespace RelationSolver { /* 关系求解 */ }
namespace SeriesIO { /* 序列化 */ }
```

---

## Common Tasks

### 添加新测试

1. 创建 `tests/test_MyFeature.cpp` 并包含 `main()` 函数
2. 添加到 `CMakeLists.txt`：
```cmake
add_executable(test_MyFeature tests/test_MyFeature.cpp ${COMMON_SOURCES})
target_link_libraries(test_MyFeature ${FIREFLY_LIBRARY} ${GMP_LIBRARY} Eigen3::Eigen)
```
3. 重新构建：`cd build && cmake --build .`

### 添加新头文件

1. 放置在 `include/MyHeader.hpp`
2. 使用包含守卫和命名空间
3. 从 `src/` 或其他头文件按需包含

### 修改核心算法

层递归算法有三层：
1. **入口**: `LayerRecursion.hpp` - `layerRecursion<T>()` 模板
2. **核心逻辑**: `LayerRecursionCore.hpp` - 函数声明
3. **实现**: `LayerRecursionCore.cpp` - 非模板实现
4. **模板实现**: `LayerRecursionCore.tpp` - 模板实现

根据修改内容编辑相应的层次。

---

## Important Notes

1. **二进制数据格式**: IBP矩阵文件使用自定义二进制格式（魔数 `"IBP1"`）。详见 `IBPMatrixLoader_Binary.hpp`。

2. **全局状态**:
   - `BINOM[MAX_VAL][MAX_VAL]` - 预计算组合数（通过 `initBinomial()` 初始化）
   - `FFInt::p` - 有限域模数（通过 `FFInt::set_new_prime()` 设置）

3. **内存管理**: `seriesCoefficient` 预分配大块连续内存。注意大的 `order` 值可能导致的内存问题。

4. **构建目录**: 不要直接编辑 `build/` 中的文件。它们由CMake生成。

5. **Mathematica脚本**: `.wl` 文件用于Mathematica/MATLAB中的数据生成。C++构建不需要它们。

---

## References

- **FireFly Library**: https://github.com/firefly-library/firefly
- **Eigen Documentation**: https://eigen.tuxfamily.org/
- **IBP Method**: Integration-by-parts identities for Feynman integral reduction (physics)
