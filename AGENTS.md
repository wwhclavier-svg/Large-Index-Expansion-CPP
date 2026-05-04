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
| FireFly | 有限域运算 | 必需（或可用本地存根替代） |
| GMP | 多精度算术（FireFly依赖） | 必需 |
| nlohmann/json | JSON解析（内置于`include/json.hpp`） | 内置 |
| Mathematica | 数据生成脚本（`.wl`、`.m`文件） | 可选 |

---

## Project Structure

```
.
├── CMakeLists.txt           # 主构建配置（使用 /root/firefly-2.0.3 的FireFly）
├── CMakeLists_test.txt      # 备用构建配置（使用本地FireFly存根）
├── CMakeLists_backup.txt    # 历史备份构建配置
├── include/                 # 头文件（模板+接口）
│   ├── firefly/
│   │   └── FFInt.hpp              # 本地最小化FFInt存根实现
│   ├── LayerRecursion.hpp         # 层递归算法入口封装
│   ├── LayerRecursionCore.hpp     # 核心递归函数声明
│   ├── LayerRecursionCore.tpp     # 模板实现
│   ├── IBPMatrixLoader_Binary.hpp # 二进制IBP矩阵加载器
│   ├── IBPMatrixLoader.hpp        # JSON格式IBP矩阵加载器
│   ├── SeriesCoefficient.hpp      # 系数存储类
│   ├── SeriesCoefficientIO.hpp    # 系数二进制序列化
│   ├── LinearSolver.hpp           # 统一求解器调度器
│   ├── LinearSolver_Eigen.hpp     # 基于Eigen的浮点求解器
│   ├── LinearSolver_FF.hpp        # 基于FireFly的有限域求解器
│   ├── RelationSolver.hpp         # 线性关系重建求解器
│   ├── IncrementalRelationSolver.hpp # 增量式多层级求解器
│   ├── RingDataLoader.hpp         # 环/代数数据加载器
│   ├── Combinatorics.hpp          # 索引工具与种子生成
│   ├── binomial.hpp               # 预计算组合数表
│   ├── binomial2.hpp              # 组合数表变体
│   ├── UnifiedStorage.hpp         # 内存管理工具
│   ├── Utilities.hpp              # 通用工具函数
│   ├── json.hpp                   # 第三方JSON库（nlohmann/json）
│   ├── RelationRecon.hpp.txt      # 存档：旧版关系重建头文件
│   └── layerRecursion.hpp.txt     # 存档：旧版层递归头文件
├── src/                     # 源文件（非模板实现）
│   ├── LayerRecursionCore.cpp     # 核心非模板函数实现（极小，22行）
│   ├── layerRecursion.cpp         # 遗留驱动实现
│   ├── main.cpp                   # 遗留主程序（未活跃使用）
│   ├── #RelationRecon.cpp.txt     # 存档：旧版关系重建实现
│   └── #layerRecursion0.cpp.txt   # 存档：旧版层递归实现
├── tests/                   # 测试可执行文件源码（每个含main()）
│   ├── test_expandFF.cpp          # 有限域展开测试
│   ├── test_relationFF.cpp        # 有限域关系求解测试
│   ├── test_IBPVerification.cpp   # IBP矩阵验证测试
│   ├── test_expand_family.cpp     # 展开族测试
│   ├── test_load_bub.cpp          # bub格式矩阵加载测试
│   ├── IBPVerification.hpp        # IBP验证辅助函数
│   └── archive/                   # 存档的历史测试文件
│       ├── test_expand.cpp
│       ├── test_recons.cpp
│       ├── test_RelationNew.cpp
│       ├── IBPMatrixLoader_Test.cpp
│       └── RingDataLoader_test.cpp
├── build/                   # 构建输出目录（自动生成）
├── docs/                    # 文档目录
│   ├── LayerRecursion_Algorithm.md    # 层递归算法详情
│   ├── RelationSolver_ComponentGuide.md    # RelationSolver组件指南
│   ├── RelationSolver_QuickReference.md    # 快速参考
│   └── RelationSolver_Documentation_Hub.md # 文档中心
├── ReconstructReductionRelation.wl          # Mathematica关系重建包
└── ReconstructReductionRelation_Documentation.md # Mathematica包文档

# 根目录中的独立测试/工具文件（未加入CMake主配置或直接在根目录）
├── test_ff_verify.cpp       # FireFly库验证测试
├── test_check_matrix.cpp    # 矩阵检查工具
├── test_matrix_dump.cpp     # 矩阵转储工具
├── test_solver.cpp          # 求解器独立测试
├── test_firefly_simple.cpp  # FireFly简单测试
├── test_ffint.cpp           # FFInt基本测试
├── test_ff_div_debug.cpp    # 除法调试工具

# 数据文件（二进制和JSON格式，必须位于工作目录/项目根目录）
├── IBPMat_*.bin             # IBP矩阵二进制数据文件
├── RingData_*.bin           # 环数据二进制文件
├── resCache_*.bin           # 展开结果缓存文件
└── ExpansionCache_*.bin     # 展开缓存文件
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

# 4. 运行测试（必须从项目根目录运行，而非build/目录）
cd ..
./build/test_expandFF
./build/test_relationFF
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

### FireFly配置

- **主配置** (`CMakeLists.txt`)：FireFly已预构建于 `/root/firefly-2.0.3`，直接链接静态库 `libfirefly.a`
- **备用配置** (`CMakeLists_test.txt`)：使用 `include/firefly/FFInt.hpp` 中的本地最小化存根，不依赖外部FireFly库，仅用于基础编译测试

### 构建选项与CMake文件说明

| 文件 | 用途 |
|------|------|
| `CMakeLists.txt` | **主配置**：使用外部FireFly（`/root/firefly-2.0.3`）+ GMP + Eigen + pthread |
| `CMakeLists_test.txt` | **备用配置**：使用本地 `include/firefly/FFInt.hpp` 存根，仅依赖Eigen + GMP |
| `CMakeLists_backup.txt` | **历史备份**：使用 `find_package` 方式查找FireFly |

如需切换配置，可覆盖主文件：
```bash
cp CMakeLists_test.txt CMakeLists.txt
cd build && cmake .. && cmake --build .
```

---

## Code Organization & Architecture

### 类型系统

项目使用**模板多态性**支持多种数值类型：

| 类型 | 用例 | 头文件 |
|------|------|--------|
| `double` | 浮点测试、数值验证 | 原生 |
| `firefly::FFInt` | 有限域计算（素数模） | `<firefly/FFInt.hpp>`（外部库或本地存根） |

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
- `i`: 解索引（0到nimax），i=0为特解，i>0为齐次解

#### `LinearSystemResult<T>` (LinearSolver_Eigen.hpp / LinearSolver_FF.hpp)
统一的线性求解器结果结构：
```cpp
template<typename T>
struct LinearSystemResult {
    bool hasSolution;
    vector<vector<T>> Mext;  // 解矩阵（特解+零空间）
    vector<int> S;           // 自由变量索引
    vector<int> pivot_cols;  // 主元列索引
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

### 层递归算法四层结构

1. **入口**: `LayerRecursion.hpp` - `layerRecursion<T>()` 模板封装
2. **核心逻辑**: `LayerRecursionCore.hpp` - 函数声明
3. **非模板实现**: `LayerRecursionCore.cpp` - 仅 `equationVariable()` 等少量非模板函数（约22行）
4. **模板实现**: `LayerRecursionCore.tpp` - 模板实现（`inhomogTerms<T>` 等核心计算）

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

### 使用模式（来自test_relationFF.cpp）

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
| `test_relationFF` | 线性关系重建 | `FFInt` | `./build/test_relationFF` |
| `test_IBPVerification` | IBP矩阵验证 | `FFInt` | `./build/test_IBPVerification` |
| `test_expand_family` | 展开族测试 | `FFInt` | `./build/test_expand_family` |
| `test_load_bub` | bub格式矩阵加载 | `FFInt` | `./build/test_load_bub` |
| `test_ff_verify` | FireFly库验证 | `FFInt` | `./build/test_ff_verify` |

**注意**: 本项目不使用`ctest`或`add_test()`。测试是独立的可执行文件。

**重要**: 测试必须从**项目根目录**运行（即 `.bin` 数据文件所在的目录），而非 `build/` 目录内。因为测试程序通过相对路径加载二进制数据文件。

### 运行测试

```bash
# 从项目根目录执行
./build/test_expandFF
./build/test_relationFF
./build/test_IBPVerification
./build/test_load_bub
./build/test_expand_family
```

### 根目录中的独立测试/工具文件

以下文件位于项目根目录，未加入 `CMakeLists.txt` 主配置，需单独编译或作为参考：

| 文件 | 用途 |
|------|------|
| `test_ff_verify.cpp` | 验证FireFly库基本功能 |
| `test_check_matrix.cpp` | 检查矩阵一致性 |
| `test_matrix_dump.cpp` | 转储矩阵内容用于调试 |
| `test_solver.cpp` | 线性求解器独立测试 |
| `test_firefly_simple.cpp` | FireFly最简功能测试 |
| `test_ffint.cpp` | 本地FFInt存根测试 |
| `test_ff_div_debug.cpp` | 有限域除法调试 |

### 存档测试文件

以下历史测试文件已移至 `tests/archive/`，不再参与主构建：
- `test_expand.cpp`、`test_recons.cpp` — 早期双精度测试
- `test_RelationNew.cpp` — 早期多regime扩展测试
- `IBPMatrixLoader_Test.cpp`、`RingDataLoader_test.cpp` — 早期数据加载器测试

### 测试数据文件

测试需要工作目录（项目根目录）中的特定二进制数据文件：
- `IBPMat_*.bin`（如 `IBPMat_DBtop.bin`、`IBPMat_DPpart_QuadriScale.bin`）- IBP矩阵
- `RingData_*.bin`（如 `RingData_DBtop.bin`、`RingData_DPpart_QuadriScale.bin`）- 环数据
- `ExpansionCache_*.bin` / `resCache_*.bin` - 系数缓存（部分自动生成）

---

## Coding Conventions

### 命名风格

| 元素 | 约定 | 示例 |
|------|------|------|
| 文件 | PascalCase（头文件） | `LayerRecursion.hpp` |
| 类 | PascalCase | `class seriesCoefficient` |
| 函数 | camelCase | `layerRecursion()`、`getIndex()` |
| 变量 | snake_case（局部）, camelCase（成员） | `int num_regs;`、`int numRegs;` |
| 常量 | UPPER_CASE | `MAX_VAL`、`BINOM` |
| 模板 | PascalCase | `typename T`、`typename Field` |
| 命名空间 | PascalCase | `namespace LayerRecursionCore` |

### 注释语言

- **主要**: 中文（简体中文）- 用于数学/算法解释
- **次要**: 英文 - 用于API文档和简短注释

### 文件组织

- 模板：头文件（`.hpp`）或`.tpp`文件
- 实现：`.cpp`文件放在`src/`中
- 每个头文件应包含包含守卫：`#ifndef FILENAME_HPP`
- 以 `#` 开头并以 `.txt` 结尾的文件（如 `#RelationRecon.cpp.txt`）为**存档的遗留代码**，不参与当前构建

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
target_link_libraries(test_MyFeature ${FIREFLY_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} Eigen3::Eigen pthread)
```
3. 重新构建：`cd build && cmake --build .`

### 添加新头文件

1. 放置在 `include/MyHeader.hpp`
2. 使用包含守卫和命名空间
3. 从 `src/` 或其他头文件按需包含

### 修改核心算法

层递归算法有四层：
1. **入口**: `LayerRecursion.hpp` - `layerRecursion<T>()` 模板
2. **核心逻辑**: `LayerRecursionCore.hpp` - 函数声明
3. **非模板实现**: `LayerRecursionCore.cpp` - 非模板实现（极小）
4. **模板实现**: `LayerRecursionCore.tpp` - 模板实现（主要逻辑）

根据修改内容编辑相应的层次。

---

## Important Notes

1. **二进制数据格式**: IBP矩阵文件使用自定义二进制格式（魔数 `"IBP1"`）。详见 `IBPMatrixLoader_Binary.hpp`。

2. **全局状态**:
   - `BINOM[MAX_VAL][MAX_VAL]` - 预计算组合数（通过 `initBinomial()` 初始化）
   - `FFInt::p` - 有限域模数（通过 `FFInt::set_new_prime()` 设置）

3. **内存管理**: `seriesCoefficient` 预分配大块连续内存。注意大的 `order` 值可能导致的内存问题。

4. **构建目录**: 不要直接编辑 `build/` 中的文件。它们由CMake生成。

5. **Mathematica脚本**: `.wl` 文件和 `.m` 文件用于Mathematica中的数据生成与结果比对。C++构建不需要它们，但它们是验证流程的一部分。

6. **测试运行目录**: 必须从项目根目录运行测试可执行文件，因为测试程序通过**相对路径**加载 `IBPMat_*.bin` 和 `RingData_*.bin` 等数据文件。

7. **多CMake配置**: 项目存在多个 `CMakeLists*.txt` 文件。主文件为 `CMakeLists.txt`，使用外部FireFly库；`CMakeLists_test.txt` 为备用配置，使用本地存根。切换前请确认数据类型需求。

8. **FFInt 类型安全陷阱**: `FFInt` 仅有 `FFInt(uint64_t)` 构造函数，没有有符号整数构造函数。`static_cast<FFInt>(-1)` 会造成 `int(-1)` 隐式转换为 `uint64_t(2^64-1)` 再 mod p，产生垃圾值（如对于 p=179424673，结果为 4099945）。编译器不会对此产生警告。正确做法是使用 `-FFInt(1)`（调用 `operator-()`）而非从负整数构造。若需添加有符号构造函数，应使用 `int64_t` 以避免与现有 `uint64_t` 重载对 `long long` 参数产生歧义。此 bug 曾导致所有积分族的展开系数错误，而 EquationVerify（自洽性检查）仍通过（因验证器复用同一错误逻辑），详见 `verify/docs/IBPVerification.md`。

---

## References

- **FireFly Library**: https://github.com/firefly-library/firefly
- **Eigen Documentation**: https://eigen.tuxfamily.org/
- **IBP Method**: Integration-by-parts identities for Feynman integral reduction (physics)
- **ComprehensiveReport.tex**: [docs/ComprehensiveReport.tex](docs/ComprehensiveReport.tex) — LIE 方法的完整理论报告（原理参考）: asymptotic-solution completeness, block-recursive structure, geometric classification of solution spaces, finite-field implementation
