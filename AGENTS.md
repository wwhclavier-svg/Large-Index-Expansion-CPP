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
├── CMakeLists.txt           # 主构建配置（使用 /home/ykm/firefly 的FireFly）
├── CMakeLists_test.txt      # 备用构建配置（使用本地FireFly存根）
├── include/                 # 头文件（模板+接口）
│   ├── firefly/
│   │   └── FFInt.hpp              # 本地最小化FFInt存根实现
│   ├── LayerRecursion.hpp         # 层递归算法入口封装
│   ├── LayerRecursionCore.hpp     # 核心递归函数声明
│   ├── LayerRecursionCore.tpp     # 模板实现
│   ├── IBPMatrixLoader_Binary.hpp # 二进制IBP矩阵加载器
│   ├── IBPMatrixLoader.hpp        # JSON格式IBP矩阵加载器
│   ├── BinaryIBPWriter.hpp        # 二进制IBP矩阵写入器
│   ├── BinaryRingWriter.hpp       # 二进制环数据写入器
│   ├── RingDataLoader.hpp         # 环/代数数据加载器
│   ├── SeriesCoefficient.hpp      # 系数存储类
│   ├── SeriesCoefficientIO.hpp    # 系数二进制序列化
│   ├── LinearSolver.hpp           # 统一求解器调度器 + LinearSystemResult
│   ├── LinearSolver_Eigen.hpp     # 基于Eigen的浮点求解器
│   ├── LinearSolver_FF.hpp        # 基于FireFly的有限域求解器
│   ├── RelationSolver.hpp         # 线性关系重建（RelationSolver组件）
│   ├── IncrementalRelationSolver.hpp # 增量式RREF维护多层级求解器
│   ├── IBPEqGenerator.hpp         # IBP方程生成器（调用Singular）
│   ├── SingularRunner.hpp         # Singular子进程管理
│   ├── RegionSolver.hpp           # 渐近区域求解器
│   ├── FamilyConfig.hpp           # 积分族配置
│   ├── IBPAnalyzer.hpp            # IBP矩阵分析
│   ├── PolyArith.hpp              # 多项式算术
│   ├── RecursionBuilder.hpp       # 递归构建器
│   ├── RingBuilder.hpp            # 环数据构建器
│   ├── Combinatorics.hpp          # 索引工具与种子生成
│   ├── binomial.hpp               # 预计算组合数表
│   ├── binomial2.hpp              # 组合数表变体
│   ├── UnifiedStorage.hpp         # 内存管理工具
│   ├── Utilities.hpp              # 通用工具函数
│   └── json.hpp                   # 第三方JSON库（nlohmann/json）
├── src/                     # 源文件（非模板实现）
│   ├── LayerRecursionCore.cpp     # 核心非模板函数实现
│   └── ... (archived .txt files in archive/)
├── tests/                   # 测试可执行文件源码（每个含main()）
│   ├── test_expandFF.cpp          # 有限域展开测试
│   ├── test_relationFF.cpp        # 有限域关系求解测试
│   ├── test_IBPVerification.cpp   # IBP矩阵验证测试
│   ├── test_compare_MMA.cpp       # MMA对比测试
│   └── ...                        # 更多测试文件
├── tools/                   # 独立工具和脚本
│   ├── family_generate.cpp        # 积分族生成工具
│   ├── diagnose_ibp.cpp           # IBP诊断工具
│   └── ...
├── data/                    # 工作二进制数据文件（IBPMat_*.bin, RingData_*.bin, ExpansionCache_*.bin）
├── relations/               # 导出的关系结果（AllRelations_*.m, RelationMeta_*.m, Relations_*.m）
├── scripts/                 # Shell脚本（compare_relation_*.sh等）
├── checkpoint/              # 会话存档和检查点文件
├── mma/                     # Mathematica管线脚本
│   ├── AllRelationsToNuFormat.wl
│   ├── ReconstructReductionRelation.wl
│   └── reference/           # MMA参考输出（ExpansionMMA_*.m, VerifyAllRelations.wl等）
├── families/                # JSON积分族定义
├── include/                 # 头文件（同上）
├── src/                     # 源文件
├── tests/                   # 测试
├── tools/                   # 工具
├── archive/                 # 存档的历史文件
├── build/                   # 构建输出目录（自动生成）
├── verify/                  # 跨验证框架（C++ vs MMA vs Kira）
│   ├── README.md
│   ├── FamilyDatabase/
│   └── VerifyUtility/
├── docs/                    # 文档目录
└── cache/                   # 调度器缓存（运行时生成）
```
**数据文件路径**（所有工作二进制数据在 `data/` 中，测试通过 `"data/IBPMat_*.bin"` 相对路径加载）：
- `data/IBPMat_*.bin` — IBP矩阵二进制数据文件
- `data/RingData_*.bin` — 环数据二进制文件
- `data/ExpansionCache_*.bin` / `data/resCache_*.bin` — 展开结果缓存文件
- `relations/AllRelations_*.m` — 导出的关系结果（MMA格式）
- `verify/<family>/` — MMA参考二进制和元数据文件

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

- **主配置** (`CMakeLists.txt`)：FireFly已预构建于 `/home/ykm/firefly`，直接链接静态库 `libfirefly.a`
- **备用配置** (`CMakeLists_test.txt`)：使用 `include/firefly/FFInt.hpp` 中的本地最小化存根，不依赖外部FireFly库，仅用于基础编译测试

### 构建选项与CMake文件说明

| 文件 | 用途 |
|------|------|
| `CMakeLists.txt` | **主配置**：使用外部FireFly（`/home/ykm/firefly`）+ GMP + Eigen + pthread |
| `CMakeLists_test.txt` | **备用配置**：使用本地 `include/firefly/FFInt.hpp` 存根，仅依赖Eigen + GMP |
| `archive/CMakeLists_backup.txt` | **历史备份**：使用 `find_package` 方式查找FireFly |

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

#### `LinearSystemResult<T>` (LinearSolver.hpp)
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
       ↓                              (RegimeData<T>, RegimeEvaluator<T>,
LayerRecursionCore.hpp                 AdaptiveEquationBuilder<T>,
       ↓                               RelationCoefficient<T>)
LayerRecursionCore.tpp                        ↓
       ↓                              IncrementalRelationSolver.hpp
LinearSolver.hpp                              ↓
  (dispatches to ↓)                  UnifiedStorage.hpp
LinearSolver_FF.hpp                  Combinatorics.hpp
LinearSolver_Eigen.hpp               Utilities.hpp
       ↓                              binomial.hpp / binomial2.hpp
SeriesCoefficient.hpp
SeriesCoefficientIO.hpp             FamilyGenerate pipeline:
IBPMatrixLoader_Binary.hpp            IBPEqGenerator.hpp
IBPMatrixLoader.hpp (JSON)            SingularRunner.hpp
RingDataLoader.hpp                    RegionSolver.hpp
                                      FamilyConfig.hpp
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
    int min_nu = 0;                    // 最小采样数（0=自动调整）
    int max_nu = 200;                  // 最大采样数（硬上限）
    double safety_factor = 1.2;        // 过定因子
    int check_interval = 1;            // 秩检查间隔
    int nullity_stable_threshold = 3;  // 零空间维度连续稳定次数
    int verification_points = 3;       // 验证用额外点数
    int plateau_size = 1;              // 稳定性确认所需额外阶数
    double tolerance = 1e-10;          // 数值容差
    bool use_special_points = true;    // 是否使用特殊点
    double random_min = 3.0;           // 随机点最小值
    double random_max = 100.0;         // 随机点最大值
    int lev_hint = 2;                  // |alpha| 上界提示
    int deg_hint = 2;                  // |beta| 上界提示
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

### 高级API：`reconstructAllRelations<T>()`（推荐）

**用途**: 多配置求解，内置 RemoveSolvedVariables 策略以减少方程冗余。

```cpp
template<typename T>
std::vector<LevDegResult<T>> reconstructAllRelations(
    const std::vector<std::vector<seriesCoefficient<T>>>& CTable,
    const std::vector<std::vector<int>>& sector,
    const std::vector<std::vector<std::vector<T>>>& A_list,
    const std::vector<std::vector<std::vector<T>>>& Ainv_list,
    int ne, int lev_max, int deg_max,
    const AdaptiveSamplingConfig& config = {});
```

### RemoveSolvedVariables 策略

嵌套循环 `lev=0..lev_max, deg=0..deg_max`，每次求解前过滤已被支配的变量：

1. **同级低次过滤**: 若 `b[α_s, β_s]` 在 `deg_s < deg` 已解出独立变量，则消除 `α == α_s, β >= β_s`（分量级）的变量
2. **跨级过滤**: 若 `b[α_s, β_s]` 在 `lev_s < lev` 已解出独立变量，则消除 `α >= α_s, β >= β_s` 的变量
3. **列压缩**: 方程矩阵压缩到活跃列后再做高斯消元，然后扩展回完整变量空间

典型效果：bub00 在 lev=2,deg=2 时变量数从 36 降到 19（约 50% 减少）。

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
   - reconstructAllRelations<T>() 嵌套循环 lev=0..lev_max, deg=0..deg_max
   - 每步用 RemoveSolvedVariables 过滤已支配变量
   - 自适应采样求解，返回 LevDegResult<T> 供下游处理
```

---

## Common Pitfalls

### 1. FFInt 类型安全：禁止强制转换负整数

`firefly::FFInt` 仅有 `FFInt(uint64_t)` 构造函数，没有有符号整数构造函数。

```cpp
// 错误: int(-1) 隐式转换为 uint64_t(2^64-1) 再 mod p，产生垃圾值
FFInt x = static_cast<FFInt>(-1);

// 正确: 对正 FFInt 使用 operator-()
FFInt x = -FFInt(1);
```

`static_cast<FFInt>(negative_int)` 始终是 bug。编译器不会对此产生警告。此问题曾影响 `Utilities.hpp` 中的 `sgn()` 及所有从负字面量构造 FFInt 的代码。

### 2. l 循环上界：始终 `incre * k`，绝非仅 `k`

`seriesCoefficient` 存储 `l ∈ [0, incre * k]` 的数据，而非 `[0, k]`。对给定 `k` 遍历 `l` 的循环必须使用 `l <= incre * k`。使用 `l <= k` 会丢弃所有 `l > k` 的系数，产生结果错误但自洽的展开结果，EquationVerify 无法检测。

受影响函数：`step2_computeG`、`step4_computeF2`（卷积）及所有遍历 `C(k, l, cid, j, i)` 中 `l` 的循环。

### 3. 自洽性检查不是独立验证

EquationVerify（`M1*C + Total == 0`）等测试复用了生成系数的同一套方程构建代码。如果 bug 出在方程构建中（如 `sgn()` 符号错误），生成器和验证器都会产生错误但一致的结果。始终将自洽性检查与独立验证方法配对——参考实现对比（MMA）或将结果代回原始（分解前）IBP 方程。

---

## Testing Strategy

### 测试可执行文件

| 测试 | 用途 | 数据类型 | 命令 |
|------|------|----------|------|
| `test_expandFF` | 展开系数计算 | `FFInt` | `./build/test_expandFF` |
| `test_relationFF` | 线性关系重建（需参数） | `FFInt` | `./build/test_relationFF [family] [order] [lev] [deg]` |
| `test_IBPVerification` | IBP矩阵验证 | `FFInt` | `./build/test_IBPVerification` |
| `test_expand_family` | 展开族测试 | `FFInt` | `./build/test_expand_family` |
| `test_load_bub` | bub格式矩阵加载 | `FFInt` | `./build/test_load_bub` |
| `test_family_generate` | 积分族生成 | `FFInt` | `./build/test_family_generate` |
| `test_family_config` | 族配置测试 | `FFInt` | `./build/test_family_config` |
| `test_ibp_analyzer` | IBP分析器测试 | `FFInt` | `./build/test_ibp_analyzer` |
| `test_poly_arith` | 多项式算术 | `FFInt` | `./build/test_poly_arith` |
| `test_recursion_builder` | 递归构建器 | `FFInt` | `./build/test_recursion_builder` |
| `test_recursion_builder_tri` | 三角族递归构建器 | `FFInt` | `./build/test_recursion_builder_tri` |
| `test_region_solver_full` | 区域求解器 | `FFInt` | `./build/test_region_solver_full` |
| `test_ring_builder` | 环构建器 | `FFInt` | `./build/test_ring_builder` |
| `verify_ring_builder` | 环构建器验证 | `FFInt` | `./build/verify_ring_builder` |
| `test_compare_MMA` | MMA对比测试 | `FFInt` | `./build/test_compare_MMA` |
| `verify_MMA_lev1_deg1` | MMA lev1/deg1验证 | `FFInt` | `./build/verify_MMA_lev1_deg1` |

**注意**: 所有二进制数据文件（`.bin`）现在位于 `data/` 目录（前身为项目根目录）。测试和工具通过相对路径 `data/IBPMat_*.bin` 加载。`verify/<family>/` 目录下的 MMA 参考 .bin 文件保持不变。

**注意**: 本项目不使用`ctest`或`add_test()`。测试是独立的可执行文件。

**重要**: 测试必须从**项目根目录**运行（即 `.bin` 数据文件所在的目录），而非 `build/` 目录内。因为测试程序通过相对路径加载二进制数据文件。

### 运行测试

```bash
# 从项目根目录执行
./build/test_expandFF
./build/test_relationFF [family] [order] [lev_range] [deg_range]
./build/test_IBPVerification
./build/test_load_bub
./build/test_expand_family
./build/test_family_generate
./build/test_region_solver_full
```

### tools/ 目录中的独立工具

| 文件 | 用途 |
|------|------|
| `family_generate.cpp` | 积分族生成主工具 |
| `diagnose_ibp.cpp` | IBP矩阵诊断 |
| `test_ff_verify.cpp` | 验证FireFly库基本功能 |
| `test_firefly_simple.cpp` | FireFly最简功能测试 |
| `test_region_solver.cpp` | 区域求解器独立测试 |
| `test_singular_runner.cpp` | Singular子进程测试 |
| `dump_mma_equations.wl` | MMA方程导出脚本 |

### 存档文件

以下历史文件已移至 `archive/`，不再参与主构建：
- 测试: `test_expand.cpp`、`test_recons.cpp`、`test_RelationNew.cpp`、`IBPMatrixLoader_Test.cpp`、`RingDataLoader_test.cpp`
- 构建: `CMakeLists_backup.txt`、`#CMakeLists_full.txt`

### 测试数据文件

测试需要项目根目录下的特定二进制数据文件，所有二进制数据已迁移至 `data/` 目录：
- `data/IBPMat_*.bin` - IBP矩阵二进制数据文件
- `data/RingData_*.bin` - 环数据二进制文件
- `data/ExpansionCache_*.bin` / `data/resCache_*.bin` - 系数缓存（部分自动生成）

**注意**: 测试从 `data/` 目录加载二进制数据。`verify/<family>/` 目录下的 MMA 参考数据文件保持不变。

---

## 验证框架（verify/）

`verify/` 目录是一个结构化跨验证框架，对比 C++ 输出与 Mathematica (MMA) 及 Kira IBP 约化结果：

```
verify/
├── README.md                    # 完整验证流程（逐步命令）
├── FamilyDatabase/
│   ├── FamilyDatabase.wl        # 统一积分族定义（19族，1L+2L+3L，按 L 和 E 排序）
│   └── README.md
├── docs/
│   ├── Test-Expand.md           # 展开验证工作流
│   ├── Test-Relation.md         # 关系验证手册
│   └── CPP-KiraVerify-Debugger.md
├── VerifyUtility/               # MMA 验证工具集脚本
│   ├── LIECoreAlgebra.wl, LIEExpand.wl, LIEReconstruct.wl, LIERegions.wl
│   ├── LIEUtility.wl, LIEFamilyDefine.wl, LIEWorkflow.wl
│   ├── SingularInterface.wl, KiraRuleLoader.wl, M2Kira.wl
│   ├── VerifyExpand-*.wl, VerifyRelation-*.wl, Compare-*.wl
│   └── ExportBinary_IBPMatrix.wl
└── results/bub00/               # 验证结果快照（.m 文件供 MMA 使用）
```

三种验证方法：
1. **CompareVerify**: C++ 输出 vs MMA 参考直接对比
2. **EquationVerify**: 自洽性检查 `M1*C + Total == 0`
3. **SeriesVerify**: 将结果代回原始 IBP 方程验证

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

8. **FFInt 类型安全陷阱**: 详见上方 [Common Pitfalls](#1-ffint-类型安全禁止强制转换负整数)。此处保留原因：此 bug 曾导致所有积分族的展开系数错误，而 EquationVerify（自洽性检查）仍通过，因验证器复用同一错误逻辑。

9. **验证框架**: `verify/` 目录提供 C++ vs MMA vs Kira 的三方交叉验证。在修改核心算法后，除自洽性检查外，应始终运行独立验证。详见 `verify/README.md`。

---

## References

- **FireFly Library**: https://github.com/firefly-library/firefly
- **Eigen Documentation**: https://eigen.tuxfamily.org/
- **IBP Method**: Integration-by-parts identities for Feynman integral reduction (physics)
- **verify/docs/**: 展开与关系重建的验证方法文档
- **docs/ReconstructAlgorithm.md**: 重建算法 MMA vs C++ 实现对比
