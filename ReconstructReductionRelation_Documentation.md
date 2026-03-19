# ReconstructReductionRelation 程序文档

## 概述

`ReconstructReductionRelation.wl` 是一个用于从大型指标展开（Large Index Expansion）中重建IBP（Integration By Parts）约化关系的Mathematica包。支持多区域、多解以及代数扩张域上的计算。

---

## 主要功能

### 1. ReconstructReductionRelation（主函数）

从大型指标展开中重建约化关系。

**语法：**
```mathematica
ReconstructReductionRelation[rankLevel, degree, order, hexpnList, aregList, topsec, vList, opts]
```

**参数：**
| 参数 | 类型 | 说明 |
|------|------|------|
| `rankLevel` | Integer | 最大秩（rank）级别 |
| `degree` | Integer | 系数多项式的最大次数 |
| `order` | Integer | 展开阶数 |
| `hexpnList` | List | 展开系数表，格式：`{{{h00,h01,...},{h10,h11,...}}, ...}`，每区域每解一个列表 |
| `aregList` | List | 区域规范列表，每个元素是Association，包含代数扩张数据 |
| `topsec` | List | 顶层sector规范 |
| `vList` | List | 变量列表 `{v1, v2, ..., vn}` |
| `opts` | Options | 算法选项 |

**选项：**
| 选项 | 默认值 | 说明 |
|------|--------|------|
| `Modulus` | `0` | 模数（0表示精确计算/Q，p表示有限域F_p） |
| `Verbose` | `False` | 调试输出 |
| `"AnsatzMode"` | `0` | 模式：0=金字塔，1=点金字塔，2=星形 |
| `"Strategy"` | `"Incremental"` | 策略："Direct"、"Incremental"、"RegionByRegion" |
| `"PlateauSize"` | `1` | 收敛稳定性平台大小 |
| `"DegreeByDegree"` | `True` | 是否逐次添加方程 |
| `"LinearSolver"` | `Automatic` | 线性求解器 |
| `"Numeric"` | `{}` | 固定数值 |

**返回值：**
- 每级别每次数的关系列表

---

## 程序结构

### 模块划分

```
ReconstructReductionRelation.wl
│
├─ Package Header（包头部）
│   ├─ 导出函数声明
│   └─ Usage定义
│
├─ Configuration and Global State（配置和全局状态）
│   ├─ $RRRVerbose - 详细输出开关
│   ├─ $LevelCount - 缩进级别计数
│   ├─ WithIndent - 缩进块
│   ├─ PrintV - 条件打印
│   └─ MyTimer - 计时器
│
├─ Options（选项定义）
│   └─ ReconstructReductionRelation选项
│
├─ Utility Functions（工具函数）
│   ├─ FirstNonZero - 首个非零索引
│   ├─ FFReduce - 模约简
│   ├─ GenerateSeeds - 生成指数向量
│   ├─ CircledMinus - 子sector结构
│   ├─ DotRankSeeds - 点-秩种子生成
│   ├─ SectorMO/SectorMOneg - Sector单序
│   └─ CoefficientEquations - 提取系数方程
│
├─ Module 1: GenerateAnsatz（生成Ansatz）
│   └─ GenerateAnsatz[mode, rank, maxDeg, limitSector, nVars, vlist]
│       模式：0=金字塔，1=点金字塔，2=星形
│
├─ Module 2: Basis Power Computation（基幂计算）
│   ├─ MonomialRulesPower - 单项式规则
│   ├─ BuildVariableMatrices - 构建变量矩阵
│   ├─ ComputePowerTable - 幂表计算（带memoization）
│   └─ ComputeBasisPower - 主入口，计算A^(-alpha)
│
├─ Module 3: Expansion Table Management（展开表管理）
│   └─ UpdateExpansionTable - 更新展开表
│
├─ Module 4: Variable Management（变量管理）
│   ├─ RemoveSolvedVariables - 移除已求解变量
│   ├─ InitializeVariables - 初始化变量（按vFixed分组）
│   ├─ SetupCoefficientTable - 设置系数表
│   ├─ DispatchLinearSolve - 分派线性求解
│   ├─ ReduceSolve - 约化求解
│   └─ ResubstituteSolutions - 回代解
│
├─ Module 5: Equation Solver Core（方程求解核心）
│   ├─ SolveDegreeEquations - 状态控制求解循环
│   ├─ SetupEquationsAll - 全区域一起设置方程
│   └─ SetupEquationsSingle - 单区域设置方程
│
└─ Main Function（主函数）
    └─ ReconstructReductionRelation - 主流程
```

---

## 关键算法流程

### 主流程（ReconstructReductionRelation）

```
1. 初始化
   └─ 设置选项、变量、vFixed

2. 生成Ansatz（Step 1）
   └─ GenerateAnsatz生成各级别ansatz

3. 准备基幂表（Step 2）
   └─ ComputeBasisPower计算每区域的A^(-alpha)

4. 初始化展开表（Step 3）
   └─ 创建空的expnTable结构

5. 主循环：遍历seed级别（Step 4-6）
   │
   ├─ 更新展开表
   │   └─ UpdateExpansionTable
   │
   └─ 循环：遍历系数次数deg
       │
       ├─ 移除已求解变量
       │   ├─ RemoveSolvedVariables
       │   └─ InitializeVariables
       │
       ├─ 设置系数表
       │   └─ SetupCoefficientTable
       │
       ├─ 求解循环
       │   └─ SolveDegreeEquations
       │       ├─ 策略分支：All regions / RegionByRegion
       │       ├─ SetupEquationsAll/SetupEquationsSingle
       │       ├─ DispatchLinearSolve
       │       └─ 稳定性检查
       │
       └─ 回代解
           └─ ResubstituteSolutions

6. 输出汇总
   └─ 返回relationList
```

### 变量管理

| 变量名 | 含义 | 生命周期 |
|--------|------|----------|
| `bVars` | 完整变量集合 | 每deg循环初始化，常量 |
| `bVarReg` | 当前未求解变量 | 动态递减 |
| `bSolAcc` | 累积的所有解 | 持续增长 |
| `bSolNew` | 本轮新求解的解 | 临时，每轮更新 |

---

## 数据结构规范

### hexpnList（展开系数表）

格式：`{{{h00, h01, ...}, {h10, h11, ...}}, ...}`

- 外层List：每区域一个元素
- 中层List：每解一个元素  
- 内层List：每阶展开一个元素（h_k为ν的多项式）

### aregList（区域规范）

每个元素是Association：
```mathematica
<|
  "VarDep" -> {A1, A2, ...},           (* 依赖变量 *)
  "VarRule" -> {A[1]->..., ...},       (* 变量替换规则 *)
  "VarIndep" -> {...},                  (* 独立变量 *)
  "MonomialBasis" -> {...},            (* 单项式基 *)
  "MonomialBasisIndex" -> {...},       (* 基索引 *)
  "MonomialBasisMatrix" -> <|...|>,    (* 基矩阵 *)
  "LimitSector" -> {...}               (* 极限sector *)
|>
```

### expnTable（展开表）

三维结构：`expnTable[[region, solution, order+1]][alpha]`

- 每个元素是Association，键为alpha，值为展开值

---

## 实现说明

### 动态方程添加机制（SetupEquationsAll，第718-749行）

代码中保留了一个被禁用的动态方程添加功能：

```mathematica
(* NOTE: Dynamic equation addition (disabled but preserved for future use)
   Purpose: When the current equation system is underdetermined 
   (Length[eqs] < Length[bVarReg]), dynamically add higher-order equations.
   
   Behavior when enabled:
   - Increment order k while k <= maxorder
   - Build relansatz for higher orders using expntable
   - Extract new equations via CoefficientEquations
   - Append to eqs until sufficient or maxorder reached
   
   Current status: Disabled via If[False, ...] wrapper.
   To enable: Change If[False, ...] to If[True, ...] or remove the wrapper.
*)
If[False,
  While[(degreebydegree == False && (k+1 <= Min[maxorder, level+deg-1])) ||
        (degreebydegree == True && (Length@eqs < Length[bVarReg]) && (k < maxorder)),
    ...
  ]
];
```

**用途**：当当前阶数产生的方程不足以求解所有待定系数时（欠定系统），自动添加更高阶的展开方程，直到：
- 方程数 ≥ 变量数，或
- 达到最大阶数 `maxorder`

**当前状态**：默认禁用（`If[False, ...]`），如需启用可改为 `If[True, ...]`。

---

## 已知限制

1. **中文字符问题**：部分注释显示为Unicode转义序列

---

## 使用示例

### 基本使用

```mathematica
(* 加载包 *)
Get["ReconstructReductionRelation.wl"]

(* 定义变量 *)
vList = {v1, v2};

(* 准备展开数据（示例）*)\nhexpnList = {
  {{
    <|{0,0} -> 1, {1,0} -> v1|>,    (* 0阶 *)
    <|{0,0} -> 0, {1,0} -> 1|>     (* 1阶 *)
  }}
};

(* 定义区域（无代数扩张）*)\naregList = {
  <|
    "VarDep" -> {A1, A2},
    "VarRule" -> {A[1]->1, A[2]->1},
    "VarIndep" -> {},
    "MonomialBasis" -> {1},
    "MonomialBasisIndex" -> {{}},
    "MonomialBasisMatrix" -> <|{} -> {{1}}|>,
    "LimitSector" -> {1, 1}
  |>
};

(* 运行 *)\nresult = ReconstructReductionRelation[
  2,        (* rankLevel *)
  1,        (* degree *)
  3,        (* order *)
  hexpnList,
  aregList,
  {1, 1},   (* topsec *)
  vList,
  Verbose -> True
];
```

---

## 版本信息

- 文件名：`ReconstructReductionRelation.wl`
- 总行数：约1004行
- 主要导出函数：`ReconstructReductionRelation`
- 辅助导出函数：`GenerateAnsatz`、`ComputeBasisPower`
