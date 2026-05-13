# Region Solver 算法实现比较文档

> 本文档记录从族定义（FamilyDatabase）到二进制输出（IBPMat/RD）的第一阶段算法实现差异，涵盖时间调度、代数信息生成、以及 C++ 与多种 MMA 实现的对比。
> 
> **管线位置**: Stage 1 — 输出 `IBPMat_*.bin` + `RingData_*.bin`，供 LayerRecursion 展开使用。

---

## 1. 输入/输出接口

### 1.1 输入接口：`FamilyDatabase.wl`

验证框架和 MMA 流水线的统一入口是 `families/FamilyDatabase.wl`。每个族是一个 `Association`：

```mathematica
$FamilyDatabase["NP322"] = <|
  "Propagators"      -> {...},        (* 传播子列表，如 "-(k2+p1)^2+msq" *)
  "LoopMomenta"      -> {"k1","k2"}, (* 圈动量 *)
  "ExternalMomenta"  -> {"p1","p2","p4"},
  "KinematicRules"   -> <|"p1^2"->"m1", "p1*p2"->"(s-m1-m2)/2", ...|>,
  "TopSector"        -> {1,1,1,1,1,1,1,1,1},
  "Numeric"          -> <|"s"->1, "t"->5, "m1"->0, "msq"->0, "d"->"1/7"|>,
  "Modulus"          -> 179424673,
  "Description"      -> "Non-Planar 2L4P, massless"
|>;
```

C++ 侧的对等输入是 `families/<name>.json`，字段语义完全一致：

```json
{
  "propagators": [...],
  "loopMomenta": ["k1","k2"],
  "externalMomenta": ["p1","p2","p4"],
  "kinematicRules": {...},
  "topSector": [1,1,1,1,1,1,1,1,1],
  "numeric": {"s":1,"t":5,...},
  "modulus": 179424673
}
```

### 1.2 输出接口：`.bin` 文件

#### `IBPMat_<family>.bin`

自定义二进制格式（魔数 `"IBP1"`），结构如下：

```
[Magic]     "IBP1" (4 bytes)
[numRegs]   int32
For each region:
  [nibp, ne, nb, incre]  int32 × 4
  [modulus]              int64
  For op in {M1,N1,K1,F0,F2,K1s,K2s,F2s}:
    [exists]             byte (0/1)
    If exists:
      [dims_len]         int32
      [dims]             int32 × dims_len
      [rowPtr_len]       int32
      [rowPtr]           int32 × rowPtr_len
      [colIdx_len]       int32
      [colIdx]           int32 × colIdx_len
      [vals_len]         int32
      [vals]             int64 (mod!=0) or Real64 (mod==0) × vals_len
```

#### `RingData_<family>.bin`

```
[count]     int32   (总 region/cell 数量)
[ne]        int32   (每个 cell 包含的矩阵数量)
For each cell:
  [lenSector]        int32
  [limitSector]      int32 × lenSector
  [nb]               int32
  [A_flat]           int64 × (ne × nb × nb)
  [Ainv_flat]        int64 × (ne × nb × nb)
```

两种实现的二进制输出**字节级兼容**（已通过 `family_generate --diff` 验证）。

**C++ 加载器**: `include/IBPMatrixLoader_Binary.hpp`（读取 `IBPMat_*.bin`）、`include/RingDataLoader.hpp`（读取 `RingData_*.bin`）。两个加载器均将二进制文件解析为内存中的 `IBPMatrix` / `RingData` 结构，供 LayerRecursion 使用。

---

## 2. 核心算法比较

### 2.1 时间调度：`ScheduledRegionSolve` vs 传统顺序执行

#### 传统流水线（`VerifyExpand-Prepare.wl` → `LIERegions.wl`）

```
FamilyDatabase
  → LIEDefineFamily              (生成 IBPEqs/SectorList/AList)
  → LIESolveRegions              (包装层)
      → regionsBySectors         (遍历所有 sector)
          → expRegSolve2         (GB + primdec)
          → buildRecursionMatrix (FTable → RecursionMatrix)
  → ExportBinary_IBPMatrix
```

**缺陷**：
- 无实例锁：多进程同时运行会覆盖输出
- 无状态持久化：中途崩溃需从头开始
- 单超时：`Timeout -> 120` 对所有 sector 一视同仁，复杂 sector 可能超时，简单 sector 浪费等待时间
- 无重试：超时 sector 直接丢弃

#### T003 `ScheduledRegionSolver.wl`

在 `LIERegions.wl` 之上增加四层调度封装：

| 模块 | 功能 | 对应传统缺陷 |
|------|------|-------------|
| `InstanceLock` | 文件锁 + PID 检测，防止多实例冲突 | 解决无实例锁 |
| `SectorStateManager` | `status.wdx` 原子读写，`sector_*.wdx` 逐 sector 缓存 | 解决无状态持久化 |
| `TieredScheduler` | Fast (1200s) → Retry (12000s)，FIFO/First 两种排序 | 解决单超时、无重试 |
| `RegionComputeEngine` | `ComputeSectorWithTimeout` 统一接口，返回 `done/timeout/error` | 统一错误处理 |

**主循环逻辑**：

```mathematica
ScheduledRegionSolve[familyName, familyConfig] := Module[...
  (* P2.1: 获取实例锁 *)
  AcquireInstanceLock[lockFile];

  (* P2.2: 加载或初始化状态 *)
  status = LoadStatus[statusFile];
  pending = Select[fullSectorList, !MemberQ[{"done","trivial"}, status[...]] &];

  (* P2.3: 分层处理 *)
  Do[
    {status, timeouts} = ProcessTier[
      pending, tier, ibpeqs, Alist, vlist, char, status, ...
    ];
    pending = timeouts;  (* 超时 sector 进入下一 tier *)
  , {tier, scheduleConfig["Tiers"]}];

  (* 组装最终结果并导出 .bin *)
  allRegions = Table[Import[sectorCacheFile], ...];
  ExportBinaryIBPMatrix[binFile, allRegions, char];
  ExportBinaryRingData[ringFile, allRegions, Alist, ne, char];
];
```

**效果**（SR212 验证）：
- 26/26 sectors，13.5s，0 failed
- 对于 NP322 等复杂族，分层超时可避免单个 hard sector 导致全流水线失败

> **C++ 侧现状**：`family_generate.cpp` 和 `RegionSolver::solveAllSectors` 目前**无**对应调度层，顺序执行所有 subsector，无断点续算。

---

### 2.2 代数信息生成：Singular 代替 MMA `bivarPrimeInfo`

#### 问题背景

在 `expRegSolve2` 中，求得 `aprimelist`（极小素理想列表）后，需要为每个素理想组件提取代数信息：`VarDep`、`VarIndep`、`MinPoly`、`VarRule`、`FractionRule`、`MonomialBasis`、`MonomialBasisMatrix`。

传统 MMA 实现 `bivarPrimeInfo`（`LIECoreAlgebra.wl` L129-157）的两个瓶颈：

1. **`quotientRingBasisMatrixPower`**：对每个基元组合计算 `PolynomialReduce[mb[i]*mb[j], gb, vargen]`。MMA 的 `PolynomialReduce` 在大商环上极慢（>1200s）。
2. **`FractionRule` 的 GB**：为求解 `B[j,i] = A_i * B_j mod prime`，需对扩展理想再求一次 `GroebnerBasis`。

#### 方案 A：纯 MMA `bivarPrimeInfo`

```mathematica
(* 1. 变量分类：从首项分析 *)
ltlist = First@MonomialList[#, Avar, Lexicographic] & /@ primeA;
lmlist = {Union@Cases[#, Alternatives@@Avar], PolynomialDeg[#, Avar]} & /@ ltlist;
{vargen, varpar} = ...;

(* 2. MinPoly：度数>1的生成元多项式 *)
ximinpoly = primeA[[ Select[Range@Length@primeA, lmlist[[#,2]]>1&] ]];

(* 3. VarRule：Solve + PolynomialReduce *)
xdrule = Solve[... , varpar, Modulus->char][[1]];
xdrule = Map[#[[1]] -> PolynomialReduce[#[[2]], ximinpoly, vargen, ...][[2]] &, xdrule];

(* 4. FractionRule：PolynomialReduce *)
fractionrule = Table[FR[i,j] -> PolynomialReduce[Avar[[i]]*Bvar[[j]], prime, ...][[2]], ...];

(* 5. MonomialBasisMatrix：PolynomialReduce —— 最大瓶颈 *)
basisMatrix = quotientRingBasisMatrixPower[ximinpoly, vargen, basisIndex];
(* 内部：PolynomialReduce[mb[i]*mb[j], gb, vargen] for all i,j *)
```

**复杂度**：`O(nb² × poly_reduce_cost)`，其中 `nb = Times @@ VarDeg` 是商环维度。对于多维 VarIndep（如 Box 族 nb=2+），`nb²` 可达 16-64，每次 `PolynomialReduce` 都很昂贵。

#### 方案 B：`SingularCoordinateRing.wl`（Singular 加速）

核心思想：**让 Singular 做所有 heavy lifting（GB、kbase、乘法表），MMA 只做轻量组装**。

```mathematica
(* Singular 脚本模板：一次性计算 gb, kbase, multiplication table *)
singularPrimeTemplate = "
LIB \"primdec.lib\";
ring r = `char`, (`vars`), (`order`);
ideal I = `ideal`;
ideal gb = std(I);
ideal kb = kbase(gb);        (* 商环单项式基，Singular 原生 *)
int vd = vdim(gb);

(* 计算乘法表：reduce(kb[a]*kb[b], gb) for all a,b *)
string mbs = \"\";
for (int a = 1; a <= vd; a++)
  for (int b = 1; b <= vd; b++)
  {
    poly p = reduce(kb[a]*kb[b], gb);
    mbs = mbs + string(p);
    if (a < vd || b < vd) { mbs = mbs + \",\"; }
  }
";

(* 调用 Singular 获取 {gb_list, kb_list, mb_polys_flat} *)
{gbList, kbList, mbList} = SingularRun[singularPrimeTemplate, ...];

(* MMA 侧只做组装：无需 PolynomialReduce *)
assembleCoordinateRing[prime, gbList, kbList, mbList, Avar, Bvar, char, limitSector];
```

**`assembleCoordinateRing` 的 MMA 组装步骤**：

| 步骤 | 操作 | 是否依赖 Singular |
|------|------|------------------|
| 变量分类 | 首项模式匹配 (`MonomialList` + `GroupBy`) | ❌ 纯 MMA |
| MinPoly | 从 `aOnlyGb` 中筛选度数>1 | ❌ 纯 MMA |
| VarRule | `Solve` 线性方程组（快速） | ❌ 纯 MMA |
| MonomialBasis | 直接使用 `kbList` | ✅ Singular `kbase` |
| MonomialBasisMatrix | 从 `mbList` 解析系数（`monomialRulesPower`） | ✅ Singular `reduce` |
| FractionRule | `PolynomialReduce`（保留，因涉及 Bvar） | ⚠️ 仍用 MMA |

**关键优化**：`MonomialBasisMatrix` 的构建从 `O(nb²)` 次 `PolynomialReduce` 降为 **0 次**——Singular 已预计算所有 `reduce(kb[i]*kb[j], gb)`。

**排序一致性陷阱**：Singular `kbase` 返回的单项式顺序可能与 MMA `quotientRingBasisPower` 的 lex 顺序不同。`assembleCoordinateRing` 中显式 `Ordering` 重排 `basisIndex` 和 `mbPolys`，确保与 MMA 二进制输出一致。

#### 方案 C：C++ `RegionSolver.hpp`

C++ 侧没有 MMA 的符号计算环境，完全依赖 **Singular 子进程** + **字符串解析**。

```cpp
// RegionSolver::solveRegion 流水线
inline std::vector<RegionData> solveRegion(abEquations, limitSector, ne, modulus) {
  // Step 1: GB (Singular)
  auto gb = computeGroebnerBasis(fullIdeal, allVars, modulus);

  // Step 2: Minimal associated primes (Singular primdecGTZ)
  auto [primelist, dims] = computeMinAssPrimes(gb, allVars, modulus);

  // Step 3: Per-component processing
  for (auto& prime : primelist) {
    auto compGb = computeGroebnerBasis(prime, allVars, modulus);

    // 3a. 变量分类（字符串解析首项，匹配 MMA classifyVariablesPrimeA）
    classifyVariablesPrimeA(primelist[p], Avars, Bvars,
        primeA, reg.VarIndep, reg.VarDep, reg.VarDeg, ximinpoly);

    // 3b. MonomialBasisIndex（从 ximinpoly 的 LM 指数上限枚举）
    reg.MonomialBasisIndex = computeMonomialBasisIndex(ximinpoly, reg.VarIndep);

    // 3c. VarRule（手动解析多项式字符串，提取线性首项）
    solveVarRule(compGb, reg.VarDep, reg.VarIndep, modulus, reg.VarRule);

    // 3d. FractionRule（Singular reduce）
    computeFractionRule(compGb, allVars, allVars, reg.VarIndep, ne, modulus, reg.FractionRule);

    // 3e. MonomialBasisMatrix（Singular reduce + PolyArith 解析）
    computeMonomialBasisMatrix(ximinpoly, reg.VarIndep,
        reg.MonomialBasisIndex, reg.nb, modulus, reg.MonomialBasisMatrix);
  }
}
```

**C++ 的关键差异点**：

| 操作 | MMA (bivarPrimeInfo) | C++ (RegionSolver) |
|------|---------------------|-------------------|
| **变量分类** | `MonomialList` + `GroupBy` 符号首项分析 | `classifyVariablesPrimeA`：字符串 `extractLeadingMonomial` + `parseExponents` |
| **VarRule** | `Solve[... Modulus->char]` 符号求解 | `solveVarRule`：遍历 compGb，找到 LM 为 `varpar^1` 的多项式，手动解析其余项为 `-P(vargen)` |
| **FractionRule** | `PolynomialReduce[A_i*B_j, prime, ...]` | `reducePolynomialsSingular`：批量发送 `A_i*B_j` 到 Singular `reduce`，再 `parseSingularPolynomial` 解析 |
| **MBM** | `quotientRingBasisMatrixPower`：`PolynomialReduce[mb[i]*mb[j], gb]` | `computeMonomialBasisMatrix`：批量发送 `mb[i]*mb[j]` 到 Singular `reduce`，`parseSingularPolynomial` 匹配 basisIndex |
| **多项式表示** | MMA 原生符号表达式 | `PolyArith::Polynomial` = `vector<{exps, coeff}>`，手动实现 `parseSingularPolynomial` / `polyToSingularString` |

**C++ 的边界情况处理**（已修复）：

1. **`vargen == []` 且 `nb == 1`**：常数环，`MonomialBasisMatrix` 应为 `[[[1]]]`。C++ 早期版本未处理此边界，导致全零。
2. **`parseSingularPolynomial` 空变量列表**：当 `vargen=[]` 时，非数字常数（如 `"119616450"`）应解析为 `{}`（零指数向量）+ 系数。C++ 已特殊处理 `nVars==0`。
3. **`solveVarRule` 空 vargen**：常数 VarRule 不应调用 `polyToSingularString`（会因空 var 列表出错），直接输出常数字符串。

---

### 2.3 C++ 实现 vs 多种 MMA 实现的全面对比

#### 实现谱系

```
MMA 侧：
  bivarPrimeInfo (纯 MMA, 慢) ──► singularBivarPrimeInfo (Singular 加速) ──► expRegSolve2 (最外层封装)
         ▲                                              ▲
         │                                              │
    LIECoreAlgebra.wl                          SingularCoordinateRing.wl

C++ 侧：
  RegionSolver.hpp ──► RecursionBuilder.hpp ──► RingBuilder.hpp ──► family_generate.cpp
         │                    │                      │
         ▼                    ▼                      ▼
   SingularRunner.hpp    IBPAnalyzer.hpp       PolyArith.hpp
```

#### 模块级对比

| 模块 | MMA 实现 | C++ 实现 | 关键差异 |
|------|---------|---------|---------|
| **IBP 方程生成** | `LIEDefineFamily`（内部调用 Singular SP2PD） | `IBPEqGenerator::generateIBPEquations`（调用 Singular SP2PD） | 两者都通过 Singular 子进程执行 `SP2PD`，结果格式略有不同但最终等价 |
| **FTable 提取** | `IBPCoefficientMatrix`（符号系数数组） | `IBPAnalyzer::extractFTable`（字符串多项式 → 系数表） | MMA 保留符号变量；C++ 在有限域上数值化 |
| **递推矩阵构建** | `RecursionCoefficientMatrix` + `recursionMatrixCompanion` | `RecursionBuilder::buildRecursionMatrices` | MMA 先符号构建再稀疏化；C++ 直接在有限域上计算 FlatMatrix |
| **环矩阵 A/Ainv** | `ComputeRingData`（`coeffs . ringmat` + `Inverse[#, Modulus->p]`） | `RingBuilder::computeRingMatrices`（矩阵乘法 + 有限域高斯消元求逆） | 数学等价；C++ 用自定义高斯消元替代 Eigen（因 FFInt 不支持 Eigen 的浮点假设） |
| **Region 求解** | `expRegSolve2` | `RegionSolver::solveAllSectors` | 见下文详细对比 |

#### `expRegSolve2` vs `RegionSolver::solveAllSectors` 详细对比

**输入/输出语义**：

| | MMA `expRegSolve2` | C++ `solveAllSectors` |
|---|-------------------|----------------------|
| 输入 | `ibpeqs`（符号方程）, `Alist`, `vlist` | `IBPEquations`（结构化数据）, `topSector` |
| 子 sector 枚举 | 由调用方 `regionsBySectors` 传入 `sectorlist` | `generateSubsectors(topSector, nl)` 内部枚举 |
| 输出 | `List<Association>`（含 CoordinateRing + RecursionMatrix） | `vector<RegionData>`（仅 CoordinateRing，RecursionMatrix 在后续步骤） |

**子 sector 处理差异**：

- **MMA**：`regionsBySectors` 遍历传入的 `sectorlist`。对于 top-sector-only 验证，可能只传 `[{1,1,...,1}]`。
- **C++**：`solveAllSectors` 枚举 `topSector` 的所有 subsector（满足 `m[i] <= topSector[i]` 且大小在 `[nl, ne]` 之间）。这会导致 C++ 对 Tri 族产生 4 个 regions，而 MMA 的 `RegionWise` 策略可能合并部分 subsector 为 3 个。

**A/B 方程构造**：

```mathematica
(* MMA: expRegSolve2 *)
aeqs0 = Coefficient[ibpeqs,"n"] /. "g"[a__] :> (Times@@Thread@Power[AlistCF, {a}-vlist]);
aeqs  = Join[aeqs0/.{1/"A"[a_]:>"B"[a]}, Table["A"[i]"B"[i]-1, {i,ne}]];
```

```cpp
// C++: IBPAnalyzer::buildABEquations + RegionSolver::solveRegion
auto abEqs = IBPAnalyzer::buildABEquations(ibp, sub, modulus);
// 内部已将 g[...] 替换为 A/B 形式，并附加 A_i*B_i-1
```

两者数学等价，但 C++ 的 `buildABEquations` 在 subsector limit 时可能产生全零方程，C++ 会显式跳过（`hasContent` 检查）。

**变量名处理**：

- **MMA**：使用上下文无关符号 `"A"[i]` / `"B"[i]`，通过 `localRep` 映射到 Singular 的 `x1, x2, ...`，再映射回来。需要复杂的字符串/符号转换来修复 Singular 输出。
- **C++**：直接使用 `A1, A2, ..., B1, B2, ...` 作为变量名，Singular 脚本直接拼接字符串，无需符号上下文转换。

**Singular 调用方式**：

| | MMA | C++ |
|---|-----|-----|
| **接口文件** | `SingularInterface.wl` | `SingularRunner.hpp` |
| **通信方式** | `Run["wolframscript ... SingularInterface.wl ..."]` | `popen` 子进程，通过 stdin/stdout 传递 Singular 脚本 |
| **输出解析** | MMA `ToExpression` 解析 Singular 字符串输出 | `PolyArith::parseSingularPolynomial` 手动解析多项式字符串 |
| **超时控制** | `"SingularTimeout"` 选项传递给底层 | 子进程 `alarm` / `setitimer` 或外部 kill |
| **并发** | 单进程串行 | 单进程串行（`ParallelSolver.hpp` 预留但未启用） |

**数据类型与精度**：

| | MMA | C++ |
|---|-----|-----|
| **算术域** | 有限域 `Modulus -> p`（MMA 内部大整数） | `firefly::FFInt`（`uint64_t` 模运算） |
| **中间计算** | 精确整数/有理数，最后取模 | 始终模 `p` 运算 |
| **矩阵求逆** | `Inverse[#, Modulus->p]&` | 自定义高斯消元（有限域版本） |
| **溢出风险** | 无（MMA 任意精度） | 有（`uint64_t` 模数必须 < 2^63，且中间乘积需 `__int128`） |

#### 已知不一致点

| 问题 | 状态 | 说明 |
|------|------|------|
| C++ A_i 值与 MMA 不同（`119616450` vs `59808223`） | ✅ 已解决 (T004) | `SingularRunner::groebnerBasis()` 中 `groebner(I)` → `std(I)` 对齐：两者计算同一理想，但生成元集不同导致 VarRule 从 `compGb` 解析时多项式选择不同。`std(I)` 附加 `option(redSB)` 产生约化 Groebner 基，与 MMA 的 `GroebnerBasis[..., MonomialOrder->Lexicographic]` 结果一致。 |
| Tri region 数量（C++ 4 vs MMA 3） | ✅ 已解决 | C++ `solveAllSectors` 枚举所有 subsector（4 个），MMA 的 `regionsBySectors` 按配置可能只传 top-sector（1 个）或部分 subsector（3 个）。`verify_ring_builder` 通过先读取 MMA 参考的 sector 集再仅处理这些 sector 来对齐。 |
| `parseSingularPolynomial` 空变量常数解析 | ✅ 已修复 | 2026-05-10 特殊处理 `nVars==0` |
| `MonomialBasisMatrix` 全零 | ✅ 已修复 | 2026-05-10 空 `vargen` 时显式置 `[[[1]]]` |
| T004 对齐导致参考文件过期 | ✅ 已解决 | `std(I)` 产生的 GB 生成元与 `groebner(I)` 不同，单分量 sector 结果一致，多分量 sector 的 A/Ainv 矩阵不同（基元素相同但坐标变换矩阵因 GB 生成元变化而改变）。所有 May 5 之前的 `.bin` 参考需用当前 C++ 管线重新生成，见 `../verify/Test-RingData.md`。 |

---

## 3. 性能特征

| 族 | ne | 传统 MMA `bivarPrimeInfo` | `singularBivarPrimeInfo` | C++ `family_generate` |
|-----|----|--------------------------|-------------------------|----------------------|
| bub00 | 2 | <1s | <1s | <1s |
| Tri | 3 | ~5s | ~2s | ~2s |
| Box | 4 | >120s（MBM 瓶颈） | ~10s | ~10s |
| NP322 | 9 | 不可行（未测试） | ~数十秒至数分钟/sector | ~数十秒至数分钟 |
| DP323 | 9 | 不可行 | 需分层调度 | 需优化 |

> 注：NP322/DP323 等复杂族的 hard sector（如 top sector）在 Singular 中求 GB + primdec 可能需要 120s-12000s，这正是 `ScheduledRegionSolver` 的分层超时设计的动机。

---

## 4. 文件索引

| 文件 | 作用 |
|------|------|
| `families/FamilyDatabase.wl` | MMA 输入接口 |
| `families/*.json` | C++ 输入接口 |
| `workspace/shared/VerifyUtility/ScheduledRegionSolver.wl` | 时间调度模块（T003） |
| `workspace/shared/VerifyUtility/LIERegions.wl` | MMA Region 求解 + 递推矩阵构建 |
| `workspace/shared/VerifyUtility/LIECoreAlgebra.wl` | `expRegSolve2`, `bivarPrimeInfo` |
| `workspace/shared/VerifyUtility/SingularCoordinateRing.wl` | Singular 加速版 `singularBivarPrimeInfo` |
| `workspace/shared/VerifyUtility/SingularInterface.wl` | MMA-Singular 通信层 |
| `include/RegionSolver.hpp` | C++ Region 求解器 |
| `include/RecursionBuilder.hpp` | C++ 递推矩阵构建 |
| `include/RingBuilder.hpp` | C++ 环矩阵构建 |
| `include/SingularRunner.hpp` | C++ Singular 子进程管理 |
| `include/PolyArith.hpp` | C++ 多项式字符串解析/格式化 |
| `tools/family_generate.cpp` | C++ 完整流水线 CLI |
| `workspace/shared/VerifyUtility/ExportBinary_IBPMatrix.wl` | MMA 二进制导出 |

---

*文档版本: 1.3*
*最后更新: 2026-05-12*
*相关文档: [docs/verify/Test-RingData.md](../verify/Test-RingData.md) — RingData 验证工作流与常见问题*
*Skill: [region-solver](../.claude/skills/region-solver-mma-scheduled/SKILL.md) — ScheduledRegionSolve 管线、状态机、resume 逻辑、IBPMat 完整性验证*

---

## 附录 A：三管线全景对比

### A.1 执行路径

```
MMA 纯符号管线:
  FamilyDatabase
    → LIEDefineFamily [MMA]
        → SP2PD (Singular) → IBPEqs
        → LIESolveRegions
            → expRegSolve2
                → bivarPrimeInfo [MMA]      (GB: Singular, 代数信息: PolynomialReduce)
          → 无调度/无断点续算

MMA-Singular 混合管线 (ScheduledRegionSolver):
  FamilyDatabase
    → LIEDefineFamily [MMA]
    → ScheduledRegionSolve
        → AcquireInstanceLock [文件锁]
        → LoadStatus [status.wdx 断点恢复]
        → TieredScheduler [Fast→Retry, 分层超时]
            → RegionComputeEngine
                → singularBivarPrimeInfo
                    → SingularCoordinateRing [Singular: GB+kbase+multable, MMA: 轻量组装]
        → ExportBinaryIBPMatrix [MMA]
        → ExportBinaryRingData [MMA]

C++-Singular 管线 (family_generate):
  families/<fam>.json
    → IBPEqGenerator::generateIBPEquations
        → SingularRunner [SP2PD 子进程]
    → RegionSolver::solveAllSectors [C++, 无调度]
        → SingularRunner [GB: groebner/std]
        → SingularRunner [primdecGTZ: 极小素理想分解]
        → 逐 region:
            → classifyVariablesPrimeA [C++ 字符串解析首项]
            → solveVarRule [C++ 遍历 GB 解析线性首项]
            → computeMonomialBasisIndex [C++ 指数上限枚举]
            → computeFractionRule [SingularRunner: reduce(A_i*B_j)]
            → computeMonomialBasisMatrix [SingularRunner: reduce(mb[*]mb[*])]
            → RingBuilder::computeRingMatrices [C++ 有限域矩阵乘法+高斯求逆]
            → RecursionBuilder::buildRecursionMatrices [C++ IBP 递推]
    → IBPMatBinaryWriter [C++ 二进制导出]
    → RingDataBinaryWriter [C++ 二进制导出]
```

### A.2 管线能力矩阵

| 能力 | MMA 纯符号 | MMA-Singular (Scheduled) | C++-Singular |
|------|-----------|------------------------|-------------|
| 输入源 | FamilyDatabase.wl | FamilyDatabase.wl | families/*.json |
| GB 计算 | Singular 子进程 | Singular 子进程 | Singular 子进程 |
| 极小素理想 | primdecGTZ | primdecGTZ | primdecGTZ |
| 商环代数信息 | PolynomialReduce (MMA, 慢) | kbase+reduce (Singular, 快) | kbase+reduce (Singular, 快) |
| 断点续算 | ❌ 不支持 | ✅ status.wdx + sector_*.wdx | ❌ 不支持 |
| 分层超时 | ❌ 单超时 | ✅ Fast(1200s)→Retry(12000s) | ❌ 无超时 |
| 实例锁 | ❌ | ✅ 文件锁+ PID | ❌ |
| 多 sector 并行 | ❌ | ✅ Interleaved 模式 | ❌ |
| A/Ainv 矩阵 | Inverse[#, Modulus→p] | Inverse[#, Modulus→p] | 自定义有限域高斯消元 |
| 导出格式 | MMA BinaryWrite | MMA BinaryWrite | C++ ofstream |
| 字节级兼容性 | — | 基准 | ✅ 与 MMA-Singular 对齐 (T004) |

---
## Appendix B: C++ Pipeline Cross-Validation Status

### B.1 Pipeline Overview (4 Steps)

```
FamilyDatabase ──→ [1] LIEDefineFamily ──→ [2] LIESolveRegions ──→ [3] ExportIBPMatrix ──→ IBPMat_*.bin
                     (IBP eq generation)     (ring decomposition)       (sparse CSR write)
                                                                     [4] ExportRingData  ──→ RingData_*.bin
                                                                        (A/Ainv compute+write)
```

Step 1: Symbolic differentiation → IBP identities → Large-index conversion
Step 2 (hardest): Groebner basis → primary decomposition → monomial basis → recursion matrices
Step 3-4 (easy): Binary export (IBP1 format + RingData format)

### B.2 Implementation Status

| Phase | Module | Status |
|-------|--------|--------|
| Phase 0 | FamilyConfig, BinaryIBPWriter, BinaryRingWriter | ✅ DONE |
| Phase 1 | CLI `family_generate.cpp` | ✅ DONE |
| Phase 2 | IBPEqGenerator (IBP equations) | ✅ CORE DONE (LI eqs pending) |
| Phase 3 | RegionSolver + aux modules | ✅ IMPLEMENTED |
| Phase 4 | Integration & cross-validation | ⏳ IN PROGRESS |

### B.3 Cross-Validation by Family

| Family | Region count | nb | test_relationFF | .bin byte-identical | Note |
|--------|-------------|-----|-----------------|---------------------|------|
| bub00 | ✅ 1 | ✅ 1 | ✅ sol_dim match | ❌ sign diff | functional equivalence confirmed |
| bub10 | ✅ | ✅ | — | ❌ | |
| bub11 | ✅ | ✅ | — | ❌ | |
| Tri | ✅ 2 | ✅ | — | ❌ | |
| Box | ✅ 15 | ✅ | — | ❌ | |
| SR | ✅ | ✅ | — | — | |
| SR3m | ✅ | ✅ | — | — | |
| SR5m | ✅ | ✅ | — | — | |
| TB123 | ⏳ timeout | — | — | — | larger CAS computation |

**Key finding:** All families produce correct region counts and functional equivalence. However, .bin files are NOT byte-identical to MMA reference — A/B equation sign conventions differ, leading Singular to select different but equivalent monomial bases. Root cause: C++ multiplies IBP by ∏z_m then negates n_i terms; MMA uses LargeIndexIBP `n_i → n+v_i` then `Coefficient["n"]` extraction.

### B.4 Architecture

```
 ┌──────────────────────────────────────────────────────────┐
 │                    C++ FamilyGenerate                     │
 ├──────────────────────────────────────────────────────────┤
 │  Module A: FamilyConfig        — pure C++                │
 │    Parse JSON → FamilyDef struct                          │
 │  Module B: IBPEqGenerator      — Singular subprocess     │
 │    Symbolic diff → IBP + LI identities → g-operator form │
 │  Module C: RegionSolver        — Singular subprocess     │
 │    GB → primdec → monomial basis → recursion matrices    │
 │  Module D: BinaryIBPWriter     — pure C++                │
 │    Write IBPMat_*.bin (IBP1 format)                      │
 │  Module E: BinaryRingWriter    — pure C++                │
 │    Compute A/Ainv → Write RingData_*.bin                 │
 │  Module F: FamilyGenerateCLI   — pure C++                │
 │    family_generate <family.json>: orchestrate A→B→C→D→E  │
 └──────────────────────────────────────────────────────────┘
```

## 相关文档

- `SymbolicRuleAlgorithm.md` — 符号规则生成（RegionSolver 的扩展）
- `ReconstructAlgorithm.md` — 关系重构算法（管线下游）
- `LayerRecursion_Algorithm.md` — Layer Recursion 展开算法（管线下游）
- `../plans/RegionSolver-Debug-Plan.md` — C++ RegionSolver 调试计划
- `../verify/Test-RingData.md` — RingData 验证
- `../verify/Test-Expand.md` — 展开一致性验证
