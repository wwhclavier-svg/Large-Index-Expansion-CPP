# 项目重组方案

目标：将项目文件分为 原始代码 / 可执行文件 / 测试工作区 三类。

## 已执行：workspace/ 按 CXXX 任务分离

每个 CXXX 任务拥有独立的工作区目录，同一 family 的不同方法（MMA vs C++-Singular）各占一个目录。

```
workspace/
├── C001-NP222var-region-solve/     219M  MMA region-solve
├── C002-NP322-region-solve/        131M  MMA region-solve
├── C003-TB123-region-solve/        3.0M  MMA region-solve
├── C004-DB313-region-solve/        4.9M  MMA region-solve
├── C005-SR212-3m-solve-test/       360M  solve + relation-test
├── C006-bub00-relation-test/       528K  relation-test
├── C007-bub10-relation-test/       544K  relation-test
├── C008-bub11-relation-test/       544K  relation-test
├── C009-Tri-relation-test/         1.6M  relation-test
├── C010-Box-relation-test/         48M   relation-test
├── C011-SR212-relation-test/       93M   relation-test
├── C012-SR212-5m-region-solve/     192K  MMA region-solve
├── C013-TB123-region-solve-cpp/    4.0K  C++-Singular (空)
├── C014-Penta1L-relation-test/     844M  relation-test
├── C015-NP222var-cpp-region-solve/ 12K   C++-Singular
├── C017-DB313-region-solve-cpp/    1.7M  C++-Singular scheduled
└── shared/                         98M   VerifyUtility + FamilyDatabase + 日志
```

### 迁移操作

| 操作 | 说明 |
|------|------|
| `mv verify/NP222var workspace/C001-NP222var-region-solve` | C001 |
| `mv verify/NP322 workspace/C002-NP322-region-solve` | C002 |
| `mv verify/TB123 workspace/C003-TB123-region-solve` | C003 |
| `mv verify/DB313 workspace/C004-DB313-region-solve` | C004 |
| `mv test/SR212-3m workspace/C005-SR212-3m-solve-test` | C005 |
| `mv test/bub00 workspace/C006-bub00-relation-test` | C006 |
| `mv test/bub10 workspace/C007-bub10-relation-test` | C007 |
| `mv test/bub11 workspace/C008-bub11-relation-test` | C008 |
| `mv test/Tri workspace/C009-Tri-relation-test` | C009 |
| `mv test/Box workspace/C010-Box-relation-test` | C010 |
| `mv test/SR212 workspace/C011-SR212-relation-test` | C011 |
| `mv verify/SR212-5m workspace/C012-SR212-5m-region-solve` | C012 |
| `mkdir workspace/C013-TB123-region-solve-cpp` | C013 (原外部路径已删) |
| `mv test/Penta1L workspace/C014-Penta1L-relation-test` | C014 |
| `mv region-solver-cppB-test/NP222var workspace/C015-NP222var-cpp-region-solve` | C015 |
| `mv cache/DB313_scheduled workspace/C017-DB313-region-solve-cpp` | C017 |
| `mv verify/VerifyUtility workspace/shared/` | 共享 MMA 库 |
| `mv verify/FamilyDatabase workspace/shared/` | 共享族数据库 |
| `mv verify/docs workspace/shared/verify-docs/` | 共享文档 |
| `mv verify/logs workspace/shared/logs/` | 共享日志 |
| `mv verify/timing workspace/shared/timing/` | 共享计时数据 |
| 重复的 .bin 文件替换为指向 data/ 的符号链接 | 节省空间 |
| 更新 task JSON/MD 中的路径引用 | 共 16 个 JSON, 14 个 MD |

### 未执行（待定）

- 源代码目录重组（`include/` `src/` `tools/` 等移入 `src/`）
- `build/` 中 `generating_cone.cpp` 源文件移出

### workspace/shared/mma-refs/（共享 MMA 参考数据，原 verify/）

```
workspace/shared/mma-refs/Box/  bub00/  bub10/  bub11/  Tri/
workspace/shared/mma-refs/DP323/  NP222/  TB123m/
workspace/shared/mma-refs/Penta1L/  SR212/  SR212-3m/
```

这些是各家族的 MMA 参考对比数据，由 T 系列任务和验证流程共用。
