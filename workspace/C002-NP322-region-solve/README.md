# NP322 Scheduled Region Solve

NP322 (Non-Planar 2L4P, 9 propagators, massless, s=1, t=5) 族 Region 计算。使用 `ScheduledRegionSolver` 进行三级超时调度。

> **Note**: 本目录遵循 `verify/DP323_recompute/` 的目录规范。当前程序正在运行中，运行结束后应整理为与 DP323_recompute 一致的结构。

## Directory Structure

```
verify/NP322/
├── run_scheduled.wl        ← ScheduledRegionSolver 调用入口
├── README.md               ← 本文档
├── cache/                  ← 计算缓存与状态
│   ├── status.wdx          ← 扇区级状态跟踪 (WDX)
│   ├── instance.lock       ← 实例锁，防重复启动
│   ├── worker_data.wdx     ← Worker 调度数据
│   ├── SectorCache/        ← 各扇区 Region 结果缓存
│   └── worker_*.log        ← Worker 日志
├── output/                 ← 输出文件目录（规划中）
├── archive/                ← 历史脚本存档（规划中）
│   ├── logs/               ← 历史运行日志
│   └── scripts/            ← 旧版测试脚本
├── IBPMat_NP322.bin        ← IBP 矩阵 (110MB, 2026-05-10)
├── RingData_NP322.bin      ← Ring 数据 (1.3MB, 2026-05-10)
├── scheduled_run.log       ← 当前运行日志 (active)
└── run.pid                 ← 当前进程 PID
```

## Current Status

| 项目 | 状态 |
|------|------|
| 进程 | 运行中 (PID `run.pid`) |
| IBPMat_NP322.bin | 已生成 (110MB @ project root: 2026-05-05, verify/NP322/: 2026-05-10) |
| RingData_NP322.bin | 已生成 (1.3MB @ verify/NP322/, 16B placeholder @ project root) |
| Sector 进度 | 90/120 cached, 10 running, 20 trivial |

## Usage

### Run

```bash
cd <project_root>
wolframscript -file verify/NP322/run_scheduled.wl
```

### Resume after interruption

`ScheduledRegionSolver` 自动将状态存入 `cache/status.wdx`。可安全 `kill` 后重启，调度器从上一检查点恢复。

### Instance lock

`cache/instance.lock` 防止多个实例并发。检测到过期锁 (dead PID) 时自动移除。

## Architecture

```
run_scheduled.wl (wrapper)
  ↓
ScheduledRegionSolve["NP322", familyConfig, scheduleConfig]
  ↓
┌──────────────────────────────────────────────┐
│  InstanceLock          — 防并发               │
│  SectorStateManager    — 原子状态持久化        │
│  TieredScheduler       — 三级超时调度           │
│  RegionComputeEngine   — Singular + 超时传递    │
└──────────────────────────────────────────────┘
```

## Schedule Tiers

| Tier | Timeout | MaxParallel | Order  | 说明 |
|------|---------|-------------|--------|------|
| Quick | 120s | 4 | First | 快速扇区初筛 |
| Fast | 1200s | 3 | First | 中等复杂度扇区 |
| Retry | 12000s | 3 | FIFO | 困难扇区重试 |

## Sector Distribution

| Props | Sectors | Status |
|-------|---------|--------|
| 1-prop | 7 | trivial (已完成) |
| 2-prop | 13 | trivial (已完成) |
| 3-prop | 9 | ✅ 全部 done |
| 4-prop | 56 | ✅ 全部 done |
| 5-prop | 27 | 25 done + 2 running |
| 6-prop | 7 | ⏳ running |
| 7-prop | 1 | ⏳ running |

Total: 120 sectors (90 done, 10 running, 20 trivial)

## Dependencies

- Mathematica (WolframScript)
- `LIEWorkflow.wl` (loads `ScheduledRegionSolver.wl`)
- `FamilyDatabase.wl` (NP322 entry)
- Singular 4.3.2

## Post-Run Cleanup

运行完成后应执行：

1. 停止后移出 `run.pid`
2. 将 `verify/NP322/IBPMat_NP322.bin` 和 `RingData_NP322.bin` 移入 `output/`
3. 将 `scheduled_run.log` 移入 `archive/logs/`
4. 将 `SingularTempFile/` 移入 `archive/temp/`
5. 将旧版测试脚本移入 `archive/scripts/`
6. 将 `IBPMat_NP322.bin` 复制到项目根目录供 C++ test_relationFF 使用
7. 检查 project root 的 `RingData_NP322.bin` 是否为完整版（当前为 16B 占位符）
