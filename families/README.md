# FamilyDatabase

IBP 积分族的集中定义数据库。所有验证脚本和测试应从此处获取族定义，确保参数一致性。

## 文件

```
families/
├── README.md              # 本文件
├── FamilyDatabase.wl      # MMA 可加载的族定义数据库 (752 行, 25 families)
└── *.json                 # C++ 族定义 (25 files, 与 WL 一致)
```

C++ JSON 与 MMA WL 中的族名、传播子、运动学规则完全一致。

## 用法 (MMA)

```mathematica
Get["families/FamilyDatabase.wl"];
fam = "bub00";
config = $FamilyDatabase[fam];
PrintFamilyInfo[fam];
PrintAllFamilies[];
```

## 已注册的积分族 (25)

### 1-Loop (L=1)

| 名称 | 描述 | E | Props |
|------|------|:---:|:---:|
| `bub00` | 1-loop bubble, msq=0 | 1 | 2 |
| `bub10` | 1-loop bubble, 1-massive | 1 | 2 |
| `bub11` | 1-loop bubble, both-massive | 1 | 2 |
| `Tri` | 1-loop triangle | 2 | 3 |
| `Box` | 1-loop box | 3 | 4 |
| `Penta1L` | 1-loop pentagon 5pt | 4 | 5 |

### 2-Loop (L=2)

| 名称 | 描述 | E | Props |
|------|------|:---:|:---:|
| `SR212` | Sunrise 2L1P, massless | 1 | 5 |
| `SR212-3m` | Sunrise 2L1P, 3-massive | 1 | 5 |
| `SR212-5m` | Sunrise 2L1P, all-massive | 1 | 5 |
| `NP222` | Non-Planar 2L2P | 2 | 7 |
| `NP222var` | Non-Planar 2L2P variant | 2 | 7 |
| `TB123` | Tri-Box 2L3P, massless | 2 | 7 |
| `TB123m` | Tri-Box 2L3P, massive | 2 | 7 |
| `DB313` | Double Box 2L4P | 3 | 9 |
| `NP322` | Non-Planar 2L4P, massless | 3 | 9 |
| `NP322m` | Non-Planar 2L4P, massive | 3 | 9 |
| `DP323` | Double Pentagon 2L5P | 4 | 11 |
| `DP2L5P` | Non-Planar DP 2L5P | 4 | 11 |
| `HB2L5P` | Non-Planar HexaBox 2L5P | 4 | 11 |
| `PB2L5P` | Planar PentaBox 2L5P | 4 | 11 |

### 3-Loop (L=3)

| 名称 | 描述 | E | Props |
|------|------|:---:|:---:|
| `BN3L` | Banana 3L0P, vacuum | 0 | 6 |
| `BN3L1P` | Banana 3L1P, massless | 1 | 9 |
| `BN3L1Pm` | Banana 3L1P, massive | 1 | 9 |
| `LA3L4Pm` | 3L Ladder A 4pt, 2-massive | 3 | 15 |
| `LB3L4Pm` | 3L Ladder B 4pt, 2-massive | 3 | 15 |

## 配置条目

| 键 | 说明 |
|------|------|
| `"Propagators"` | 传播子列表（符号形式） |
| `"LoopMomenta"` | 圈动量列表 |
| `"ExternalMomenta"` | 独立外动量列表 |
| `"KinematicRules"` | 标量积替换规则 |
| `"TopSector"` | 顶角扇区 (1=主传播子, 0=辅助 ISP) |
| `"Numeric"` | 数值替换规则 (含 `"d"` = 维度参数) |
| `"Modulus"` | 有限域特征 p (默认 Prime[10000000]) |

## 添加新积分族

1. 在 `FamilyDatabase.wl` 的 `$FamilyDatabase` 中添加 Association 条目
2. 在 `families/` 中创建对应的 `<name>.json`（字段语义一致）
3. 更新本 README 的表格
4. 运行 `diff` 验证两侧族名集合一致

## 文献参考

部分族定义来自：
- arXiv:2009.07803 — 2L 5-point 拓扑 (PB2L5P, HB2L5P, DP2L5P)
- arXiv:2410.15431 — 3L 梯子族 (LA3L4Pm, LB3L4Pm)
- TP412&NP322_Eigenvalue.nb — 3L Banana 族 (BN3L, BN3L1P, BN3L1Pm)
