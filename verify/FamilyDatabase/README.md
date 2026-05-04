# FamilyDatabase

IBP 积分族的集中定义数据库。所有验证脚本和测试应从此处获取族定义，确保参数一致性。

## 文件

```
verify/FamilyDatabase/
├── README.md              # 本文件
└── FamilyDatabase.wl      # MMA 可加载的族定义数据库
```

## 用法

```mathematica
(* 加载数据库 *)
Get["verify/FamilyDatabase/FamilyDatabase.wl"];

(* 获取积分族配置 *)
fam = "bub00";
config = $FamilyDatabase[fam];

(* 查看配置 *)
PrintFamilyInfo[fam];

(* 列出所有注册的族 *)
PrintAllFamilies[];

(* 在 LIEDefineFamily 中使用 *)
data = LIEDefineFamily[
    config["Propagators"], config["LoopMomenta"],
    config["ExternalMomenta"], config["KinematicRules"],
    config["TopSector"],
    "Numeric" -> config["Numeric"],
    Modulus -> config["Modulus"]
];
```

## 已注册的积分族 (19)

### 1-Loop (L=1)

| 名称 | 描述 | E | Props | 质量 |
|------|------|:---:|:---:|------|
| `bub00` | 1-loop bubble, msq=0 | 1 | 2 | 无 |
| `bub10` | 1-loop bubble, 1-massive | 1 | 2 | msq=1 |
| `bub11` | 1-loop bubble, both-massive | 1 | 2 | msq=1 |
| `Tri` | 1-loop triangle | 2 | 3 | msq=0 |
| `Box` | 1-loop box | 3 | 4 | msq=0 |

### 2-Loop (L=2)

| 名称 | 描述 | E | Props | 质量 |
|------|------|:---:|:---:|------|
| `SR212` | Sunrise 2L1P, massless | 1 | 5 | 无 |
| `SR212-3m` | Sunrise 2L1P, 3-massive | 1 | 5 | 3/5 massive |
| `SR212-5m` | Sunrise 2L1P, all-massive | 1 | 5 | 全 massive |
| `NP222` | Non-Planar 2L2P | 2 | 7 | 参数化 (m=0) |
| `TB123` | Tri-Box 2L3P, massless | 2 | 7 (6+1) | 无 |
| `TB123m` | Tri-Box 2L3P, massive | 2 | 7 (6+1) | msq=1 |
| `DB313` | Double Box 2L4P | 3 | 9 (7+2) | 无 |
| `NP322` | Non-Planar 2L4P, massless | 3 | 9 | 无 |
| `NP322m` | Non-Planar 2L4P, massive | 3 | 9 | msq=1 |
| `DP323` | Double Pentagon 2L5P | 4 | 11 (8+3) | msq=3 |

### 3-Loop (L=3, 新提取)

| 名称 | 描述 | E | Props | 质量 |
|------|------|:---:|:---:|------|
| `BN3L` | Banana 3L0P, vacuum | 0 | 6 | msq=1 |
| `BN3L1P` | Banana 3L1P, massless | 1 | 9 (8+1) | 无 |
| `BN3L1Pm` | Banana 3L1P, massive | 1 | 9 (8+1) | msq=1 |

## 配置条目说明

每个族的 Association 包含以下键：

| 键 | 类型 | 说明 |
|------|------|------|
| `"Description"` | String | 人类可读描述 |
| `"Propagators"` | List[Expr] | 传播子列表（符号形式，numeric 未代入） |
| `"LoopMomenta"` | List[Symbol] | 圈动量列表 |
| `"ExternalMomenta"` | List[Symbol] | 独立外动量列表 |
| `"KinematicRules"` | List[Rule] | 标量积替换规则（符号形式） |
| `"TopSector"` | List[0\|1] | 顶角扇区（1=主传播子, 0=辅助） |
| `"Numeric"` | List[Rule] | 数值替换规则（含 "d" 维度参数） |
| `"Modulus"` | Integer | 有限域特征 p |

## 添加新积分族

1. 在 `FamilyDatabase.wl` 的 `$FamilyDatabase` 中添加条目
2. 在 MMA-Mini 中创建对应的 `Compare-FamilyGenerate-<name>.wl`（引用数据库）
3. 运行生成 `.bin` 文件
4. 更新本 README 的积分族表格

## C++ 侧的对应关系

C++ 测试通过 CLI 参数 `famname` 选择积分族，加载 `IBPMat_<famname>.bin` 和 `RingData_<famname>.bin`。族名必须与此数据库中的键名一致。
