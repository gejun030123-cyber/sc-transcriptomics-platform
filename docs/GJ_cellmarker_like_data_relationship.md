# GJ 单细胞数据的 CellMarker-like 网站数据关系

## 1. 文档范围

本文档只整理当前 `/home/oelab/AnaData/GJ/01_analysis_results` 中已有的 CSV
结果以及样本名称中明确出现的信息。

- 不考虑 treatment。
- 不根据样本名称推断实验组、对照组、基因型或处理方式。
- 不补充外部 CellMarker、文献、Cell Ontology 或 Uberon 数据。
- 不进行细胞类型注释。
- `condition` 当前等于 `sample_id`，网站中不将其解释为处理分组。
- 无法根据名称确定类器官类型或培养天数时，统一记录为“未识别”。

因此，这里所说的 CellMarker-like，指采用类似的“项目—样本—组织/类器官—
分群—基因”浏览和检索方式，不表示当前数据中的基因已经是经过验证的
CellMarker。

## 2. 当前数据概况

| 项目 | 样本数 | QC 后细胞数 | Leiden cluster 数 | 导出基因数 |
|---|---:|---:|---:|---:|
| CXY_megtacolon | 20 | 124,614 | 33 | 26,007 |
| WYQ_VLO | 10 | 199,332 | 30 | 26,593 |
| LHS | 13 | 219,512 | 33 | 27,010 |
| WYQ_rett | 10 | 248,347 | 34 | 26,249 |
| **合计** | **53** | **791,805** | **130 个项目内 cluster** | — |

Cluster 编号只在所属项目内有效。例如，`CXY_megtacolon` 的 cluster 0 与
`LHS` 的 cluster 0 不是同一个 cluster。网站内部应使用组合编号：

```text
CXY_megtacolon__cluster_0
WYQ_VLO__cluster_0
LHS__cluster_0
WYQ_rett__cluster_0
```

## 3. 样本名称解析规则

解析时不区分大小写，但保留原始 `sample_id` 不变。

### 3.1 类器官类型

| 名称编码 | 网站标准名称 | 示例 |
|---|---|---|
| `IO`，包括 `HIO`、`VIO` | 肠类器官 | `CTHIOD20`、`VIO-D20` |
| `LO`，包括 `VLO`、`VasLO`、`eLO` | 肝类器官 | `H9-VLO`、`eVasLO-D30` |
| `BO`，包括 `VBO` | 脑类器官 | `H9-VBO-E`、`RTT-BO-d60` |
| `HO` | 心脏类器官 | `IPSC-HOd20`、`RTT-HO-day-30` |
| 不包含以上编码 | 未识别 | `AD-E`、`RTT-D30` |

网站建议保存两个字段：

```text
sample_id             # 原始样本名，不修改
organoid_type         # 肠类器官/肝类器官/脑类器官/心脏类器官/未识别
```

### 3.2 培养天数

只解析明确的 `D数字`、`d数字` 或 `day-数字`：

| 名称形式 | day 值 | 示例 |
|---|---:|---|
| `D7`、`d7` | 7 | `LO-D7` |
| `D20`、`d20` | 20 | `CTLOD20`、`IPSC-HOd20` |
| `D30`、`day-30` | 30 | `eLO-D30`、`RTT-HO-day-30` |
| `D60`、`d60` | 60 | `RTT-BO-ctrl-d60` |
| 没有明确 D/day 编码 | 未识别 | `LO0301`、`Dilation0507` |

`0301`、`0507`、`0716`、`250507`、`251015` 等没有明确 D/day 前缀的数字
不解释为培养天数。

## 4. 类器官与天数汇总

### 4.1 类器官类型

| 类器官类型 | 样本数 |
|---|---:|
| 肠类器官 | 6 |
| 肝类器官 | 14 |
| 脑类器官 | 4 |
| 心脏类器官 | 4 |
| 未识别 | 25 |
| **合计** | **53** |

### 4.2 培养天数

| day | 样本数 |
|---:|---:|
| 7 | 2 |
| 20 | 10 |
| 30 | 10 |
| 60 | 2 |
| 未识别 | 29 |
| **合计** | **53** |

### 4.3 类器官—天数关系

| 类器官类型 | D7 | D20 | D30 | D60 | 天数未识别 | 合计 |
|---|---:|---:|---:|---:|---:|---:|
| 肠类器官 | 0 | 4 | 2 | 0 | 0 | 6 |
| 肝类器官 | 2 | 4 | 4 | 0 | 4 | 14 |
| 脑类器官 | 0 | 0 | 0 | 2 | 2 | 4 |
| 心脏类器官 | 0 | 2 | 2 | 0 | 0 | 4 |
| 未识别 | 0 | 0 | 2 | 0 | 23 | 25 |
| **合计** | **2** | **10** | **10** | **2** | **29** | **53** |

## 5. 样本整理结果

### 5.1 CXY_megtacolon

该项目的样本名称未出现 IO、LO、BO 或 HO 编码，也没有明确的 D/day 天数编码。

| sample_id | 类器官类型 | day |
|---|---|---:|
| Dilation | 未识别 | 未识别 |
| Dilation0507 | 未识别 | 未识别 |
| Dilation0716 | 未识别 | 未识别 |
| IBD | 未识别 | 未识别 |
| IN | 未识别 | 未识别 |
| IP3 | 未识别 | 未识别 |
| NEC_251015 | 未识别 | 未识别 |
| NORMAL3 | 未识别 | 未识别 |
| Normal1 | 未识别 | 未识别 |
| Normal2 | 未识别 | 未识别 |
| Normal3 | 未识别 | 未识别 |
| Stenosis | 未识别 | 未识别 |
| Stenosis0507 | 未识别 | 未识别 |
| Stenosis0716 | 未识别 | 未识别 |
| Transition | 未识别 | 未识别 |
| X_250507 | 未识别 | 未识别 |
| X_251015 | 未识别 | 未识别 |
| k_250507 | 未识别 | 未识别 |
| k_251015 | 未识别 | 未识别 |
| k_25618 | 未识别 | 未识别 |

### 5.2 WYQ_VLO

| sample_id | 类器官类型 | day | 识别依据 |
|---|---|---:|---|
| CTHIOD20-0203 | 肠类器官 | 20 | HIO、D20 |
| LO-D7 | 肝类器官 | 7 | LO、D7 |
| LO-N-D7 | 肝类器官 | 7 | LO、D7 |
| LO0301 | 肝类器官 | 未识别 | LO；0301不解释为day |
| VIO-D20 | 肠类器官 | 20 | VIO、D20 |
| eLO-D30 | 肝类器官 | 30 | LO、D30 |
| eVasLO-D30 | 肝类器官 | 30 | VasLO、D30 |
| lo | 肝类器官 | 未识别 | LO |
| vLO0302 | 肝类器官 | 未识别 | VLO；0302不解释为day |
| vaslo | 肝类器官 | 未识别 | VasLO |

### 5.3 LHS

| sample_id | 类器官类型 | day | 识别依据 |
|---|---|---:|---|
| AD-A-K | 未识别 | 未识别 | 无类器官和day编码 |
| AD-E | 未识别 | 未识别 | 无类器官和day编码 |
| AD-K | 未识别 | 未识别 | 无类器官和day编码 |
| H9-VBO-E | 脑类器官 | 未识别 | VBO |
| H9-VBO-K | 脑类器官 | 未识别 | VBO |
| IPSC-HOd20 | 心脏类器官 | 20 | HO、d20 |
| RTT-BO-ctrl-d60 | 脑类器官 | 60 | BO、d60；不使用ctrl作为分组 |
| RTT-BO-d60 | 脑类器官 | 60 | BO、d60 |
| RTT-CTRL-D30 | 未识别 | 30 | D30；不使用CTRL作为分组 |
| RTT-D30 | 未识别 | 30 | D30 |
| RTT-HO-day-30 | 心脏类器官 | 30 | HO、day-30 |
| RTT-HOd20 | 心脏类器官 | 20 | HO、d20 |
| iPSC-HO-day-30 | 心脏类器官 | 30 | HO、day-30 |

### 5.4 WYQ_rett

| sample_id | 类器官类型 | day | 识别依据 |
|---|---|---:|---|
| CT-HIO-D30 | 肠类器官 | 30 | HIO、D30 |
| CT-LO-D30 | 肝类器官 | 30 | LO、D30 |
| CTHIOD20 | 肠类器官 | 20 | HIO、D20 |
| CTLO-D20 | 肝类器官 | 20 | LO、D20 |
| CTLOD20 | 肝类器官 | 20 | LO、D20 |
| REHHIOD20 | 肠类器官 | 20 | HIO、D20 |
| REHLOD20 | 肝类器官 | 20 | LO、D20 |
| RETT-HIO-D30 | 肠类器官 | 30 | HIO、D30 |
| RETT-LO-D30 | 肝类器官 | 30 | LO、D30 |
| RETTLO-D20 | 肝类器官 | 20 | LO、D20 |

## 6. 当前 CSV 的数据关系

| 实体 | 唯一键 | 现有字段/来源 | 关系 |
|---|---|---|---|
| Project | `project_id` | `experiment_id` | 一个项目包含多个样本和多个cluster |
| Sample | `sample_id` | sample metadata/inventory | 一个样本属于一个项目并包含多个细胞 |
| SampleOrganoid | `sample_id` | 从样本名按本文规则解析 | 一个样本对应一个类器官类型和一个day值 |
| Cell | `cell_id` | `sc_batch_cell_metadata.csv` | 一个细胞属于一个样本和一个cluster |
| Cluster | `project_id + leiden` | `leiden` | 一个cluster属于一个项目并包含多个细胞 |
| SampleCluster | `sample_id + cluster_id` | sample cluster proportions | 连接样本与cluster，保存细胞数和比例 |
| Gene | `GeneID` | pseudobulk CSV | 一个基因可出现在多个项目和样本中 |
| SampleExpression | `sample_id + GeneID` | counts/log2CPM | 连接样本与基因表达量 |

关系结构：

```text
Project
├── Sample
│   ├── SampleOrganoid（organoid_type、day）
│   ├── Cell
│   │   └── Cluster
│   └── SampleExpression ── Gene
├── Cluster
└── SampleCluster ───────── Sample + Cluster
```

## 7. 与 CellMarker-like 页面概念的对应关系

这里只建立页面概念对应，不引入外部 CellMarker 数据。

| 当前数据 | 网站中的 CellMarker-like 概念 | 当前可解释范围 |
|---|---|---|
| `project_id` | Dataset/Project | 可直接展示 |
| `organoid_type` | Tissue/Organoid 分类 | 只按样本名编码展示 |
| `day` | Development day | 只展示明确 D/day 编码 |
| `sample_id` | Sample | 可直接展示 |
| `leiden` | Cell cluster | 只能称为分群，不能称为细胞类型 |
| `cell_id`、UMAP | Single-cell view | 可绘制UMAP并按样本/cluster着色 |
| `GeneID/GeneName` | Gene | 可进行基因检索 |
| counts/log2CPM | Sample gene expression | 可展示样本层面表达量 |
| cluster proportion | Cell composition | 可展示样本内cluster组成 |

当前不能建立以下关系：

```text
Cluster -> Cell type
Cluster -> Validated marker
Gene -> Literature evidence
Gene -> Standard CellMarker entry
```

原因是当前导出结果没有细胞类型注释、cluster差异marker或文献证据。这些关系
在网站中应显示为“暂无数据”，而不是根据名称或表达量自动补充。

## 8. 建议的网站浏览层级

在不增加新数据的前提下，网站可以提供以下页面：

1. **项目页**：展示项目、样本数、细胞数、cluster数和可识别的类器官类型。
2. **类器官页**：按肠、肝、脑、心脏和未识别分类浏览样本。
3. **天数页**：按 D7、D20、D30、D60 和未识别浏览样本。
4. **样本页**：展示所属项目、类器官类型、day、QC后细胞数及cluster比例。
5. **Cluster页**：展示项目内cluster编号、细胞数、样本构成及UMAP位置。
6. **基因页**：展示该基因在各项目和样本中的counts及log2CPM。
7. **UMAP页**：按project、sample、organoid type、day或cluster着色。

建议网站筛选顺序：

```text
Project -> Organoid type -> Day -> Sample -> Cluster -> Gene
```

## 9. 网站导入时建议增加的派生字段

不修改原始 CSV，只在导入数据库时增加：

```text
project_id
sample_id
organoid_type
day
cluster_id
cluster_display_name
```

其中：

```text
cluster_id = project_id + "__cluster_" + leiden
cluster_display_name = "Cluster " + leiden
```

`organoid_type` 与 `day` 必须保留“未识别”状态，禁止把缺失信息自动补成某个
类器官或时间点。
