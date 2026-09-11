# 托管基因集注册表

`functional_state` 的标准通路评分使用受管理员控制的本地快照，不在分析任务中联网，也不接受网页传入的文件路径。干净克隆默认使用随仓库提供的只读基线：

```text
resources/functional_state_resources/gene_sets/
├── gene_set_registry.json
├── hallmark/
├── reactome/
├── wikipathways/
└── go/
```

`gene_set_registry.json` 保存每个 GMT 的来源 URL、版本、许可证、SHA-256 和 term 数。随仓库分发的是
运行所需 GMT 与注册表；管理员同步过程产生的原始 Reactome ZIP、GO OBO 和 Human GAF 校验对象只可保留在受控
管理员目录，不能提交到仓库。分析开始时会重新计算所选 GMT 的 SHA-256，校验失败即停止，不会悄悄使用被替换的基因集。

## 管理员同步

随库基线已覆盖 Hallmark、Reactome、WikiPathways Human、GO 与 CollecTRI v2.0，不必联网同步即可运行默认功能。
其中 WikiPathways 使用 Enrichr 发布的 `WikiPathways_2024_Human` gene-symbol 快照，保证可与常见 Human
symbol DEG 直接匹配；其上游 WikiPathways 内容为 CC0，仍应在结果或方法中引用 WikiPathways 和 Enrichr。若管理员需要替换或更新 Hallmark、Reactome 或 GO 通路快照，先将资源根设置为独立的受控目录，再执行：

```bash
export FUNCTIONAL_STATE_RESOURCE_DIR=/srv/sc-platform/references/functional_state
python scripts/sync_managed_gene_sets.py --resource-dir "$FUNCTIONAL_STATE_RESOURCE_DIR"
```

未设置 `FUNCTIONAL_STATE_RESOURCE_DIR` 时，平台会继续使用随库基线；该命令只访问固定的官方来源并写入管理员资源目录：

- MSigDB Hallmark Human 2026.1.Hs；
- Reactome `ReactomePathways.gmt`；
- Gene Ontology `go-basic.obo` 与 `HUMAN-uniprot.gaf.gz`，生成 direct-annotation 的 Human BP/MF/CC GMT。

GO 转换不进行祖先 term 传播；该策略与 GO release、GAF 生成日期共同写入注册表。同步前应确认 MSigDB 的适用许可；Reactome 数据为 CC0，GO 数据为 CC BY 4.0，许可证链接随运行 manifest 输出。替换 TF 网络时，还须同时放入 `collectri_human.tsv` 与同名 `.metadata.json`，并核对 SHA-256、来源和许可。随库的第三方资源说明见 [`resources/THIRD_PARTY_DATA_NOTICES.md`](../resources/THIRD_PARTY_DATA_NOTICES.md)。

## 分析时的 QC

平台只加载本次明确选择的 term，最多 60 个，并在 `gene_set_coverage.csv` 中记录原始基因数、实际检测基因数、缺失基因与覆盖度：

- ≥70%：正常评分；
- 40–70%：保留评分并标记警告；
- <40% 的标准通路：不计算评分；
- 托管标准通路另要求至少检测到 10 个基因。

内置 FAO/炎症小面板仍保留为透明的研究假设签名：它们与 Hallmark/Reactome 标准通路分开标识，允许在至少 2 个检测基因时探索性评分，不能替代标准库的结果。

选择 `IBD 类器官上皮主线`预设时，平台从同一冻结 Hallmark 快照读取有限的炎症、应激、代谢、WNT 与细胞周期 term；不会在任务中联网或下载基因集。若管理员尚未执行同步，预设会明确停止并提示初始化资源，而不会降级到未经版本控制的在线库。
