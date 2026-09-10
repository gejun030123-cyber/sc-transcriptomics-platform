# 托管基因集注册表

`functional_state` 的标准通路评分使用受管理员控制的本地快照，不在分析任务中联网，也不接受网页传入的文件路径。资源默认位于：

```text
data/functional_state_resources/gene_sets/
├── gene_set_registry.json
├── hallmark/
├── reactome/
├── go/
└── sources/
```

`gene_set_registry.json` 保存每个 GMT 的来源 URL、版本、许可证、SHA-256 和 term 数；`sources/` 保留 Reactome ZIP、GO OBO 与 Human GAF 的原始校验对象。分析开始时会重新计算所选 GMT 的 SHA-256，校验失败即停止，不会悄悄使用被替换的基因集。

## 管理员同步

在部署环境中执行：

```bash
python scripts/sync_managed_gene_sets.py --resource-dir "$FUNCTIONAL_STATE_RESOURCE_DIR"
```

未设置 `FUNCTIONAL_STATE_RESOURCE_DIR` 时，使用 `data/functional_state_resources`。该命令只访问固定的官方来源并写入管理员资源目录：

- MSigDB Hallmark Human 2026.1.Hs；
- Reactome `ReactomePathways.gmt`；
- Gene Ontology `go-basic.obo` 与 `HUMAN-uniprot.gaf.gz`，生成 direct-annotation 的 Human BP/MF/CC GMT。

GO 转换不进行祖先 term 传播；该策略与 GO release、GAF 生成日期共同写入注册表。同步前应确认 MSigDB 的适用许可；Reactome 数据为 CC0，GO 数据为 CC BY 4.0，许可证链接随运行 manifest 输出。

## 分析时的 QC

平台只加载本次明确选择的 term，最多 60 个，并在 `gene_set_coverage.csv` 中记录原始基因数、实际检测基因数、缺失基因与覆盖度：

- ≥70%：正常评分；
- 40–70%：保留评分并标记警告；
- <40% 的标准通路：不计算评分；
- 托管标准通路另要求至少检测到 10 个基因。

内置 FAO/炎症小面板仍保留为透明的研究假设签名：它们与 Hallmark/Reactome 标准通路分开标识，允许在至少 2 个检测基因时探索性评分，不能替代标准库的结果。

选择 `IBD 类器官上皮主线`预设时，平台从同一冻结 Hallmark 快照读取有限的炎症、应激、代谢、WNT 与细胞周期 term；不会在任务中联网或下载基因集。若管理员尚未执行同步，预设会明确停止并提示初始化资源，而不会降级到未经版本控制的在线库。
