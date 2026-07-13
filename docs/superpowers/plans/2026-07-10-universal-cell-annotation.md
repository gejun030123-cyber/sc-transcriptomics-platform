# 通用细胞注释补充方案

## 原则

不把单一 marker 表冒充为通用注释器。平台采用“先广后细、证据可见、低置信度不强判”的机制：未知组织先做大谱系初注释；已知样本背景再用 TME、免疫、血液或 PBMC marker 集，或用户自定义 marker，做细分复核。

## 已实现

1. 新增默认 `Universal` marker 集：上皮、内皮、成纤维、周细胞/平滑肌、髓系、T、NK、B、浆细胞、肥大和增殖细胞等大谱系。
2. 人/鼠符号大小写不一致时按不区分大小写匹配 marker。
3. 按当前数据的 marker 覆盖度过滤候选类型；命中 marker 不足的类型不参与打分。
4. 输出 marker 覆盖度图，并将逐类型覆盖统计写入任务摘要。
5. 可用最低 marker score 将缺乏证据的细胞标为 `Unknown`；现有 entropy/score-margin 置信度仍可叠加使用。

## 推荐操作

未知来源数据：`Universal` -> 查看 marker 覆盖度、score 热图和 UMAP -> 对可疑大谱系运行子簇精细分析 -> 使用对应场景 marker 集或自定义 marker 精注释。CellTypist 只作为具备匹配人类参考模型时的交叉证据，不替代人工审阅。

## 人类多证据模式（已开始实施）

`multi_evidence` 将规则标签作为最终可解释候选，记录 marker Top 1/Top 2、得分、cluster 多数一致率和可用的 CellTypist 标签；不部署 R，不接入 SingleR/Azimuth。Universal 标签同时写入 `cell_type_l1`、`cell_type_l2`、`cell_type_l3`、`cell_ontology_id`、`final_annotation` 和 `annotation_status`，使后续人工审核可以在 h5ad 中追溯自动判断。
