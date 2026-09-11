# 随仓库分发的公开单细胞资源

本目录只包含运行离线单细胞 ORA/GSEA 与 `functional_state` 所需的公开、版本固定资源，
不包含任何项目数据、患者信息、CellTypist 模型或用户上传文件。资源文件和注册表均为只读；
运行时会校验所选 GMT/TF 网络的 SHA-256，网页不能指定下载地址或任意服务器路径。

| 资源 | 随库版本 / 文件 | 许可与署名要求 | 官方来源 |
| --- | --- | --- | --- |
| MSigDB Hallmark Human | `2026.1.Hs`; `gene_sets/hallmark/h.all.v2026.1.Hs.symbols.gmt` | [CC BY 4.0](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp)；保留 MSigDB 的版权与每个 gene set 的附加条款 | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) |
| Reactome Pathways | `2026-06-21`; `gene_sets/reactome/reactome_human_current.gmt` | Reactome 数据为 [CC0 1.0](https://reactome.org/license)；建议署名 Reactome | [Reactome download](https://reactome.org/download-data) |
| Gene Ontology + HUMAN-uniprot GAF | ontology release `2026-07-26`; `gene_sets/go/` | [CC BY 4.0](https://geneontology.org/docs/go-citation-policy/)；公开使用或再分发时须保留 GO 版本、来源及署名 | [GO current releases](https://current.geneontology.org/) |
| CollecTRI Human regulons | v2.0; `functional_state_resources/collectri_human.tsv` | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)；须署名 Müller-Dott *et al.*，并链接 DOI | [Zenodo 10.5281/zenodo.8192729](https://doi.org/10.5281/zenodo.8192729) |

使用 CollecTRI 时请同时引用其数据集与相应论文；使用 GO 时请在报告中写明上述 release。`gene_set_registry.json`
与 `collectri_human.metadata.json` 保存每次随库快照的来源、版本、校验值和生成信息。更新资源必须
通过管理员维护流程完成，并先复核对应许可证；不要把 KEGG、WikiPathways 或任何来源不明的本地
`genesets/` 文件直接加入仓库。
