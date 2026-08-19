# Bulk RNA secondary figure audit

本次迁移覆盖标准化、QC、辅助 PCA、DEG 箱线图/整合比较和 Time-course。所有新图
使用 `FigureSpec + NatureStyle + FigureValidator`，默认输出 SVG、PDF、PNG，并在
最终 89/183 mm 画布上检查字体、标签重叠、边界、背景和 DPI。

## 已修复

- 标准化分布图不再给每个样本绘制拥挤的长标签；最终尺寸最多显示 8 个代表性标签。
- QC 总览会自动省略恒定指标（例如连续表达输入下全为 0 的 MT%），不会留下空面板。
- QC PCA、PCA PC1/PC3、PC2/PC3 在颜色组过多时改用单色 + marker，避免彩虹图例和
  红绿冲突。
- QC pairs/violin 会省略恒定指标、压缩组标签，并统一使用可辨识的 marker/色板。
- PCA 方差/载荷图、DEG 箱线图、DEG Jaccard/方向/log2FC/相关性矩阵和筛选器图已
  不再调用 legacy `save_matplotlib_figure`。
- Time-course Q-Q、cluster trajectory、gene×time 和 pairwise heatmap 已接入固定模板；
  热图展示行数受最终字号约束，完整结果仍保留在 CSV。
- PDF/TIFF 现被结果文件合同、结果页和富集结果 API 识别；SVG 保留可编辑文字。

## 真实项目回归

使用仓库 36 样本项目 `bec48c50-bc9` 的上传表达矩阵/中间 H5AD，在独立运行目录生成
before/after 结果。代表性 readiness：

| 模块 | 代表性结果 | readiness |
|---|---|---:|
| 标准化 | `bulk_norm_boxplot_compare`、`bulk_norm_pca_compare` | 100 |
| QC | overview、pairs、elbow、violin、correlation | 96–100 |
| PCA | 主图、PC1/3、PC2/3、variance、loadings | 92–100 |
| DEG boxplot | `bulk_deg_box_g1` fixture | 100 |

所有报告均确认物理画布符合 89/183 mm，PNG 为 600 DPI，SVG/PDF 文件有效且没有
真实文字越界。92 分的 PCA 报告仅保留“12 个自动分组使用 marker”语义提示，不是
布局失败；96 分的 QC overview 仅提示已省略恒定 MT 指标。

回归输出目录（可直接查看 PNG、SVG、PDF 与 readiness JSON）：

- `data/runtime_tmp/bulk_norm_new_48djjy0b/`
- `data/runtime_tmp/bulk_qc_final_mhxxouji/`
- `data/runtime_tmp/bulk_pca_final_4n9lduim/`

完整测试：`84 passed, 2 warnings`。

## 尚未迁移

PPI、Forest plot、Deconvolution 等尚未有固定 Bulk 模板；这些图仍应在新增模板后
再进入 Nature Portfolio 门禁，不能把 legacy 图直接标记为投稿就绪。
