# 单细胞 & Bulk RNA-seq 生信分析平台

基于 Flask + Scanpy + OmicVerse 的 Web 端生物信息学分析平台，支持单细胞转录组和 Bulk RNA-seq 全流程分析。

## 功能特性

### 单细胞转录组分析（10 个模块）

| 模块 | 功能 | 核心工具 |
|------|------|----------|
| 质控 (QC) | MT/ribo/hb 过滤、Scrublet 双细胞检测、QC 散点图 | scanpy, omicverse |
| 预处理 | 标准化、高变异基因选择（支持批次校正） | scanpy |
| 降维分析 | PCA + UMAP/t-SNE，方差比例图 | scanpy, sklearn |
| 批次校正 | Harmony / ComBat / SysVI | omicverse, scvi |
| 聚类分析 | 多分辨率 Leiden 聚类、UMAP 比较、marker dotplot | scanpy |
| QC 重新评估 | 聚类后 doublet/MT 检测、低质量簇标记 | scanpy |
| 细胞注释 | TME/Immune marker 基因集、自定义 marker、dotplot 验证 | scanpy |
| 差异表达 | Wilcoxon/t-test/logreg、火山图、dotplot、基因 UMAP | scanpy |
| 轨迹分析 | Diffusion Map + DPT 拟时序分析 | scanpy |
| 比例分析 | 细胞比例差异（卡方检验）、堆叠柱状图/饼图 | scipy |

### Bulk RNA-seq 分析（7 个模块）

| 模块 | 功能 | 核心工具 |
|------|------|----------|
| 数据质控 | 文库大小、基因检测、离群值过滤、样本相关性热图 | scanpy |
| 数据标准化 | DESeq2 size factors / CPM / 分位数标准化 | scipy |
| 差异表达分析 | t-test / Mann-Whitney / DESeq2（基于 OmicVerse）、火山图、MA 图、基因箱线图 | omicverse |
| PCA / UMAP | 降维可视化、载荷图、肘部图 | scanpy, sklearn |
| 热图分析 | Top 差异基因热图、样本相关性热图、分组注释条 | scipy |
| 通路富集 | ORA / GSEA（GO/KEGG/WikiPathways/Reactome） | omicverse, gseapy |
| 时序分析 | 多时间点差异基因（spline F-test）、模糊 c-means 轨迹聚类 | patsy, statsmodels |

### 通用功能

- 参数提示系统：每个参数附带中文说明和建议范围
- 实时进度追踪：带时间戳的运行日志
- 分析流程引导：每步完成后自动推荐下一步
- 交互式图表：Plotly.js 支持缩放、平移、导出
- 多格式输入：支持 .h5ad, .csv, .txt, .xlsx, .xls, .mtx.gz
- 自动分组检测：从样本名中自动提取实验分组

## 技术栈

- **后端**: Flask, SQLite (WAL), ThreadPoolExecutor
- **分析**: scanpy, omicverse, scipy, statsmodels, patsy, scikit-learn
- **前端**: Bootstrap 5, Plotly.js, Jinja2
- **数据库**: SQLite（projects, analysis_tasks, result_files）

## 安装

```bash
# 克隆仓库
git clone https://github.com/gejun030123-cyber/sc-transcriptomics-platform.git
cd sc-transcriptomics-platform

# 安装依赖
pip install flask flask-cors scanpy omicverse plotly psutil patsy statsmodels gseapy pydeseq2 scikit-learn

# 启动服务
python app.py
```

服务将在 `http://localhost:5000` 启动。

## 使用流程

1. 创建项目 → 上传数据文件
2. 选择分析类型（单细胞 / Bulk RNA-seq）
3. 按流程顺序执行各分析模块
4. 每步完成后查看结果图表和数据表格
5. 下载 h5ad 中间文件和 CSV 结果

### 单细胞典型流程

```
质控 → 预处理 → 降维 → 聚类 → QC 重新评估 → 细胞注释 → 差异表达 → 轨迹分析 → 比例分析
```

### Bulk RNA-seq 典型流程

```
数据质控 → 标准化 → 差异表达 → PCA/UMAP → 热图 → 通路富集 / 时序分析
```

## 项目结构

```
├── app.py                 # Flask 应用工厂
├── config.py              # 全局配置
├── database.py            # SQLite 初始化
├── models.py              # ORM 模型
├── worker.py              # 异步任务执行器
├── routes/
│   ├── main.py            # 首页控制台
│   ├── projects.py        # 项目 CRUD
│   ├── upload.py          # 文件上传
│   ├── analysis.py        # 分析模块配置
│   ├── results.py         # 结果展示
│   └── api.py             # REST API
├── modules/
│   ├── base.py            # 分析基类
│   ├── io_utils.py        # 数据读取
│   ├── visualization.py   # Plotly 可视化工具
│   ├── qc.py              # 单细胞质控
│   ├── preprocess.py      # 预处理
│   ├── dimred.py          # 降维
│   ├── batch_correct.py   # 批次校正
│   ├── clustering.py      # 聚类
│   ├── qc_reassess.py     # QC 重新评估
│   ├── annotation.py      # 细胞注释
│   ├── deg.py             # 差异表达
│   ├── trajectory.py      # 轨迹分析
│   ├── proportion.py      # 比例分析
│   ├── bulk_qc.py         # Bulk 质控
│   ├── bulk_normalize.py  # Bulk 标准化
│   ├── bulk_deg.py        # Bulk DEG
│   ├── bulk_pca.py        # Bulk PCA
│   ├── bulk_heatmap.py    # Bulk 热图
│   ├── bulk_enrichment.py # 通路富集
│   └── bulk_timecourse.py # 时序分析
└── templates/
    ├── base.html
    ├── index.html
    ├── project_new.html
    ├── project_detail.html
    ├── upload.html
    ├── sc_analysis.html
    ├── bulk_analysis.html
    ├── analysis_select.html
    ├── analysis_result.html
    └── results_gallery.html
```

## 配置

编辑 `config.py` 修改：

- `DATA_DIR`: 数据存储目录
- `DB_PATH`: SQLite 数据库路径
- `MAX_WORKERS`: 并行任务数
- `CUDA_DEVICES`: GPU 设备号（SysVI 批次校正用）

## License

MIT
