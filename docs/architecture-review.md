# 架构审查报告：`/data/GJ/platform`

**日期：** 2026-06-23
**范围：** 全项目架构深度分析

---

## 候选 1：BaseAnalysis 深化 — 消除 13 个模块的重复样板代码

**推荐强度：`Strong`**

**涉及文件：**
- `modules/base.py`（176 行）
- 13 个单细胞模块：`qc.py`, `normalize.py`, `hvg.py`, `dimred.py`, `batch_correct.py`, `clustering.py`, `deg.py`, `trajectory.py`, `proportion.py`, `cell_communication.py`, `annotation.py`, `qc_reassess.py`
- 8 个 Bulk 模块（类似模式）

**问题：**

每个模块的 `run()` 方法都复制粘贴了相同的 10-15 行样板代码：

```python
# 出现在 13 个模块中，完全相同
adata = sc.read_h5ad(input_path)
adata = remap_var_names(adata)
# ... 分析逻辑 ...
plots_dir = os.path.join(self.project_dir, 'plots')
os.makedirs(plots_dir, exist_ok=True)
result_files = []
# ... 绘图 ...
intermediate_dir = os.path.join(self.project_dir, 'intermediate')
os.makedirs(intermediate_dir, exist_ok=True)
output_path = os.path.join(intermediate_dir, '<name>_output.h5ad')
adata.write_h5ad(output_path)
return {'output_adata': output_path, 'result_files': result_files, 'summary': {...}}
```

同时，`base.py` 提供了 `get_plotly_layout()`（第 87 行）、`export_results()`（第 134 行）、`apply_filters()`（第 42 行），但 **几乎没有任何模块调用它们**。唯一例外是 `proportion.py` 调用了 `get_plotly_layout()`。base 类的接口是**声明性的装饰**，而非实际的深度抽象。

**此外：**
- `validate_input()` 在 15/21 个模块中是 `return None` 的死代码
- `PARAM_SCHEMA = {}` 从未被任何模块填充
- `INPUT_REQUIRES` 声明了但从未被框架检查

**Before / After：**

```
  Before (当前)                          After (深化后)
+---------------------+            +-------------------------+
|   BaseAnalysis       |            |   BaseAnalysis           |
|  +---------------+   |            |  +-------------------+  |
|  | get_plotly_   |   |            |  | run(input_path)   |  |  <-- 模板方法
|  |  layout()     |   | <-- 无人调用 |  |  load_data()      |  |  <-- 新增
|  | export_       |   |            |  |  save_output()    |  |  <-- 新增
|  |  results()    |   |            |  |  build_layout()   |  |  <-- 统一
|  | apply_        |   |            |  |  validate_input() |  |  <-- 默认实现
|  |  filters()    |   |            |  +-------------------+  |
|  +---------------+   |            +-------------------------+
+---------------------+                      |
| 13 个子类各含:        |            +-------+--------------+
|  253行 run() 方法     |            | 13 个子类只需实现:     |
|  重复样板代码 x13     |            |  analyze(adata)      |  <-- 纯分析逻辑
|  手写 plotly layout   |            |  build_plots()       |  <-- 可选
+---------------------+            +----------------------+
```

**收益：**
- **Leverage：** 一个接口（`run()` 模板方法），13 个调用点受益
- **Locality：** 数据加载/保存逻辑集中在一个模块，改一处修全部
- 测试可以只测 `analyze()` 纯逻辑，无需构造文件系统

---

## 候选 2：PARAM_SCHEMAS 从路由文件提取到模块类

**推荐强度：`Strong`**

**涉及文件：**
- `routes/analysis.py`（499 行，其中 338 行是 schema 数据）
- 21 个模块类

**问题：**

`routes/analysis.py` 第 50-388 行是一个 338 行的 `PARAM_SCHEMAS` 字典，为每个分析模块定义参数配置。这是**领域数据**，不是路由逻辑。它让 `analysis.py` 成为项目中维护成本最高的文件之一。

同时，`analysis.py` 的 `analyze()` 函数（第 411-498 行）是一个**胖处理器**，包含：
- 参数解析和类型强制转换（第 429-436 行）
- JSON 反序列化（第 438-460 行）
- 任务创建和数据库持久化（第 461-467 行）
- Worker 调度（第 468 行）
- 文件系统扫描（第 473-483 行）
- `bulk_enrichment` 的特殊业务逻辑（第 486-492 行）

**Before / After：**

```
  Before                              After
+----------------------+        +----------------------+
| routes/analysis.py    |        | routes/analysis.py    |
|  PARAM_SCHEMAS (338行)| <-- 膨胀  |  analyze() (~30行)   | <-- 瘦路由
|  analyze() (~90行)    |        |  只做: 调用服务层      |
|  参数解析             |        +----------+-----------+
|  任务创建             |                   |
|  文件扫描             |        +----------+-----------+
|  特殊业务逻辑         |        | modules/base.py       |
+----------------------+        |  PARAM_SCHEMA = {...}  | <-- 每个模块自描述
                                |  (或单独 schema 文件)   |
+----------------------+        +-----------------------+
| modules/qc.py         |        | modules/qc.py          |
|  PARAM_SCHEMA = {}    | <-- 空    |  PARAM_SCHEMA = {      | <-- 自描述
|  run() 253行          |        |   'n_genes': {...},    |
+----------------------+        |   'pct_mito': {...}   |
                                |  }                     |
                                +-----------------------+
```

**收益：**
- **Leverage：** 模块自描述，路由层只需 `registry[key].PARAM_SCHEMA`
- **Locality：** 参数定义和使用它的模块在同一位置
- `analysis.py` 从 499 行降至 ~50 行

---

## 候选 3：统一数据库连接管理

**推荐强度：`Strong`**

**涉及文件：**
- `database.py`（64 行）— 连接工厂 1：`get_db()` 单例
- `models.py`（225 行）— 连接工厂 2：`_get_conn()` 每次调用新建
- `worker.py`（91 行）— 连接工厂 3：任务级独立连接

**问题：**

项目中存在 **3 个独立的 SQLite 连接策略**，各自设置相同的 PRAGMA，互不共享：

| 来源 | 连接方式 | 使用者 |
|------|----------|--------|
| `database.py` `get_db()` | 单例 | 仅 `init_db()` |
| `models.py` `_get_conn()` | 每次方法调用新建 | 所有路由处理器 |
| `worker.py` 第 26-29 行 | 每个任务新建 | 后台线程 |

`worker.py` 完全绕过 `AnalysisTask` 模型，直接用原始 SQL 更新任务状态、进度和结果文件（第 36-86 行）。这意味着任何 schema 变更都必须手动同步 `models.py` 和 `worker.py`。

**Before / After：**

```
  Before                              After
+-------------+                  +---------------------+
|database.py  | <-- 连接1           | database.py          |
| get_db()    |                   |  get_conn()          | <-- 唯一连接工厂
+-------------+                   |  (线程安全, 带PRAGMA) |
+-------------+                   +----------+----------+
| models.py   | <-- 连接2                      |
| _get_conn() |                   +----------+----------+
+-------------+                   | models.py            |
+-------------+                   |  使用 get_conn()     | <-- 统一
| worker.py   | <-- 连接3 + 原始SQL  +---------------------+
| 独立连接     |                   | worker.py            |
| 绕过ORM     |                   |  使用 AnalysisTask   | <-- 通过模型
+-------------+                   |  .save() 更新状态     |
                                  +---------------------+
```

**收益：**
- **Leverage：** 一个连接工厂，3 个调用点
- **Locality：** schema 变更只需改一处
- worker 通过模型更新状态，消除原始 SQL 同步负担

---

## 候选 4：提取路径构建工具 + 消除跨文件重复

**推荐强度：`Worth exploring`**

**涉及文件：**
- `routes/analysis.py`（第 473、487 行）
- `routes/projects.py`（第 34 行）
- `routes/upload.py`（第 87、103、117 行）
- `routes/api.py`（第 131、275、323 行）

**问题：**

`os.path.join(Config.DATA_DIR, 'projects', pid, 'uploads')` 和 `os.path.join(Config.DATA_DIR, 'projects', pid, 'results')` 在 5+ 个文件中重复出现。没有工具函数计算项目子目录路径。

同时，`analysis.py`（第 473-483 行）和 `projects.py`（第 34-40 行）有重复的上传文件扫描逻辑。

**收益：**
- **Leverage：** 一个路径工具，N 个调用点
- 改变存储结构时只改一处

---

## 候选 5：api.py 域逻辑提取

**推荐强度：`Worth exploring`**

**涉及文件：**
- `routes/api.py`（665 行，项目最大路由文件）

**问题：**

`api.py` 包含 ~200 行与 HTTP 路由无关的域逻辑：

- `adata_info()`（第 124-226 行）：~100 行 anndata 内省逻辑
- `data_info_full()`（第 496-576 行）：与 `adata_info()` 重复 ~80 行近乎相同的代码
- `validate_filter_expression()`（第 310-410 行）：~100 行分析逻辑（加载 CSV、构建矩阵、调用表达式解析器）
- `deg_comparisons()`（第 273-307 行）：文件系统扫描 + 数据库查询内联

**收益：**
- **Locality：** anndata 内省逻辑集中，消除 80 行重复
- 路由文件从 665 行降至 ~300 行

---

## 候选 6：死代码清理 + 共享常量提取

**推荐强度：`Worth exploring`**

**涉及文件：**
- `modules/preprocess.py`（76 行，未注册，完全死代码）
- `modules/qc.py`（第 4-20 行）和 `modules/hvg.py`（第 5-21 行）：`S_GENES` 和 `G2M_GENES` 重复定义
- `routes/analysis.py`（第 42-48 行）：`STATUS_MAP` 被其他路由文件当共享常量导入

**收益：**
- 删除 `preprocess.py` 消除维护歧义
- 基因列表提取到共享位置，编辑一次同步全部
- `STATUS_MAP` 移到 `models.py` 或独立常量模块

---

## 首选推荐

**候选 1（BaseAnalysis 深化）** 是最值得首先处理的。

理由：它影响 13 个模块，是最大的 leverage 机会。深化 base 类后，每个模块只需实现纯分析逻辑，样板代码从每个模块 10-15 行降至 0 行。测试可以通过 `analyze()` 方法直接测试，无需构造文件系统。同时为候选 2（PARAM_SCHEMAS 提取）铺路——模块自描述后，路由层自然变瘦。
