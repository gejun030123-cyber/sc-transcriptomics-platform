# 平台分析流程可自定义化增强设计

**日期**: 2026-06-08
**状态**: 已完成（Phase 1-4 全部实施，2026-06-18）
**方案**: 方案 C — 底层统一 schema + 表层渐进
**用户群体**: 有经验的生信用户（了解分析原理，会调参数，不想写代码）

---

## 背景与目标

当前平台的分析流程存在以下限制：
- 模块执行顺序由 `PIPELINE_ORDER` 硬编码，用户无法调整
- 参数无法保存/复用，每次分析需手动填写
- 图表样式由各模块硬编码，无法统一控制
- 数据过滤阈值固定，无法灵活组合条件

本设计引入统一的 JSON 配置 schema，分 4 个阶段逐步增强自定义能力，每个阶段独立交付、向后兼容。

---

## 第 1 部分：配置 Schema 设计

所有自定义能力统一存储为 JSON 格式，每个顶层 section 对应一个实现阶段：

```json
{
  "name": "PBMC 3k 标准流程",
  "description": "单细胞标准 QC → 聚类 → 注释流程",
  "version": 1,
  "analysis_type": "sc",
  "created_at": "2026-06-08T10:00:00",

  "pipeline": {
    "modules": ["qc", "preprocess", "dimred", "clustering", "annotation"],
    "skip_validation": false
  },

  "params": {
    "qc": {
      "min_genes": 200,
      "max_genes": 5000,
      "max_pct_mt": 20,
      "scrublet": true
    },
    "clustering": {
      "resolution": 0.8
    }
  },

  "filters": {
    "qc": [
      {"column": "total_counts", "op": ">=", "value": 1000, "label": "最小文库大小"},
      {"column": "pct_counts_mt", "op": "<=", "value": 15, "label": "最大线粒体比例"}
    ]
  },

  "visualization": {
    "theme": "default",
    "color_palette": ["#e53935", "#1a237e", "#43a047", "#fb8c00"],
    "figure_width": 800,
    "figure_height": 500,
    "font_family": "Arial",
    "font_size": 12,
    "bg_color": "white",
    "grid": false,
    "export_formats": ["plotly_json"]
  }
}
```

### 设计要点

- 每个顶层 section 可独立实现，未设置的 section 使用平台默认值
- `pipeline.modules` 列表顺序即执行顺序，受依赖约束
- `params` 中的 key 与现有 `PARAM_SCHEMAS` 完全对齐，无需额外映射
- `version` 字段用于未来 schema 升级时的兼容处理
- **config 传递机制**: `routes/analysis.py` 在构建模块参数时，将 `visualization` 和 `filters` 以 `_visualization` 和 `_filters` 键注入 `self.params`，模块内部通过 `self.params.get('_visualization', {})` 访问。`pipeline` 部分由路由层直接消费，不传递给模块。

---

## 第 2 部分：Phase 1 — 参数预设系统

最高频、最独立的功能。用户可保存/加载常用参数配置。

### 存储设计

```
projects/{project_id}/presets/
├── preset_001.json
├── preset_002.json
└── _global/
    ├── pbmc_standard.json
    └── bulk_deseq2.json
```

- **项目级预设**: 存储在项目目录内，仅当前项目可用
- **全局预设**: 存储在 `_global/` 目录，所有项目可选
- 每个预设文件包含 `name`、`description`、`analysis_type`、`params` 字段

### API 端点

| 方法 | 路径 | 说明 |
|------|------|------|
| `GET` | `/api/presets?project_id=X&type=sc` | 列出项目级 + 全局预设 |
| `POST` | `/api/presets` | 保存当前参数为预设 |
| `DELETE` | `/api/presets/<id>` | 删除预设 |
| `GET` | `/api/presets/<id>` | 获取预设详情 |

### 前端交互

在 `analysis_select.html` 参数表单上方增加工具栏：

```
[保存预设 ▾] [加载预设 ▾]
```

- **保存**: 弹出模态框输入名称和描述 → 调用 `POST /api/presets`
- **加载**: 下拉列表显示所有预设（项目级 + 全局） → 选择后自动填充当前模块的表单字段
- 加载时只覆盖当前模块的参数，不影响其他模块
- 预设文件缺失的字段使用 `PARAM_SCHEMAS` 中的 `default` 值

### 向后兼容

无预设时行为完全不变。现有参数表单的默认值逻辑不受影响。

---

## 第 3 部分：Phase 2 — 模块选择与流程编排

让用户自由选择执行哪些模块、以什么顺序执行。

### 当前限制

`PIPELINE_ORDER` 硬编码模块顺序，用户只能勾选/取消，不能调整顺序。

### 依赖约束

部分模块有硬性前置依赖，前端实时校验：

```
qc → preprocess → dimred → clustering → annotation
qc → preprocess → dimred → clustering → deg
qc → preprocess → dimred → clustering → trajectory
qc → preprocess → dimred → clustering → qc_reassess
qc → preprocess → dimred → clustering → proportion
preprocess → batch_correct → dimred（batch_correct 可选插入）
```

规则：
- `preprocess` 必须在 `dimred` 之前
- `dimred` 必须在 `clustering` 之前
- `clustering` 必须在 `annotation`、`deg`、`trajectory`、`qc_reassess`、`proportion` 之前
- `batch_correct` 必须在 `dimred` 之前（如使用）
- 用户拖拽排序时，违反约束的排列显示红色警告，阻止提交

### 前端交互

`analysis_select.html` 模块选择区域改造为可拖拽列表：

```
┌──────────────────────────────────────────────────┐
│  ☑ QC                 [≡ 拖拽]  ⚙ 参数展开        │
│  ☑ 预处理              [≡ 拖拽]  ⚙ 参数展开        │
│  ☐ 批次校正            [≡ 拖拽]  ⚙ 参数展开        │
│  ☑ 降维                [≡ 拖拽]  ⚙ 参数展开        │
│  ☑ 聚类                [≡ 拖拽]  ⚙ 参数展开        │
│  ☑ 注释                [≡ 拖拽]  ⚙ 参数展开        │
│  ────────────────────────────────────────────── │
│  [保存为流程模板]  [加载流程模板]                    │
└──────────────────────────────────────────────────┘
```

- 每行: 复选框 + 模块名 + 拖拽手柄 + 参数折叠/展开按钮
- 依赖模块被取消勾选时，自动取消其下游模块
- 参数区域按模块折叠/展开（现有行为保留）
- 流程模板与参数预设共享存储目录，文件包含 `pipeline` + `params` 部分

### API 变更

`POST /api/analyze` 请求体新增 `pipeline` 字段：

```json
{
  "project_id": "...",
  "file_path": "...",
  "pipeline": {
    "modules": ["qc", "preprocess", "dimred", "clustering"],
    "skip_validation": false
  },
  "params": { ... }
}
```

后端 `routes/analysis.py` 根据 `pipeline.modules` 列表替代硬编码的 `PIPELINE_ORDER`。`skip_validation` 为 `true` 时跳过依赖检查（高级用户用于实验性流程）。

---

## 第 4 部分：Phase 3 — 可视化样式自定义

统一控制所有图表的视觉风格。

### 样式参数

```json
{
  "visualization": {
    "theme": "default",
    "color_palette": ["#e53935", "#1a237e", "#43a047"],
    "figure_width": 800,
    "figure_height": 500,
    "font_family": "Arial",
    "font_size": 12,
    "bg_color": "white",
    "grid": false,
    "export_formats": ["plotly_json"]
  }
}
```

### 预设主题

| 主题 | 背景 | 配色风格 | 适用场景 |
|------|------|----------|----------|
| `default` | 白色 | 平台默认色板 | 通用 |
| `nature` | 白色 | Nature/Science 论文配色 | 投稿 |
| `dark` | #1a1a2e | 高对比度 | 演示/海报 |

### 实现方式

在 `modules/base.py` 的 `BaseAnalysis` 中新增 `get_plotly_layout()` 方法：

```python
def get_plotly_layout(self, title='', **overrides):
    viz = self.params.get('_visualization', {})
    layout = {
        'title': title,
        'width': viz.get('figure_width', 800),
        'height': viz.get('figure_height', 500),
        'plot_bgcolor': viz.get('bg_color', 'white'),
        'font': {
            'family': viz.get('font_family', 'Arial'),
            'size': viz.get('font_size', 12),
        },
    }
    layout.update(overrides)
    return layout
```

各模块的 `fig.update_layout(...)` 渐进迁移为 `fig.update_layout(self.get_plotly_layout(title='...'))`，不需要一次性改完所有模块。

### 配色策略

- `color_palette` 作为分类变量的默认色板
- 各模块内部的语义颜色映射（如 volcano 的 Up/Down/NS）优先保留
- `color_palette` 仅在无显式映射时生效，不破坏现有图表语义

### 静态导出

当 `export_formats` 包含 `"svg"` 或 `"png"` 时：
- 后端调用 `plotly.io.write_image()` 生成静态文件
- 依赖 `kaleido` 库（`config.py` 中已有相关配置）
- 静态文件与 Plotly JSON 并行输出，不替代交互式图表

### 前端交互

在项目设置页面或分析页面顶部增加"可视化设置"折叠面板：

```
▾ 可视化设置
  主题:  [default ▾]      配色:  [默认 ▾]
  宽度:  [800]            高度:  [500]
  字号:  [12]             字体:  [Arial ▾]
  导出:  ☑ plotly_json    ☐ svg    ☐ png
```

---

## 第 5 部分：Phase 4 — 数据过滤条件自定义

用声明式规则定义过滤条件，替代硬编码阈值。

### 过滤规则格式

```json
{
  "filters": {
    "qc": [
      {"column": "total_counts", "op": ">=", "value": 1000, "label": "最小文库大小"},
      {"column": "pct_counts_mt", "op": "<=", "value": 15, "label": "最大线粒体比例"},
      {"column": "n_genes_by_counts", "op": ">=", "value": 200, "label": "最小基因数"}
    ],
    "annotation": [
      {"column": "celltype", "op": "in", "value": ["T cell", "B cell", "NK cell"], "label": "仅保留免疫细胞"}
    ]
  }
}
```

### 支持的操作符

| 操作符 | 含义 | value 类型 | 示例 |
|--------|------|-----------|------|
| `>=` | 大于等于 | number | `{"column": "total_counts", "op": ">=", "value": 1000}` |
| `<=` | 小于等于 | number | `{"column": "pct_counts_mt", "op": "<=", "value": 15}` |
| `==` | 等于 | string/number | `{"column": "batch", "op": "==", "value": "control"}` |
| `!=` | 不等于 | string/number | `{"column": "leiden", "op": "!=", "value": "0"}` |
| `in` | 在列表中 | array | `{"column": "celltype", "op": "in", "value": ["A","B"]}` |
| `not_in` | 不在列表中 | array | `{"column": "batch", "op": "not_in", "value": ["bad"]}` |
| `between` | 范围 | [min, max] | `{"column": "total_counts", "op": "between", "value": [500, 50000]}` |

### 执行逻辑

在 `BaseAnalysis` 中新增 `apply_filters(adata, module_name)` 方法：

```python
def apply_filters(self, adata, module_name):
    filters = self.params.get('_filters', {}).get(module_name, [])
    for rule in filters:
        col, op, val = rule['column'], rule['op'], rule['value']
        if col not in adata.obs.columns:
            continue
        if op == '>=':       mask = adata.obs[col] >= val
        elif op == '<=':     mask = adata.obs[col] <= val
        elif op == '==':     mask = adata.obs[col] == val
        elif op == '!=':     mask = adata.obs[col] != val
        elif op == 'in':     mask = adata.obs[col].isin(val)
        elif op == 'not_in': mask = ~adata.obs[col].isin(val)
        elif op == 'between': mask = adata.obs[col].between(val[0], val[1])
        else: continue
        n_before = adata.n_obs
        adata = adata[mask].copy()
        n_filtered = n_before - adata.n_obs
        if n_filtered > 0:
            self.progress(-1, f"过滤 {col} {op} {val}: 移除 {n_filtered} 个样本")
    return adata
```

各模块在 `run()` 开头调用 `adata = self.apply_filters(adata, 'module_name')` 即可激活过滤。

### 与现有过滤的关系

- 现有 `PARAM_SCHEMAS` 中的阈值参数（如 `min_genes`、`max_pct_mt`）保留作为默认值
- `_filters` 中的规则在默认阈值之后执行，作为"二次过滤"
- 不破坏现有参数的语义，提供更灵活的补充

### 前端交互

在每个模块的参数区域底部增加"自定义过滤规则"折叠区：

```
▾ 自定义过滤规则（可选）
  列名: [total_counts ▾]  操作: [>= ▾]  值: [1000]  [+ 添加]
  列名: [pct_counts_mt]   操作: [<=]     值: [15]    [× 删除]
```

- 列名下拉自动填充 `adata.obs.columns`（通过 `/api/obs-columns` 获取）
- 操作符下拉为固定列表
- 可添加多条规则，规则之间为 AND 关系
- 规则随参数预设一起保存/加载

---

## 实施阶段总览

| 阶段 | 功能 | 优先级 | 预计时间 | 依赖 |
|------|------|--------|----------|------|
| Phase 1 | 参数预设系统 | 最高 | 2h | 无 |
| Phase 2 | 模块选择与流程编排 | 高 | 3h | Phase 1（共享存储） |
| Phase 3 | 可视化样式自定义 | 中 | 2h | Phase 1 |
| Phase 4 | 数据过滤条件自定义 | 中 | 2h | Phase 1 |

Phase 1 和 Phase 2 可并行开发（代码无直接依赖，存储目录独立）。Phase 3 和 Phase 4 可并行开发。

---

## 改动文件汇总

| 文件 | Phase | 改动类型 | 说明 |
|------|-------|----------|------|
| `modules/base.py` | 1,3,4 | 修改 | 新增 `apply_filters()`、`get_plotly_layout()` |
| `routes/api.py` | 1 | 修改 | 新增 `/api/presets` 端点 |
| `routes/analysis.py` | 2 | 修改 | 支持 `pipeline.modules` 替代 `PIPELINE_ORDER` |
| `templates/analysis_select.html` | 1,2 | 修改 | 预设工具栏 + 可拖拽模块列表 |
| `modules/__init__.py` | 2 | 修改 | 依赖约束校验逻辑 |
| 各分析模块 | 3 | 修改 | 渐进迁移 `fig.update_layout` |

---

## 错误处理

| 场景 | 处理方式 |
|------|----------|
| 预设文件格式错误（JSON 解析失败） | 忽略该预设，在列表中标记"损坏" |
| 预设中包含不存在的参数 key | 忽略未知 key，使用 `PARAM_SCHEMAS` 默认值 |
| 流程编排违反依赖约束 | 前端阻止提交 + 后端二次校验返回错误 |
| 过滤规则引用不存在的列 | 跳过该规则，在 summary 中报告 |
| `kaleido` 未安装时请求 SVG 导出 | 跳过静态导出，仅生成 Plotly JSON，在 summary 中提示 |
| 全局预设目录不存在 | 自动创建，不影响功能 |

---

## 不做的事情（YAGNI）

- 不支持 YAML/TOML 配置格式（JSON 足够，前端原生支持）
- 不引入完整的 workflow DSL（如 Snakemake/Nextflow）
- 不支持跨项目的参数继承链
- 不实现实时协作编辑配置
- 不添加配置版本控制系统（Git 管理即可）

---

## 测试与验收

### Phase 1 验收

- 保存当前参数为预设 → 重新加载后表单自动填充
- 项目级预设仅在当前项目可见
- 全局预设在所有项目可见
- 预设缺失的字段回退到默认值

### Phase 2 验收

- 取消勾选模块 → 该模块不执行
- 拖拽调整顺序 → 执行顺序变更
- 违反依赖约束 → 前端显示警告，阻止提交
- 保存/加载流程模板 → 模块选择和顺序恢复

### Phase 3 验收

- 切换主题 → 所有新生成的图表样式变更
- 修改尺寸/字号 → 图表按新值渲染
- SVG 导出 → 生成 `.svg` 文件

### Phase 4 验收

- 添加过滤规则 → 对应模块的数据被过滤
- 多条规则 → AND 关系生效
- 引用不存在的列 → 该规则被跳过，不报错

### 验证命令

```bash
# 验证配置 schema 解析
python -c "
import json
with open('test_config.json') as f:
    config = json.load(f)
assert 'pipeline' in config
assert 'params' in config
assert 'filters' in config
assert 'visualization' in config
print('Schema OK')
"

# 验证依赖约束
python -c "
from modules import PIPELINE_DEPS
assert 'dimred' in PIPELINE_DEPS.get('clustering', [])
print('Dependencies OK')
"
```
