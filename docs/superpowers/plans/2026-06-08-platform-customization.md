# 平台分析流程可自定义化增强 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为现有生信分析平台添加统一的 JSON 配置系统，分 4 个阶段实现参数预设、流程编排、可视化自定义和数据过滤自定义。

**Architecture:** 在 `modules/base.py` 添加共享方法（`get_plotly_layout`、`apply_filters`），在 `routes/api.py` 添加预设 CRUD 端点，在 `modules/__init__.py` 添加依赖约束定义，在 `templates/analysis_select.html` 添加预设工具栏和可拖拽模块列表。配置以 JSON 文件存储在 `data/projects/{id}/presets/` 和 `data/presets/_global/`。

**Tech Stack:** Flask, Jinja2, Plotly, SQLite (元数据), JSON (配置存储)

**Spec:** `docs/superpowers/specs/2026-06-08-platform-customization-design.md`

---

## 文件结构总览

| 文件 | 操作 | Phase | 职责 |
|------|------|-------|------|
| `modules/base.py` | 修改 | 1,3,4 | `apply_filters()`、`get_plotly_layout()` |
| `routes/api.py` | 修改 | 1 | `/api/presets` CRUD 端点 |
| `modules/__init__.py` | 修改 | 2 | `PIPELINE_DEPS` 依赖约束定义 |
| `routes/analysis.py` | 修改 | 2 | 支持 `pipeline.modules` 替代硬编码顺序 |
| `templates/analysis_select.html` | 修改 | 1,2 | 预设工具栏 + 可拖拽模块列表 |

---

## Task 1: 基础设施 — base.py 新增 `apply_filters()` 和 `get_plotly_layout()`

**Files:**
- Modify: `modules/base.py`

### Step 1: 添加 `apply_filters()` 方法

在 `BaseAnalysis` 类中，在 `run()` 抽象方法之前添加：

```python
def apply_filters(self, adata, module_name):
    """根据 self.params['_filters'][module_name] 中的声明式规则过滤 adata。"""
    filters = self.params.get('_filters', {}).get(module_name, [])
    for rule in filters:
        col = rule.get('column', '')
        op = rule.get('op', '')
        val = rule.get('value')
        if col not in adata.obs.columns:
            continue
        if op == '>=':
            mask = adata.obs[col] >= val
        elif op == '<=':
            mask = adata.obs[col] <= val
        elif op == '==':
            mask = adata.obs[col] == val
        elif op == '!=':
            mask = adata.obs[col] != val
        elif op == 'in':
            mask = adata.obs[col].isin(val)
        elif op == 'not_in':
            mask = ~adata.obs[col].isin(val)
        elif op == 'between':
            mask = adata.obs[col].between(val[0], val[1])
        else:
            continue
        n_before = adata.n_obs
        adata = adata[mask].copy()
        n_filtered = n_before - adata.n_obs
        if n_filtered > 0:
            self.progress(-1, f"过滤 {col} {op} {val}: 移除 {n_filtered} 个样本")
    return adata
```

### Step 2: 添加 `get_plotly_layout()` 方法

紧接 `apply_filters()` 之后添加：

```python
VISUALIZATION_THEMES = {
    'default': {
        'bg_color': 'white',
        'color_palette': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'],
        'font_family': 'Arial',
    },
    'nature': {
        'bg_color': 'white',
        'color_palette': ['#E64B35', '#4DBBD5', '#00A087', '#3C5488',
                          '#F39B7F', '#8491B4', '#91D1C2', '#DC0000'],
        'font_family': 'Helvetica',
    },
    'dark': {
        'bg_color': '#1a1a2e',
        'color_palette': ['#e94560', '#0f3460', '#16213e', '#533483',
                          '#e94560', '#00b4d8', '#48cae4', '#90e0ef'],
        'font_family': 'Arial',
    },
}

def get_plotly_layout(self, title='', **overrides):
    """根据 self.params['_visualization'] 生成统一的 Plotly layout dict。"""
    viz = self.params.get('_visualization', {})
    theme_name = viz.get('theme', 'default')
    theme = VISUALIZATION_THEMES.get(theme_name, VISUALIZATION_THEMES['default'])
    layout = {
        'title': title,
        'width': viz.get('figure_width', 800),
        'height': viz.get('figure_height', 500),
        'plot_bgcolor': viz.get('bg_color', theme['bg_color']),
        'font': {
            'family': viz.get('font_family', theme['font_family']),
            'size': viz.get('font_size', 12),
        },
    }
    layout.update(overrides)
    return layout
```

### Step 3: 验证

```bash
cd /data/GJ/platform && python -c "
from modules.base import BaseAnalysis, VISUALIZATION_THEMES
assert 'default' in VISUALIZATION_THEMES
assert 'nature' in VISUALIZATION_THEMES
assert 'dark' in VISUALIZATION_THEMES
assert hasattr(BaseAnalysis, 'apply_filters')
assert hasattr(BaseAnalysis, 'get_plotly_layout')
print('base.py OK')
"
```

### Step 4: Commit

```bash
git add modules/base.py
git commit -m "feat(base): add apply_filters() and get_plotly_layout() for customization"
```

---

## Task 2: Phase 1 — 预设 CRUD API 端点

**Files:**
- Modify: `routes/api.py`

### Step 1: 添加预设目录常量和辅助函数

在 `routes/api.py` 顶部 import 区域之后、现有路由之前添加：

```python
import uuid
from config import Config

PRESETS_GLOBAL_DIR = os.path.join(Config.DATA_DIR, 'presets', '_global')

def _get_project_presets_dir(project_id):
    return os.path.join(Config.DATA_DIR, 'projects', project_id, 'presets')

def _ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def _load_preset(filepath):
    """加载单个预设文件，返回 dict 或 None（解析失败时）。"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        data['_filepath'] = filepath
        data['_id'] = os.path.splitext(os.path.basename(filepath))[0]
        return data
    except (json.JSONDecodeError, IOError):
        return None

def _list_presets_in_dir(directory, scope):
    """列出目录下所有预设文件，返回 list[dict]。"""
    if not os.path.isdir(directory):
        return []
    presets = []
    for fname in sorted(os.listdir(directory)):
        if not fname.endswith('.json'):
            continue
        p = _load_preset(os.path.join(directory, fname))
        if p:
            p['_scope'] = scope
            presets.append(p)
    return presets
```

### Step 2: 添加 GET /api/presets 列表端点

```python
@api_bp.route('/presets')
def list_presets():
    """列出项目级 + 全局预设。可选 ?project_id=X&type=sc 过滤。"""
    project_id = request.args.get('project_id', '')
    filter_type = request.args.get('type', '')

    presets = []

    # 全局预设
    _ensure_dir(PRESETS_GLOBAL_DIR)
    presets.extend(_list_presets_in_dir(PRESETS_GLOBAL_DIR, 'global'))

    # 项目级预设
    if project_id:
        proj_dir = _get_project_presets_dir(project_id)
        presets.extend(_list_presets_in_dir(proj_dir, 'project'))

    # 按 type 过滤
    if filter_type:
        presets = [p for p in presets if p.get('analysis_type') == filter_type]

    # 清理内部字段，返回给前端
    result = []
    for p in presets:
        result.append({
            'id': p['_id'],
            'name': p.get('name', ''),
            'description': p.get('description', ''),
            'analysis_type': p.get('analysis_type', ''),
            'scope': p['_scope'],
        })

    return jsonify({'presets': result})
```

### Step 3: 添加 GET /api/presets/<preset_id> 详情端点

```python
@api_bp.route('/presets/<preset_id>')
def get_preset(preset_id):
    """获取预设详情。先查项目级，再查全局。"""
    project_id = request.args.get('project_id', '')

    # 先查项目级
    if project_id:
        fpath = os.path.join(_get_project_presets_dir(project_id), f'{preset_id}.json')
        p = _load_preset(fpath)
        if p:
            return jsonify({'preset': p})

    # 再查全局
    fpath = os.path.join(PRESETS_GLOBAL_DIR, f'{preset_id}.json')
    p = _load_preset(fpath)
    if p:
        return jsonify({'preset': p})

    return jsonify({'error': '预设不存在'}), 404
```

### Step 4: 添加 POST /api/presets 保存端点

```python
@api_bp.route('/presets', methods=['POST'])
def save_preset():
    """保存参数为预设。body: {name, description, analysis_type, params, project_id?, scope?}"""
    data = request.get_json()
    if not data or not data.get('name'):
        return jsonify({'error': '预设名称不能为空'}), 400

    preset_id = str(uuid.uuid4())[:8]
    preset = {
        'name': data['name'],
        'description': data.get('description', ''),
        'analysis_type': data.get('analysis_type', ''),
        'params': data.get('params', {}),
        'pipeline': data.get('pipeline', {}),
        'filters': data.get('filters', {}),
        'visualization': data.get('visualization', {}),
    }

    scope = data.get('scope', 'project')
    project_id = data.get('project_id', '')

    if scope == 'global':
        target_dir = PRESETS_GLOBAL_DIR
    else:
        if not project_id:
            return jsonify({'error': '项目级预设需要 project_id'}), 400
        target_dir = _get_project_presets_dir(project_id)

    _ensure_dir(target_dir)
    fpath = os.path.join(target_dir, f'{preset_id}.json')
    with open(fpath, 'w', encoding='utf-8') as f:
        json.dump(preset, f, ensure_ascii=False, indent=2)

    return jsonify({'id': preset_id, 'message': '预设已保存'})
```

### Step 5: 添加 DELETE /api/presets/<preset_id> 删除端点

```python
@api_bp.route('/presets/<preset_id>', methods=['DELETE'])
def delete_preset(preset_id):
    """删除预设。"""
    project_id = request.args.get('project_id', '')

    # 先查项目级
    if project_id:
        fpath = os.path.join(_get_project_presets_dir(project_id), f'{preset_id}.json')
        if os.path.isfile(fpath):
            os.remove(fpath)
            return jsonify({'message': '预设已删除'})

    # 再查全局
    fpath = os.path.join(PRESETS_GLOBAL_DIR, f'{preset_id}.json')
    if os.path.isfile(fpath):
        os.remove(fpath)
        return jsonify({'message': '预设已删除'})

    return jsonify({'error': '预设不存在'}), 404
```

### Step 6: 验证

```bash
cd /data/GJ/platform && python -c "
import ast
ast.parse(open('routes/api.py').read())
print('api.py syntax OK')
"
```

### Step 7: Commit

```bash
git add routes/api.py
git commit -m "feat(api): add /api/presets CRUD endpoints for parameter presets"
```

---

## Task 3: Phase 1 — 预设前端 UI（保存/加载工具栏）

**Files:**
- Modify: `templates/analysis_select.html`

### Step 1: 在参数表单上方添加预设工具栏

在 `analysis_select.html` 的 `<form method="post" class="param-form">` 之后、`<div class="mb-3">` (输入数据) 之前，插入：

```html
<!-- 预设工具栏 -->
<div class="d-flex gap-2 mb-3 align-items-center" id="preset-toolbar">
    <div class="dropdown">
        <button class="btn btn-outline-secondary btn-sm dropdown-toggle" type="button"
                data-bs-toggle="dropdown" id="btn-load-preset">加载预设</button>
        <ul class="dropdown-menu" id="preset-list">
            <li><span class="dropdown-item text-muted small">加载中...</span></li>
        </ul>
    </div>
    <button type="button" class="btn btn-outline-secondary btn-sm" id="btn-save-preset"
            data-bs-toggle="modal" data-bs-target="#savePresetModal">保存预设</button>
</div>

<!-- 保存预设模态框 -->
<div class="modal fade" id="savePresetModal" tabindex="-1">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <h5 class="modal-title">保存参数预设</h5>
                <button type="button" class="btn-close" data-bs-dismiss="modal"></button>
            </div>
            <div class="modal-body">
                <div class="mb-3">
                    <label class="form-label">预设名称</label>
                    <input type="text" class="form-control" id="preset-name" placeholder="如: PBMC 标准 QC 参数">
                </div>
                <div class="mb-3">
                    <label class="form-label">描述（可选）</label>
                    <input type="text" class="form-control" id="preset-desc" placeholder="简要说明用途">
                </div>
                <div class="mb-3">
                    <label class="form-label">保存范围</label>
                    <select class="form-select" id="preset-scope">
                        <option value="project">当前项目</option>
                        <option value="global">全局（所有项目可用）</option>
                    </select>
                </div>
            </div>
            <div class="modal-footer">
                <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">取消</button>
                <button type="button" class="btn btn-primary" id="btn-confirm-save">保存</button>
            </div>
        </div>
    </div>
</div>
```

### Step 2: 在 `{% block extra_js %}` 的 `<script>` 中添加预设交互逻辑

在现有 `<script>` 标签的开头（`const inputSelect` 之前）插入：

```javascript
// === 预设系统 ===
const PROJECT_ID = '{{ project.id }}';
const MODULE_NAME = '{{ module.name }}';
const ANALYSIS_TYPE = MODULE_NAME.startsWith('bulk_') ? 'bulk' : 'sc';

function collectCurrentParams() {
    const params = {};
    document.querySelectorAll('.param-form input, .param-form select, .param-form textarea').forEach(el => {
        const name = el.name;
        if (!name || name === 'input_path') return;
        if (el.type === 'checkbox') {
            params[name] = el.checked;
        } else {
            params[name] = el.value;
        }
    });
    return params;
}

function applyPresetToForm(presetParams) {
    if (!presetParams) return;
    Object.entries(presetParams).forEach(([key, value]) => {
        const el = document.querySelector(`[name="${key}"]`);
        if (!el) return;
        if (el.type === 'checkbox') {
            el.checked = !!value;
        } else {
            el.value = value;
        }
    });
}

async function loadPresets() {
    try {
        const resp = await fetch(`/api/presets?project_id=${PROJECT_ID}&type=${ANALYSIS_TYPE}`);
        const data = await resp.json();
        const list = document.getElementById('preset-list');
        list.innerHTML = '';
        if (!data.presets || data.presets.length === 0) {
            list.innerHTML = '<li><span class="dropdown-item text-muted small">暂无预设</span></li>';
            return;
        }
        data.presets.forEach(p => {
            const li = document.createElement('li');
            const scopeLabel = p.scope === 'global' ? ' [全局]' : '';
            li.innerHTML = `<a class="dropdown-item" href="#" data-preset-id="${p.id}"
                              data-scope="${p.scope}">${p.name}${scopeLabel}<br>
                              <small class="text-muted">${p.description || ''}</small></a>`;
            li.querySelector('a').addEventListener('click', async (e) => {
                e.preventDefault();
                const pid = e.currentTarget.dataset.presetId;
                const scope = e.currentTarget.dataset.scope;
                const scopeParam = scope === 'global' ? '' : `&project_id=${PROJECT_ID}`;
                const r = await fetch(`/api/presets/${pid}?project_id=${PROJECT_ID}`);
                const d = await r.json();
                if (d.preset && d.preset.params) {
                    applyPresetToForm(d.preset.params[MODULE_NAME] || {});
                }
            });
            list.appendChild(li);
        });
    } catch (e) {
        console.error('加载预设失败:', e);
    }
}

document.getElementById('btn-confirm-save').addEventListener('click', async () => {
    const name = document.getElementById('preset-name').value.trim();
    if (!name) { alert('请输入预设名称'); return; }
    const desc = document.getElementById('preset-desc').value.trim();
    const scope = document.getElementById('preset-scope').value;
    const params = collectCurrentParams();

    const body = {
        name, description: desc, scope, analysis_type: ANALYSIS_TYPE,
        project_id: PROJECT_ID,
        params: { [MODULE_NAME]: params },
    };

    try {
        const resp = await fetch('/api/presets', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify(body),
        });
        const data = await resp.json();
        if (resp.ok) {
            bootstrap.Modal.getInstance(document.getElementById('savePresetModal')).hide();
            document.getElementById('preset-name').value = '';
            document.getElementById('preset-desc').value = '';
            loadPresets();
        } else {
            alert(data.error || '保存失败');
        }
    } catch (e) {
        alert('保存失败: ' + e.message);
    }
});

// 页面加载时获取预设列表
loadPresets();

// 表单提交时注入 _visualization 配置
document.querySelector('.param-form').addEventListener('submit', function() {
    // 如果已加载预设包含 visualization 配置，注入到表单
    let vizHidden = this.querySelector('input[name="_visualization"]');
    if (!vizHidden) {
        vizHidden = document.createElement('input');
        vizHidden.type = 'hidden';
        vizHidden.name = '_visualization';
        this.appendChild(vizHidden);
    }
    // 默认值，后续 Phase 3 前端面板会覆盖此值
    if (!vizHidden.value) {
        vizHidden.value = JSON.stringify({});
    }
});
```

### Step 3: 验证

手动在浏览器中访问任意模块分析页面，确认：
- "加载预设" 和 "保存预设" 按钮出现在参数表单上方
- 点击"保存预设"弹出模态框
- 保存后下拉列表刷新

### Step 4: Commit

```bash
git add templates/analysis_select.html
git commit -m "feat(ui): add preset save/load toolbar to analysis form"
```

---

## Task 4: Phase 2 — 依赖约束定义

**Files:**
- Modify: `modules/__init__.py`

### Step 1: 添加 PIPELINE_DEPS 字典

在 `PIPELINE_ORDER` 之后添加：

```python
# 模块依赖约束：value 中的模块必须在 key 之前执行
PIPELINE_DEPS = {
    'preprocess': ['qc'],
    'dimred': ['preprocess'],
    'batch_correct': ['preprocess'],
    'clustering': ['dimred'],
    'qc_reassess': ['clustering'],
    'annotation': ['clustering'],
    'deg': ['clustering'],
    'trajectory': ['clustering'],
    'proportion': ['clustering'],
}

def validate_pipeline_order(modules):
    """校验模块执行顺序是否满足依赖约束。
    返回 (is_valid: bool, errors: list[str])。
    """
    errors = []
    seen = set()
    for mod in modules:
        deps = PIPELINE_DEPS.get(mod, [])
        for dep in deps:
            if dep not in seen and dep in modules:
                errors.append(f"'{mod}' 依赖 '{dep}'，但 '{dep}' 未在其之前执行")
        seen.add(mod)
    return len(errors) == 0, errors
```

### Step 2: 验证

```bash
cd /data/GJ/platform && python -c "
from modules import validate_pipeline_order, PIPELINE_DEPS

# 合法顺序
ok, errs = validate_pipeline_order(['qc', 'preprocess', 'dimred', 'clustering', 'annotation'])
assert ok, f'Expected valid, got {errs}'

# 非法顺序：dimred 在 preprocess 之前
ok, errs = validate_pipeline_order(['qc', 'dimred', 'preprocess', 'clustering'])
assert not ok
assert any('dimred' in e and 'preprocess' in e for e in errs)

# 合法：跳过 batch_correct
ok, errs = validate_pipeline_order(['qc', 'preprocess', 'dimred', 'clustering'])
assert ok

print('validate_pipeline_order OK')
"
```

### Step 3: Commit

```bash
git add modules/__init__.py
git commit -m "feat(modules): add PIPELINE_DEPS and validate_pipeline_order()"
```

---

## Task 5: Phase 2 — 路由层支持 pipeline 参数注入

**Files:**
- Modify: `routes/analysis.py`

> **Scope note:** 本 Task 仅处理 `_visualization` 和 `_filters` 的注入机制。Spec 中的 `POST /api/analyze` 批量执行端点（一次提交整个 pipeline）需要 worker 层改造，将在后续独立计划中实现。

### Step 1: 导入依赖校验函数

在 `routes/analysis.py` 顶部 import 区添加：

```python
from modules import MODULE_REGISTRY, validate_pipeline_order
```

### Step 2: 修改 analyze 路由支持 pipeline 参数

找到 `if request.method == 'POST':` 块，在获取 `input_path` 之后、创建 `task` 之前，添加 pipeline 校验逻辑。将现有代码：

```python
        input_path = request.form.get('input_path', '')
        if not input_path:
            flash('请选择输入文件', 'danger')
            return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))
        task = AnalysisTask(project_id=pid, module_name=module_name,
                           params_json=json.dumps(params))
```

改为：

```python
        input_path = request.form.get('input_path', '')
        if not input_path:
            flash('请选择输入文件', 'danger')
            return redirect(url_for('analysis.analyze', pid=pid, module_name=module_name))

        # 注入 _visualization 和 _filters 到 params
        viz_json = request.form.get('_visualization', '')
        filters_json = request.form.get('_filters', '')
        if viz_json:
            try:
                params['_visualization'] = json.loads(viz_json)
            except json.JSONDecodeError:
                pass
        if filters_json:
            try:
                params['_filters'] = json.loads(filters_json)
            except json.JSONDecodeError:
                pass

        task = AnalysisTask(project_id=pid, module_name=module_name,
                           params_json=json.dumps(params))
```

### Step 3: 验证

```bash
cd /data/GJ/platform && python -c "
import ast
ast.parse(open('routes/analysis.py').read())
print('analysis.py syntax OK')
"
```

### Step 4: Commit

```bash
git add routes/analysis.py
git commit -m "feat(analysis): inject _visualization and _filters into module params"
```

---

## Task 6: Phase 2 — 可拖拽模块列表前端（模板保存/加载）

**Files:**
- Modify: `templates/analysis_select.html`

> **Note:** 当前阶段的可拖拽列表用于**保存/加载流程模板**（记录模块顺序偏好），不改变单模块提交的行为。批量执行（一次提交多个模块）需要额外的 worker 层改造，不在本计划范围内。

### Step 1: 替换左侧模块列表为可拖拽版本

将现有侧边栏的模块列表部分（`<a href="/projects/...">` 循环）替换为可拖拽列表。在 `analysis_select.html` 的 `<div class="card p-3">` 中，将模块列表循环改为：

```html
<div class="card p-3">
    <h6 class="mb-2">
        {% if module.name.startswith('bulk_') %}Bulk RNA-seq 流程{% else %}单细胞分析流程{% endif %}
    </h6>
    <div id="module-list" class="module-drag-list">
        {% for m in all_modules %}
        <div class="module-item d-flex align-items-center py-1 {% if m.name == module.name %}active{% endif %}"
             data-module="{{ m.name }}" draggable="true">
            <span class="drag-handle me-2" style="cursor:grab; color:#999;">☰</span>
            <a href="/projects/{{ project.id }}/analyze/{{ m.name }}"
               class="text-decoration-none small flex-grow-1">{{ m.display }}</a>
            {% set completed = completed_tasks|selectattr('module_name', 'equalto', m.name)|list %}
            {% if completed %}<span class="text-success small">&#10003;</span>{% endif %}
        </div>
        {% endfor %}
    </div>
    <div class="mt-2 border-top pt-2">
        <button type="button" class="btn btn-outline-secondary btn-sm w-100" id="btn-save-pipeline"
                data-bs-toggle="modal" data-bs-target="#savePipelineModal">保存为流程模板</button>
    </div>
</div>
```

### Step 2: 添加拖拽 CSS

在 `{% block extra_css %}` 或 `<style>` 标签中添加：

```css
.module-drag-list .module-item {
    padding: 4px 8px;
    border-radius: 4px;
    transition: background 0.15s;
}
.module-drag-list .module-item:hover {
    background: #f0f0f0;
}
.module-drag-list .module-item.active {
    background: #e7f1ff;
    font-weight: 600;
}
.module-drag-list .module-item.dragging {
    opacity: 0.4;
}
.module-drag-list .module-item.drag-over {
    border-top: 2px solid #0d6efd;
}
```

### Step 3: 添加拖拽 JS 逻辑

在 `<script>` 中添加：

```javascript
// === 模块拖拽排序 ===
(function() {
    const list = document.getElementById('module-list');
    if (!list) return;
    let dragItem = null;

    list.querySelectorAll('.module-item').forEach(item => {
        item.addEventListener('dragstart', (e) => {
            dragItem = item;
            item.classList.add('dragging');
            e.dataTransfer.effectAllowed = 'move';
        });
        item.addEventListener('dragend', () => {
            if (dragItem) dragItem.classList.remove('dragging');
            list.querySelectorAll('.drag-over').forEach(el => el.classList.remove('drag-over'));
            dragItem = null;
        });
        item.addEventListener('dragover', (e) => {
            e.preventDefault();
            e.dataTransfer.dropEffect = 'move';
            const target = item;
            if (target !== dragItem) {
                target.classList.add('drag-over');
            }
        });
        item.addEventListener('dragleave', () => {
            item.classList.remove('drag-over');
        });
        item.addEventListener('drop', (e) => {
            e.preventDefault();
            item.classList.remove('drag-over');
            if (dragItem && dragItem !== item) {
                const allItems = [...list.querySelectorAll('.module-item')];
                const fromIdx = allItems.indexOf(dragItem);
                const toIdx = allItems.indexOf(item);
                if (fromIdx < toIdx) {
                    item.after(dragItem);
                } else {
                    item.before(dragItem);
                }
            }
        });
    });
})();
```

### Step 4: 添加流程模板保存/加载模态框

在 `savePresetModal` 模态框之后添加：

```html
<!-- 保存流程模板模态框 -->
<div class="modal fade" id="savePipelineModal" tabindex="-1">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <h5 class="modal-title">保存流程模板</h5>
                <button type="button" class="btn-close" data-bs-dismiss="modal"></button>
            </div>
            <div class="modal-body">
                <div class="mb-3">
                    <label class="form-label">模板名称</label>
                    <input type="text" class="form-control" id="pipeline-name" placeholder="如: 标准单细胞流程">
                </div>
                <div class="mb-3">
                    <label class="form-label">描述（可选）</label>
                    <input type="text" class="form-control" id="pipeline-desc">
                </div>
            </div>
            <div class="modal-footer">
                <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">取消</button>
                <button type="button" class="btn btn-primary" id="btn-confirm-pipeline-save">保存</button>
            </div>
        </div>
    </div>
</div>
```

### Step 5: 流程模板保存 JS

在 `<script>` 中添加：

```javascript
document.getElementById('btn-confirm-pipeline-save').addEventListener('click', async () => {
    const name = document.getElementById('pipeline-name').value.trim();
    if (!name) { alert('请输入模板名称'); return; }
    const desc = document.getElementById('pipeline-desc').value.trim();

    // 按当前 DOM 顺序收集模块名
    const modules = [...document.querySelectorAll('#module-list .module-item')]
        .map(el => el.dataset.module);

    const body = {
        name, description: desc,
        analysis_type: ANALYSIS_TYPE,
        project_id: PROJECT_ID,
        scope: 'project',
        pipeline: { modules, skip_validation: false },
        params: collectCurrentParams(),
    };

    try {
        const resp = await fetch('/api/presets', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify(body),
        });
        if (resp.ok) {
            bootstrap.Modal.getInstance(document.getElementById('savePipelineModal')).hide();
            document.getElementById('pipeline-name').value = '';
            document.getElementById('pipeline-desc').value = '';
        }
    } catch (e) {
        alert('保存失败: ' + e.message);
    }
});
```

### Step 6: 验证

手动在浏览器中确认：
- 模块列表显示拖拽手柄 ☰
- 拖拽可调整模块顺序
- "保存为流程模板"按钮可用

### Step 7: Commit

```bash
git add templates/analysis_select.html
git commit -m "feat(ui): add draggable module list and pipeline template save"
```

---

## Task 7: Phase 3 — 迁移 2 个模块使用 get_plotly_layout()

**Files:**
- Modify: `modules/bulk_deg.py`
- Modify: `modules/proportion.py`

### Step 1: bulk_deg.py — 迁移火山图和 MA 图布局

在 `_run_single_comparison()` 函数中，找到火山图的 `fig_vol.update_layout(...)` 调用（约 line 128），替换为：

```python
    fig_vol.update_layout(
        **self.get_plotly_layout(title=f'火山图 ({group1} vs {group2})',
                                 xaxis_title='log2(Fold Change)', yaxis_title='-log10(padj)')
    )
```

但 `_run_single_comparison` 是模块级函数，不是类方法，没有 `self`。需要将 `get_plotly_layout` 作为参数传入，或者在函数外调用。

更简洁的方案：在 `_run_single_comparison` 的签名中添加 `viz_params=None` 参数，在 `run()` 调用时传入 `self.params.get('_visualization', {})`。

修改 `_run_single_comparison` 签名（line 39）：

```python
def _run_single_comparison(adata, counts, group1_samples, group2_samples, group1, group2,
                           method, fc_threshold, pval_threshold, top_n, gene_id_to_name,
                           plots_dir, results_dir, suffix='', viz_params=None):
```

在函数开头添加：

```python
    viz = viz_params or {}
    fig_width = viz.get('figure_width', 700)
    fig_height = viz.get('figure_height', 500)
    bg_color = viz.get('bg_color', 'white')
    font_family = viz.get('font_family', 'Arial')
    font_size = viz.get('font_size', 12)
```

将火山图的 `update_layout` 改为：

```python
    fig_vol.update_layout(
        title=f'火山图 ({group1} vs {group2})',
        xaxis_title='log2(Fold Change)', yaxis_title='-log10(padj)',
        plot_bgcolor=bg_color, width=fig_width, height=fig_height,
        font=dict(family=font_family, size=font_size),
    )
```

同样修改 MA 图的 `update_layout`。

在 `run()` 中调用 `_run_single_comparison` 时传入 `viz_params`：

```python
deg_df, files, _, _ = _run_single_comparison(
    adata, counts, g1_samples, g2_samples, g1, g2,
    method, fc_threshold, pval_threshold, top_n,
    gene_id_to_name, plots_dir, results_dir, suffix=str(idx),
    viz_params=self.params.get('_visualization', {}))
```

### Step 2: proportion.py — 迁移柱状图布局

在 `ProportionAnalysis.run()` 中，找到主柱状图的 `fig.update_layout(...)` 调用，替换为：

```python
        fig.update_layout(**self.get_plotly_layout(
            title='Cell Proportions by Group',
            barmode='stack',
            xaxis_title=batch_key, yaxis_title='Proportion',
        ))
```

同样修改子比较图的 `fig_sub.update_layout(...)`：

```python
                fig_sub.update_layout(**self.get_plotly_layout(
                    title=f'Cell Proportions: {group_a} vs {group_b} (p={pval_sub:.4f})',
                    barmode='stack',
                    xaxis_title=batch_key, yaxis_title='Proportion',
                ))
```

### Step 3: 验证

```bash
cd /data/GJ/platform && python -c "
import ast
ast.parse(open('modules/bulk_deg.py').read())
ast.parse(open('modules/proportion.py').read())
print('bulk_deg.py + proportion.py syntax OK')
"
```

### Step 4: Commit

```bash
git add modules/bulk_deg.py modules/proportion.py
git commit -m "feat(viz): migrate bulk_deg and proportion to use get_plotly_layout()"
```

---

## Task 8: Phase 4 — 在 2 个模块中集成 apply_filters()

**Files:**
- Modify: `modules/qc.py`
- Modify: `modules/bulk_qc.py`

### Step 1: qc.py — 在数据加载后调用 apply_filters

在 `QCAnalysis.run()` 方法中，找到数据加载完成、开始 QC 过滤之前的位置（`self.progress(15, ...)` 之前），添加：

```python
        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'qc')
```

### Step 2: bulk_qc.py — 同样集成

在 `BulkQCAnalysis.run()` 方法中，数据加载完成后、QC 过滤之前添加：

```python
        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'bulk_qc')
```

### Step 3: 前端添加过滤规则 UI

在 `analysis_select.html` 的参数循环 `{% endfor %}` 之后、提交按钮之前，添加过滤规则折叠区：

```html
<!-- 自定义过滤规则 -->
<div class="card p-3 mb-3" id="filter-rules-section">
    <div class="d-flex justify-content-between align-items-center" data-bs-toggle="collapse"
         data-bs-target="#filter-rules-body" style="cursor:pointer;">
        <h6 class="mb-0">自定义过滤规则（可选）</h6>
        <span class="text-muted small">▾</span>
    </div>
    <div class="collapse mt-2" id="filter-rules-body">
        <div id="filter-rules-list"></div>
        <button type="button" class="btn btn-outline-primary btn-sm mt-2" id="btn-add-filter">+ 添加规则</button>
    </div>
</div>
```

在 `<script>` 中添加过滤规则交互逻辑：

```javascript
// === 自定义过滤规则 ===
const FILTER_OPS = [
    {value: '>=', label: '>='},
    {value: '<=', label: '<='},
    {value: '==', label: '=='},
    {value: '!=', label: '!='},
    {value: 'in', label: 'in (列表)'},
    {value: 'not_in', label: 'not in (列表)'},
    {value: 'between', label: 'between (范围)'},
];

function addFilterRule(column = '', op = '>=', value = '') {
    const list = document.getElementById('filter-rules-list');
    const row = document.createElement('div');
    row.className = 'd-flex gap-2 align-items-center mb-2 filter-rule';

    const colInput = document.createElement('input');
    colInput.type = 'text';
    colInput.className = 'form-control form-control-sm';
    colInput.style.width = '160px';
    colInput.placeholder = '列名';
    colInput.value = column;
    colInput.name = '_filter_col';

    const opSelect = document.createElement('select');
    opSelect.className = 'form-select form-select-sm';
    opSelect.style.width = '130px';
    opSelect.name = '_filter_op';
    FILTER_OPS.forEach(o => {
        const opt = document.createElement('option');
        opt.value = o.value;
        opt.textContent = o.label;
        if (o.value === op) opt.selected = true;
        opSelect.appendChild(opt);
    });

    const valInput = document.createElement('input');
    valInput.type = 'text';
    valInput.className = 'form-control form-control-sm';
    valInput.style.width = '140px';
    valInput.placeholder = '值';
    valInput.value = typeof value === 'object' ? JSON.stringify(value) : value;
    valInput.name = '_filter_val';

    const delBtn = document.createElement('button');
    delBtn.type = 'button';
    delBtn.className = 'btn btn-outline-danger btn-sm';
    delBtn.textContent = '×';
    delBtn.addEventListener('click', () => row.remove());

    row.append(colInput, opSelect, valInput, delBtn);
    list.appendChild(row);
}

document.getElementById('btn-add-filter').addEventListener('click', () => addFilterRule());

function collectFilterRules() {
    const rules = [];
    document.querySelectorAll('.filter-rule').forEach(row => {
        const col = row.querySelector('[name="_filter_col"]').value.trim();
        const op = row.querySelector('[name="_filter_op"]').value;
        let valRaw = row.querySelector('[name="_filter_val"]').value.trim();
        if (!col) return;
        let val;
        if (op === 'in' || op === 'not_in') {
            val = valRaw.split(',').map(s => s.trim()).filter(Boolean);
        } else if (op === 'between') {
            val = valRaw.split(',').map(s => parseFloat(s.trim()));
        } else {
            val = isNaN(Number(valRaw)) ? valRaw : Number(valRaw);
        }
        rules.push({column: col, op, value: val});
    });
    return rules;
}

// 在表单提交时注入 filters
document.querySelector('.param-form').addEventListener('submit', function() {
    const filters = collectFilterRules();
    if (filters.length > 0) {
        let hidden = this.querySelector('input[name="_filters"]');
        if (!hidden) {
            hidden = document.createElement('input');
            hidden.type = 'hidden';
            hidden.name = '_filters';
            this.appendChild(hidden);
        }
        hidden.value = JSON.stringify({[MODULE_NAME]: filters});
    }
});
```

### Step 4: 验证

```bash
cd /data/GJ/platform && python -c "
import ast
ast.parse(open('modules/qc.py').read())
ast.parse(open('modules/bulk_qc.py').read())
print('qc.py + bulk_qc.py syntax OK')
"
```

### Step 5: Commit

```bash
git add modules/qc.py modules/bulk_qc.py templates/analysis_select.html
git commit -m "feat(filters): add apply_filters() to qc modules + filter rule UI"
```

---

## Task 9: 最终集成验证

### Step 1: 全量语法检查

```bash
cd /data/GJ/platform && python -c "
import ast, os
for root, dirs, files in os.walk('.'):
    if 'node_modules' in root or '__pycache__' in root:
        continue
    for f in files:
        if f.endswith('.py'):
            path = os.path.join(root, f)
            try:
                ast.parse(open(path).read())
            except SyntaxError as e:
                print(f'SYNTAX ERROR: {path}: {e}')
print('All .py files syntax OK')
"
```

### Step 2: 验证依赖约束

```bash
cd /data/GJ/platform && python -c "
from modules import validate_pipeline_order, PIPELINE_DEPS, MODULE_REGISTRY

# 验证所有注册模块都有依赖定义或无依赖
for mod in MODULE_REGISTRY:
    deps = PIPELINE_DEPS.get(mod, [])
    for d in deps:
        assert d in MODULE_REGISTRY, f'{mod} depends on unknown module {d}'

# 标准单细胞流程
ok, _ = validate_pipeline_order(['qc', 'preprocess', 'dimred', 'clustering', 'annotation'])
assert ok

# 标准 bulk 流程
ok, _ = validate_pipeline_order(['bulk_qc', 'bulk_normalize', 'bulk_deg', 'bulk_pca', 'bulk_heatmap', 'bulk_enrichment'])
assert ok

print('Integration validation OK')
"
```

### Step 3: 验证预设 API 端点存在

```bash
cd /data/GJ/platform && python -c "
from app import app
rules = [rule.rule for rule in app.url_map.iter_rules()]
assert '/api/presets' in rules or any('presets' in r for r in rules), 'presets endpoint not found'
print('API endpoints OK')
"
```

### Step 4: 最终 Commit

```bash
git add -A
git commit -m "feat: platform customization — presets, pipeline, viz, filters (Phase 1-4)"
```
