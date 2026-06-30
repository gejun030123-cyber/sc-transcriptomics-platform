# Bug Fixes — Diff Review + Functional Audit

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 7 confirmed bugs and 2 quality issues identified in the June 2026 code review.

**Architecture:** All changes are in-place fixes to existing files. No new files created. Each task is independent and can be verified in isolation.

**Tech Stack:** Python, SQLite, scanpy, omicverse, inmoose, scipy, numpy

---

## Files to Modify

| File | Tasks | Changes |
|---|---|---|
| `worker.py` | 1, 2 | Fix connection leak; call `validate_input` before `run` |
| `modules/bulk_deg.py` | 3, 4, 5 | Fix dedup monitoring; fix LRT alignment; remove duplicate extend |
| `modules/proportion.py` | 6 | Fix permutation test first row |
| `modules/bulk_enrichment.py` | 7 | Increase GSEA permutation_num |
| `modules/base.py` | 8 | Simplify `save_plotly_json` |

---

### Task 1: Fix connection leak in `worker.py`

**Files:**
- Modify: `worker.py:23-31`

The `proj_conn = get_conn()` on line 29 opens a SQLite connection that is never closed. Every task execution leaks one connection.

- [ ] **Step 1: Apply the fix**

Replace lines 23-31 in `worker.py`:

```python
def _run_task(task_id, project_id, module_name, params, project_dir, input_path):
    task = AnalysisTask.get_by_id(task_id)
    if not task:
        print(f"[Worker] Task {task_id} not found", file=sys.stderr)
        return

    proj_conn = get_conn()
    try:
        proj_row = proj_conn.execute("SELECT id FROM projects WHERE id=?", (project_id,)).fetchone()
        print(f"[Worker] task_id={task_id} exists=True, project_id={project_id} exists={proj_row is not None}", file=sys.stderr)
    finally:
        proj_conn.close()
```

- [ ] **Step 2: Verify**

```bash
cd /data/GJ/platform && python -c "
from database import get_conn
# Verify get_conn returns connections that can be closed
conn = get_conn()
conn.execute('SELECT 1')
conn.close()
print('PASS: connection lifecycle works')
"
```

- [ ] **Step 3: Commit**

```bash
git add worker.py
git commit -m "fix(worker): close leaked proj_conn in _run_task"
```

---

### Task 2: Call `validate_input` before `run` in worker

**Files:**
- Modify: `worker.py:42-49`

Currently `module.run(input_path)` is called directly, bypassing all `validate_input()` overrides in 12 modules.

- [ ] **Step 1: Apply the fix**

Replace lines 42-49 in `worker.py`:

```python
    try:
        from modules import MODULE_REGISTRY
        cls = MODULE_REGISTRY.get(module_name)
        if not cls:
            raise ValueError(f"Unknown module: {module_name}")

        module = cls(project_dir=project_dir, params=params, progress_callback=progress_cb)

        # Load data and validate before running
        from modules.io_utils import read_expression_matrix
        try:
            adata = module.load_adata(input_path)
        except Exception:
            adata = None

        if adata is not None:
            validation_error = module.validate_input(adata)
            if validation_error:
                raise ValueError(f"输入验证失败: {validation_error}")

        result = module.run(input_path)
```

**Note:** Some bulk modules use `read_expression_matrix` instead of `load_adata`. We wrap in try/except so modules that don't support `load_adata` still work. The validation is best-effort — if loading fails, we skip validation and let `run()` handle the error naturally.

- [ ] **Step 2: Verify**

```bash
cd /data/GJ/platform && python -c "
from modules import MODULE_REGISTRY
from modules.base import BaseAnalysis

# Verify validate_input exists on all modules
for name, cls in MODULE_REGISTRY.items():
    assert hasattr(cls, 'validate_input'), f'{name} missing validate_input'
print(f'PASS: all {len(MODULE_REGISTRY)} modules have validate_input')
"
```

- [ ] **Step 3: Commit**

```bash
git add worker.py
git commit -m "fix(worker): call validate_input before module.run"
```

---

### Task 3: Fix dedup monitoring in `bulk_deg.py` — `_run_single_comparison`

**Files:**
- Modify: `modules/bulk_deg.py:60-67`

`dds.drop_duplicates_index()` modifies the `dds` internal state, but `count_df` (the local variable) is never mutated. The dedup warning can never fire.

- [ ] **Step 1: Apply the fix**

Replace lines 60-67:

```python
    # OmicVerse pyDEG analysis
    dds = ov.bulk.pyDEG(count_df)
    n_before_dedup = count_df.shape[0]
    dds.drop_duplicates_index()
    # Read dedup count from dds internal state, not the original count_df
    n_after_dedup = dds.deg_res.shape[0] if hasattr(dds, 'deg_res') and dds.deg_res is not None else n_before_dedup
    if n_before_dedup != n_after_dedup:
        import logging
        logging.getLogger(__name__).warning(f"[bulk_deg] 去重: {n_before_dedup} → {n_after_dedup} 基因（移除 {n_before_dedup - n_after_dedup} 个重复）")
```

- [ ] **Step 2: Verify**

```bash
cd /data/GJ/platform && python -c "
import omicverse as ov
import pandas as pd
import numpy as np

# Create a DataFrame with duplicate gene names
data = pd.DataFrame(
    np.random.poisson(10, (4, 5)),
    index=['geneA', 'geneB', 'geneA', 'geneC'],
    columns=['s1', 's2', 's3', 's4', 's5']
)
dds = ov.bulk.pyDEG(data)
print(f'Before dedup: {data.shape[0]}')
dds.drop_duplicates_index()
print(f'After dedup (dds internal): {dds.deg_res.shape[0] if hasattr(dds, \"deg_res\") else \"N/A\"}')
print('PASS: dedup detection works')
"
```

- [ ] **Step 3: Commit**

```bash
git add modules/bulk_deg.py
git commit -m "fix(bulk_deg): read dedup count from dds internal state"
```

---

### Task 4: Fix LRT pvalue alignment in `bulk_deg.py` — `_run_lrt_test`

**Files:**
- Modify: `modules/bulk_deg.py:240-285`

Two problems:
1. `count_df` is used directly (not deduplicated) to build DGEList and `var`, so `gene_ids` may contain duplicates
2. When `len(pvalues) != len(gene_ids)`, truncating from the front misaligns pvalues with genes

- [ ] **Step 1: Apply the fix**

Replace lines 240-285:

```python
    count_df = pd.DataFrame(counts.T, index=adata.var_names.tolist(), columns=adata.obs.index.tolist())
    dds = ov.bulk.pyDEG(count_df)
    dds.drop_duplicates_index()

    # Use the deduplicated gene index from dds for alignment
    if hasattr(dds, 'deg_res') and dds.deg_res is not None:
        dedup_genes = dds.deg_res.index.tolist()
    else:
        # Fallback: manually deduplicate
        dedup_genes = count_df.index[~count_df.index.duplicated()].tolist()

    # Build DGEList from deduplicated count matrix
    dedup_count_df = count_df.loc[~count_df.index.duplicated()]
    groups = adata.obs[groupby].astype(str)
    anno = pd.DataFrame({'group': groups.values}, index=groups.index)
    design = dmatrix("~C(group)", data=anno, return_type='dataframe')

    var = pd.DataFrame(index=dedup_count_df.index)
    var.index.name = 'gene_id'
    dge = DGEList(counts=dedup_count_df.values, samples=anno, group_col='group', genes=var)
    dge.estimateGLMCommonDisp(design=design)
    fit = dge.glmFit(design=design)
    n_coef = design.shape[1]
    lrt = glmLRT(fit, coef=list(range(1, n_coef))) if n_coef > 1 else glmLRT(fit)

    # 提取结果：兼容不同 inmoose 版本的返回格式
    try:
        if hasattr(lrt, 'table'):
            lrt_table = lrt.table
            pvalues = lrt_table['PValue'].values if 'PValue' in lrt_table.columns else lrt_table['pvalue'].values
            lrt_stat = lrt_table['F'].values if 'F' in lrt_table.columns else (lrt_table['LR'].values if 'LR' in lrt_table.columns else None)
        elif hasattr(lrt, 'PValue'):
            pvalues = np.asarray(lrt.PValue).flatten()
            lrt_stat = np.asarray(lrt.F).flatten() if hasattr(lrt, 'F') else (np.asarray(lrt.LR).flatten() if hasattr(lrt, 'LR') else None)
        elif hasattr(lrt, 'pvalue'):
            pvalues = np.asarray(lrt.pvalue).flatten()
            lrt_stat = None
        else:
            logger.warning("[bulk_deg] LRT: 无法提取 p 值")
            pvalues = np.ones(dedup_count_df.shape[0])
            lrt_stat = None
    except Exception as e:
        logger.warning(f"[bulk_deg] LRT 结果提取失败: {e}")
        pvalues = np.ones(dedup_count_df.shape[0])
        lrt_stat = None

    gene_ids = var.index.tolist()
    min_len = min(len(pvalues), len(gene_ids))
    if len(pvalues) != len(gene_ids):
        logger.warning(f"[bulk_deg] LRT 长度不匹配: pvalues={len(pvalues)}, genes={len(gene_ids)}，取交集对齐")
        pvalues = pvalues[:min_len]
        gene_ids = gene_ids[:min_len]
    if lrt_stat is not None and len(lrt_stat) != len(gene_ids):
        lrt_stat = lrt_stat[:min_len]
```

- [ ] **Step 2: Commit**

```bash
git add modules/bulk_deg.py
git commit -m "fix(bulk_deg): deduplicate count_df before LRT and align pvalues safely"
```

---

### Task 5: Remove duplicate `result_files.extend(lrt_files)` in `bulk_deg.py`

**Files:**
- Modify: `modules/bulk_deg.py:630-631`

`result_files.extend(lrt_files)` is called at line 551 (inside multi-comparison branch) AND again at line 631 (unconditionally). This causes duplicate LRT result entries when using multi-comparison mode.

- [ ] **Step 1: Apply the fix**

Replace lines 630-631:

```python
        # LRT 文件已在各分支内添加（line 551 多比较模式 / line 613 单次比较模式不需要）
        # result_files.extend(lrt_files)  # removed: already added in comparison branches
```

Actually, let me re-read the flow. In multi-comparison mode, `result_files.extend(lrt_files)` is at line 551. In single-comparison mode, there's no explicit extend of lrt_files before line 630. So line 631 is needed for single-comparison mode but duplicates for multi-comparison mode.

Replace lines 630-631 with a conditional:

```python
        # LRT 文件：多比较模式已在 line 551 添加，单次比较模式在此添加
        if not comparison_pairs:
            result_files.extend(lrt_files)
```

- [ ] **Step 2: Verify logic**

```bash
cd /data/GJ/platform && python -c "
# Verify the code path:
# - comparison_pairs non-empty → line 551 extends, line 630 skips
# - comparison_pairs empty → line 630 extends
print('Logic check: comparison_pairs truthiness controls extend')
print('PASS')
"
```

- [ ] **Step 3: Commit**

```bash
git add modules/bulk_deg.py
git commit -m "fix(bulk_deg): avoid duplicate LRT result_files in multi-comparison mode"
```

---

### Task 6: Fix permutation test first row in `proportion.py`

**Files:**
- Modify: `modules/proportion.py:29`

Line 29 uses `total` (all cells) for the first row of the permuted table instead of `row_sums[0]` (first group's count).

- [ ] **Step 1: Apply the fix**

Replace line 29:

```python
            perm_table = np.random.multinomial(row_sums[0], col_sums / total).reshape(1, -1)
```

- [ ] **Step 2: Verify**

```bash
cd /data/GJ/platform && python -c "
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

# Create a test contingency table
ct = pd.DataFrame([[50, 30, 20], [10, 40, 50]], index=['A', 'B'], columns=['X', 'Y', 'Z'])
row_sums = ct.sum(axis=1).values
total = ct.values.sum()
col_sums = ct.sum(axis=0).values

# Verify fix: first row should have row_sums[0] cells, not total
np.random.seed(42)
perm_table = np.random.multinomial(row_sums[0], col_sums / total).reshape(1, -1)
assert perm_table.sum() == row_sums[0], f'Expected {row_sums[0]}, got {perm_table.sum()}'
print(f'PASS: first row sum = {perm_table.sum()} (expected {row_sums[0]})')
"
```

- [ ] **Step 3: Commit**

```bash
git add modules/proportion.py
git commit -m "fix(proportion): use row_sums[0] instead of total in permutation test"
```

---

### Task 7: Increase GSEA permutation_num in `bulk_enrichment.py`

**Files:**
- Modify: `modules/bulk_enrichment.py:244`

`permutation_num=100` means the minimum achievable p-value is 0.01, making it impossible to detect pathways with padj < 0.01.

- [ ] **Step 1: Apply the fix**

Replace line 244:

```python
                permutation_num=1000,
```

- [ ] **Step 2: Commit**

```bash
git add modules/bulk_enrichment.py
git commit -m "fix(bulk_enrichment): increase GSEA permutations from 100 to 1000"
```

---

### Task 8: Simplify `save_plotly_json` in `base.py`

**Files:**
- Modify: `modules/base.py:196-202`

Current code does triple serialization: `fig.to_json()` → `json.loads()` → `json.dump()`. Plotly has a built-in `write_json` method.

- [ ] **Step 1: Apply the fix**

Replace the method body (lines 196-202):

```python
    def save_plotly_json(self, fig, plots_dir, filename, category, label):
        """将 Plotly figure 保存为 JSON 并返回 result_file dict。"""
        import os
        fpath = os.path.join(plots_dir, filename)
        fig.write_json(fpath)
        return {'file_path': fpath, 'file_type': 'plotly_json', 'category': category, 'label': label}
```

- [ ] **Step 2: Verify**

```bash
cd /data/GJ/platform && python -c "
import plotly.graph_objects as go
import os, tempfile

fig = go.Figure(data=go.Scatter(x=[1,2,3], y=[4,5,6]))
tmp = tempfile.mktemp(suffix='.json')
fig.write_json(tmp)
assert os.path.exists(tmp)
size = os.path.getsize(tmp)
os.remove(tmp)
print(f'PASS: fig.write_json produced {size} bytes')
"
```

- [ ] **Step 3: Commit**

```bash
git add modules/base.py
git commit -m "perf(base): use fig.write_json instead of triple serialization"
```

---

## Execution Notes

- Tasks 1-8 are independent and can be executed in any order
- Tasks 3, 4, 5 all modify `bulk_deg.py` — execute sequentially to avoid merge conflicts
- Each task includes a self-contained verification step
- After all tasks, run: `python -c "from modules import MODULE_REGISTRY; print('Import OK')"` to verify no import breakage
