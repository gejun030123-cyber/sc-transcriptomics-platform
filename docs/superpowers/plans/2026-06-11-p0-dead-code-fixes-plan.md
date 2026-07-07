
# P0: 死代码修复 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 4 dead-code parameters in SC modules — wire in use_mde, annotation method, DE reference, and trajectory method

**Architecture:** Minimal surgical changes to 4 module files + 1 schema cleanup. No new files, no interface changes.

**Tech Stack:** scanpy, omicverse, pymde (optional for MDE)

---

### Task 1: Wire `use_mde` in dimred.py

**Files:**
- Modify: `modules/dimred.py:34-35`
- Modify: `routes/analysis.py` (dimred schema, remove help text about unimplemented)

- [ ] **Step 1: Read current dimred.py and implement use_mde branch**

In `modules/dimred.py`, after line 23 add:
```python
use_mde = self.params.get('use_mde', False)
```

Replace lines 34-35:
```python
self.progress(70, "Computing UMAP...")
sc.tl.umap(adata)
```
With:
```python
self.progress(70, "Computing 2D embedding...")
if use_mde:
    try:
        import pymde
        mde = pymde.preserve_neighbors(adata.obsm['X_pca'], embedding_dim=2, device='cpu')
        embedding = mde.embed(verbose=False)
        adata.obsm['X_mde'] = embedding.cpu().numpy() if hasattr(embedding, 'cpu') else embedding.numpy()
        adata.obsm['X_umap'] = adata.obsm['X_mde']
    except ImportError:
        self.progress(71, "pymde not installed, falling back to UMAP...")
        sc.tl.umap(adata)
else:
    sc.tl.umap(adata)
```

Update summary to include embedding method used:
```python
'embedding_method': 'mde' if use_mde else 'umap',
```

- [ ] **Step 2: Commit**

```bash
git add modules/dimred.py
git commit -m "feat(dimred): wire use_mde parameter with pymde fallback"
```

---

### Task 2: Wire `method` in annotation.py (auto_marker vs manual)

**Files:**
- Modify: `modules/annotation.py:21-167` (run method)

- [ ] **Step 1: Read annotation.py fully and implement method branching**

At line ~76 (after loading data), add:
```python
method = self.params.get('method', 'auto_marker')
```

The existing scoring logic (lines ~84-115) stays as the `auto_marker` path. Wrap it in:
```python
if method == 'manual':
    # Manual mode: user provides ClusterID:CellType mapping in custom_markers
    manual_mapping = {}
    for line in custom_markers_str.split('\n'):
        line = line.strip()
        if not line or ':' not in line:
            continue
        cluster_id, cell_type = line.split(':', 1)
        manual_mapping[cluster_id.strip()] = cell_type.strip()
    
    if not manual_mapping:
        self.progress(50, "No manual mapping provided, falling back to auto_marker...")
        method = 'auto_marker'

if method == 'auto_marker':
    # ... existing scoring + assignment logic ...
```

For `manual` mode, after parsing the mapping:
```python
if method == 'manual':
    self.progress(50, "Applying manual cell type mapping...")
    adata.obs['celltype'] = adata.obs[leiden_key].astype(str).map(manual_mapping)
    adata.obs['celltype'] = adata.obs['celltype'].fillna('Unknown')
```

- [ ] **Step 2: Commit**

```bash
git add modules/annotation.py
git commit -m "feat(annotation): wire method param for manual cluster-to-celltype mapping"
```

---

### Task 3: Wire `reference` in deg.py

**Files:**
- Modify: `modules/deg.py:32`

- [ ] **Step 1: Wire reference into rank_genes_groups call**

At line ~27 (after extracting params), add:
```python
reference = self.params.get('reference', 'rest')
```

Replace line 32:
```python
sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=100)
```
With:
```python
ref_kwarg = {} if reference == 'rest' else {'reference': str(reference)}
sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=100, **ref_kwarg)
```

- [ ] **Step 2: Commit**

```bash
git add modules/deg.py
git commit -m "feat(deg): wire reference param for vs-specific-group DE"
```

---

### Task 4: Clean trajectory.py — remove slingshot option, wire cluster_key

**Files:**
- Modify: `modules/trajectory.py:25-32`
- Modify: `routes/analysis.py` (trajectory schema)

- [ ] **Step 1: Clean trajectory.py**

The `method` param is read but only `diffusion_map` code exists. The `cluster_key` is read but never used. Changes:

1. Remove `method` extraction (line 25) — keep only diffusion_map
2. Use `cluster_key` for coloring plots (currently hardcoded to 'leiden' in plot sections)

Replace lines 25-26:
```python
method = self.params.get('method', 'diffusion_map')
cluster_key = self.params.get('cluster_key', 'leiden')
```
With:
```python
cluster_key = self.params.get('cluster_key', 'leiden')
```

In the plot generation sections, replace any hardcoded `'leiden'` color references with `cluster_key`.

Update progress message at line 28 to remove method reference:
```python
self.progress(20, "Computing diffusion map...")
```

- [ ] **Step 2: Update trajectory schema in analysis.py**

In `routes/analysis.py`, update the trajectory schema to remove `slingshot` option:
```python
'trajectory': [
    {'key': 'cluster_key', 'label': '聚类列名', 'type': 'text', 'default': 'leiden', 'help': '用于轨迹推断的聚类列名。'},
    {'key': 'plot_genes', 'label': '拟时序基因表达（可选）', 'type': 'textarea', 'default': '',
     'help': '手动输入基因名，逗号或换行分隔。生成这些基因沿拟时序的表达曲线图。最多 10 个基因。'},
],
```

- [ ] **Step 3: Commit**

```bash
git add modules/trajectory.py routes/analysis.py
git commit -m "fix(trajectory): remove dead slingshot option, wire cluster_key for plots"
```

---

### Task 5: Verify all fixes work

- [ ] **Step 1: Run existing tests**

```bash
cd /data/GJ/platform && python -m pytest tests/ -v
```

- [ ] **Step 2: Manual smoke test — start app and verify each module's param form renders**

```bash
cd /data/GJ/platform && python app.py &
# Visit each module config page and verify params render correctly
```

- [ ] **Step 3: Final commit if any test fixes needed**
