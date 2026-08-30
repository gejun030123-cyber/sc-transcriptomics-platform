# tests/test_theme_focus_tools.py
"""主题驱动富集（P0）测试：theme_lexicon、AI 只读工具、sc_cell_go focus 回归。

对应 docs/ai-theme-driven-analysis-plan.md 第 5、6 节验收标准。
"""
import json

from config import Config


def _write_test_gene_sets(data_dir, library="Test_Gene_Set_2023"):
    """在平台管理的 go_gene_sets 目录下伪造一个最小 GMT 库。"""
    go_dir = data_dir / "go_gene_sets"
    go_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "fatty acid metabolic process (GO:0035336)\t\tFADS1\tACOX1\tCPT1A\tACADL\tHADHB",
        "lipid transport (GO:0006869)\t\tAPOE\tAPOB\tAPOC2\tLCAT",
        "inflammatory response (GO:0006954)\t\tIL1B\tIL6\tTNF\tPTGS2\tNFKB1",
        "chromosome segregation (GO:0007059)\t\tCENPE\tBUB1\tMAD2L1",
    ]
    path = go_dir / f"{library}.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return library, path


# ---------------------------------------------------------------------------
# theme_lexicon
# ---------------------------------------------------------------------------

class TestThemeLexicon:

    def test_list_themes_contains_expected_keys(self):
        from modules import theme_lexicon

        themes = {item["key"] for item in theme_lexicon.list_themes()}
        assert {"lipid_metabolism", "inflammation", "immune_response", "emt"} <= themes
        for item in theme_lexicon.list_themes():
            assert item["label"]
            assert item["keywords"]
            assert item["keywords_en"]

    def test_resolve_themes_chinese_mixed(self):
        from modules import theme_lexicon

        assert theme_lexicon.resolve_themes("脂代谢和炎症") == [
            "lipid_metabolism", "inflammation",
        ]
        assert theme_lexicon.resolve_themes("hypoxia related") == ["hypoxia"]
        assert theme_lexicon.resolve_themes("无关内容") == []

    def test_search_terms_matches_and_sorts(self, test_project):
        from modules import theme_lexicon

        library, _ = _write_test_gene_sets(
            _data_dir_for(test_project), "Test_Gene_Set_2023",
        )
        outcome = theme_lexicon.search_terms(
            "脂代谢", libraries=[library], limit=10,
        )
        assert outcome["themes"] == ["lipid_metabolism"]
        terms = outcome["terms"]
        matched_names = {item["term"] for item in terms}
        assert "fatty acid metabolic process (GO:0035336)" in matched_names
        assert "lipid transport (GO:0006869)" in matched_names
        assert all("chromosome segregation" not in name for name in matched_names)
        first = terms[0]
        assert first["library"] == library
        assert first["n_genes"] >= 4
        assert first["example_genes"]
        assert first["matched_keywords"]

    def test_load_library_terms_uses_managed_snapshot(self, test_project):
        from modules import theme_lexicon

        library, _ = _write_test_gene_sets(
            _data_dir_for(test_project), "Test_Gene_Set_2023",
        )
        terms = theme_lexicon.load_library_terms(library)
        assert "lipid transport (GO:0006869)" in terms

    def test_search_terms_two_themes_combined(self, test_project):
        from modules import theme_lexicon

        library, _ = _write_test_gene_sets(
            _data_dir_for(test_project), "Test_Gene_Set_2023",
        )
        outcome = theme_lexicon.search_terms(
            "脂代谢和炎症", libraries=[library],
        )
        assert set(outcome["themes"]) == {"lipid_metabolism", "inflammation"}
        matched_names = {item["term"] for item in outcome["terms"]}
        assert "inflammatory response (GO:0006954)" in matched_names

    def test_search_terms_missing_library_recorded_not_fatal(self, test_project):
        from modules import theme_lexicon

        _write_test_gene_sets(_data_dir_for(test_project), "Test_Gene_Set_2023")
        outcome = theme_lexicon.search_terms(
            "inflammation", libraries=["Test_Gene_Set_2023", "Absent_Library"],
        )
        assert "Absent_Library" in outcome["library_errors"]
        assert outcome["terms"]

    def test_search_terms_invalid_library_name_rejected(self, test_project):
        import pytest
        from modules import theme_lexicon

        with pytest.raises(ValueError):
            theme_lexicon.search_terms("inflammation", libraries=["../etc"])

    def test_search_terms_no_keywords_returns_empty(self, test_project):
        from modules import theme_lexicon

        _write_test_gene_sets(_data_dir_for(test_project), "Test_Gene_Set_2023")
        outcome = theme_lexicon.search_terms("完全无关的主题词")
        assert outcome["n_terms"] == 0
        assert outcome["terms"] == []


# ---------------------------------------------------------------------------
# AI 工具：search_pathway_terms / read_task_table
# ---------------------------------------------------------------------------

class TestSearchPathwayTermsTool:

    def test_list_themes_mode(self):
        from modules.ai_tools import execute_tool

        result = execute_tool("search_pathway_terms", {"list_themes": True}, "any")
        keys = {item["key"] for item in result["themes"]}
        assert "lipid_metabolism" in keys

    def test_requires_query(self):
        from modules.ai_tools import execute_tool

        result = execute_tool("search_pathway_terms", {}, "any")
        assert "error" in result

    def test_search_via_tool(self, test_project):
        from modules.ai_tools import execute_tool

        library, _ = _write_test_gene_sets(
            _data_dir_for(test_project), "Test_Gene_Set_2023",
        )
        result = execute_tool(
            "search_pathway_terms",
            {"query": "脂代谢和炎症", "libraries": [library]},
            test_project,
        )
        assert result["themes"]
        assert result["terms"]
        assert "focus_terms" in result["usage_note"]


class TestReadTaskTableTool:

    def _make_result_file(self, test_project, rows=None):
        import pandas as pd
        from models import AnalysisTask, ResultFile

        task = AnalysisTask(
            project_id=test_project, module_name="sc_cell_go", params_json="{}",
        )
        task.save()
        frame = pd.DataFrame(rows if rows is not None else [
            {"Term": "lipid transport (GO:0006869)", "Adjusted P-value": 0.001,
             "Genes": "APOE;APOB"},
            {"Term": "inflammatory response (GO:0006954)", "Adjusted P-value": 0.04,
             "Genes": "IL1B;IL6"},
            {"Term": "chromosome segregation (GO:0007059)", "Adjusted P-value": 0.6,
             "Genes": "CENPE"},
        ])
        path = Config.results_dir(test_project) + "/fake_enrichment.csv"
        import os
        os.makedirs(Config.results_dir(test_project), exist_ok=True)
        frame.to_csv(path, index=False)
        record = ResultFile.create(
            task_id=task.id, project_id=test_project, file_type="csv",
            category="table", label="fake enrichment", file_path=path,
        )
        return task, record

    def test_basic_read(self, test_project):
        from modules.ai_tools import execute_tool

        _, record = self._make_result_file(test_project)
        result = execute_tool(
            "read_task_table", {"file_id": record.id}, test_project,
        )
        assert result["total_rows"] == 3
        assert result["returned_rows"] == 3
        assert {"Term", "Adjusted P-value"} <= set(result["columns"])

    def test_contains_and_fdr_filters(self, test_project):
        from modules.ai_tools import execute_tool

        _, record = self._make_result_file(test_project)
        result = execute_tool(
            "read_task_table",
            {"file_id": record.id, "contains": "lipid", "fdr_max": 0.05},
            test_project,
        )
        assert result["filtered_rows"] == 1
        assert "lipid transport" in result["rows"][0]["Term"]

    def test_fdr_filter_applies_before_column_selection(self, test_project):
        from modules.ai_tools import execute_tool

        _, record = self._make_result_file(test_project)
        result = execute_tool(
            "read_task_table",
            {"file_id": record.id, "columns": ["Term"], "fdr_max": 0.05},
            test_project,
        )
        assert result["filtered_rows"] == 2
        assert result["columns"] == ["Term"]
        assert all(set(row) == {"Term"} for row in result["rows"])

    def test_limit_truncation(self, test_project):
        from modules.ai_tools import execute_tool

        rows = [
            {"Term": f"term {i}", "Adjusted P-value": 0.01} for i in range(10)
        ]
        _, record = self._make_result_file(test_project, rows=rows)
        result = execute_tool(
            "read_task_table", {"file_id": record.id, "limit": 4}, test_project,
        )
        assert result["returned_rows"] == 4
        assert result["truncated"] is True

    def test_rejects_foreign_project_file(self, test_project):
        from modules.ai_tools import execute_tool

        _, record = self._make_result_file(test_project)
        result = execute_tool(
            "read_task_table", {"file_id": record.id}, "other_project",
        )
        assert "error" in result

    def test_rejects_unknown_file(self, test_project):
        from modules.ai_tools import execute_tool

        result = execute_tool(
            "read_task_table", {"file_id": "missing"}, test_project,
        )
        assert "error" in result

    def test_rejects_outside_project_path(self, test_project, tmp_path):
        from modules.ai_tools import execute_tool

        _, record = self._make_result_file(test_project)
        outside_path = tmp_path / "outside.csv"
        outside_path.write_text("Term,Adjusted P-value\noutside,0.01\n", encoding="utf-8")
        record.file_path = str(outside_path)
        record.save()
        result = execute_tool("read_task_table", {"file_id": record.id}, test_project)
        assert "error" in result

    def test_rejects_symlinked_result_directory(self, test_project):
        import os
        from pathlib import Path

        from models import ResultFile
        from modules.ai_tools import execute_tool

        task, _ = self._make_result_file(test_project)
        results_dir = Path(Config.results_dir(test_project))
        real_dir = results_dir / "real"
        real_dir.mkdir(parents=True, exist_ok=True)
        target_path = real_dir / "table.csv"
        target_path.write_text("Term,Adjusted P-value\nlipid,0.01\n", encoding="utf-8")
        symlink_dir = results_dir / "linked"
        os.symlink(real_dir, symlink_dir)
        record = ResultFile.create(
            task_id=task.id, project_id=test_project, file_type="csv",
            category="table", label="symlinked table",
            file_path=str(symlink_dir / "table.csv"),
        )
        result = execute_tool("read_task_table", {"file_id": record.id}, test_project)
        assert "error" in result


# ---------------------------------------------------------------------------
# sc_cell_go focus 回归（helper 级；完整 run 级见 test_sc_de_contract.py）
# ---------------------------------------------------------------------------

class TestFocusTermHelpers:

    def test_parse_accepts_list_text_and_stringified_list(self):
        from modules.sc_cell_go import _parse_focus_terms

        assert _parse_focus_terms(["a", "b"]) == ["a", "b"]
        assert _parse_focus_terms("a;b, c\nd") == ["a", "b", "c", "d"]
        # _validate_analysis_params 会把 list 强制 str() 化
        assert _parse_focus_terms(str(["lipid transport (GO:0006869)", "x"])) == [
            "lipid transport (GO:0006869)", "x",
        ]
        assert _parse_focus_terms(str(["term, with comma", "x"])) == [
            "term, with comma", "x",
        ]
        assert _parse_focus_terms(None) == []
        assert _parse_focus_terms("") == []

    def test_parse_rejects_too_many_terms(self):
        import pytest
        from modules.sc_cell_go import _parse_focus_terms

        with pytest.raises(ValueError):
            _parse_focus_terms([f"term{i}" for i in range(61)])

    def test_mask_matches_full_and_display_names(self):
        import pandas as pd
        from modules.sc_cell_go import _focus_term_mask

        frame = pd.DataFrame({
            "Term": [
                "lipid transport (GO:0006869)",
                "inflammatory response (GO:0006954)",
                "chromosome segregation (GO:0007059)",
            ],
        })
        mask = _focus_term_mask(frame, [
            "Lipid Transport (GO:0006869)",       # 大小写不敏感全名
            "inflammatory response",               # 去掉 ID 的显示名
        ])
        assert mask.tolist() == [True, True, False]
        # 不是子串匹配：短关键词不应误聚焦
        mask2 = _focus_term_mask(frame, ["lipid"])
        assert mask2.tolist() == [False, False, False]

    def test_go_focus_rows_keep_only_full_library_fdr_matches(self):
        import pandas as pd
        from modules.sc_cell_go import _focus_significant_rows

        frame = pd.DataFrame({
            "status": ["completed", "completed", "completed", "completed"],
            "Significant": [True, True, False, True],
            "gene_set": [
                "GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                "GO_Molecular_Function_2023", "KEGG_2021_Human",
            ],
            "Term": ["lipid transport", "lipoprotein particle", "cholesterol binding", "lipid pathway"],
            "method": ["ORA"] * 4,
            "direction": ["Up", "Up", "Down", "Up"],
        })
        selected = _focus_significant_rows(
            frame, ["lipid transport", "lipoprotein particle", "cholesterol binding", "lipid pathway"],
        )
        assert selected["Term"].tolist() == ["lipid transport", "lipoprotein particle"]
        assert selected["Database"].tolist() == ["GO_BP", "GO_CC"]

    def test_focus_params_registered_in_schema(self):
        from modules.schemas import PARAM_SCHEMAS

        keys = {item["key"] for item in PARAM_SCHEMAS["sc_cell_go"]}
        assert "focus_terms" in keys

    def test_validate_analysis_params_keeps_focus_terms(self):
        from modules.ai_tools import _validate_analysis_params

        cleaned, error = _validate_analysis_params(
            "sc_cell_go",
            {"focus_terms": ["lipid transport (GO:0006869)"], "method": "ORA"},
        )
        assert error is None
        assert "focus_terms" in cleaned
        assert "method" in cleaned


def _data_dir_for(pid):
    """返回 test_project fixture 隔离后的 DATA_DIR path 对象。"""
    from pathlib import Path

    return Path(Config.DATA_DIR)
