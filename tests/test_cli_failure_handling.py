"""
Unit and regression tests for CLI failure handling in tf_explorer.cli.
Verifies that complete analysis failure produces a non-zero exit status,
while valid batch processing with partial successes continues and exits cleanly.
"""

import sys
import os
import pytest
from unittest.mock import patch, MagicMock

from tf_explorer import cli


def test_cli_complete_failure_single_gene_exits_nonzero(monkeypatch, tmp_path):
    """When a single requested gene analysis raises an exception, CLI must exit with status 1."""
    test_args = ["tf-explorer", "--gene", "NONEXISTENT_GENE_FAIL", "--out", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", test_args)

    with patch("tf_explorer.analysis.analyze_gene", side_effect=RuntimeError("Gene not found")):
        with pytest.raises(SystemExit) as exc_info:
            cli.main()
        assert exc_info.value.code == 1


def test_cli_complete_failure_batch_genes_exits_nonzero(monkeypatch, tmp_path):
    """When all genes in a batch file fail analysis, CLI must exit with status 1."""
    batch_file = tmp_path / "genes.txt"
    batch_file.write_text("GENE1\nGENE2\nGENE3\n", encoding="utf-8")

    test_args = ["tf-explorer", "--batch-genes", str(batch_file), "--out", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", test_args)

    with patch("tf_explorer.analysis.analyze_gene", side_effect=RuntimeError("Database down")):
        with pytest.raises(SystemExit) as exc_info:
            cli.main()
        assert exc_info.value.code == 1


def test_cli_partial_failure_batch_genes_preserves_batch_and_exits_zero(monkeypatch, tmp_path):
    """When some genes succeed and some fail in batch mode, CLI continues and exits 0 with warnings."""
    batch_file = tmp_path / "genes.txt"
    batch_file.write_text("GOOD_GENE1\nBAD_GENE\nGOOD_GENE2\n", encoding="utf-8")

    test_args = ["tf-explorer", "--batch-genes", str(batch_file), "--out", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", test_args)

    def mock_analyze(gene_name, **kwargs):
        if gene_name == "BAD_GENE":
            raise ValueError("Invalid gene")
        return None

    with patch("tf_explorer.analysis.analyze_gene", side_effect=mock_analyze) as mock_fn:
        # Should not raise SystemExit
        cli.main()
        assert mock_fn.call_count == 3


def test_cli_success_single_gene_exits_zero(monkeypatch, tmp_path):
    """When gene analysis succeeds, CLI finishes with status 0."""
    test_args = ["tf-explorer", "--gene", "TP53", "--out", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", test_args)

    with patch("tf_explorer.analysis.analyze_gene", return_value=None) as mock_fn:
        cli.main()
        assert mock_fn.call_count == 1


def test_cli_missing_batch_file_exits_nonzero(monkeypatch, tmp_path):
    """When the specified batch file does not exist, CLI exits with status 1."""
    missing_file = tmp_path / "does_not_exist.txt"
    test_args = ["tf-explorer", "--batch-genes", str(missing_file), "--out", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", test_args)

    with pytest.raises(SystemExit) as exc_info:
        cli.main()
    assert exc_info.value.code == 1
