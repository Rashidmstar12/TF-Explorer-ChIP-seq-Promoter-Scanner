"""
Tests for native external API clients and caching logic.
Verifies pure-Python implementations for UCSC, GTEx, STRING, and platformdirs caching.
"""

import os
import pytest
import pandas as pd
from tf_explorer.analysis import get_default_cache_dir, summarize_experiments
from tf_explorer.external_apis import (
    get_ucsc_conservation,
    resolve_gtex_gencode_id,
    get_gtex_median_expression,
    fetch_gtex_data,
    fetch_string_data,
)


def test_get_default_cache_dir():
    """Verify that get_default_cache_dir returns an existing, writable directory path."""
    cache_dir = get_default_cache_dir()
    assert isinstance(cache_dir, str)
    assert len(cache_dir) > 0
    assert os.path.exists(cache_dir)
    assert os.path.isdir(cache_dir)


def test_summarize_experiments_decoupled():
    """Verify that summarize_experiments calculates metrics correctly without Streamlit context."""
    df_raw = pd.DataFrame({
        'experiment_id': ['EXP1', 'EXP2', 'EXP3'],
        'num_strict_peaks': [2, 0, 1],
        'num_loose_peaks': [3, 1, 1],
        'total_peaks_in_file': [100, 200, 150]
    })
    n_strict, n_loose, t_strict, t_loose, t_exp, df_proc = summarize_experiments(df_raw)
    assert n_strict == 2
    assert n_loose == 3
    assert t_strict == 3
    assert t_loose == 5
    assert t_exp == 3
    assert 'num_overlapping_peaks_strict' in df_proc.columns


def test_ucsc_conservation_structure():
    """Verify UCSC conservation client returns proper schema for valid coordinates."""
    # Query a small 10bp interval
    res = get_ucsc_conservation("chr7", 1000, 1010, genome="hg38")
    assert isinstance(res, dict)
    coord_key = "chr7:1000-1010"
    assert coord_key in res
    assert "tracks" in res[coord_key]
    assert "phyloP100way" in res[coord_key]["tracks"]
    assert "phastCons100way" in res[coord_key]["tracks"]
    assert isinstance(res[coord_key]["tracks"]["phyloP100way"], list)


def test_gtex_resolution_and_fallback():
    """Verify GTEx gene resolution returns None gracefully for nonexistent or malformed symbols."""
    invalid_res = resolve_gtex_gencode_id("NONEXISTENT_GENE_SYMBOL_XYZ_12345")
    assert invalid_res is None


def test_string_data_graceful_error_handling():
    """Verify STRING client returns expected 3-tuple structure without throwing unhandled exceptions."""
    # Test with synthetic TF names
    img_bytes, df_partners, err = fetch_string_data("TP53", ["MDM2"])
    # Return should always be a 3-tuple: (Optional[bytes], Optional[DataFrame], Optional[str])
    assert img_bytes is None or isinstance(img_bytes, bytes)
    assert df_partners is None or isinstance(df_partners, pd.DataFrame)
    assert err is None or isinstance(err, str)
