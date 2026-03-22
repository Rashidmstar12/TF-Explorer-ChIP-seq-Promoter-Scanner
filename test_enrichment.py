"""Unit tests for tf_explorer.enrichment module."""

import sys
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for CI

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from tf_explorer.enrichment import (
    calc_gc_profile,
    find_cpg_islands,
    calc_tf_cobinding,
    plot_gc_cpg_profile,
    plot_cobinding_heatmap,
)


# ---------------------------------------------------------------------------
# calc_gc_profile
# ---------------------------------------------------------------------------

def test_calc_gc_profile_basic():
    """GC profile returns sensible values for a known sequence."""
    seq = "GCGCGCGCGCGCGCGCGCGC"  # 100 % GC
    positions, gc_fracs = calc_gc_profile(seq, window=4)
    assert len(positions) == len(gc_fracs)
    assert len(positions) > 0
    assert all(f == 1.0 for f in gc_fracs), "Expected 100 % GC"


def test_calc_gc_profile_at_only():
    """Pure AT sequence → 0 % GC in every window."""
    seq = "ATATATATAT" * 10
    positions, gc_fracs = calc_gc_profile(seq, window=10)
    assert all(f == 0.0 for f in gc_fracs), "Expected 0 % GC"


def test_calc_gc_profile_mixed():
    """50 % GC sequence → all windows near 0.5."""
    # "GCTA" repeats: each 4-bp unit is exactly 50 % GC, so every window
    # of size 4 aligned on the repeat boundary is exactly 0.5.
    seq = "GCTA" * 100             # 400 bp, exactly 50 % GC per window=4
    positions, gc_fracs = calc_gc_profile(seq, window=4)
    for f in gc_fracs:
        assert abs(f - 0.5) < 1e-9, f"Expected exactly 50 % GC, got {f:.4f}"


def test_calc_gc_profile_window_larger_than_seq():
    """When window > sequence length, result should be empty."""
    seq = "ACGT"
    positions, gc_fracs = calc_gc_profile(seq, window=100)
    assert len(positions) == 0
    assert len(gc_fracs) == 0


# ---------------------------------------------------------------------------
# find_cpg_islands
# ---------------------------------------------------------------------------

def test_find_cpg_islands_synthetic():
    """A synthetic CG-dense region embedded in AT flanks is detected.

    The island scanner uses a 200-bp sliding window, so qualifying windows
    can start as early as the point where the window first becomes ≥50% GC
    (i.e., ~100bp before the core). We therefore check that the detected
    island *covers* the CG-rich core (end position > core start) rather
    than asserting on its exact start.
    """
    flank = "AT" * 150                        # 300 bp, 0 % GC, 0 CpG
    core = "CGCGCGCGCG" * 40                  # 400 bp, 100 % GC, many CpG
    seq = flank + core + flank
    islands = find_cpg_islands(seq)
    assert len(islands) >= 1, "Expected at least one CpG island in CG-rich core"
    # At least one island must extend past the start of the core (pos 300)
    assert any(isl["end"] > 300 for isl in islands), (
        "Expected an island covering the CG-rich core region"
    )


def test_find_cpg_islands_no_island():
    """Pure AT sequence should have no CpG islands."""
    seq = "AT" * 500               # 1000 bp, 0 % GC
    islands = find_cpg_islands(seq)
    assert islands == [], "Expected no CpG islands in AT-only sequence"


def test_find_cpg_islands_short_seq():
    """Sequence shorter than window returns empty list."""
    seq = "CGCGCGCG"              # 8 bp < 200 bp window
    islands = find_cpg_islands(seq)
    assert islands == [], "Expected no islands when seq < window"


def test_find_cpg_islands_fields():
    """Returned dicts contain required keys."""
    core = "CGCG" * 100           # 400 bp, high GC & CpG OE
    seq = "AT" * 100 + core + "AT" * 100
    islands = find_cpg_islands(seq)
    if islands:
        required = {"start", "end", "length", "gc_fraction", "cpg_oe", "cpg_count"}
        for isl in islands:
            assert required.issubset(isl.keys()), f"Missing keys: {required - isl.keys()}"
            assert isl["end"] > isl["start"]
            assert isl["length"] == isl["end"] - isl["start"]
            assert 0.0 <= isl["gc_fraction"] <= 1.0


# ---------------------------------------------------------------------------
# calc_tf_cobinding
# ---------------------------------------------------------------------------

def _make_hits_df():
    """Helper: three TFs, three biosamples with partial overlap."""
    data = {
        "tf":        ["TFA", "TFB", "TFA", "TFC", "TFB", "TFC"],
        "biosample": ["CL1", "CL1", "CL2", "CL2", "CL3", "CL3"],
        "overlap":   [True,  True,  True,  True,  True,  True],
        "experiment": ["E1", "E2", "E3", "E4", "E5", "E6"],
        "peak_start": [100, 200, 150, 250, 100, 300],
        "peak_end":   [150, 250, 200, 300, 150, 350],
    }
    return pd.DataFrame(data)


def test_calc_tf_cobinding_shape():
    """Matrix has correct shape (n_tfs × n_tfs)."""
    df = _make_hits_df()
    matrix = calc_tf_cobinding(df)
    assert matrix.shape == (3, 3), f"Expected 3×3, got {matrix.shape}"


def test_calc_tf_cobinding_symmetric():
    """Matrix is symmetric."""
    df = _make_hits_df()
    matrix = calc_tf_cobinding(df)
    for tf1 in matrix.index:
        for tf2 in matrix.columns:
            assert matrix.loc[tf1, tf2] == matrix.loc[tf2, tf1], "Matrix not symmetric"


def test_calc_tf_cobinding_diagonal():
    """Diagonal equals the number of biosamples that TF binds in."""
    df = _make_hits_df()
    matrix = calc_tf_cobinding(df)
    # TFA binds in CL1 and CL2 → diagonal[TFA] == 2
    assert matrix.loc["TFA", "TFA"] == 2
    # TFB binds in CL1 and CL3 → diagonal[TFB] == 2
    assert matrix.loc["TFB", "TFB"] == 2


def test_calc_tf_cobinding_single_tf_returns_empty():
    """Fewer than two TFs → empty DataFrame."""
    df = pd.DataFrame(
        {"tf": ["TFA", "TFA"], "biosample": ["CL1", "CL2"], "overlap": [True, True]}
    )
    matrix = calc_tf_cobinding(df)
    assert matrix.empty, "Expected empty matrix for single TF"


def test_calc_tf_cobinding_empty_input():
    """Empty input → empty DataFrame."""
    assert calc_tf_cobinding(pd.DataFrame()).empty


# ---------------------------------------------------------------------------
# Plotting smoke tests (just verify no exception is raised)
# ---------------------------------------------------------------------------

def test_plot_gc_cpg_profile_runs():
    seq = ("CGCGCGCGCG" * 40) + ("AT" * 100) + ("CGCGCGCGCG" * 40)
    fig = plot_gc_cpg_profile(seq, promoter_up=500, promoter_down=200, gene_name="TESTGENE")
    assert fig is not None


def test_plot_cobinding_heatmap_runs():
    df = _make_hits_df()
    from tf_explorer.enrichment import calc_tf_cobinding
    matrix = calc_tf_cobinding(df)
    fig = plot_cobinding_heatmap(matrix, "TESTGENE")
    assert fig is not None


if __name__ == "__main__":
    test_calc_gc_profile_basic()
    test_calc_gc_profile_at_only()
    test_calc_gc_profile_mixed()
    test_calc_gc_profile_window_larger_than_seq()
    test_find_cpg_islands_synthetic()
    test_find_cpg_islands_no_island()
    test_find_cpg_islands_short_seq()
    test_find_cpg_islands_fields()
    test_calc_tf_cobinding_shape()
    test_calc_tf_cobinding_symmetric()
    test_calc_tf_cobinding_diagonal()
    test_calc_tf_cobinding_single_tf_returns_empty()
    test_calc_tf_cobinding_empty_input()
    test_plot_gc_cpg_profile_runs()
    test_plot_cobinding_heatmap_runs()
    print("All tests passed!")
