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


# ---------------------------------------------------------------------------
# find_core_promoter_elements
# ---------------------------------------------------------------------------

from tf_explorer.enrichment import find_core_promoter_elements, _reverse_complement


def _make_promoter_with_tata(promoter_up=2000, promoter_down=500):
    """Return a promoter sequence with a clear TATA box at -28 bp (canonical)."""
    total_len = promoter_up + promoter_down
    seq = list("A" * total_len)
    # Place TATAAAA at position (promoter_up - 28) which is TSS-rel = -28
    # "TATAAAA" matches the pattern TATA[AT]A[AT]
    tata_pos = promoter_up - 28
    for i, c in enumerate("TATAAAA"):
        if tata_pos + i < total_len:
            seq[tata_pos + i] = c
    return "".join(seq)


def test_find_core_elements_tata_detected():
    """TATA box inserted at -28 bp is detected in canonical position."""
    seq = _make_promoter_with_tata()
    hits = find_core_promoter_elements(seq, promoter_up=2000)
    tata_hits = [h for h in hits if h["name"] == "TATA box"]
    assert len(tata_hits) >= 1, "Expected at least one TATA box hit"
    canonical = [h for h in tata_hits if h["canonical_position"]]
    assert len(canonical) >= 1, "Expected at least one canonical TATA hit at -28"


def test_find_core_elements_required_keys():
    """Every returned hit contains all required keys."""
    seq = _make_promoter_with_tata()
    hits = find_core_promoter_elements(seq, promoter_up=2000)
    required = {
        "name", "start", "end", "start_tss_rel", "end_tss_rel",
        "matched_seq", "strand", "canonical_position", "description",
    }
    for h in hits:
        missing = required - h.keys()
        assert not missing, f"Hit is missing keys: {missing}"
        assert h["end"] > h["start"], "end must be > start"
        assert isinstance(h["canonical_position"], bool)


def test_find_core_elements_empty_seq():
    """Short sequence returns empty list without error."""
    hits = find_core_promoter_elements("ACGT", promoter_up=2)
    assert hits == []


def test_find_core_elements_gc_box():
    """Synthetic GC box in promoter body is detected (possibly non-canonical)."""
    # Place GGGCGG 100 bp upstream of TSS
    seq = list("A" * 2500)
    pos = 2000 - 100
    for i, c in enumerate("GGGCGG"):
        seq[pos + i] = c
    hits = find_core_promoter_elements("".join(seq), promoter_up=2000)
    gc_hits = [h for h in hits if "GC box" in h["name"]]
    assert len(gc_hits) >= 1, "Expected at least one GC box hit"


def test_reverse_complement():
    """_reverse_complement returns correct reverse complement."""
    assert _reverse_complement("ATCG") == "CGAT"
    assert _reverse_complement("AAAA") == "TTTT"
    assert _reverse_complement("GCGC") == "GCGC"
    assert _reverse_complement("") == ""


def test_plot_core_promoter_elements_runs():
    """plot_core_promoter_elements does not raise with real data."""
    seq = _make_promoter_with_tata()
    from tf_explorer.enrichment import find_core_promoter_elements, plot_core_promoter_elements
    hits = find_core_promoter_elements(seq, promoter_up=2000)
    fig = plot_core_promoter_elements(hits, 2000, 500, "TESTGENE")
    assert fig is not None


def test_plot_core_promoter_elements_empty():
    """plot_core_promoter_elements works with empty hit list."""
    from tf_explorer.enrichment import plot_core_promoter_elements
    fig = plot_core_promoter_elements([], 2000, 500, "TESTGENE")
    assert fig is not None


# ---------------------------------------------------------------------------
# calc_consensus_peaks
# ---------------------------------------------------------------------------

from tf_explorer.enrichment import calc_consensus_peaks


def _make_consensus_df():
    """Helper: two experiments each having a peak in the same region."""
    return pd.DataFrame(
        {
            "peak_start":    [1000, 1010, 5000],
            "peak_end":      [1100, 1110, 5100],
            "experiment":    ["E1", "E2", "E1"],
            "tf":            ["TFA", "TFA", "TFA"],
            "biosample":     ["CL1", "CL2", "CL1"],
            "overlap":       [True,  True,  True],
            "signal":        [20.0,  30.0,  10.0],
            "distance_to_tss": [-50, -40, 1000],
        }
    )


def test_calc_consensus_peaks_basic():
    """Two overlapping peaks from different experiments form one consensus peak."""
    df = _make_consensus_df()
    result = calc_consensus_peaks(df, min_experiments=2, merge_distance=50)
    assert len(result) >= 1, "Expected at least one consensus peak"
    assert result.iloc[0]["supporting_experiments"] >= 2


def test_calc_consensus_peaks_min_experiments_filter():
    """Increasing min_experiments to 3 excludes regions with only 2 supporting exps."""
    df = _make_consensus_df()
    result3 = calc_consensus_peaks(df, min_experiments=3, merge_distance=50)
    # Only 2 unique experiments in our mock data, so nothing should pass
    assert result3.empty, "Expected no consensus peaks when min_experiments=3 > available=2"


def test_calc_consensus_peaks_empty():
    """Empty DataFrame returns empty result."""
    assert calc_consensus_peaks(pd.DataFrame()).empty


def test_calc_consensus_peaks_required_columns():
    """Returned DataFrame has expected columns."""
    df = _make_consensus_df()
    result = calc_consensus_peaks(df, min_experiments=2)
    if not result.empty:
        for col in ["start", "end", "width", "supporting_experiments"]:
            assert col in result.columns, f"Missing column: {col}"


# ---------------------------------------------------------------------------
# plot_signal_distribution
# ---------------------------------------------------------------------------

from tf_explorer.enrichment import plot_signal_distribution


def _make_signal_df():
    return pd.DataFrame(
        {
            "biosample": ["CL1", "CL1", "CL2", "CL2", "CL3"],
            "tf":        ["TFA", "TFA", "TFB", "TFB", "TFA"],
            "signal":    [10.0,  20.0,  15.0,  25.0,  30.0],
            "score":     [100,   200,   150,   250,   300],
            "overlap":   [True,  True,  True,  True,  True],
        }
    )


def test_plot_signal_distribution_by_biosample():
    """Returns a figure when grouping by biosample."""
    df = _make_signal_df()
    fig = plot_signal_distribution(df, "TESTGENE", group_col="biosample")
    assert fig is not None


def test_plot_signal_distribution_by_tf():
    """Returns a figure when grouping by tf."""
    df = _make_signal_df()
    fig = plot_signal_distribution(df, "TESTGENE", group_col="tf")
    assert fig is not None


def test_plot_signal_distribution_empty():
    """Returns None for empty DataFrame."""
    result = plot_signal_distribution(pd.DataFrame(), "TESTGENE")
    assert result is None


def test_plot_signal_distribution_all_zero_signal():
    """Returns None when all signal values are zero."""
    df = pd.DataFrame({"signal": [0, 0, 0], "biosample": ["CL1", "CL1", "CL2"]})
    result = plot_signal_distribution(df, "TESTGENE")
    assert result is None


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
    test_find_core_elements_tata_detected()
    test_find_core_elements_required_keys()
    test_find_core_elements_empty_seq()
    test_find_core_elements_gc_box()
    test_reverse_complement()
    test_plot_core_promoter_elements_runs()
    test_plot_core_promoter_elements_empty()
    test_calc_consensus_peaks_basic()
    test_calc_consensus_peaks_min_experiments_filter()
    test_calc_consensus_peaks_empty()
    test_calc_consensus_peaks_required_columns()
    test_plot_signal_distribution_by_biosample()
    test_plot_signal_distribution_by_tf()
    test_plot_signal_distribution_empty()
    test_plot_signal_distribution_all_zero_signal()
    print("All tests passed!")
