"""Promoter enrichment analysis utilities.

Provides:
  - GC-content sliding-window profiling
  - CpG island detection  (Gardiner-Garden & Frommer 1987 algorithm)
  - TF co-binding frequency matrix
  - Ready-to-use matplotlib figures for both analyses
"""

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GC Content & CpG Island Detection
# ---------------------------------------------------------------------------

def calc_gc_profile(seq: str, window: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Compute sliding-window GC content across a DNA sequence.

    Parameters
    ----------
    seq    : DNA string (A/T/G/C; case-insensitive).
    window : Window size in bp (default 100).

    Returns
    -------
    positions    : Centre of each window (bp, 0-indexed from sequence start).
    gc_fractions : GC fraction in [0, 1] for each window.
    """
    seq = seq.upper()
    n = len(seq)
    step = max(1, window // 4)   # 25% step → 75% overlap
    positions, gc_fractions = [], []
    for i in range(0, n - window + 1, step):
        w = seq[i : i + window]
        gc = sum(1 for c in w if c in "GC")
        positions.append(i + window // 2)
        gc_fractions.append(gc / window)
    return np.array(positions, dtype=float), np.array(gc_fractions, dtype=float)


def find_cpg_islands(
    seq: str,
    window: int = 200,
    min_length: int = 200,
    gc_threshold: float = 0.50,
    cpg_oe_threshold: float = 0.60,
) -> List[Dict]:
    """Detect CpG islands using the Gardiner-Garden & Frommer (1987) algorithm.

    A region qualifies when a sliding window satisfies **all** of:
      * GC content ≥ ``gc_threshold``
      * CpG observed / expected ≥ ``cpg_oe_threshold``
        where  expected CpG = (count_C × count_G) / window

    Contiguous qualifying windows are merged; the merged region must be
    ≥ ``min_length`` bp to be reported.

    Parameters
    ----------
    seq              : DNA string (case-insensitive).
    window           : Sliding-window size (bp). Default 200.
    min_length       : Minimum island length to report (bp). Default 200.
    gc_threshold     : Minimum GC fraction. Default 0.50.
    cpg_oe_threshold : Minimum CpG O/E ratio. Default 0.60.

    Returns
    -------
    List of dicts with keys: start, end, length, gc_fraction, cpg_oe, cpg_count.
    Coordinates are 0-based, half-open [start, end).
    """
    seq = seq.upper()
    n = len(seq)
    if n < window:
        return []

    island_mask = np.zeros(n, dtype=bool)
    for i in range(n - window + 1):
        w = seq[i : i + window]
        c_count = w.count("C")
        g_count = w.count("G")
        gc = (c_count + g_count) / window
        if gc < gc_threshold:
            continue
        cpg_exp = (c_count * g_count) / window
        if cpg_exp == 0:
            continue
        cpg_oe = w.count("CG") / cpg_exp
        if cpg_oe >= cpg_oe_threshold:
            island_mask[i : i + window] = True

    # Merge contiguous True regions into islands
    islands: List[Dict] = []
    in_island = False
    start = 0
    for i in range(n + 1):
        currently_in = (i < n) and island_mask[i]
        if currently_in and not in_island:
            in_island = True
            start = i
        elif not currently_in and in_island:
            in_island = False
            length = i - start
            if length >= min_length:
                region = seq[start:i]
                c_count = region.count("C")
                g_count = region.count("G")
                gc = (c_count + g_count) / length
                cpg_obs = region.count("CG")
                cpg_exp = (c_count * g_count) / length
                cpg_oe = cpg_obs / cpg_exp if cpg_exp > 0 else 0.0
                islands.append(
                    {
                        "start": start,
                        "end": i,
                        "length": length,
                        "gc_fraction": round(gc, 3),
                        "cpg_oe": round(cpg_oe, 3),
                        "cpg_count": cpg_obs,
                    }
                )
    return islands


# ---------------------------------------------------------------------------
# TF Co-binding Analysis
# ---------------------------------------------------------------------------

def calc_tf_cobinding(df_encode: pd.DataFrame) -> pd.DataFrame:
    """Build a symmetric TF co-binding frequency matrix.

    For each pair of TFs, count the number of unique biosamples where *both*
    TFs have at least one overlapping peak (``overlap == True``, when that
    column is present) at the analysed promoter.

    Parameters
    ----------
    df_encode : DataFrame from ``analysis.analyze_gene`` (encode_hits.csv).
                Must contain columns ``tf`` and ``biosample``.
                If an ``overlap`` column is present, only True rows are used.

    Returns
    -------
    Square DataFrame indexed/columned by TF name with co-binding counts.
    An empty DataFrame is returned when fewer than two TFs are present.
    """
    if df_encode.empty or "tf" not in df_encode.columns:
        return pd.DataFrame()

    # Use only confirmed promoter-overlapping peaks when column is present
    if "overlap" in df_encode.columns:
        df = df_encode[df_encode["overlap"] == True].copy()
    else:
        df = df_encode.copy()

    tfs = sorted(df["tf"].unique())
    if len(tfs) < 2:
        return pd.DataFrame()

    cobind_matrix = pd.DataFrame(0, index=tfs, columns=tfs, dtype=int)

    if "biosample" not in df.columns:
        return cobind_matrix

    for _biosample, grp in df.groupby("biosample"):
        tfs_in_sample = grp["tf"].unique()
        for i, tf1 in enumerate(tfs_in_sample):
            for tf2 in tfs_in_sample[i:]:
                cobind_matrix.loc[tf1, tf2] += 1
                if tf1 != tf2:
                    cobind_matrix.loc[tf2, tf1] += 1

    return cobind_matrix


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_gc_cpg_profile(
    seq: str,
    promoter_up: int,
    promoter_down: int,
    gene_name: str,
    window: int = 100,
) -> plt.Figure:
    """Two-panel figure: GC-content sliding-window profile (top) + CpG island
    annotations (bottom).

    The x-axis shows distance from the TSS (negative = upstream).
    """
    positions, gc_fractions = calc_gc_profile(seq, window=window)
    islands = find_cpg_islands(seq)

    # Convert 0-based sequence positions → TSS-relative coords
    pos_tss = positions - promoter_up

    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(12, 5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    # Top panel: GC content
    ax1.fill_between(pos_tss, gc_fractions * 100, alpha=0.3, color="steelblue")
    ax1.plot(pos_tss, gc_fractions * 100, color="steelblue", linewidth=1.5)
    ax1.axhline(50, color="gray", linestyle="--", linewidth=0.8, label="50 % GC")
    ax1.axvline(0, color="black", linestyle="--", linewidth=1.0, label="TSS")
    ax1.set_ylabel("GC Content (%)")
    ax1.set_title(f"Promoter GC Profile & CpG Islands – {gene_name}")
    ax1.set_ylim(0, 100)
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True, alpha=0.2)

    # Bottom panel: CpG islands
    ax2.axvline(0, color="black", linestyle="--", linewidth=1.0)
    ax2.set_xlim(-promoter_up, promoter_down)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.set_xlabel("Distance to TSS (bp)")
    ax2.set_ylabel("CpG\nIslands", fontsize=8)

    if islands:
        for isl in islands:
            rel_s = isl["start"] - promoter_up
            rel_e = isl["end"] - promoter_up
            ax2.barh(
                0.5, rel_e - rel_s, left=rel_s,
                height=0.5, color="forestgreen", alpha=0.8,
            )
            if (rel_e - rel_s) > 80:
                ax2.text(
                    (rel_s + rel_e) / 2, 0.5, "CGI",
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold",
                )
    else:
        ax2.text(
            0, 0.5, "No CpG Islands detected",
            ha="center", va="center", fontsize=9,
            color="gray", style="italic",
        )

    plt.tight_layout()
    return fig


def plot_cobinding_heatmap(cobind_matrix: pd.DataFrame, gene_name: str) -> plt.Figure:
    """Annotated heatmap of TF co-binding frequencies (shared biosamples)."""
    n = len(cobind_matrix)
    fig, ax = plt.subplots(figsize=(max(5, n * 1.2), max(4, n * 1.0)))
    sns.heatmap(
        cobind_matrix,
        annot=True,
        fmt="d",
        cmap="Blues",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Shared biosamples"},
    )
    ax.set_title(
        f"TF Co-binding at {gene_name} Promoter\n"
        "(number of biosamples where both TFs bind)"
    )
    plt.tight_layout()
    return fig
