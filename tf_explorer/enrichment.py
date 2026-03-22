"""Promoter enrichment analysis utilities.

Provides:
  - GC-content sliding-window profiling
  - CpG island detection  (Gardiner-Garden & Frommer 1987 algorithm)
  - TF co-binding frequency matrix
  - Core promoter element detection (TATA, CCAAT, GC-box, Inr)
  - Consensus / high-confidence peak identification
  - ChIP-seq signal intensity distribution plots
  - Ready-to-use matplotlib figures for all analyses
"""

import logging
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple

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


# ---------------------------------------------------------------------------
# Core Promoter Element Detection
# ---------------------------------------------------------------------------

#: Canonical core promoter elements with IUPAC regex patterns and expected
#: TSS-relative position ranges (upstream = negative).
_CORE_ELEMENTS = [
    {
        "name": "TATA box",
        "pattern": r"TATA[AT]A[AT]",           # TATAAA consensus (relaxed)
        "pos_min": -40,
        "pos_max": -20,
        "description": "Core promoter element ~25-30 bp upstream of TSS; recruits TBP/TFIID",
    },
    {
        "name": "CCAAT box",
        "pattern": r"[CG]CAAT",                # NF-Y binding site
        "pos_min": -100,
        "pos_max": -60,
        "description": "Recruits NF-Y and CBF; found at ~-80 bp in many TATA-less promoters",
    },
    {
        "name": "GC box (SP1)",
        "pattern": r"GGGCGG|CCGCCC",           # SP1/KLF consensus
        "pos_min": -200,
        "pos_max": 50,
        "description": "SP1/KLF family binding site (GC-rich); common in CpG island promoters",
    },
    {
        "name": "Initiator (Inr)",
        "pattern": r"[CT][CT]A[ACGT][AT][CT][CT]",  # YYANWYY
        "pos_min": -5,
        "pos_max": 5,
        "description": "Overlaps the TSS; mediates basal transcription in TATA-less promoters",
    },
    {
        "name": "DPE (Downstream Promoter Element)",
        "pattern": r"[AG]G[AT][CT][GT]",       # consensus RG[AT][CT][GT]
        "pos_min": 28,
        "pos_max": 34,
        "description": "Downstream core element at +28-34; co-operates with Inr",
    },
]


def find_core_promoter_elements(seq: str, promoter_up: int) -> List[Dict]:
    """Scan the promoter sequence for classical core promoter elements.

    Searches for TATA box, CCAAT box, GC box (SP1 site), Initiator element,
    and Downstream Promoter Element (DPE) using IUPAC-based regex patterns.
    Hits outside the canonical position window for each element are still
    reported but flagged with ``canonical_position = False``.

    Parameters
    ----------
    seq         : Promoter DNA sequence (case-insensitive).
                  Position 0 in the string corresponds to ``-promoter_up`` bp
                  relative to the TSS.
    promoter_up : Number of bp upstream of the TSS in ``seq``.

    Returns
    -------
    List of dicts with keys:
      name, pattern, start, end, start_tss_rel, end_tss_rel,
      matched_seq, strand, canonical_position, description.
    Coordinates are 0-based, half-open [start, end).
    """
    seq_upper = seq.upper()
    seq_rc = _reverse_complement(seq_upper)
    hits: List[Dict] = []

    for elem in _CORE_ELEMENTS:
        for strand, s in (("+", seq_upper), ("-", seq_rc)):
            for m in re.finditer(elem["pattern"], s):
                if strand == "+":
                    start = m.start()
                    end = m.end()
                else:
                    # Map reverse-complement position back to forward strand
                    end = len(seq) - m.start()
                    start = len(seq) - m.end()

                start_tss = start - promoter_up
                end_tss = end - promoter_up
                canonical = (
                    elem["pos_min"] <= start_tss <= elem["pos_max"]
                    or elem["pos_min"] <= end_tss - 1 <= elem["pos_max"]
                )

                hits.append(
                    {
                        "name": elem["name"],
                        "pattern": elem["pattern"],
                        "start": start,
                        "end": end,
                        "start_tss_rel": start_tss,
                        "end_tss_rel": end_tss,
                        "matched_seq": seq_upper[start:end] if strand == "+" else seq_rc[m.start():m.end()],
                        "strand": strand,
                        "canonical_position": canonical,
                        "description": elem["description"],
                    }
                )

    return hits


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA string (uppercase only)."""
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def plot_core_promoter_elements(
    hits: List[Dict],
    promoter_up: int,
    promoter_down: int,
    gene_name: str,
) -> plt.Figure:
    """Draw a linear map of detected core promoter elements relative to TSS.

    Canonical-position hits are shown as filled boxes; non-canonical as
    hatched boxes. Each element type gets a distinct colour.
    """
    if not hits:
        fig, ax = plt.subplots(figsize=(12, 2))
        ax.text(0, 0.5, "No core promoter elements detected.",
                ha="center", va="center", fontsize=11, color="gray", style="italic")
        ax.set_xlim(-promoter_up, promoter_down)
        ax.axis("off")
        ax.set_title(f"Core Promoter Elements – {gene_name}")
        plt.tight_layout()
        return fig

    names = sorted({h["name"] for h in hits})
    palette = sns.color_palette("tab10", len(names))
    color_map = dict(zip(names, palette))

    # Separate forward and reverse
    fwd = [h for h in hits if h["strand"] == "+"]
    rev = [h for h in hits if h["strand"] == "-"]

    n_rows = 2  # forward + reverse
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 3), sharex=True,
                              gridspec_kw={"height_ratios": [1, 1]})
    fig.suptitle(f"Core Promoter Elements – {gene_name}", fontsize=12, y=1.02)

    for ax, panel_hits, strand_label, y_track in zip(axes, [fwd, rev], ["Forward (+)", "Reverse (-)"], [0.5, 0.5]):
        ax.axvline(0, color="black", linestyle="--", linewidth=1.0, label="TSS")
        ax.set_xlim(-promoter_up, promoter_down)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(strand_label, fontsize=9)
        ax.grid(True, axis="x", alpha=0.15)

        for h in panel_hits:
            s = h["start_tss_rel"]
            e = h["end_tss_rel"]
            w = max(e - s, 1)
            color = color_map[h["name"]]
            hatch = "" if h["canonical_position"] else "///"
            rect = plt.Rectangle((s, 0.2), w, 0.6,
                                  facecolor=color, edgecolor="black",
                                  linewidth=0.8, hatch=hatch, alpha=0.8)
            ax.add_patch(rect)
            if w > 4:
                ax.text(s + w / 2, 0.5, h["name"].split(" ")[0],
                        ha="center", va="center", fontsize=6.5, color="white",
                        fontweight="bold", clip_on=True)

    axes[-1].set_xlabel("Distance to TSS (bp)")

    # Legend
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=color_map[n], label=n, edgecolor="black") for n in names]
    legend_handles.append(Patch(facecolor="white", hatch="///", label="Non-canonical position", edgecolor="black"))
    fig.legend(handles=legend_handles, loc="lower center", ncol=min(len(legend_handles), 4),
               bbox_to_anchor=(0.5, -0.25), fontsize=8, frameon=True)
    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Consensus Peak Identification
# ---------------------------------------------------------------------------

def calc_consensus_peaks(
    df_encode: pd.DataFrame,
    min_experiments: int = 2,
    merge_distance: int = 50,
) -> pd.DataFrame:
    """Identify high-confidence (consensus) peaks found in multiple experiments.

    A consensus peak must overlap or be within ``merge_distance`` bp of another
    peak from a *different* experiment. The returned DataFrame contains each
    unique peak region with the count of supporting experiments.

    Parameters
    ----------
    df_encode        : encode_hits.csv DataFrame (must have ``peak_start``,
                       ``peak_end``, ``experiment`` columns).
    min_experiments  : Minimum number of experiments a peak must appear in
                       (default 2).
    merge_distance   : Max gap (bp) between peaks to merge them (default 50).

    Returns
    -------
    DataFrame with columns: chrom, start, end, width, supporting_experiments,
    experiment_ids, tfs, biosamples, mean_signal.
    """
    required = {"peak_start", "peak_end", "experiment"}
    if df_encode.empty or not required.issubset(df_encode.columns):
        return pd.DataFrame()

    # Use overlapping peaks only when column present
    df = df_encode.copy()
    if "overlap" in df.columns:
        df = df[df["overlap"] == True]

    if df.empty:
        return pd.DataFrame()

    # Sort peaks
    df = df.sort_values(["peak_start", "peak_end"]).reset_index(drop=True)

    # Build a normalised working frame with stable column names
    _tf_src  = "tf"       if "tf"       in df.columns else "experiment"
    _bio_src = "biosample" if "biosample" in df.columns else "experiment"
    _sig_src = "signal"   if ("signal"   in df.columns and df["signal"].max() > 0) else (
               "score"    if "score"    in df.columns else "peak_start")

    records = df[["peak_start", "peak_end", "experiment", _tf_src, _bio_src, _sig_src]].copy()
    records.columns = ["start", "end", "exp", "tf", "biosample", "signal"]
    records = records.sort_values("start").reset_index(drop=True)

    clusters: List[Dict] = []
    cur_start = int(records.loc[0, "start"])
    cur_end   = int(records.loc[0, "end"])
    cur_rows  = [0]

    for i in range(1, len(records)):
        row_start = int(records.loc[i, "start"])
        row_end   = int(records.loc[i, "end"])
        if row_start <= cur_end + merge_distance:
            cur_end = max(cur_end, row_end)
            cur_rows.append(i)
        else:
            clusters.append({"start": cur_start, "end": cur_end, "rows": cur_rows})
            cur_start = row_start
            cur_end   = row_end
            cur_rows  = [i]
    clusters.append({"start": cur_start, "end": cur_end, "rows": cur_rows})

    # Filter clusters with enough unique experiments
    results = []
    for cl in clusters:
        sub = records.iloc[cl["rows"]]
        unique_exps = sub["exp"].nunique()
        if unique_exps >= min_experiments:
            mean_sig = sub["signal"].replace(0, np.nan).mean()
            results.append(
                {
                    "start": cl["start"],
                    "end": cl["end"],
                    "width": cl["end"] - cl["start"],
                    "supporting_experiments": unique_exps,
                    "experiment_ids": "; ".join(sorted(sub["exp"].unique())),
                    "tfs": "; ".join(sorted(sub["tf"].unique())),
                    "biosamples": "; ".join(sorted(sub["biosample"].unique())),
                    "mean_signal": round(float(mean_sig), 2) if not np.isnan(mean_sig) else 0.0,
                }
            )

    if not results:
        return pd.DataFrame()

    return pd.DataFrame(results).sort_values("supporting_experiments", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Signal Intensity Distribution
# ---------------------------------------------------------------------------

def plot_signal_distribution(
    df_encode: pd.DataFrame,
    gene_name: str,
    signal_col: str = "signal",
    group_col: str = "biosample",
) -> Optional[plt.Figure]:
    """Violin + strip plot of ChIP-seq signal intensity per group.

    Returns None when there is insufficient data.

    Parameters
    ----------
    df_encode  : encode_hits.csv DataFrame.
    gene_name  : Gene name for the plot title.
    signal_col : Column containing signal values (default ``"signal"``).
    group_col  : Grouping column (default ``"biosample"``).
    """
    if df_encode.empty:
        return None

    # Choose best available signal column
    for col in [signal_col, "score"]:
        if col in df_encode.columns and df_encode[col].max() > 0:
            signal_col = col
            break
    else:
        return None

    # Use only overlapping peaks when column present
    df = df_encode.copy()
    if "overlap" in df.columns:
        df = df[df["overlap"] == True]

    if df.empty or signal_col not in df.columns:
        return None

    # Drop zero / null signals
    df = df[df[signal_col] > 0].dropna(subset=[signal_col])
    if df.empty:
        return None

    group_col_used = group_col if group_col in df.columns else (
        "tf" if "tf" in df.columns else None
    )

    fig, ax = plt.subplots(figsize=(max(8, len(df[group_col_used].unique()) * 1.2 if group_col_used else 6), 5))

    if group_col_used and df[group_col_used].nunique() >= 2:
        order = (
            df.groupby(group_col_used)[signal_col]
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        sns.violinplot(data=df, x=group_col_used, y=signal_col, order=order,
                       hue=group_col_used, legend=False,
                       ax=ax, palette="Set2", inner="box", cut=0)
        sns.stripplot(data=df, x=group_col_used, y=signal_col, order=order,
                      hue=group_col_used, legend=False,
                      ax=ax, palette="dark:black", alpha=0.3, size=3, jitter=True)
        ax.set_xlabel(group_col_used.capitalize())
        plt.xticks(rotation=45, ha="right")
    else:
        # Single group or no grouping: use histogram
        ax.hist(df[signal_col].dropna(), bins=30, color="steelblue", edgecolor="white", alpha=0.85)
        ax.set_xlabel(signal_col.capitalize())
        ax.set_ylabel("Count")

    ax.set_ylabel(f"ChIP-seq {signal_col.capitalize()}")
    ax.set_title(f"ChIP-seq Signal Distribution – {gene_name}")
    ax.grid(True, axis="y", alpha=0.2)
    plt.tight_layout()
    return fig
