import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import yaml
import time
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any
import random

from . import encode_client, genome_client, motifs

# Configure logging
logger = logging.getLogger(__name__)

def setup_logging(output_dir: str):
    """Sets up file logging to the output directory."""
    log_file = os.path.join(output_dir, "analysis_log.txt")
    file_handler = logging.FileHandler(log_file)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Add to root logger
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    root_logger.setLevel(logging.INFO)

def save_config(output_dir: str, config: Dict):
    """Saves the configuration used for the run."""
    with open(os.path.join(output_dir, "config_used.yaml"), "w") as f:
        yaml.dump(config, f)

from .external_apis import get_ucsc_conservation

def calculate_promoter_coordinates(tss: int, strand: str, upstream: int = 2000, downstream: int = 500) -> Tuple[int, int]:
    """
    Calculate strand-aware promoter genomic coordinates [start, end).
    - Positive strand (+): [tss - upstream, tss + downstream)
    - Negative strand (-): [tss - downstream, tss + upstream)
    Ensures start < end on the genomic coordinate axis.
    """
    if strand == "+":
        return tss - upstream, tss + downstream
    else:
        return tss - downstream, tss + upstream

def get_default_cache_dir() -> str:
    """
    Returns a standard, cross-platform cache directory for TF-Explorer data.
    Uses platformdirs to adhere to OS-specific conventions:
    - Linux: ~/.cache/tf_explorer
    - macOS: ~/Library/Caches/tf_explorer
    - Windows: %LOCALAPPDATA%/tf_explorer/Cache
    Falls back to a local '.tf_explorer_cache' folder if permissions issues arise.
    """
    try:
        import platformdirs
        cache_dir = platformdirs.user_cache_dir("tf_explorer")
        os.makedirs(cache_dir, exist_ok=True)
        return cache_dir
    except Exception as e:
        logger.warning(f"Failed to initialize user cache dir via platformdirs: {e}. Falling back to local cache.")
        local_cache = os.path.join(os.getcwd(), ".tf_explorer_cache")
        os.makedirs(local_cache, exist_ok=True)
        return local_cache

def summarize_experiments(df_stats: pd.DataFrame):
    """
    Standardizes experiment metrics calculation decoupled from UI.
    Input: df_stats (raw DataFrame from experiment_stats.csv or analysis pipeline)
    Returns: (n_strict, n_loose, total_strict_peaks, total_loose_peaks, total_experiments, processed_df)
    """
    col_map = {
        'num_strict_peaks': 'num_overlapping_peaks_strict',
        'num_loose_peaks': 'num_overlapping_peaks_loose',
        'total_peaks_in_file': 'total_peaks'
    }
    df = df_stats.rename(columns=col_map)
    
    for col in ['num_overlapping_peaks_strict', 'num_overlapping_peaks_loose']:
        if col not in df.columns:
            df[col] = 0
            
    n_strict = int((df['num_overlapping_peaks_strict'] > 0).sum())
    n_loose = int((df['num_overlapping_peaks_loose'] > 0).sum())
            
    total_strict_peaks = int(df['num_overlapping_peaks_strict'].sum())
    total_loose_peaks = int(df['num_overlapping_peaks_loose'].sum())
    total_experiments = len(df)
    
    return n_strict, n_loose, total_strict_peaks, total_loose_peaks, total_experiments, df

def plot_cell_line_comparison(summary_df: pd.DataFrame):
    """
    Generates a bar chart comparing binding frequencies across cell lines.
    Deduplicates technical replicate files by experiment accession to prevent
    pseudoreplication skew (DEF-015).
    """
    if summary_df.empty or 'biosample' not in summary_df.columns:
        return None
        
    col = 'num_overlapping_peaks_strict' if 'num_overlapping_peaks_strict' in summary_df.columns else 'num_strict_peaks'
    if col not in summary_df.columns:
        return None
        
    # Check if experiment_id is available for deduplication
    has_experiment = 'experiment_id' in summary_df.columns
    
    if has_experiment:
        # Group by biosample and experiment_id first to get peak presence per experiment
        exp_level = summary_df.groupby(['biosample', 'experiment_id'])[col].max().reset_index()
        stats = exp_level.groupby('biosample').agg(
            total_units=('experiment_id', 'nunique'),
            units_with_peaks=(col, lambda x: (x > 0).sum())
        ).reset_index()
        metric_label = "Occupancy Rate (% of Experiments)"
    else:
        stats = summary_df.groupby('biosample').agg(
            total_units=('file_accession', 'count') if 'file_accession' in summary_df.columns else ('biosample', 'count'),
            units_with_peaks=(col, lambda x: (x > 0).sum())
        ).reset_index()
        metric_label = "Binding Rate (% of Files)"
    
    stats['binding_rate'] = (stats['units_with_peaks'] / stats['total_units']) * 100
    stats = stats.sort_values('binding_rate', ascending=False)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=stats, x='biosample', y='binding_rate', hue='biosample', legend=False, ax=ax, palette='viridis')
    
    ax.set_title(f"TF Binding Rate by Cell Line (Strict Window) - {metric_label}", fontsize=13)
    ax.set_ylabel(metric_label, fontsize=12)
    ax.set_xlabel("Cell Line", fontsize=12)
    plt.tight_layout()
    
    return fig

def scan_motif_synergy(df_motifs, max_distance=100) -> pd.DataFrame:
    """

    Scans predicted motifs to find clusters of different TFs within max_distance bp of each other.
    """
    if df_motifs.empty or len(df_motifs) < 2:
        return pd.DataFrame(columns=["start", "end", "tfs", "motifs_count", "description"])
        
    # Sort motifs by start coordinate
    df_sorted = df_motifs.sort_values("start").copy()
    
    synergy_regions = []
    n = len(df_sorted)
    
    # We want to find contiguous groups where the distance between adjacent motifs is <= max_distance
    # and there are at least 2 different TFs involved.
    i = 0
    while i < n:
        cluster = [df_sorted.iloc[i]]
        j = i + 1
        while j < n:
            curr_start = df_sorted.iloc[j]["start"]
            if curr_start - cluster[0]["start"] <= max_distance:
                cluster.append(df_sorted.iloc[j])
                j += 1
            else:
                break
                
        # Check if the cluster has at least 2 different TFs
        unique_tfs = set(m["tf_name"] for m in cluster)
        if len(unique_tfs) >= 2:
            cluster_start = min(m["start"] for m in cluster)
            cluster_end = max(m["end"] for m in cluster)
            tfs_str = ", ".join(sorted(unique_tfs))
            synergy_regions.append({
                "start": int(cluster_start),
                "end": int(cluster_end),
                "tfs": tfs_str,
                "motifs_count": len(cluster),
                "description": f"Synergy hotspot: {tfs_str} co-occurring within {max_distance}bp"
            })
            i = j
        else:
            i += 1
            
    return pd.DataFrame(synergy_regions)

def validate_gene_symbol(gene_name: str) -> str:
    """
    Validates gene symbol to prevent path traversal, injection, or malformed inputs.
    """
    if not isinstance(gene_name, str) or not gene_name.strip():
        raise ValueError("Gene symbol must be a non-empty string.")
    cleaned = gene_name.strip()
    forbidden = ["/", "\\", "..", "\0", "\n", "\r", ";", "&", "|", "`", "$"]
    if any(ch in cleaned for ch in forbidden):
        raise ValueError(f"Invalid gene symbol '{cleaned}': contains forbidden characters.")
    return cleaned

def validate_genome_assembly(genome: str) -> str:
    """
    Validates and standardizes genome assembly identifier.
    Supported: hg38, GRCh38, hg19, GRCh37.
    """
    if not isinstance(genome, str) or not genome.strip():
        return "hg38"
    cleaned = genome.strip().lower()
    if cleaned in ("hg38", "grch38"):
        return "hg38"
    if cleaned in ("hg19", "grch37"):
        return "hg19"
    raise ValueError(f"Unsupported genome assembly '{genome}'. Supported: hg38, hg19 (GRCh38, GRCh37).")

def analyze_gene(
    gene_name: str,
    tf_list: List[str],
    jaspar_ids: List[str] = None,
    experiments_list: List[Dict] = None,
    genome: str = "hg38",
    promoter_up: int = 2000,
    promoter_down: int = 500,
    loose_promoter_up: int = 5000,
    loose_promoter_down: int = 5000,
    pwm_threshold: float = 8.0,
    output_dir: str = "results",
    cache_dir: str = None,
    bed_output: bool = False,
    plot_track: bool = False,
    force_download: bool = False,
    random_seed: int = 42,
    progress_callback=None,
    explicit_tss: int = None,
    transcript_id: str = None
):
    """
    Main analysis pipeline.
    """
    gene_name = validate_gene_symbol(gene_name)
    genome = validate_genome_assembly(genome)

    os.makedirs(output_dir, exist_ok=True)
    setup_logging(output_dir)
    
    # Save config
    config = {
        "gene_name": gene_name,
        "tf_list": tf_list,
        "genome": genome,

        "promoter_up": promoter_up,
        "promoter_down": promoter_down,
        "loose_promoter_up": loose_promoter_up,
        "loose_promoter_down": loose_promoter_down,
        "pwm_threshold": pwm_threshold,
        "timestamp": datetime.now().isoformat()
    }
    save_config(output_dir, config)
    
    if experiments_list is None:
        experiments_list = []
        
    logger.info(f"Starting analysis for {gene_name}...")
    
    # 1. Get Gene Coordinates
    if explicit_tss is not None:
        # Use provided TSS
        # We still need chrom and strand, so we fetch coords but override TSS
        coords = genome_client.get_gene_coordinates(gene_name, genome)
        if not coords:
             msg = f"Could not find gene coordinates for '{gene_name}'."
             logger.error(msg)
             raise ValueError(msg)
        chrom, start, end, strand = coords
        tss = explicit_tss
        logger.info(f"Using explicit TSS for {gene_name} (Transcript {transcript_id}): {chrom}:{tss} ({strand})")
    else:
        coords = genome_client.get_gene_coordinates(gene_name, genome)
        if not coords:
             msg = f"Could not find gene coordinates for '{gene_name}'. Please check the gene symbol (e.g., 'DNMT3B' instead of 'dnmtb')."
             logger.error(msg)
             raise ValueError(msg)
        chrom, start, end, strand = coords
        
        # Calculate TSS
        tss = start if strand == "+" else end
        logger.info(f"TSS for {gene_name} is at {chrom}:{tss} ({strand})")
    
    # 2. Get Promoter Sequence
    promoter_seq = genome_client.get_promoter_sequence(chrom, tss, strand, promoter_up, promoter_down, genome)
    if not promoter_seq:
        logger.error("Could not fetch promoter sequence. Aborting.")
        return
        
    # 3. ENCODE Search & Overlap
    encode_hits = []
    nearest_nonoverlapping_peaks = []
    
    # Define promoter genomic region for overlap check
    if strand == "+":
        promoter_genomic_start = tss - promoter_up
        promoter_genomic_end = tss + promoter_down
        loose_genomic_start = tss - loose_promoter_up
        loose_genomic_end = tss + loose_promoter_down
    else:
        promoter_genomic_start = tss - promoter_down
        promoter_genomic_end = tss + promoter_up
        loose_genomic_start = tss - loose_promoter_down
        loose_genomic_end = tss + loose_promoter_up
        
    logger.info(f"Promoter genomic region: chr{chrom}:{promoter_genomic_start}-{promoter_genomic_end}")
    logger.info(f"Loose genomic region: chr{chrom}:{loose_genomic_start}-{loose_genomic_end}")
    
    # Collect all experiments first to know total count for progress
    all_experiments = []
    
    if experiments_list:
        all_experiments = experiments_list
        logger.info(f"Using {len(all_experiments)} provided experiments.")
    else:
        for tf in tf_list:
            exps = encode_client.search_encode_tf_chipseq(tf, "Homo sapiens", assembly=genome)
            for exp in exps:
                exp['tf_name'] = tf # Keep track of TF name
                all_experiments.append(exp)
    total_files = len(all_experiments)
    logger.info(f"Found total {total_files} peak files to process.")
    
    # Use cache_dir if provided, else use system default cache directory
    if cache_dir:
        download_dir = cache_dir
    else:
        download_dir = get_default_cache_dir()
    
    os.makedirs(download_dir, exist_ok=True)
    errors = []
    
    # Detailed stats for UI
    experiment_stats = []
    
    for i, exp in enumerate(all_experiments):
        tf = exp['tf_name']
        fname = f"{exp['file_accession']}.bed.gz"
        local_path = os.path.join(download_dir, fname)
        
        # Force download if requested
        if force_download and os.path.exists(local_path):
            try:
                os.remove(local_path)
            except OSError:
                pass
        
        # Update progress
        if progress_callback:
            progress_callback(i, total_files, f"Processing {tf}: {fname}")
            
        # Small delay to be nice to ENCODE API/Server
        if i > 0:
            time.sleep(2.0) 
            
        # Check if file exists and is valid
        if os.path.exists(local_path):
            try:
                import gzip
                with gzip.open(local_path, 'rt') as f:
                    f.read(1) 
            except Exception:
                logger.warning(f"Corrupt file found: {local_path}. Deleting and re-downloading.")
                os.remove(local_path)

        # Download if not exists
        if not os.path.exists(local_path):
            encode_client.download_peak_file(exp['download_url'], local_path)
        
        # Check overlap using helper function
        try:
            hits_strict, hits_loose, total_peaks, closest = find_overlaps_for_experiment(
                local_path, chrom, 
                promoter_genomic_start, promoter_genomic_end, 
                loose_genomic_start, loose_genomic_end,
                tss, strand, exp, tf
            )
            
            # Add to main hits list (using strict hits for main CSV)
            if hits_strict:
                encode_hits.extend(hits_strict)
            elif closest:
                nearest_nonoverlapping_peaks.append(closest)
                
            # Collect stats
            experiment_stats.append({
                "experiment_id": exp['dataset_accession'],
                "biosample": exp.get('biosample', 'Unknown'),
                "file_accession": exp['file_accession'],
                "binding_in_promoter": len(hits_strict) > 0,
                "binding_in_loose": len(hits_loose) > 0,
                "num_strict_peaks": len(hits_strict),
                "num_loose_peaks": len(hits_loose),
                "total_peaks_in_file": total_peaks
            })
            
        except Exception as e:
            logger.warning(f"Failed to process {fname}: {e}")
            errors.append(f"{tf} ({exp['dataset_accession']}): {str(e)}")
            
    if progress_callback:
        progress_callback(total_files, total_files, "Processing complete.")
                
    # Save ENCODE hits (strict overlapping peaks only)
    if not encode_hits:
        logger.info("No ENCODE strict hits found. Creating empty DataFrame with columns.")
        df_encode = pd.DataFrame(columns=["tf", "experiment", "biosample", "file", "peak_chrom", "peak_start", "peak_end", "distance_to_tss", "score", "signal", "overlap"])
    else:
        df_encode = pd.DataFrame(encode_hits)
    logger.info(f"Writing ENCODE strict hits to {os.path.join(output_dir, f'{gene_name}_encode_hits.csv')}")
    df_encode.to_csv(os.path.join(output_dir, f"{gene_name}_encode_hits.csv"), index=False)
    
    # Save nearest non-overlapping peaks separately
    if not nearest_nonoverlapping_peaks:
        df_nearest = pd.DataFrame(columns=["tf", "experiment", "biosample", "file", "peak_chrom", "peak_start", "peak_end", "distance_to_tss", "score", "signal", "overlap"])
    else:
        df_nearest = pd.DataFrame(nearest_nonoverlapping_peaks)
    logger.info(f"Writing nearest non-overlapping peaks to {os.path.join(output_dir, f'{gene_name}_nearest_nonoverlapping_peaks.csv')}")
    df_nearest.to_csv(os.path.join(output_dir, f"{gene_name}_nearest_nonoverlapping_peaks.csv"), index=False)
    
    # Save Detailed Stats
    df_stats = pd.DataFrame(experiment_stats)
    df_stats.to_csv(os.path.join(output_dir, f"{gene_name}_experiment_stats.csv"), index=False)
    
    # 4. Motif Scanning
    motif_hits = motifs.scan_promoter_with_jaspar(promoter_seq, jaspar_ids, threshold=pwm_threshold)
    if not motif_hits:
        df_motifs = pd.DataFrame(columns=["tf_name", "jaspar_id", "strand", "start", "end", "score", "sequence"])
    else:
        df_motifs = pd.DataFrame(motif_hits)
    df_motifs.to_csv(os.path.join(output_dir, f"{gene_name}_motif_predictions.csv"), index=False)
    
    # Fetch UCSC Conservation
    conservation_data = {}
    try:
        conservation_res = get_ucsc_conservation(f"chr{chrom}", promoter_genomic_start, promoter_genomic_end, genome)
        if conservation_res:
            key = list(conservation_res.keys())[0]
            conservation_data = conservation_res[key].get("tracks", {})
    except Exception as e:
        logger.warning(f"Error fetching conservation data: {e}")
        
    # Process base-by-base conservation data
    # Missing data is represented as np.nan to distinguish from true 0.0 scores or fast evolution (negative phyloP)
    base_data = {}
    for pos in range(promoter_genomic_start, promoter_genomic_end + 1):
        base_data[pos] = {"phyloP": np.nan, "phastCons": np.nan}
        
    phylo_list = conservation_data.get("phyloP100way", []) if conservation_data else []
    for entry in phylo_list:
        start_pos = entry.get("start")
        val = entry.get("value")
        if start_pos in base_data and val is not None:
            base_data[start_pos]["phyloP"] = float(val)
            
    phast_list = conservation_data.get("phastCons100way", []) if conservation_data else []
    for entry in phast_list:
        start_pos = entry.get("start")
        val = entry.get("value")
        if start_pos in base_data and val is not None:
            base_data[start_pos]["phastCons"] = float(val)
            
    rows = []
    for pos, vals in base_data.items():
        rel_pos = pos - tss if strand == "+" else tss - pos
        rows.append({
            "genomic_pos": pos,
            "distance_to_tss": rel_pos,
            "phyloP": vals["phyloP"],
            "phastCons": vals["phastCons"]
        })
    df_cons = pd.DataFrame(rows).sort_values("distance_to_tss")
    df_cons.to_csv(os.path.join(output_dir, f"{gene_name}_conservation.csv"), index=False)
    
    # High-Confidence Peak Identification (post-processing)
    if not df_encode.empty and not df_cons.empty:
        max_cons_list = []
        high_conf_list = []
        for _, row in df_encode.iterrows():
            try:
                p_start = int(row["peak_start"])
                p_end = int(row["peak_end"])
                peak_scores = df_cons[(df_cons["genomic_pos"] >= p_start) & (df_cons["genomic_pos"] <= p_end)]["phastCons"].dropna()
                if not peak_scores.empty:
                    max_c = float(peak_scores.max())
                    is_high_conf = bool(max_c > 0.8)
                else:
                    max_c = np.nan
                    is_high_conf = False
                max_cons_list.append(max_c)
                high_conf_list.append(is_high_conf)
            except Exception:
                max_cons_list.append(np.nan)
                high_conf_list.append(False)
        df_encode["max_conservation"] = max_cons_list
        df_encode["peak_overlaps_highly_conserved_sequence"] = high_conf_list
        df_encode["high_confidence_site"] = high_conf_list
    else:
        df_encode["max_conservation"] = np.nan
        df_encode["peak_overlaps_highly_conserved_sequence"] = False
        df_encode["high_confidence_site"] = False
        
    df_encode.to_csv(os.path.join(output_dir, f"{gene_name}_encode_hits.csv"), index=False)
    
    # Motif Synergy Scan
    df_synergy = scan_motif_synergy(df_motifs, max_distance=100)
    df_synergy.to_csv(os.path.join(output_dir, f"{gene_name}_synergy_hotspots.csv"), index=False)
    
    # 5. Combined Summary
    # Just a simple count summary for now
    biosamples = df_encode['biosample'].unique().tolist() if not df_encode.empty and 'biosample' in df_encode.columns else []
    biosamples_str = ", ".join([str(b) for b in biosamples if str(b) != "Unknown"])
    
    summary = {
        "gene": gene_name,
        "promoter_length": len(promoter_seq),
        "encode_peaks_found": len(df_encode),
        "experiments_with_peaks": df_encode['experiment'].nunique() if not df_encode.empty else 0,
        "biosamples_with_peaks": biosamples_str,
        "motif_predictions": len(df_motifs),
        "coords": f"{chrom}:{start}-{end} ({strand})",
        "promoter_window": f"chr{chrom}:{promoter_genomic_start}-{promoter_genomic_end}",
        "loose_window": f"chr{chrom}:{loose_genomic_start}-{loose_genomic_end}",
        "errors": "; ".join(errors) if errors else ""
    }
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(os.path.join(output_dir, f"{gene_name}_combined_summary.csv"), index=False)
    
    # 6. Visualization
    plot_binding_summary(df_encode, df_motifs, gene_name, output_dir)
    
    if plot_track:
        plot_promoter_track(df_encode, df_motifs, promoter_up, promoter_down, gene_name, output_dir)
        
    # 7. BED Output
    if bed_output:
        generate_bed_output(df_encode, df_motifs, chrom, promoter_genomic_start, gene_name, output_dir, strand=strand, promoter_end=promoter_genomic_end)

    logger.info("Analysis complete.")
    return summary_df

def plot_binding_summary(df_encode, df_motifs, gene_name, output_dir):
    """Generates a bar plot of binding counts."""
    plt.figure(figsize=(10, 6))
    
    data = []
    if not df_encode.empty:
        # Filter for overlapping peaks ONLY
        df_overlapping = df_encode[df_encode['overlap'] == True]
        
        if not df_overlapping.empty:
            encode_counts = df_overlapping['tf'].value_counts().reset_index()
            encode_counts.columns = ['TF', 'Count']
            encode_counts['Source'] = 'ENCODE ChIP-seq'
            data.append(encode_counts)
        
    if not df_motifs.empty:
        motif_counts = df_motifs['tf_name'].value_counts().reset_index()
        motif_counts.columns = ['TF', 'Count']
        motif_counts['Source'] = 'JASPAR PWM'
        data.append(motif_counts)
        
    if data:
        df_plot = pd.concat(data)
        sns.barplot(data=df_plot, x='TF', y='Count', hue='Source')
        plt.title(f"TF Binding Summary for {gene_name}")
        plt.savefig(os.path.join(output_dir, f"{gene_name}_tf_binding_plot.png"))
    plt.close()


def create_promoter_track_fig(df_encode, df_motifs, promoter_up, promoter_down, gene_name, high_confidence_threshold=None, top_n=None, df_cons=None, df_synergy=None):
    """
    Creates a multi-track matplotlib figure for the promoter region including ChIP peaks, motifs,
    shaded synergy hotspots, and base-by-base evolutionary conservation tracks.
    """
    import numpy as np
    from matplotlib.lines import Line2D
    
    has_cons = df_cons is not None and not df_cons.empty
    
    if has_cons:
        fig, (ax, ax_cons) = plt.subplots(
            2, 1, 
            figsize=(12, 8), 
            sharex=True, 
            gridspec_kw={'height_ratios': [2, 1]}
        )
    else:
        fig, ax = plt.subplots(figsize=(12, 5))
        ax_cons = None
        
    # Draw TSS on both tracks
    ax.axvline(x=0, color='black', linestyle='--', label='TSS')
    if ax_cons is not None:
        ax_cons.axvline(x=0, color='black', linestyle='--')
        
    ax.set_xlim(-promoter_up, promoter_down)
    ax.set_ylim(0, 2)
    ax.set_ylabel("Epigenetic / Binding Tracks")
    ax.set_yticks([])
    
    # Plot Shaded Synergy Hotspots
    if df_synergy is not None and not df_synergy.empty:
        for idx, row in df_synergy.iterrows():
            start_rel = row['start'] - promoter_up
            end_rel = row['end'] - promoter_up
            ax.axvspan(start_rel, end_rel, color='orange', alpha=0.15, label='Synergy Hotspot' if idx == 0 else "")
            ax.text((start_rel + end_rel) / 2, 0.2, "Synergy", color='darkorange', fontsize=7, ha='center', fontweight='bold')
            
    # Plot Motifs
    if not df_motifs.empty:
        for _, row in df_motifs.iterrows():
            start_rel = row['start'] - promoter_up
            end_rel = row['end'] - promoter_up
            ax.plot([start_rel, end_rel], [1.0, 1.0], linewidth=5, color='blue', alpha=0.5)
            ax.text((start_rel+end_rel)/2, 1.1, row['tf_name'], ha='center', fontsize=8, color='blue', fontweight='semibold')
            
    # Plot ENCODE peaks
    if not df_encode.empty:
        window_peaks = df_encode[
            (df_encode['distance_to_tss'] >= -promoter_up) & 
            (df_encode['distance_to_tss'] <= promoter_down)
        ].copy()
        
        if not window_peaks.empty:
            for _, row in window_peaks.iterrows():
                dist = row['distance_to_tss']
                # Plot base peak
                ax.plot(dist, 1.3, marker='v', color='lightcoral', markersize=6, linestyle='None', alpha=0.6)
            
            # High-Confidence Peaks
            high_conf_peaks = pd.DataFrame()
            if 'high_confidence_site' in window_peaks.columns:
                high_conf_peaks = window_peaks[window_peaks['high_confidence_site'] == True]
            elif top_n is not None:
                if 'signal' in window_peaks.columns and window_peaks['signal'].max() > 0:
                    high_conf_peaks = window_peaks.sort_values('signal', ascending=False).head(top_n)
                elif 'score' in window_peaks.columns:
                    high_conf_peaks = window_peaks.sort_values('score', ascending=False).head(top_n)
            
            if not high_conf_peaks.empty:
                for _, row in high_conf_peaks.iterrows():
                    dist = row['distance_to_tss']
                    ax.plot(dist, 1.35, marker='v', color='red', markersize=10, linestyle='None')
                    val = row.get('signal', row.get('score', ''))
                    if val:
                        ax.text(dist, 1.45, f"{val:.1f}", ha='center', fontsize=7, color='red', rotation=90)
                        
            # Build legend
            legend_elements = [
                Line2D([0], [0], color='black', linestyle='--', label='TSS'),
                Line2D([0], [0], marker='v', color='lightcoral', label='ENCODE Peak', linestyle='None'),
                Line2D([0], [0], marker='v', color='red', label='High-Conf Peak (Conserved)', linestyle='None', markersize=10)
            ]
            if not df_motifs.empty:
                legend_elements.append(Line2D([0], [0], color='blue', lw=4, label='Motif', alpha=0.5))
            if df_synergy is not None and not df_synergy.empty:
                from matplotlib.patches import Patch
                legend_elements.append(Patch(facecolor='orange', edgecolor='none', alpha=0.2, label='Synergy Hotspot'))
                
            ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
            
    # Plot Evolutionary Conservation
    if ax_cons is not None:
        # Sort by distance
        df_cons_sorted = df_cons.sort_values("distance_to_tss")
        x_vals = df_cons_sorted["distance_to_tss"].values
        phast_vals = df_cons_sorted["phastCons"].values
        phylo_vals = df_cons_sorted["phyloP"].values
        
        # Plot phastCons as filled area
        ax_cons.fill_between(x_vals, phast_vals, color='teal', alpha=0.3, label='phastCons (100-way)')
        ax_cons.plot(x_vals, phast_vals, color='teal', lw=1.5)
        
        ax_cons.set_ylabel("phastCons score")
        ax_cons.set_ylim(0, 1.1)
        ax_cons.set_xlabel("Distance to TSS (bp)")
        
        # High conservation shading threshold line
        ax_cons.axhline(y=0.8, color='darkred', linestyle=':', alpha=0.5, label='High Cons. Threshold (0.8)')
        ax_cons.legend(loc='upper right', fontsize=8)
        ax_cons.set_title("Evolutionary Conservation (UCSC)", fontsize=10)
        
    ax.set_title(f"TF Binding Sites and Epigenetic Tracks on {gene_name} Promoter")
    plt.tight_layout()
    return fig

def plot_promoter_track(df_encode, df_motifs, promoter_up, promoter_down, gene_name, output_dir):
    """Generates and saves the promoter track plot."""
    cons_path = os.path.join(output_dir, f"{gene_name}_conservation.csv")
    synergy_path = os.path.join(output_dir, f"{gene_name}_synergy_hotspots.csv")
    
    df_cons = None
    if os.path.exists(cons_path):
        df_cons = pd.read_csv(cons_path)
        
    df_synergy = None
    if os.path.exists(synergy_path):
        df_synergy = pd.read_csv(synergy_path)
        
    fig = create_promoter_track_fig(df_encode, df_motifs, promoter_up, promoter_down, gene_name, df_cons=df_cons, df_synergy=df_synergy)
    fig.savefig(os.path.join(output_dir, f"{gene_name}_track_plot.png"))
    plt.close(fig)

def generate_bed_output(df_encode, df_motifs, chrom, promoter_start, gene_name, output_dir, strand="+", promoter_end=None):
    """Generates a BED file of all binding sites."""
    bed_path = os.path.join(output_dir, f"{gene_name}_overlaps.bed")
    
    with open(bed_path, 'w') as f:
        # Header
        f.write(f"track name='{gene_name}_TF_Sites' description='TF binding sites for {gene_name}' visibility=2 itemRgb='On'\n")
        
        # ENCODE Peaks
        if not df_encode.empty:
            for _, row in df_encode.iterrows():
                # BED format: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
                # Overlapping peaks in Red, others in Gray
                color = "255,0,0" if row.get('overlap', False) else "128,128,128"
                name = f"{row['tf']}_{row['experiment']}"
                
                # Ensure coordinates are integers
                try:
                    start = int(row['peak_start'])
                    end = int(row['peak_end'])
                    f.write(f"chr{row['peak_chrom']}\t{start}\t{end}\t{name}\t0\t.\t{start}\t{end}\t{color}\n")
                except (ValueError, KeyError):
                    continue

        # Motif Predictions
        if not df_motifs.empty:
            for _, row in df_motifs.iterrows():
                color = "0,0,255"  # Blue for motifs
                name = f"{row['tf_name']}_{row['jaspar_id']}"
                try:
                    if strand == "-" and promoter_end is not None:
                        # For negative-strand genes, offset 0 is promoter_end (5' end of transcript)
                        m_start = int(promoter_end) - int(row['end'])
                        m_end = int(promoter_end) - int(row['start'])
                        # A motif on transcript's '+' strand maps to genomic '-' strand
                        m_strand = '-' if str(row.get('strand', '+')) == '+' else '+'
                    else:
                        m_start = int(promoter_start) + int(row['start'])
                        m_end = int(promoter_start) + int(row['end'])
                        m_strand = str(row.get('strand', '+'))
                    score_val = float(row.get('score', 0.0))
                    bed_score = max(0, min(1000, int(score_val * 10)))
                    clean_chr = str(chrom).replace('chr', '')
                    f.write(f"chr{clean_chr}\t{m_start}\t{m_end}\t{name}\t{bed_score}\t{m_strand}\t{m_start}\t{m_end}\t{color}\n")
                except Exception:
                    continue

def calculate_distance_to_tss(peak_center: int, tss: int, strand: str) -> int:
    """Calculates signed distance from peak center to TSS (upstream is negative)."""
    dist = peak_center - tss
    if strand == "-":
        dist = -dist
    return dist

def find_overlaps_for_experiment(file_path, chrom, p_start_strict, p_end_strict, p_start_loose, p_end_loose, tss, strand, exp, tf):
    """
    Parses a BED file and finds peaks overlapping the promoter region (strict and loose).
    Returns: (hits_strict, hits_loose, total_peaks, closest_peak)
    """
    import gzip
    
    hits_strict = []
    hits_loose = []
    closest_peak = None
    min_abs_dist = float('inf')
    total_peaks = 0
    
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            
            # Basic parsing
            p_chrom = parts[0].replace("chr", "")
            try:
                peak_start = int(parts[1])
                peak_end = int(parts[2])
                
                # Parse score (col 5) and signalValue (col 7) if available
                # BED is 0-indexed for columns in python split
                score = 0
                signal = 0.0
                
                if len(parts) > 4:
                    try:
                        score = int(parts[4])
                    except ValueError:
                        pass
                
                if len(parts) > 6:
                    try:
                        signal = float(parts[6])
                    except ValueError:
                        pass
                        
            except ValueError:
                continue
                
            # Check chrom match
            if p_chrom != chrom.replace("chr", ""):
                continue
                
            total_peaks += 1
            
            # Calculate signed distance from peak center to TSS
            peak_center = (peak_start + peak_end) // 2
            dist_to_tss = calculate_distance_to_tss(peak_center, tss, strand)
            
            abs_dist = abs(dist_to_tss)

            # Track closest peak
            if abs_dist < min_abs_dist:
                min_abs_dist = abs_dist
                closest_peak = {
                    "tf": tf,
                    "experiment": exp['dataset_accession'],
                    "biosample": exp.get('biosample', 'Unknown'),
                    "file": exp['file_accession'],
                    "peak_chrom": p_chrom,
                    "peak_start": peak_start,
                    "peak_end": peak_end,
                    "distance_to_tss": dist_to_tss,
                    "score": score,
                    "signal": signal,
                    "overlap": False
                }

            # Check Strict Overlap
            if max(p_start_strict, peak_start) < min(p_end_strict, peak_end):
                hit_data = {
                    "tf": tf,
                    "experiment": exp['dataset_accession'],
                    "biosample": exp.get('biosample', 'Unknown'),
                    "file": exp['file_accession'],
                    "peak_chrom": p_chrom,
                    "peak_start": peak_start,
                    "peak_end": peak_end,
                    "distance_to_tss": dist_to_tss,
                    "score": score,
                    "signal": signal,
                    "overlap": True
                }
                hits_strict.append(hit_data)
                
            # Check Loose Overlap
            if max(p_start_loose, peak_start) < min(p_end_loose, peak_end):
                hit_data = {
                    "tf": tf,
                    "experiment": exp['dataset_accession'],
                    "biosample": exp.get('biosample', 'Unknown'),
                    "file": exp['file_accession'],
                    "peak_chrom": p_chrom,
                    "peak_start": peak_start,
                    "peak_end": peak_end,
                    "distance_to_tss": dist_to_tss,
                    "score": score,
                    "signal": signal,
                    "overlap": True # It overlaps the loose window
                }
                hits_loose.append(hit_data)
                
    return hits_strict, hits_loose, total_peaks, closest_peak
