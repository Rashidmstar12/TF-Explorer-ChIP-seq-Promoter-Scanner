import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import yaml
import time
from datetime import datetime
from typing import List, Dict
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
    progress_callback=None
):
    """
    Main analysis pipeline.
    """
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
            exps = encode_client.search_encode_tf_chipseq(tf, "Homo sapiens")
            for exp in exps:
                exp['tf_name'] = tf # Keep track of TF name
                all_experiments.append(exp)
    total_files = len(all_experiments)
    logger.info(f"Found total {total_files} peak files to process.")
    
    # Use cache_dir if provided, else use temp_downloads inside output_dir
    if cache_dir:
        download_dir = cache_dir
    else:
        download_dir = os.path.join(output_dir, "temp_downloads")
    
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
                encode_hits.append(closest)
                
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
                
    # Save ENCODE hits
    if not encode_hits:
        logger.info("No ENCODE hits found. Creating empty DataFrame with columns.")
        df_encode = pd.DataFrame(columns=["tf", "experiment", "biosample", "file", "peak_chrom", "peak_start", "peak_end", "distance_to_tss", "overlap"])
    else:
        df_encode = pd.DataFrame(encode_hits)
    logger.info(f"Writing ENCODE hits to {os.path.join(output_dir, f'{gene_name}_encode_hits.csv')}")
    df_encode.to_csv(os.path.join(output_dir, f"{gene_name}_encode_hits.csv"), index=False)
    
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
        generate_bed_output(df_encode, df_motifs, chrom, promoter_genomic_start, gene_name, output_dir)

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


def create_promoter_track_fig(df_encode, df_motifs, promoter_up, promoter_down, gene_name, high_confidence_threshold=None, top_n=None):
    """Creates a matplotlib figure for the promoter track."""
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Draw TSS
    ax.axvline(x=0, color='black', linestyle='--', label='TSS')
    ax.set_xlim(-promoter_up, promoter_down)
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_yticks([])
    ax.set_ylim(0, 2)
    
    # Plot Motifs
    if not df_motifs.empty:
        for _, row in df_motifs.iterrows():
            start_rel = row['start'] - promoter_up
            end_rel = row['end'] - promoter_up
            ax.plot([start_rel, end_rel], [1, 1], linewidth=4, color='blue', alpha=0.5)
            ax.text((start_rel+end_rel)/2, 1.1, row['tf_name'], ha='center', fontsize=8, color='blue')
            
    # Plot ENCODE peaks
    if not df_encode.empty:
        # Filter peaks within window
        window_peaks = df_encode[
            (df_encode['distance_to_tss'] >= -promoter_up) & 
            (df_encode['distance_to_tss'] <= promoter_down)
        ].copy()
        
        if not window_peaks.empty:
            # Base track (all peaks)
            # Plot as light red triangles
            for _, row in window_peaks.iterrows():
                dist = row['distance_to_tss']
                ax.plot(dist, 1.2, marker='v', color='lightcoral', markersize=6, linestyle='None', alpha=0.6)
            
            # High-Confidence Track
            high_conf_peaks = pd.DataFrame()
            
            # Determine high confidence
            if top_n is not None:
                # Sort by signal (col 7) or score (col 5)
                # We stored 'signal' and 'score' in df_encode
                if 'signal' in window_peaks.columns and window_peaks['signal'].max() > 0:
                    high_conf_peaks = window_peaks.sort_values('signal', ascending=False).head(top_n)
                elif 'score' in window_peaks.columns:
                    high_conf_peaks = window_peaks.sort_values('score', ascending=False).head(top_n)
                else:
                    high_conf_peaks = window_peaks.head(top_n) # Fallback
                    
            elif high_confidence_threshold is not None:
                if 'signal' in window_peaks.columns and window_peaks['signal'].max() > 0:
                     high_conf_peaks = window_peaks[window_peaks['signal'] >= high_confidence_threshold]
                elif 'score' in window_peaks.columns:
                     high_conf_peaks = window_peaks[window_peaks['score'] >= high_confidence_threshold]
            
            # Plot High Confidence
            if not high_conf_peaks.empty:
                for _, row in high_conf_peaks.iterrows():
                    dist = row['distance_to_tss']
                    ax.plot(dist, 1.25, marker='v', color='red', markersize=10, linestyle='None', label='High-Conf Peak')
                    # Optional: Add label for signal
                    val = row.get('signal', row.get('score', ''))
                    if val:
                        ax.text(dist, 1.35, f"{val:.1f}", ha='center', fontsize=7, rotation=90)

            # Legend
            # Create custom legend handles
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color='black', linestyle='--', label='TSS'),
                Line2D([0], [0], marker='v', color='lightcoral', label='ENCODE Peak', linestyle='None'),
                Line2D([0], [0], marker='v', color='red', label='High-Conf Peak', linestyle='None', markersize=10)
            ]
            if not df_motifs.empty:
                legend_elements.append(Line2D([0], [0], color='blue', lw=4, label='Motif', alpha=0.5))
                
            ax.legend(handles=legend_elements, loc='upper right')

    ax.set_title(f"TF Binding Sites on {gene_name} Promoter")
    plt.tight_layout()
    return fig

def plot_promoter_track(df_encode, df_motifs, promoter_up, promoter_down, gene_name, output_dir):
    """Generates and saves the promoter track plot."""
    fig = create_promoter_track_fig(df_encode, df_motifs, promoter_up, promoter_down, gene_name)
    fig.savefig(os.path.join(output_dir, f"{gene_name}_track_plot.png"))
    plt.close(fig)

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
            dist_to_tss = peak_center - tss
            if strand == "-":
                dist_to_tss = -dist_to_tss
            
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
