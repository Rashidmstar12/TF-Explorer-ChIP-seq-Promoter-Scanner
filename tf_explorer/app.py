import streamlit as st
import pandas as pd
import os
import sys
import tempfile
import shutil

# Ensure we import the local tf_explorer package, not the installed one
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tf_explorer import analysis, motifs, comparative
import importlib
importlib.reload(motifs)
importlib.reload(analysis)
importlib.reload(comparative)

# Set page configuration
st.set_page_config(
    page_title="TF-Explorer v1.1",
    page_icon="🧬",
    layout="wide"
)


def summarize_experiments(df_stats):
    """
    Standardizes experiment metrics calculation.
    Input: df_stats (raw DataFrame from experiment_stats.csv)
    Returns: (n_strict, n_loose, total_strict_peaks, total_loose_peaks, total_experiments, processed_df)
    """
    # 1. Rename columns to standardized names if needed
    # Map CSV columns to internal names
    col_map = {
        'num_strict_peaks': 'num_overlapping_peaks_strict',
        'num_loose_peaks': 'num_overlapping_peaks_loose',
        'total_peaks_in_file': 'total_peaks'
    }
    df = df_stats.rename(columns=col_map)
    
    # Ensure required columns exist
    for col in ['num_overlapping_peaks_strict', 'num_overlapping_peaks_loose']:
        if col not in df.columns:
            df[col] = 0
            
    # 2. Count files directly (no grouping)
    # User wants to count every file/sample as an "experiment"
    
    n_strict = (df['num_overlapping_peaks_strict'] > 0).sum()
    n_loose = (df['num_overlapping_peaks_loose'] > 0).sum()
            
    total_strict_peaks = df['num_overlapping_peaks_strict'].sum()
    total_loose_peaks = df['num_overlapping_peaks_loose'].sum()
    total_experiments = len(df) # Total files
    
    return n_strict, n_loose, total_strict_peaks, total_loose_peaks, total_experiments, df

def plot_cell_line_comparison(summary_df):
    """
    Generates a bar chart comparing binding frequencies across cell lines.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if summary_df.empty or 'biosample' not in summary_df.columns:
        return None
        
    # Group by biosample
    # Calculate % of files with strict peaks
    stats = summary_df.groupby('biosample').agg(
        total_files=('file_accession', 'count'),
        files_with_peaks=('num_strict_peaks', lambda x: (x > 0).sum())
    ).reset_index()
    
    stats['binding_rate'] = (stats['files_with_peaks'] / stats['total_files']) * 100
    
    # Sort by binding rate
    stats = stats.sort_values('binding_rate', ascending=False)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=stats, x='biosample', y='binding_rate', ax=ax, palette='viridis')
    
    ax.set_title("TF Binding Rate by Cell Line (Strict Window)", fontsize=14)
    ax.set_ylabel("Binding Rate (% of Files)", fontsize=12)
    ax.set_xlabel("Cell Line", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    return fig

def main():
    st.title("🧬 TF-Explorer v1.2 – ChIP-seq Promoter Scanner")
    st.markdown("""
    Analyze transcription factor binding sites for **ANY GENE**.
    1. **Search** for ENCODE ChIP-seq experiments.
    2. **Select** the experiments you want to analyze.
    3. **Run** the analysis (with optional JASPAR motif prediction).
    """)

    # --- Sidebar ---
    st.sidebar.header("Configuration")

    # Gene Input
    gene_name = st.sidebar.text_input("Gene Symbol", value="", help="Enter the official gene symbol (e.g., TP53, MYC).")

    # TF Input
    default_tfs = ""
    tf_input = st.sidebar.text_area("Transcription Factors", value=default_tfs, help="Comma-separated list of TFs to search in ENCODE.")
    
    # JASPAR Input
    default_jaspar = ""
    jaspar_input = st.sidebar.text_area("JASPAR Matrix IDs", value=default_jaspar, help="Comma-separated list of JASPAR IDs for motif scanning.")
    
    include_jaspar = st.sidebar.checkbox("Include JASPAR Motif Predictions", value=False, help="Check to scan promoter for JASPAR motifs.")

    # Advanced Settings
    with st.sidebar.expander("Advanced Settings"):
        promoter_up = st.number_input("Promoter Upstream (bp)", value=2000, step=100)
        promoter_down = st.number_input("Promoter Downstream (bp)", value=500, step=100)
        pwm_threshold = st.slider("PWM Threshold", min_value=0.0, max_value=20.0, value=8.0, step=0.1)
        genome = st.selectbox("Genome Assembly", ["hg38"], index=0)
        force_download = st.checkbox("Force Re-download Files", value=False, help="Check to delete cached files and re-download them. Use if you suspect corruption.")
    
    # --- Session State ---
    if "encode_results" not in st.session_state:
        st.session_state.encode_results = []
    if "search_performed" not in st.session_state:
        st.session_state.search_performed = False
    if "editor_key" not in st.session_state:
        st.session_state.editor_key = 0

    # Step 1: Search
    if st.sidebar.button("Search ENCODE", type="primary"):
        if not tf_input:
            st.error("Please enter at least one Transcription Factor.")
        else:
            tf_list = [tf.strip() for tf in tf_input.split(",") if tf.strip()]
            st.session_state.encode_results = []
            
            with st.spinner("Searching ENCODE..."):
                from tf_explorer import encode_client
                all_results = []
                for tf in tf_list:
                    # Map genome to organism (simple mapping for now)
                    organism = "Homo sapiens" 
                    results = encode_client.search_encode_tf_chipseq(tf, organism)
                    for res in results:
                        # Filter by selected genome assembly
                        # Handle GRCh38 as equivalent to hg38
                        res_assembly = res.get('assembly', '')
                        if res_assembly == genome or (genome == 'hg38' and res_assembly == 'GRCh38'):
                            res['tf_name'] = tf
                            res['Select'] = False # Default selection state
                            all_results.append(res)
                
                st.session_state.encode_results = all_results
                st.session_state.search_performed = True
                
    # Step 2: Select & Analyze
    if st.session_state.search_performed:
        st.subheader("Found Experiments")
        
        if not st.session_state.encode_results:
            st.warning("No experiments found for the specified TFs.")
        else:
            # Filter by Cell Line
            all_biosamples = sorted(list(set(item['biosample'] for item in st.session_state.encode_results)))
            selected_biosamples = st.multiselect("Filter by Cell Line", all_biosamples, default=all_biosamples)
            
            # Filter results based on selection
            filtered_results = [item for item in st.session_state.encode_results if item['biosample'] in selected_biosamples]
            
            # Select All Button
            if st.button("Select All Filtered Experiments"):
                for item in st.session_state.encode_results:
                    if item['biosample'] in selected_biosamples:
                        item['Select'] = True
                    else:
                        item['Select'] = False
                st.session_state.editor_key += 1
                st.rerun()

            # Display data editor for selection
            df = pd.DataFrame(filtered_results)
            
            # Reorder columns for better visibility
            cols = ['Select', 'tf_name', 'biosample', 'description', 'dataset_accession', 'file_accession', 'assembly']
            # Ensure all cols exist
            cols = [c for c in cols if c in df.columns]
            
            edited_df = st.data_editor(
                df[cols],
                key=f"editor_{st.session_state.editor_key}",
                column_config={
                    "Select": st.column_config.CheckboxColumn(
                        "Select",
                        help="Select experiments to analyze",
                        default=False,
                    )
                },
                disabled=["tf_name", "biosample", "dataset_accession", "file_accession", "assembly"],
                hide_index=True,
                use_container_width=True
            )
            
            # Get selected experiments
            selected_experiments = []
            for idx, row in edited_df.iterrows():
                if row['Select']:
                    acc = row['file_accession']
                    original_obj = next((item for item in st.session_state.encode_results if item['file_accession'] == acc), None)
                    if original_obj:
                        selected_experiments.append(original_obj)
            
            st.write(f"Selected **{len(selected_experiments)}** experiments.")
            
            if st.button("Run Analysis on Selected"):
                if not gene_name:
                    st.error("Please enter a gene symbol.")
                elif not selected_experiments:
                    st.error("Please select at least one experiment.")
                else:
                    # Prepare JASPAR IDs
                    jaspar_ids = []
                    if include_jaspar:
                        raw_jaspar = [jid.strip() for jid in jaspar_input.split(",") if jid.strip()]
                        from tf_explorer import motifs
                        for item in raw_jaspar:
                            if item.upper() in motifs.COMMON_TFS:
                                jaspar_ids.append(motifs.COMMON_TFS[item.upper()])
                            else:
                                jaspar_ids.append(item)
                    
                    # Run Analysis
                    cache_dir = os.path.join(os.getcwd(), "encode_cache")
                    os.makedirs(cache_dir, exist_ok=True)
                    out_dir = tempfile.mkdtemp(prefix="tf_explorer_")
                    
                    status_container = st.status("Running Analysis...", expanded=True)
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    def progress_callback(current, total, message):
                        progress = current / total if total > 0 else 0
                        progress_bar.progress(progress)
                        status_text.text(f"{message} ({current}/{total})")
                    
                    try:
                        with status_container:
                            st.write(f"Analyzing **{gene_name}** with {len(selected_experiments)} experiments...")
                            
                            summary_df = analysis.analyze_gene(
                                gene_name=gene_name,
                                tf_list=[tf_input],
                                experiments_list=selected_experiments,
                                genome=genome,
                                promoter_up=promoter_up,
                                promoter_down=promoter_down,
                                loose_promoter_up=5000,
                                loose_promoter_down=5000,
                                pwm_threshold=pwm_threshold,
                                output_dir=out_dir,
                                cache_dir=cache_dir,
                                plot_track=True,
                                bed_output=True,
                                force_download=force_download,
                                progress_callback=progress_callback
                            )
                    except ValueError as e:
                        st.error(str(e))
                        status_container.update(label="Analysis Failed", state="error", expanded=True)
                        st.stop()
                    except Exception as e:
                        import traceback
                        st.error(f"An unexpected error occurred: {e}")
                        st.subheader("Traceback")
                        st.code(traceback.format_exc())
                        status_container.update(label="Analysis Failed", state="error", expanded=True)
                        st.stop()

                    if summary_df is None:
                        st.error("Analysis failed. Please check logs.")
                        status_container.update(label="Analysis Failed", state="error", expanded=True)
                        log_path = os.path.join(out_dir, "analysis_log.txt")
                        if os.path.exists(log_path):
                            with open(log_path, "r") as f:
                                st.code(f.read())
                        st.stop()

                    else:
                        # Store results in session state
                        st.session_state['analysis_results'] = {
                            'summary_df': summary_df,
                            'out_dir': out_dir,
                            'gene_name': gene_name,
                            'tf_input': tf_input,
                            'promoter_up': promoter_up,
                            'promoter_down': promoter_down,
                            'total_checked': len(selected_experiments)
                        }
                        
                        st.write("Analysis complete!")
                        progress_bar.progress(1.0)
                        status_text.text("Done!")
                        status_container.update(label="Analysis Complete!", state="complete", expanded=False)
            
            # Display Results (Persistent)
            if 'analysis_results' in st.session_state:
                results = st.session_state['analysis_results']
                summary_df = results['summary_df']
                out_dir = results['out_dir']
                gene_name = results['gene_name']
                tf_input = results['tf_input']
                promoter_up = results['promoter_up']
                promoter_down = results['promoter_down']
                total_checked = results['total_checked']
                
                # Always show logs for debugging
                log_path = os.path.join(out_dir, "analysis_log.txt")
                if os.path.exists(log_path):
                    with open(log_path, "r") as f:
                        st.sidebar.download_button("Download Analysis Log", f, file_name="analysis_log.txt")

                st.subheader("Analysis Results")
                
                # Display Summary Metrics
                # Load stats if available
                stats_path = os.path.join(out_dir, "{}_experiment_stats.csv".format(gene_name))
                
                n_strict = 0
                n_loose = 0
                total_strict_peaks = 0
                total_loose_peaks = 0
                
                if os.path.exists(stats_path):
                    stats_df_raw = pd.read_csv(stats_path)
                    # Use shared helper function
                    n_strict, n_loose, total_strict_peaks, total_loose_peaks, total_experiments, stats_df = summarize_experiments(stats_df_raw)
                else:
                    # Fallback if stats file missing (shouldn't happen if analysis succeeded)
                    stats_df = pd.DataFrame()
                    total_experiments = total_checked
                
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("Files Analyzed", "{}".format(total_experiments))
                col2.metric("Files with Promoter Peaks (strict)", "{} / {}".format(n_strict, total_experiments))
                col3.metric("Files with Loose Peaks (±5kb)", "{} / {}".format(n_loose, total_experiments))
                col4.metric("Overlapping Peak Regions (All Exp)", "{}".format(total_strict_peaks))

                row1 = summary_df.iloc[0]

                # Display Window Coordinates
                st.markdown("### Window Coordinates")
                st.markdown(f"**Promoter Window (Strict):** `{row1.get('promoter_window', 'N/A')}`")
                st.markdown(f"**Loose Window (±5kb):** `{row1.get('loose_window', 'N/A')}`")

                # Generate Report Text
                biosamples_str = row1['biosamples_with_peaks']
                coords_str = row1['coords']
                
                if n_strict > 0:
                    st.success("**Positive Result:** Found binding evidence in **{}** experiments (Strict Mode).".format(n_strict))
                    
                    report_text = f"""Analysis Report for {gene_name}
--------------------------------------------------
Target: {gene_name}
Transcription Factors: {tf_input}
Genome: {genome}
Promoter Region: -{promoter_up} bp to +{promoter_down} bp relative to TSS

Conclusion:
Positive result. Binding evidence found in {n_strict} experiments.

Biological Context:
Binding was observed in the following biosamples (cell lines/tissues):
{biosamples_str}

Evidence:
- {n_strict} independent ChIP-seq experiments confirm this interaction.
- The binding sites are located at {coords_str}.
"""
                    st.markdown(f"""
                    This analysis provides **supporting evidence** that **{tf_input}** binds to the **{gene_name}** promoter.
                    - **{n_strict}** independent ChIP-seq experiments confirm this interaction.
                    - **Biosamples:** {biosamples_str}
                    - The binding sites are located at **{coords_str}**.
                    """)
                    
                    # Detailed Per-Experiment Table
                    if total_checked > 1 and os.path.exists(stats_path):
                        st.subheader("Detailed Experiment Statistics")
                        st.markdown("Check **num_overlapping_peaks_loose** to see if peaks were found nearby (±5kb).")
                        
                        display_df = stats_df.rename(columns={
                            'num_strict_peaks': 'num_overlapping_peaks_strict',
                            'num_loose_peaks': 'num_overlapping_peaks_loose',
                            'total_peaks_in_file': 'total_peaks'
                        })
                        
                        cols_to_show = ['experiment_id', 'biosample', 'file_accession', 'num_overlapping_peaks_strict', 'num_overlapping_peaks_loose', 'total_peaks']
                        cols_to_show = [c for c in cols_to_show if c in display_df.columns]
                        
                        st.dataframe(display_df[cols_to_show])
                        with open(stats_path, "rb") as f:
                            st.download_button("Download Stats CSV", f, file_name=f"{gene_name}_experiment_stats.csv", mime="text/csv")
                    
                    if total_checked == 1 and os.path.exists(stats_path):
                         with open(stats_path, "rb") as f:
                            st.download_button("Download Stats CSV", f, file_name=f"{gene_name}_experiment_stats.csv", mime="text/csv")


                else:
                    st.warning(f"**Negative Result:** No binding sites found in any of the **{total_checked}** analyzed experiments.")
                    report_text = f"""Analysis Report for {gene_name}
--------------------------------------------------
Target: {gene_name}
Transcription Factors: {tf_input}
Genome: {genome}
Promoter Region: -{promoter_up} bp to +{promoter_down} bp relative to TSS

Conclusion:
Negative result. No binding sites found in any of the {total_checked} analyzed experiments.

Recommendation:
- Consider checking different cell lines or experimental conditions.
- Verify that the promoter region ({promoter_up}bp upstream) covers the expected binding site.
"""
                    st.markdown(f"""
                    This suggests that **{tf_input}** may not bind to the **{gene_name}** promoter under the conditions tested in these experiments.
                    - Consider checking different cell lines or experimental conditions.
                    - Verify that the promoter region ({promoter_up}bp upstream) covers the expected binding site.
                    """)
                
                st.download_button("Download Report", report_text, file_name=f"{gene_name}_analysis_report.txt")

                # Plots
                st.subheader("Visualizations")
                
                tab1, tab2, tab3, tab4, tab5 = st.tabs(["Promoter Track", "Cell Line Comparison", "TF Comparison", "Biosample Distribution", "TF Binding Counts"])
                
                with tab1:
                    # ... (Existing Promoter Track Code) ...
                    # Interactive Promoter Track
                    hits_path = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    motifs_path = os.path.join(out_dir, f"{gene_name}_motif_predictions.csv")
                    
                    df_encode_plot = pd.DataFrame()
                    df_motifs_plot = pd.DataFrame()
                    
                    if os.path.exists(hits_path):
                        df_encode_plot = pd.read_csv(hits_path)
                    if os.path.exists(motifs_path):
                        df_motifs_plot = pd.read_csv(motifs_path)
                        
                    # UI Controls for High-Peak Track
                    st.markdown("### High-Peak Visualization")
                    st.markdown("Overlay high-confidence peaks on the promoter track to see where binding is strongest.")
                    
                    col_ctrl1, col_ctrl2 = st.columns(2)
                    with col_ctrl1:
                        show_high_conf = st.checkbox("Highlight High-Confidence Peaks", value=False)
                    
                    top_n = None
                    threshold = None
                    
                    if show_high_conf:
                        with col_ctrl2:
                            filter_mode = st.radio("Filter By:", ["Top N Peaks", "Signal/Score Threshold"], horizontal=True)
                            if filter_mode == "Top N Peaks":
                                top_n = st.slider("Number of Top Peaks", min_value=1, max_value=50, value=10)
                            else:
                                max_sig = 100.0
                                if not df_encode_plot.empty and 'signal' in df_encode_plot.columns:
                                    max_sig = df_encode_plot['signal'].max()
                                    if pd.isna(max_sig) or max_sig == 0:
                                        max_sig = 100.0
                                threshold = st.number_input("Minimum Signal/Score", min_value=0.0, max_value=float(max_sig), value=float(max_sig/2))

                    # Generate Plot
                    if not df_encode_plot.empty or not df_motifs_plot.empty:
                        fig = analysis.create_promoter_track_fig(
                            df_encode_plot, df_motifs_plot, 
                            promoter_up, promoter_down, gene_name,
                            high_confidence_threshold=threshold,
                            top_n=top_n
                        )
                        st.pyplot(fig)
                        
                        if show_high_conf:
                            st.caption("High-Confidence Peaks Table")
                            window_peaks = df_encode_plot[
                                (df_encode_plot['distance_to_tss'] >= -promoter_up) & 
                                (df_encode_plot['distance_to_tss'] <= promoter_down)
                            ].copy()
                            
                            high_conf_table = pd.DataFrame()
                            if top_n is not None:
                                if 'signal' in window_peaks.columns and window_peaks['signal'].max() > 0:
                                    high_conf_table = window_peaks.sort_values('signal', ascending=False).head(top_n)
                                elif 'score' in window_peaks.columns:
                                    high_conf_table = window_peaks.sort_values('score', ascending=False).head(top_n)
                            elif threshold is not None:
                                if 'signal' in window_peaks.columns and window_peaks['signal'].max() > 0:
                                    high_conf_table = window_peaks[window_peaks['signal'] >= threshold]
                                elif 'score' in window_peaks.columns:
                                    high_conf_table = window_peaks[window_peaks['score'] >= threshold]
                                    
                            if not high_conf_table.empty:
                                cols = ['peak_chrom', 'peak_start', 'peak_end', 'experiment', 'biosample', 'distance_to_tss']
                                if 'signal' in high_conf_table.columns: cols.append('signal')
                                if 'score' in high_conf_table.columns: cols.append('score')
                                
                                st.dataframe(high_conf_table[cols].sort_values('distance_to_tss'))
                            else:
                                st.info("No peaks meet the selected criteria.")
                                
                    else:
                        st.warning("No data to plot.")

                with tab2:
                    st.markdown("### Cell Line Comparison")
                    
                    # Load detailed stats to get cell lines
                    stats_path_comp = os.path.join(out_dir, f"{gene_name}_experiment_stats.csv")
                    stats_df_comp = pd.DataFrame()
                    if os.path.exists(stats_path_comp):
                        stats_df_comp = pd.read_csv(stats_path_comp)
                    
                    # 1. Multi-Cell Line Overview (Existing)
                    st.markdown("#### Overview: Binding Rates by Cell Line")
                    # Use detailed stats instead of summary_df
                    fig_comp = plot_cell_line_comparison(stats_df_comp)
                    if fig_comp:
                        st.pyplot(fig_comp)
                    else:
                        st.info("Insufficient data for overview.")
                        
                    st.divider()
                    
                    # 2. Pairwise Comparison (New v1.2)
                    st.markdown("#### Pairwise Comparison (Advanced)")
                    st.markdown("Select two or more cell lines to compare binding patterns in detail.")
                    
                    if not stats_df_comp.empty and 'biosample' in stats_df_comp.columns:
                        available_cls = sorted([x for x in stats_df_comp['biosample'].unique() if x != 'Unknown'])
                        
                        if len(available_cls) >= 2:
                            selected_cls = st.multiselect("Select Cell Lines to Compare", available_cls, default=available_cls[:2])
                            
                            if st.button("Compare Selected Cell Lines"):
                                if len(selected_cls) < 2:
                                    st.error("Please select at least two cell lines.")
                                else:
                                    # Run Comparative Analysis
                                    # Load the hits dataframe first
                                    hits_path_comp = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                                    if os.path.exists(hits_path_comp):
                                        hits_df_comp = pd.read_csv(hits_path_comp)
                                        
                                        # Instantiate with list of cell lines
                                        comp_analysis = comparative.ComparativeAnalysis(hits_df_comp, selected_cls, gene_name, group_by='biosample')
                                        
                                        metrics = comp_analysis.calculate_metrics()
                                        interpretation = comp_analysis.generate_interpretation(metrics)
                                        
                                        st.markdown("##### Results")
                                        
                                        # 1. Multi-Track Plot
                                        st.markdown("**Comparative Binding Tracks**")
                                        fig_tracks = comp_analysis.plot_comparison(promoter_up, promoter_down)
                                        st.pyplot(fig_tracks)
                                        
                                        # 2. Jaccard Heatmap
                                        st.markdown("**Similarity Heatmap (Jaccard Index)**")
                                        fig_heatmap = comp_analysis.plot_jaccard_heatmap(metrics)
                                        st.pyplot(fig_heatmap)
                                        
                                        st.markdown("##### Interpretation")
                                        st.info(interpretation)
                                        
                                        st.markdown("##### Metrics")
                                        with st.expander("Jaccard Similarity Matrix", expanded=True):
                                            st.dataframe(metrics['jaccard_matrix'].style.format("{:.2f}"))
                                            
                                        with st.expander("Unique Binding Sites (bp)", expanded=True):
                                            st.caption("Number of base pairs covered uniquely by each cell line (not shared with others).")
                                            unique_df = pd.DataFrame(list(metrics['unique_bases_counts'].items()), columns=['Cell Line', 'Unique BP'])
                                            st.bar_chart(unique_df.set_index('Cell Line'))
                                            
                                        with st.expander("Full Metrics JSON"):
                                            safe_metrics = {k: v for k, v in metrics.items() if not isinstance(v, pd.DataFrame)}
                                            st.json(safe_metrics)
                                    else:
                                        st.error("Could not load hits data for comparison.")
                        else:
                            st.warning(f"Need at least 2 different cell lines in the results to perform comparison. Found: {available_cls}")
                            st.info("Tip: Make sure you have selected experiments from at least two different cell lines in the 'Select Experiments' table and clicked 'Run Analysis'.")
                    else:
                        st.warning("No biosample data available in the analysis results.")


                with tab3:
                    st.markdown("### TF Comparison")
                    st.markdown("Compare binding patterns of different Transcription Factors on this gene.")
                    
                    hits_path_tf = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    if os.path.exists(hits_path_tf):
                        hits_df_tf = pd.read_csv(hits_path_tf)
                        if not hits_df_tf.empty and 'tf' in hits_df_tf.columns:
                            available_tfs = sorted([x for x in hits_df_tf['tf'].unique()])
                            
                            if len(available_tfs) >= 2:
                                selected_tfs = st.multiselect("Select TFs to Compare", available_tfs, default=available_tfs[:2])
                                
                                if st.button("Compare Selected TFs"):
                                    if len(selected_tfs) < 2:
                                        st.error("Please select at least two TFs.")
                                    else:
                                        # Run Comparative Analysis with group_by='tf'
                                        comp_analysis_tf = comparative.ComparativeAnalysis(hits_df_tf, selected_tfs, gene_name, group_by='tf')
                                        
                                        metrics_tf = comp_analysis_tf.calculate_metrics()
                                        interpretation_tf = comp_analysis_tf.generate_interpretation(metrics_tf)
                                        
                                        st.markdown("##### Results")
                                        
                                        # 1. Multi-Track Plot
                                        st.markdown("**Comparative Binding Tracks**")
                                        fig_tracks_tf = comp_analysis_tf.plot_comparison(promoter_up, promoter_down)
                                        st.pyplot(fig_tracks_tf)
                                        
                                        # 2. Jaccard Heatmap
                                        st.markdown("**Similarity Heatmap (Jaccard Index)**")
                                        fig_heatmap_tf = comp_analysis_tf.plot_jaccard_heatmap(metrics_tf)
                                        st.pyplot(fig_heatmap_tf)
                                        
                                        st.markdown("##### Interpretation")
                                        st.info(interpretation_tf)
                                        
                                        st.markdown("##### Metrics")
                                        with st.expander("Jaccard Similarity Matrix", expanded=True):
                                            st.dataframe(metrics_tf['jaccard_matrix'].style.format("{:.2f}"))
                                            
                                        with st.expander("Unique Binding Sites (bp)", expanded=True):
                                            st.caption("Number of base pairs covered uniquely by each TF (not shared with others).")
                                            unique_tf_df = pd.DataFrame(list(metrics_tf['unique_bases_counts'].items()), columns=['TF', 'Unique BP'])
                                            st.bar_chart(unique_tf_df.set_index('TF'))
                            else:
                                st.warning(f"Need at least 2 different TFs in the results to perform comparison. Found: {available_tfs}")
                        else:
                            st.warning("No TF data available.")
                    else:
                        st.warning("No analysis results found.")

                with tab4:
                    hits_path = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    if os.path.exists(hits_path):
                        hits_df = pd.read_csv(hits_path)
                        if not hits_df.empty and 'biosample' in hits_df.columns:
                            bio_counts = hits_df['biosample'].value_counts()
                            st.bar_chart(bio_counts)
                            st.caption("Number of binding events per biosample (cell line/tissue).")
                        else:
                            st.info("No biosample data available.")

                with tab5:
                    plot_path = os.path.join(out_dir, f"{gene_name}_tf_binding_plot.png")
                    if os.path.exists(plot_path):
                        st.image(plot_path, caption=f"TF Binding Summary for {gene_name}", use_container_width=True)
                    else:
                        st.warning("Plot not generated.")

                # Data Tables
                st.subheader("Detailed Data")
                
                hits_path = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                motifs_path = os.path.join(out_dir, f"{gene_name}_motif_predictions.csv")
                
                if os.path.exists(hits_path):
                    hits_df = pd.read_csv(hits_path)
                    with st.expander("ENCODE ChIP-seq Hits"):
                        st.dataframe(hits_df)
                        with open(hits_path, "rb") as f:
                            st.download_button("Download CSV", f, file_name=f"{gene_name}_encode_hits.csv", mime="text/csv")

                if os.path.exists(motifs_path):
                    motifs_df = pd.read_csv(motifs_path)
                    with st.expander("Motif Predictions"):
                        st.dataframe(motifs_df)
                        with open(motifs_path, "rb") as f:
                            st.download_button("Download CSV", f, file_name=f"{gene_name}_motif_predictions.csv", mime="text/csv")

                # BED File
                bed_path = os.path.join(out_dir, f"{gene_name}_overlaps.bed")
                if os.path.exists(bed_path):
                    with open(bed_path, "rb") as f:
                        st.download_button("Download BED File", f, file_name=f"{gene_name}_overlaps.bed", mime="text/plain")

if __name__ == "__main__":
    main()
