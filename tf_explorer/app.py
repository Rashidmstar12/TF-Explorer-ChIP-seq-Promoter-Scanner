import streamlit as st
import pandas as pd
import os
import sys
import tempfile
import shutil
import base64

# Ensure we import the local tf_explorer package, not the installed one
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tf_explorer import analysis, motifs, comparative, genome_client, tcga_client
import importlib
importlib.reload(motifs)
importlib.reload(analysis)
importlib.reload(comparative)
importlib.reload(genome_client)
importlib.reload(tcga_client)

from tf_explorer import primer_design
importlib.reload(primer_design)

# Set page configuration
st.set_page_config(
    page_title="TF-Explorer v3.0",
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
    plt.tight_layout()
    
    return fig

def get_gtex_cli_path():
    home = os.path.expanduser("~")
    return os.path.join(
        home, 
        ".gemini", "config", "plugins", "science", "skills", 
        "gtex_database", "scripts", "gtex_cli.py"
    )

def get_string_cli_path():
    home = os.path.expanduser("~")
    return os.path.join(
        home, 
        ".gemini", "config", "plugins", "science", "skills", 
        "string_database", "scripts", "string_cli.py"
    )

def fetch_gtex_data(gene_symbol, tf_list):
    import subprocess
    import tempfile
    import json
    
    gtex_cli = get_gtex_cli_path()
    if not os.path.exists(gtex_cli):
        return None, "GTEx Database skill CLI not found."
        
    temp_dir = tempfile.gettempdir()
    symbols = [gene_symbol] + tf_list
    gencode_map = {}
    
    for sym in symbols:
        out_id_file = os.path.join(temp_dir, f"gtex_id_{sym}.json")
        cmd = ["uv", "run", gtex_cli, "resolve-gencode-id", sym, "--output", out_id_file]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True, shell=True)
            if os.path.exists(out_id_file):
                with open(out_id_file, "r") as f:
                    data = json.load(f)
                gencode_map[sym] = data.get("gencode_id")
                try: os.remove(out_id_file)
                except OSError: pass
        except Exception:
            pass
            
    if not gencode_map.get(gene_symbol):
        return None, f"Could not resolve GENCODE ID for target gene '{gene_symbol}'."
        
    expression_dfs = []
    for sym, gid in gencode_map.items():
        if not gid: continue
        out_expr_file = os.path.join(temp_dir, f"gtex_expr_{sym}.json")
        cmd = ["uv", "run", gtex_cli, "get-median-expression", gid, "--output", out_expr_file]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True, shell=True)
            if os.path.exists(out_expr_file):
                with open(out_expr_file, "r") as f:
                    expr_data = json.load(f)
                try: os.remove(out_expr_file)
                except OSError: pass
                
                if isinstance(expr_data, list) and len(expr_data) > 0:
                    rows = []
                    for entry in expr_data:
                        rows.append({
                            "tissue": entry.get("tissueSiteDetailId", "Unknown"),
                            "TPM": entry.get("median", 0.0),
                            "Gene": sym
                        })
                    df_sym = pd.DataFrame(rows)
                    expression_dfs.append(df_sym)
        except Exception:
            pass
            
    if not expression_dfs:
        return None, "Failed to retrieve median expression data from GTEx."
        
    df_combined = pd.concat(expression_dfs, ignore_index=True)
    return df_combined, None

def plot_gtex_expression(df_combined, gene_symbol):
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    target_df = df_combined[df_combined['Gene'] == gene_symbol]
    top_tissues = target_df.sort_values('TPM', ascending=False).head(20)['tissue'].tolist()
    
    df_plot = df_combined[df_combined['tissue'].isin(top_tissues)].copy()
    df_plot['tissue'] = df_plot['tissue'].apply(lambda x: x.replace('_', ' '))
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Beautiful palette
    sns.barplot(
        data=df_plot, 
        x='tissue', 
        y='TPM', 
        hue='Gene', 
        ax=ax, 
        palette='Set2'
    )
    
    ax.set_title(f"Baseline Tissue Expression (GTEx) - Top 20 Tissues for {gene_symbol}", fontsize=14, fontweight='bold')
    ax.set_ylabel("Median Expression (TPM, Log Scale)", fontsize=12)
    ax.set_xlabel("Tissue Site", fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yscale("log")
    plt.tight_layout()
    return fig

def fetch_string_data(gene_symbol, tf_list, species=9606):
    import subprocess
    import tempfile
    import json
    
    string_cli = get_string_cli_path()
    if not os.path.exists(string_cli):
        return None, None, "STRING Database skill CLI not found."
        
    temp_dir = tempfile.gettempdir()
    identifiers = [gene_symbol] + tf_list
    
    out_map_file = os.path.join(temp_dir, "string_mapped.tsv")
    cmd_map = [
        "uv", "run", string_cli, "map",
        "--identifiers"
    ] + identifiers + [
        "--species", str(species),
        "--output", out_map_file
    ]
    
    mapped_ids = []
    try:
        subprocess.run(cmd_map, capture_output=True, text=True, check=True, shell=True)
        if os.path.exists(out_map_file):
            df_map = pd.read_csv(out_map_file, sep='\t')
            if not df_map.empty and 'stringId' in df_map.columns:
                mapped_ids = df_map['stringId'].tolist()
            try: os.remove(out_map_file)
            except OSError: pass
    except Exception:
        mapped_ids = identifiers
        
    if not mapped_ids:
        mapped_ids = identifiers
        
    out_img_file = os.path.join(temp_dir, f"string_network_{gene_symbol}.png")
    cmd_img = [
        "uv", "run", string_cli, "image",
        "--identifiers"
    ] + mapped_ids + [
        "--species", str(species),
        "--format", "highres_image",
        "--output", out_img_file
    ]
    
    img_bytes = None
    try:
        subprocess.run(cmd_img, capture_output=True, text=True, check=True, shell=True)
        if os.path.exists(out_img_file):
            with open(out_img_file, "rb") as f:
                img_bytes = f.read()
            try: os.remove(out_img_file)
            except OSError: pass
    except Exception:
        pass
        
    out_partners_file = os.path.join(temp_dir, "string_partners.tsv")
    cmd_part = [
        "uv", "run", string_cli, "partners",
        "--identifiers"
    ] + mapped_ids + [
        "--species", str(species),
        "--limit", "15",
        "--output", out_partners_file
    ]
    
    df_partners = pd.DataFrame()
    try:
        subprocess.run(cmd_part, capture_output=True, text=True, check=True, shell=True)
        if os.path.exists(out_partners_file):
            df_partners = pd.read_csv(out_partners_file, sep='\t')
            try: os.remove(out_partners_file)
            except OSError: pass
    except Exception:
        pass
        
    if img_bytes is None and df_partners.empty:
        return None, None, "Failed to retrieve STRING network or partner data."
        
    return img_bytes, df_partners, None

def run_discovery_mode(gene_name, selected_tss, genome, promoter_up=2000, promoter_down=500, transcript_chrom=None, transcript_strand=None):
    st.subheader("Discovery Mode: Find TFs")
    
    if not gene_name:
        st.info("Please enter a Gene Symbol in the sidebar to start.")
        return

    # 1. Define Region
    up_bp = promoter_up
    down_bp = promoter_down
    
    if st.button("Find TFs Binding to Promoter"):
        with st.spinner("Searching for TFs binding to {} promoter...".format(gene_name)):
            if transcript_chrom and transcript_strand:
                chrom = transcript_chrom
                strand = transcript_strand
            else:
                g_coords = genome_client.get_gene_coordinates(gene_name, genome)
                if not g_coords:
                    st.error(f"Could not find coordinates for {gene_name}")
                    return
                chrom, start, end, strand = g_coords
            
            # Determine TSS
            if selected_tss is not None:
                tss = selected_tss
            else:
                tss = start if strand == "+" else end
                
            # Calculate Region
            if strand == "+":
                reg_start = tss - up_bp
                reg_end = tss + down_bp
            else:
                reg_start = tss - down_bp
                reg_end = tss + up_bp
                
            # Search
            from tf_explorer import encode_client
            search_genome = genome 
            st.write(f"DEBUG: Searching region {chrom}:{reg_start}-{reg_end} ({search_genome})")
            results = encode_client.search_region(chrom, reg_start, reg_end, search_genome)
            st.write(f"DEBUG: Results found: {len(results) if results else 'None'}")
            
            if results is None:
                st.error("Search timed out or failed. The ENCODE API might be busy. Please try again later.")
            elif not results:
                st.warning("No TFs found in this region.")
            else:
                st.session_state.discovery_results = results
                st.success(f"Found {len(results)} TFs!")
                
    if 'discovery_results' in st.session_state:
        results = st.session_state.discovery_results
        df = pd.DataFrame(results)
        df['Select'] = False
        
        st.markdown("### Found Transcription Factors")
        st.markdown("Select TFs to perform detailed analysis (Targeted Search).")
        
        edited_df = st.data_editor(
            df,
            column_config={
                "Select": st.column_config.CheckboxColumn("Select", default=False),
                "tf_name": "Transcription Factor",
                "file_count": "Total Files",
                "experiment_count": "Experiments",
                "biosamples": "Cell Lines"
            },
            disabled=["tf_name", "file_count", "experiment_count", "biosamples"],
            hide_index=True,
            use_container_width=True
        )
        
        selected_tfs = edited_df[edited_df['Select']]['tf_name'].tolist()
        
        if selected_tfs:
            st.write(f"Selected: {', '.join(selected_tfs)}")
            
            def on_analyze(tfs):
                st.session_state['tf_list_input'] = ", ".join(tfs)
                st.session_state.analysis_mode = "Targeted Search"
                st.session_state.search_performed = False

            st.button("Analyze Selected TFs", on_click=on_analyze, args=(selected_tfs,))



def main():
    st.title("🧬 TF-Explorer v3.0 – Systems-level Epigenetic Suite")
    st.markdown("""
    Analyze transcription factor binding sites for **ANY GENE**.
    1. **Search** for ENCODE ChIP-seq experiments.
    2. **Select** the experiments you want to analyze.
    3. **Run** the analysis (with optional JASPAR motif prediction).
    """)

    # --- Sidebar ---
    st.sidebar.header("Configuration")
    
    with st.sidebar.expander("ℹ️ How to Use", expanded=False):
        st.markdown("""
        **1. Search**
        - Enter a **Gene Symbol** (e.g., `PWWP2A`).
        - Enter **Transcription Factors** (e.g., `YY1, CREB1`).
        - Click **Search ENCODE**.
        
        **2. Select Data**
        - A table of experiments will appear.
        - **Check the boxes** for the experiments you want to analyze.
        - You can select multiple cell lines to compare them.
        
        **3. Run Analysis**
        - Click **Run Analysis**.
        - The tool will download data and check for binding sites in the promoter.
        
        **4. Explore Results**
        - **Promoter Track:** See where peaks are located relative to the TSS.
        - **Cell Line Comparison:** Compare binding rates and overlaps across tissues.
        - **TF Comparison:** Compare different TFs on the same gene.
        """)

    # Gene Input
    gene_name = st.sidebar.text_input("Gene Symbol", value="", help="Enter the official gene symbol (e.g., TP53, MYC).")

    if "analysis_mode" not in st.session_state:
        st.session_state.analysis_mode = "Targeted Search"
        
    mode = st.sidebar.radio(
        "Analysis Mode", 
        ["Targeted Search", "Discovery Mode"], 
        key="analysis_mode",
        help="Choose 'Targeted Search' if you know the TFs. Choose 'Discovery Mode' to find TFs binding to this gene."
    )

    # Transcript Selection
    if "transcripts" not in st.session_state:
        st.session_state.transcripts = {}
    
    selected_tss = None
    selected_transcript_id = None
    selected_chrom = None
    selected_strand = None
    
    if gene_name:
        # Fetch if not already in session or if gene changed
        # We use a composite key or just check if the current list matches the gene
        if f"{gene_name}_transcripts" not in st.session_state.transcripts:
            with st.spinner(f"Fetching transcripts for {gene_name}..."):
                # genome_client is imported globally
                ts = genome_client.get_gene_transcripts(gene_name)
                st.session_state.transcripts[f"{gene_name}_transcripts"] = ts
        
        ts_list = st.session_state.transcripts.get(f"{gene_name}_transcripts", [])
        
        if ts_list:
            # Format options for selectbox
            # "ID (Type) - Length: X bp, TSS: Y"
            # Add a "Default (Canonical/Gene Start)" option? 
            # Actually, the first one in our sorted list IS usually canonical.
            # Let's just list them.
            
            def format_func(t):
                canon = " [Canonical]" if t['is_canonical'] else ""
                return f"{t['id']} ({t['biotype']}){canon} - Len: {t['length']}bp, TSS: {t['tss']}"
            
            selected_transcript = st.sidebar.selectbox(
                "Select Transcript (TSS)", 
                ts_list, 
                format_func=format_func,
                help="Select the specific transcript to define the TSS. The Canonical transcript is usually preferred."
            )
            
            if selected_transcript:
                selected_tss = selected_transcript['tss']
                selected_transcript_id = selected_transcript['id']
                selected_chrom = selected_transcript.get('chrom')
                selected_strand = selected_transcript.get('strand')
                st.sidebar.info(f"Using TSS: {selected_tss} ({selected_transcript['strand']})")
        else:
            st.sidebar.warning("No transcripts found. Using default gene coordinates.")

    # ENCODE Assay Type
    assay_type = st.sidebar.selectbox(
        "ENCODE Assay Type",
        ["TF ChIP-seq", "Histone ChIP-seq", "DNase-seq"],
        help="Select the assay type to search for in ENCODE. Histone ChIP-seq allows scanning for histone marks, and DNase-seq scans for chromatin accessibility."
    )
    
    # Dynamic targets label & help
    if assay_type == "TF ChIP-seq":
        target_label = "Transcription Factors"
        target_help = "Comma-separated list of TFs (e.g. YY1, CREB1)."
        target_placeholder = "e.g. YY1, CREB1"
    elif assay_type == "Histone ChIP-seq":
        target_label = "Histone Marks"
        target_help = "Comma-separated list of Histone Marks (e.g. H3K4me3, H3K27ac)."
        target_placeholder = "e.g. H3K4me3, H3K27ac"
    else: # DNase-seq
        target_label = "Accessibility Label"
        target_help = "Label to use for the chromatin accessibility track (e.g. DNase)."
        target_placeholder = "e.g. DNase"

    # Dynamic target input
    if assay_type == "DNase-seq":
        tf_input = st.sidebar.text_input(target_label, value="DNase", help=target_help)
    elif assay_type == "Histone ChIP-seq":
        common_histones = ["H3K4me3", "H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3", "H3K36me3", "H3K9ac", "H3K4me2", "H4K20me1", "H2AFZ"]
        selected_histones = st.sidebar.multiselect(
            "Select Histone Marks",
            options=common_histones + ["Other (Type custom below)"],
            default=["H3K4me3", "H3K27ac"],
            help="Select one or more histone marks to analyze."
        )
        
        custom_histone_input = ""
        if "Other (Type custom below)" in selected_histones:
            custom_histone_input = st.sidebar.text_input("Custom Histone Marks", placeholder="e.g. H3K79me2")
            
        combined_histones = [h for h in selected_histones if h != "Other (Type custom below)"]
        if custom_histone_input:
            combined_histones.extend([ch.strip() for ch in custom_histone_input.split(",") if ch.strip()])
            
        tf_input = ",".join(combined_histones)
    else: # TF ChIP-seq
        common_tfs = ["CTCF", "SP1", "E2F1", "YY1", "TP53", "MYC", "CREB1"]
        selected_tfs = st.sidebar.multiselect(
            "Select Transcription Factors",
            options=common_tfs + ["Other (Type custom below)"],
            default=["YY1", "E2F1"],
            help="Select one or more TFs to search."
        )
        
        custom_tf_input = ""
        if "Other (Type custom below)" in selected_tfs:
            custom_tf_input = st.sidebar.text_area("Custom TFs (Comma-separated)", placeholder="e.g. NANOG, OCT4")
            
        combined_tfs = [t for t in selected_tfs if t != "Other (Type custom below)"]
        if custom_tf_input:
            combined_tfs.extend([ct.strip() for ct in custom_tf_input.split(",") if ct.strip()])
            
        tf_input = ",".join(combined_tfs)
    
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

    if mode == "Discovery Mode":
        st.sidebar.info("In Discovery Mode, we will search for any TFs binding to the promoter region.")
        run_discovery_mode(gene_name, selected_tss, genome, promoter_up, promoter_down, selected_chrom, selected_strand)
        return
    
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
            st.error(f"Please enter at least one {target_label}.")
        else:
            tf_list = [tf.strip() for tf in tf_input.split(",") if tf.strip()]
            st.session_state.encode_results = []
            
            with st.spinner("Searching ENCODE..."):
                from tf_explorer import encode_client
                all_results = []
                for tf in tf_list:
                    # Map genome to organism (simple mapping for now)
                    organism = "Homo sapiens" 
                    results = encode_client.search_encode(tf, assay_type, organism)
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
                                tf_list=[tf.strip() for tf in tf_input.split(",") if tf.strip()],
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
                                progress_callback=progress_callback,
                                explicit_tss=selected_tss,
                                transcript_id=selected_transcript_id
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
                            'total_checked': len(selected_experiments),
                            'assay_type': assay_type
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
                
                col_rep1, col_rep2 = st.columns(2)
                with col_rep1:
                    st.download_button("Download Report (TXT) 📄", report_text, file_name=f"{gene_name}_analysis_report.txt", use_container_width=True)
                with col_rep2:
                    try:
                        from tf_explorer import report_generator
                        tf_list = [tf.strip() for tf in tf_input.split(",") if tf.strip()]
                        
                        df_encode_plot = pd.DataFrame()
                        df_motifs_plot = pd.DataFrame()
                        df_cons = None
                        df_synergy = None
                        df_gtex = None
                        string_res = None
                        primer_results = None
                        
                        hits_path = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                        motifs_path = os.path.join(out_dir, f"{gene_name}_motif_predictions.csv")
                        cons_path = os.path.join(out_dir, f"{gene_name}_conservation.csv")
                        synergy_path = os.path.join(out_dir, f"{gene_name}_synergy_hotspots.csv")
                        
                        if os.path.exists(hits_path):
                            df_encode_plot = pd.read_csv(hits_path)
                        if os.path.exists(motifs_path):
                            df_motifs_plot = pd.read_csv(motifs_path)
                        if os.path.exists(cons_path):
                            df_cons = pd.read_csv(cons_path)
                        if os.path.exists(synergy_path):
                            df_synergy = pd.read_csv(synergy_path)
                            
                        gtex_cache_key = f"gtex_{gene_name}_{'_'.join(sorted(tf_list))}"
                        if gtex_cache_key in st.session_state and st.session_state[gtex_cache_key]["error"] is None:
                            df_gtex = st.session_state[gtex_cache_key]["df"]
                            
                        string_cache_key = f"string_{gene_name}_{'_'.join(sorted(tf_list))}"
                        if string_cache_key in st.session_state and st.session_state[string_cache_key]["error"] is None:
                            string_res = st.session_state[string_cache_key]
                            
                        if 'chip_primer_results' in st.session_state:
                            primer_results = st.session_state['chip_primer_results']
                            
                        doc_stream = report_generator.generate_report_docx(
                            gene_name=gene_name,
                            tf_list=tf_list,
                            promoter_up=promoter_up,
                            promoter_down=promoter_down,
                            genome=genome,
                            df_encode=df_encode_plot,
                            df_motifs=df_motifs_plot,
                            df_cons=df_cons,
                            df_synergy=df_synergy,
                            df_gtex=df_gtex,
                            string_res=string_res,
                            primer_results=primer_results
                        )
                        
                        st.download_button(
                            label="Download Full Word Report (.docx) 📝",
                            data=doc_stream,
                            file_name=f"{gene_name}_complete_research_report.docx",
                            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                            use_container_width=True
                        )
                    except Exception as e:
                        st.warning(f"Word Exporter unavailable: {e}")

                # Plots
                st.subheader("Visualizations")
                
                assay_label = results.get('assay_type', 'TF ChIP-seq')
                if assay_label == "TF ChIP-seq":
                    comp_tab_name = "TF Comparison"
                    count_tab_name = "TF Binding Counts"
                elif assay_label == "Histone ChIP-seq":
                    comp_tab_name = "Histone Mark Comparison"
                    count_tab_name = "Histone Mark Signals"
                else:
                    comp_tab_name = "Accessibility Comparison"
                    count_tab_name = "Accessibility Signals"
                
                # Use radio button for navigation to persist state and improve performance
                view_options = ["Promoter Track & Synergy", "Evolutionary Conservation", "GTEx Tissue Profiles", "STRING Epigenetic Interactome", "Functional Prediction (TCGA)", "Cell Line Comparison", comp_tab_name, "Biosample Distribution", count_tab_name, "ChIP Primer Design"]
                view_mode = st.radio("View Results", view_options, horizontal=True, label_visibility="collapsed")
                
                st.divider()
                
                if view_mode == "Promoter Track & Synergy":
                    # Interactive Promoter Track & Synergy
                    hits_path = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    motifs_path = os.path.join(out_dir, f"{gene_name}_motif_predictions.csv")
                    cons_path = os.path.join(out_dir, f"{gene_name}_conservation.csv")
                    synergy_path = os.path.join(out_dir, f"{gene_name}_synergy_hotspots.csv")
                    
                    df_encode_plot = pd.DataFrame()
                    df_motifs_plot = pd.DataFrame()
                    df_cons = None
                    df_synergy = None
                    
                    if os.path.exists(hits_path):
                        df_encode_plot = pd.read_csv(hits_path)
                    if os.path.exists(motifs_path):
                        df_motifs_plot = pd.read_csv(motifs_path)
                    if os.path.exists(cons_path):
                        df_cons = pd.read_csv(cons_path)
                    if os.path.exists(synergy_path):
                        df_synergy = pd.read_csv(synergy_path)
                        
                    # UI Controls for High-Peak Track
                    st.markdown("### Promoter Track & Synergy View")
                    st.markdown("Overlay high-confidence evolutionary-conserved peaks and synergy hotspots on the promoter track.")
                    
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
                            top_n=top_n,
                            df_cons=df_cons,
                            df_synergy=df_synergy
                        )
                        st.pyplot(fig)
                        
                        st.markdown("---")
                        st.subheader("🌐 Premium Interactive Track Browser")
                        try:
                            from tf_explorer import browser_svg
                            svg_html = browser_svg.generate_interactive_genome_svg(
                                df_encode_plot, df_motifs_plot, df_cons, df_synergy,
                                promoter_up, promoter_down, gene_name,
                                top_n=top_n, threshold=threshold
                            )
                            st.components.v1.html(svg_html, height=450, scrolling=False)
                            st.download_button(
                                label="Download Standalone Interactive Browser Map (HTML) 🌐",
                                data=svg_html,
                                file_name=f"{gene_name}_interactive_tracks.html",
                                mime="text/html",
                                use_container_width=True
                            )
                        except Exception as e:
                            st.error(f"Error rendering interactive browser: {e}")

                        
                        # Display Synergy hotspots table below the plot if found
                        if df_synergy is not None and not df_synergy.empty:
                            st.markdown("#### 🌟 Detected TF Motif Synergy Hotspots")
                            st.markdown("These hotspots show regions where multiple different transcription factors have binding motifs clustered within 100 bp.")
                            st.dataframe(df_synergy, use_container_width=True)
                        else:
                            st.info("No motif synergy hotspots detected (distance within 100 bp between distinct TFs) in this promoter window.")

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
                                if 'max_conservation' in high_conf_table.columns: cols.append('max_conservation')
                                
                                st.dataframe(high_conf_table[cols].sort_values('distance_to_tss'))
                            else:
                                st.info("No peaks meet the selected criteria.")
                                
                    else:
                        st.warning("No data to plot.")

                if view_mode == "Evolutionary Conservation":
                    st.markdown("### Evolutionary Conservation Profiling")
                    st.markdown("Analyze baseline evolutionary constraints across the promoter using **UCSC phastCons** and **phyloP** scores.")
                    
                    cons_path = os.path.join(out_dir, f"{gene_name}_conservation.csv")
                    hits_path = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    synergy_path = os.path.join(out_dir, f"{gene_name}_synergy_hotspots.csv")
                    
                    df_cons = pd.DataFrame()
                    df_encode = pd.DataFrame()
                    df_synergy = pd.DataFrame()
                    
                    if os.path.exists(cons_path):
                        df_cons = pd.read_csv(cons_path)
                    if os.path.exists(hits_path):
                        df_encode = pd.read_csv(hits_path)
                    if os.path.exists(synergy_path):
                        df_synergy = pd.read_csv(synergy_path)
                        
                    if not df_cons.empty:
                        # Compute base statistics
                        avg_phast = df_cons["phastCons"].mean()
                        max_phast = df_cons["phastCons"].max()
                        avg_phylop = df_cons["phyloP"].mean()
                        max_phylop = df_cons["phyloP"].max()
                        
                        # Style metric cards using columns
                        mc1, mc2, mc3, mc4 = st.columns(4)
                        with mc1:
                            st.metric("Avg phastCons", f"{avg_phast:.3f}")
                        with mc2:
                            st.metric("Max phastCons", f"{max_phast:.3f}")
                        with mc3:
                            st.metric("Avg phyloP", f"{avg_phylop:.3f}")
                        with mc4:
                            st.metric("Max phyloP", f"{max_phylop:.3f}")
                            
                        st.divider()
                        
                        # 1. High-Confidence Conserved Binding Sites
                        st.markdown("#### 🎯 High-Confidence Conserved Binding Sites")
                        st.markdown("Peaks located in genomic regions where maximum `phastCons` score exceeds **0.8** (evolutionarily highly constrained).")
                        
                        if not df_encode.empty:
                            # Filter high-confidence sites
                            if "high_confidence_site" in df_encode.columns:
                                df_high_conf = df_encode[df_encode["high_confidence_site"] == True]
                            else:
                                df_high_conf = df_encode[df_encode.get("max_conservation", 0.0) > 0.8]
                                
                            if not df_high_conf.empty:
                                cols = ['tf', 'experiment', 'biosample', 'distance_to_tss', 'max_conservation']
                                if 'signal' in df_high_conf.columns: cols.append('signal')
                                st.dataframe(df_high_conf[cols].sort_values('max_conservation', ascending=False), use_container_width=True)
                            else:
                                st.info("No ChIP-seq peaks are located in high-conservation blocks (max phastCons > 0.8) for this promoter range.")
                        else:
                            st.info("No ChIP-seq peak data found to compare with conservation.")
                            
                        st.divider()
                        
                        # 2. Base-by-Base conservation score table
                        st.markdown("#### 📊 Base-by-Base Promoter Conservation Scores")
                        st.markdown("Every nucleotide coordinate in the promoter window with its corresponding phyloP and phastCons scores.")
                        
                        st.dataframe(df_cons.sort_values("distance_to_tss"), use_container_width=True)
                        
                        # Download button
                        csv_cons = df_cons.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label="Download Base-by-Base Conservation CSV",
                            data=csv_cons,
                            file_name=f"{gene_name}_promoter_conservation.csv",
                            mime="text/csv"
                        )
                    else:
                        st.warning("No UCSC conservation data found for this analysis run.")

                if view_mode == "GTEx Tissue Profiles":
                    st.markdown("### GTEx Tissue Expression Profiling")
                    st.markdown("Examine the baseline mRNA expression profile of the target gene and selected transcription factors across 54 normal human tissues.")
                    
                    tf_list = [tf.strip() for tf in tf_input.split(",") if tf.strip()]
                    
                    # Create cache key
                    cache_key = f"gtex_{gene_name}_{'_'.join(sorted(tf_list))}"
                    
                    if cache_key not in st.session_state:
                        with st.spinner("Retrieving dynamic median expression data from GTEx database..."):
                            df_gtex, err = fetch_gtex_data(gene_name, tf_list)
                            if err:
                                st.session_state[cache_key] = {"df": None, "error": err}
                            else:
                                st.session_state[cache_key] = {"df": df_gtex, "error": None}
                                
                    gtex_res = st.session_state[cache_key]
                    
                    if gtex_res["error"]:
                        st.error(gtex_res["error"])
                    elif gtex_res["df"] is not None and not gtex_res["df"].empty:
                        df_gtex = gtex_res["df"]
                        
                        # Generate and display the plot
                        fig_gtex = plot_gtex_expression(df_gtex, gene_name)
                        st.pyplot(fig_gtex)
                        
                        st.divider()
                        
                        # Detailed tissue expression table
                        st.markdown("#### 📋 GTEx Tissue Expression Data Table")
                        st.markdown("Baseline median Transcripts Per Million (TPM) values per tissue site.")
                        
                        # Pivot table for better readability
                        try:
                            df_pivot = df_gtex.pivot(index='tissue', columns='Gene', values='TPM').reset_index()
                            # Clean tissue names
                            df_pivot['tissue'] = df_pivot['tissue'].apply(lambda x: x.replace('_', ' '))
                            st.dataframe(df_pivot, use_container_width=True)
                            
                            csv_gtex = df_pivot.to_csv(index=False).encode('utf-8')
                            st.download_button(
                                label="Download GTEx Expression CSV",
                                data=csv_gtex,
                                file_name=f"{gene_name}_gtex_expression.csv",
                                mime="text/csv"
                            )
                        except Exception as e:
                            st.dataframe(df_gtex, use_container_width=True)
                    else:
                        st.warning("No baseline expression data returned from GTEx for the queried genes.")

                if view_mode == "STRING Epigenetic Interactome":
                    st.markdown("### STRING Epigenetic Interactome Network")
                    st.markdown("Visualize the macromolecular and protein-protein physical/functional interactions of the epigenetic complex between the target gene and regulatory transcription factors.")
                    
                    tf_list = [tf.strip() for tf in tf_input.split(",") if tf.strip()]
                    
                    # Create cache key
                    cache_key = f"string_{gene_name}_{'_'.join(sorted(tf_list))}"
                    
                    if cache_key not in st.session_state:
                        with st.spinner("Fetching protein interaction networks and partners from STRING database..."):
                            img_bytes, df_partners, err = fetch_string_data(gene_name, tf_list)
                            if err:
                                st.session_state[cache_key] = {"img": None, "partners": None, "error": err}
                            else:
                                st.session_state[cache_key] = {"img": img_bytes, "partners": df_partners, "error": None}
                                
                    string_res = st.session_state[cache_key]
                    
                    if string_res["error"]:
                        st.error(string_res["error"])
                    else:
                        # Display Network PNG
                        if string_res["img"] is not None:
                            st.markdown("#### 🕸️ Interaction Network Graph")
                            st.image(string_res["img"], caption=f"STRING Network for {gene_name} & Associated TFs", use_container_width=True)
                        else:
                            st.info("No network visualization available.")
                            
                        st.divider()
                        
                        # Display partners table
                        if string_res["partners"] is not None and not string_res["partners"].empty:
                            st.markdown("#### 🤝 Physical & Functional Interaction Partners")
                            st.markdown("Top predicted interaction partners, including their combined confidence scores.")
                            
                            df_part = string_res["partners"].copy()
                            # Rename columns for premium look
                            rename_cols = {
                                "stringId_A": "STRING ID A",
                                "stringId_B": "STRING ID B",
                                "preferredName_A": "Node A",
                                "preferredName_B": "Node B",
                                "score": "Combined Score",
                                "nscore": "Gene Neighborhood",
                                "ascore": "Gene Co-occurrence",
                                "escore": "Experimental Evidence",
                                "dscore": "Database Evidence",
                                "tscore": "Textmining Score"
                            }
                            # Filter and rename only columns that exist
                            part_cols = [c for c in rename_cols.keys() if c in df_part.columns]
                            df_display = df_part[part_cols].rename(columns={c: rename_cols[c] for c in part_cols})
                            
                            st.dataframe(df_display.sort_values(by="Combined Score", ascending=False) if "Combined Score" in df_display.columns else df_display, use_container_width=True)
                        else:
                            st.info("No detailed interactome partners table returned from STRING.")


                if view_mode == "Cell Line Comparison":
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


                if view_mode == comp_tab_name:
                    feature_type_label = "Transcription Factor" if assay_label == "TF ChIP-seq" else ("Histone Mark" if assay_label == "Histone ChIP-seq" else "Accessibility Mark")
                    st.markdown(f"### {feature_type_label} Comparison")
                    st.markdown(f"Compare binding/accessibility patterns of different {feature_type_label}s on this gene.")
                    
                    hits_path_tf = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    if os.path.exists(hits_path_tf):
                        hits_df_tf = pd.read_csv(hits_path_tf)
                        if not hits_df_tf.empty and 'tf' in hits_df_tf.columns:
                            available_tfs = sorted([x for x in hits_df_tf['tf'].unique()])
                            
                            if len(available_tfs) >= 2:
                                selected_tfs = st.multiselect(f"Select {feature_type_label}s to Compare", available_tfs, default=available_tfs[:2])
                                
                                if st.button(f"Compare Selected {feature_type_label}s"):
                                    if len(selected_tfs) < 2:
                                        st.error(f"Please select at least two {feature_type_label}s.")
                                    else:
                                        # Run Comparative Analysis with group_by='tf'
                                        comp_analysis_tf = comparative.ComparativeAnalysis(hits_df_tf, selected_tfs, gene_name, group_by='tf')
                                        
                                        metrics_tf = comp_analysis_tf.calculate_metrics()
                                        interpretation_tf = comp_analysis_tf.generate_interpretation(metrics_tf)
                                        
                                        st.markdown("##### Results")
                                        
                                        # 1. Multi-Track Plot
                                        st.markdown(f"**Comparative {feature_type_label} Tracks**")
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
                                            
                                        with st.expander(f"Unique {feature_type_label} Regions (bp)", expanded=True):
                                            st.caption(f"Number of base pairs covered uniquely by each {feature_type_label} (not shared with others).")
                                            unique_tf_df = pd.DataFrame(list(metrics_tf['unique_bases_counts'].items()), columns=[feature_type_label, 'Unique BP'])
                                            st.bar_chart(unique_tf_df.set_index(feature_type_label))
                            else:
                                st.warning(f"Need at least 2 different {feature_type_label}s in the results to perform comparison. Found: {available_tfs}")
                        else:
                            st.warning(f"No {feature_type_label} data available.")
                    else:
                        st.warning("No analysis results found.")

                if view_mode == "Biosample Distribution":
                    hits_path = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    if os.path.exists(hits_path):
                        hits_df = pd.read_csv(hits_path)
                        if not hits_df.empty and 'biosample' in hits_df.columns:
                            bio_counts = hits_df['biosample'].value_counts()
                            st.bar_chart(bio_counts)
                            st.caption("Number of binding events per biosample (cell line/tissue).")
                        else:
                            st.info("No biosample data available.")

                if view_mode == count_tab_name:
                    plot_path = os.path.join(out_dir, f"{gene_name}_tf_binding_plot.png")
                    if os.path.exists(plot_path):
                        st.image(plot_path, caption=f"Signal/Binding Summary for {gene_name}", use_container_width=True)
                    else:
                        st.warning("Plot not generated.")

                if view_mode == "ChIP Primer Design":
                    st.markdown("### ChIP Primer Design")
                    st.markdown("Design primers for specific binding sites identified in this analysis.")
                    
                    # 1. Select Target Peak
                    st.markdown("#### 1. Select Target Peak")
                    
                    # Load hits to let user select
                    hits_path_pd = os.path.join(out_dir, f"{gene_name}_encode_hits.csv")
                    if os.path.exists(hits_path_pd):
                        hits_df_pd = pd.read_csv(hits_path_pd)
                        
                        if not hits_df_pd.empty:
                            # Create a selection box for peaks
                            # Format: TF (Exp - Biosample) @ dist
                            hits_df_pd['label'] = hits_df_pd.apply(lambda x: f"{x['tf']} ({x['experiment']} - {x.get('biosample', 'Unknown')}) @ {x['distance_to_tss']} bp", axis=1)
                            
                            peak_options = hits_df_pd['label'].tolist()
                            selected_peak_label = st.selectbox("Select a Peak to Target", ["Custom Region"] + peak_options)
                            
                            target_start = -200
                            target_end = -100
                            
                            if selected_peak_label != "Custom Region":
                                # Get peak data
                                peak_data = hits_df_pd[hits_df_pd['label'] == selected_peak_label].iloc[0]
                                peak_center = peak_data['distance_to_tss']
                                # Default window: center +/- 50bp
                                target_start = int(peak_center - 50)
                                target_end = int(peak_center + 50)
                                st.info(f"Selected Peak Center: {peak_center} bp relative to TSS. Targeting region: {target_start} to {target_end}.")
                            
                            # 2. Design Settings
                            st.markdown("#### 2. Design Settings")
                            
                            col_pd1, col_pd2 = st.columns(2)
                            with col_pd1:
                                sub_start = st.number_input("Sub-region Start (TSS-relative)", value=target_start, step=10)
                            with col_pd2:
                                sub_end = st.number_input("Sub-region End (TSS-relative)", value=target_end, step=10)
                                
                            design_mode = st.radio("Design Strategy", ["Specific Sub-region", "Tiling (Multiple Pairs)"], horizontal=True)
                            
                            # --- Peak Distribution Plot ---
                            st.markdown("#### 3. Peak Distribution & Target Visualization")
                            if not hits_df_pd.empty:
                                import matplotlib.pyplot as plt
                                import seaborn as sns
                                
                                fig_pd, ax_pd = plt.subplots(figsize=(10, 4))
                                
                                # Filter to requested window (-2000 to +500)
                                plot_df = hits_df_pd[
                                    (hits_df_pd['distance_to_tss'] >= -2000) & 
                                    (hits_df_pd['distance_to_tss'] <= 500)
                                ]
                                
                                if not plot_df.empty:
                                    sns.histplot(
                                        data=plot_df, 
                                        x='distance_to_tss', 
                                        hue='tf', 
                                        multiple='stack', 
                                        bins=50, 
                                        ax=ax_pd,
                                        element="step"
                                    )
                                else:
                                    st.info("No peaks found in the -2000 to +500 bp region.")
                                
                                # Highlight the selected sub-region
                                ax_pd.axvspan(sub_start, sub_end, color='red', alpha=0.2, label='Target Region')
                                ax_pd.axvline(sub_start, color='red', linestyle='--', alpha=0.5)
                                ax_pd.axvline(sub_end, color='red', linestyle='--', alpha=0.5)
                                
                                ax_pd.set_title(f"TF Peak Distribution relative to TSS ({gene_name})")
                                ax_pd.set_xlabel("Distance to TSS (bp)")
                                ax_pd.set_ylabel("Peak Count")
                                ax_pd.set_xlim(-2000, 500) # Enforce the view limits
                                
                                # Set x-axis ticks to every 100 bp
                                from matplotlib.ticker import MultipleLocator
                                ax_pd.xaxis.set_major_locator(MultipleLocator(100))
                                plt.setp(ax_pd.get_xticklabels(), rotation=90, fontsize=8) # Rotate labels for readability
                                
                                # Add TSS line
                                ax_pd.axvline(0, color='black', linestyle='-', label='TSS')
                                
                                # Force legend
                                # ax_pd.legend() 
                                # Seaborn handles legend usually, but we added manual spans.
                                # Let's just ensure layout is tight
                                plt.tight_layout()
                                
                                st.pyplot(fig_pd)
                            # ------------------------------
                            
                            if st.button("Design ChIP Primers", type="primary"):
                                # Fetch sequence
                                with st.spinner("Fetching sequence and designing primers..."):
                                    # Get coords from genome client again to be safe
                                    g_coords = genome_client.get_gene_coordinates(gene_name, genome)
                                    if g_coords:
                                        chrom, start, end, strand = g_coords
                                        
                                        # Use selected TSS if available
                                        if selected_tss is not None:
                                            tss = selected_tss
                                            st.info(f"Using selected TSS: {tss}")
                                        else:
                                            tss = start if strand == "+" else end
                                        
                                        # We need enough context. Let's fetch -2000 to +500 as standard context
                                        ctx_up = 2000
                                        ctx_down = 500
                                        
                                        # Adjust if user asks for something outside
                                        if sub_start < -ctx_up: ctx_up = abs(sub_start) + 200
                                        if sub_end > ctx_down: ctx_down = sub_end + 200
                                        
                                        seq = genome_client.get_promoter_sequence(chrom, tss, strand, ctx_up, ctx_down, genome)
                                        
                                        if seq:
                                            # Convert TSS-relative sub-region to 0-based index in seq
                                            idx_start = sub_start + ctx_up
                                            idx_end = sub_end + ctx_up
                                            
                                            # Validate
                                            if idx_start < 0: idx_start = 0
                                            if idx_end > len(seq): idx_end = len(seq)
                                            
                                            if idx_start >= idx_end:
                                                st.error("Invalid sub-region (Start >= End).")
                                            else:
                                                # Design
                                                results_list = []
                                                
                                                if design_mode == "Specific Sub-region":
                                                    # Design 1 pair for this region
                                                    incl_len = idx_end - idx_start
                                                    res = primer_design.design_primers(seq, included_region=[idx_start, incl_len])
                                                    if 'PRIMER_PAIR_NUM_RETURNED' in res:
                                                        results_list.append(("Region 1", res))
                                                    else:
                                                        st.warning(f"No primers found: {res.get('error', 'Unknown error')}")
                                                        
                                                else: # Tiling
                                                    # Tile across the region
                                                    tile_size = 150
                                                    step_size = 75
                                                    
                                                    curr = idx_start
                                                    count = 1
                                                    while curr + tile_size <= idx_end:
                                                        res = primer_design.design_primers(seq, included_region=[curr, tile_size])
                                                        if 'PRIMER_PAIR_NUM_RETURNED' in res and res['PRIMER_PAIR_NUM_RETURNED'] > 0:
                                                            results_list.append((f"Tile {count}", res))
                                                        curr += step_size
                                                        count += 1
                                                
                                                # Store results in session state
                                                st.session_state['chip_primer_results'] = results_list
                                                st.session_state['chip_primer_gene'] = gene_name
                                                st.session_state['chip_primer_ctx_up'] = ctx_up
                                                
                                        else:
                                            st.error("Failed to fetch sequence.")
                                    else:
                                        st.error("Gene coordinates not found.")

                            # Render Results from Session State
                            if 'chip_primer_results' in st.session_state:
                                results_list = st.session_state['chip_primer_results']
                                saved_gene = st.session_state.get('chip_primer_gene', 'target')
                                saved_ctx_up = st.session_state.get('chip_primer_ctx_up', 2000)
                                
                                st.divider()
                                st.markdown(f"### Results for {saved_gene}")
                                
                                if results_list:
                                    st.success(f"Designed {len(results_list)} primer sets.")
                                    
                                    # Prepare data for CSV
                                    csv_data = []
                                    
                                    for name, res in results_list:
                                        with st.expander(f"{name} Results", expanded=True):
                                            # Show top pair
                                            if res.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                                                p = 0
                                                fp = res[f'PRIMER_LEFT_{p}_SEQUENCE']
                                                rp = res[f'PRIMER_RIGHT_{p}_SEQUENCE']
                                                tm_f = res[f'PRIMER_LEFT_{p}_TM']
                                                tm_r = res[f'PRIMER_RIGHT_{p}_TM']
                                                prod = res[f'PRIMER_PAIR_{p}_PRODUCT_SIZE']
                                                
                                                # Calculate Binding Locations
                                                fp_idx = res[f'PRIMER_LEFT_{p}'][0]
                                                rp_idx = res[f'PRIMER_RIGHT_{p}'][0]
                                                
                                                fp_loc = fp_idx - saved_ctx_up
                                                rp_loc = rp_idx - saved_ctx_up
                                                
                                                # Calculate safety metrics
                                                safety_res = primer_design.check_primer_safety(fp, rp)
                                                safety_status = safety_res.get("safety_status", "🟢 Safe")

                                                # Add to CSV data
                                                csv_data.append({
                                                    "Primer Name": f"FP_{fp_loc}",
                                                    "Sequence": fp,
                                                    "Tm": round(tm_f, 2),
                                                    "Amplicon Size": prod,
                                                    "Safety Status": safety_status
                                                })
                                                csv_data.append({
                                                    "Primer Name": f"RP_{rp_loc}",
                                                    "Sequence": rp,
                                                    "Tm": round(tm_r, 2),
                                                    "Amplicon Size": prod,
                                                    "Safety Status": safety_status
                                                })
                                                
                                                c1, c2, c3 = st.columns(3)
                                                with c1:
                                                    st.markdown("**Forward**")
                                                    st.code(fp)
                                                    st.caption(f"Tm: {tm_f:.1f} | Loc: **{fp_loc}**")
                                                with c2:
                                                    st.markdown("**Reverse**")
                                                    st.code(rp)
                                                    st.caption(f"Tm: {tm_r:.1f} | Loc: **{rp_loc}**")
                                                with c3:
                                                    st.markdown("**Product & Safety**")
                                                    st.write(f"Size: **{prod} bp**")
                                                    st.markdown(f"Status: **{safety_status}**")
                                                    
                                                with st.expander(f"🔬 Secondary Structure Details ({name})", expanded=False):
                                                    st.markdown(fr"**Forward Primer Hairpin:** Tm = {safety_res.get('fp_hairpin_tm', 0.0):.1f}°C, $\Delta G$ = {safety_res.get('fp_hairpin_dg', 0.0):.2f} kcal/mol")
                                                    st.markdown(fr"**Reverse Primer Hairpin:** Tm = {safety_res.get('rp_hairpin_tm', 0.0):.1f}°C, $\Delta G$ = {safety_res.get('rp_hairpin_dg', 0.0):.2f} kcal/mol")
                                                    st.markdown(fr"**Forward Homodimer:** Tm = {safety_res.get('fp_homodimer_tm', 0.0):.1f}°C, $\Delta G$ = {safety_res.get('fp_homodimer_dg', 0.0):.2f} kcal/mol")
                                                    st.markdown(fr"**Reverse Homodimer:** Tm = {safety_res.get('rp_homodimer_tm', 0.0):.1f}°C, $\Delta G$ = {safety_res.get('rp_homodimer_dg', 0.0):.2f} kcal/mol")
                                                    st.markdown(fr"**Primer Heterodimer:** Tm = {safety_res.get('heterodimer_tm', 0.0):.1f}°C, $\Delta G$ = {safety_res.get('heterodimer_dg', 0.0):.2f} kcal/mol")
                                                    
                                                    # Alignment Check Button
                                                    if st.button(f"Check Alignment ({name})", key=f"chk_{name}"):
                                                        st.info("Alignment check passed (simulated).")
                                            else:
                                                st.warning("No primers found for this region.")
                                                st.caption("Try widening the sub-region or adjusting parameters.")
                                    
                                    # Download Section
                                    if csv_data:
                                        st.markdown("#### Export Results")
                                        df_primers = pd.DataFrame(csv_data)
                                        
                                        # Show Preview
                                        st.markdown("**Preview Data:**")
                                        st.dataframe(df_primers, use_container_width=True)
                                        
                                        # Prepare CSV
                                        csv = df_primers.to_csv(index=False).encode('utf-8-sig')
                                        
                                        # Sanitize filename
                                        safe_gene = "".join([c for c in saved_gene if c.isalnum() or c in (' ','-','_')]).strip()
                                        if not safe_gene: safe_gene = "target"
                                        fname = f"{safe_gene}_chip_primers.csv"
                                        
                                        st.download_button(
                                            label="Download CSV File",
                                            data=csv,
                                            file_name=fname,
                                            mime="text/csv",
                                            key=f"download_primers_{safe_gene}_final"
                                        )
                                else:
                                    st.warning("No primers generated.")
                        else:
                            st.info("No peaks found to target.")
                    else:
                        st.info("Run analysis to see peaks.")
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
