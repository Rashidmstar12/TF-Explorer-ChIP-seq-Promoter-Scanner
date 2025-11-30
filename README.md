# TF-Explorer v1.2 – ChIP-seq Promoter Scanner

A general Python command-line tool and **Streamlit Web App** to analyze transcription factor (TF) binding sites for **ANY GENE**.

## Features

- **Multi-Cell Line Comparison**: Compare binding patterns across multiple cell lines with Jaccard heatmaps and multi-track plots.
- **Multi-TF Comparison**: Compare different Transcription Factors on the same gene to identify co-binding.
- **Advanced Visualizations**: Interactive Jaccard Similarity Heatmaps, Stacked Signal Tracks, and Unique Binding Site counts.
- **ENCODE Integration**: Query and download TF ChIP-seq peaks.
- **Gene Agnostic**: Works for any human gene (e.g., PWWP2A, TP53, MYC).
- **Interactive GUI**: User-friendly Streamlit interface with flexible multi-select controls.
- **High-Peak Visualization**: Overlay high-confidence peaks on the promoter track.
- **Strict & Loose Windows**: Analyze binding in both strict promoter regions and wider (±5kb) windows.
- **File-Based Counting**: Accurate reporting of total files/samples analyzed.
- **Motif Analysis**: Predict binding sites using JASPAR PWMs.

## Installation

```bash
pip install -r requirements.txt
pip install .
```

## How to Run Locally

To quickly verify the installation and see the tool in action, you can run the provided demo script:

1.  Open a terminal in the project directory.
2.  Run the `run_demo.bat` script:
    ```cmd
    run_demo.bat
    ```
    Or run the command directly:
    ```bash
    py -m tf_explorer.cli --gene PWWP2A --tf-list "E2F1,YY1" --jaspar-ids "MA0024.3,MA0095.2" --bed-output --plot-track --out demo_results
    ```

This will analyze the **PWWP2A** gene for **E2F1** and **YY1** binding sites, generating results in the `demo_results` folder.

## GUI Usage

To use the interactive web interface:

1.  Run the `run_gui.bat` script:
    ```cmd
    run_gui.bat
    ```
    Or run via command line:
    ```bash
    py -m streamlit run tf_explorer/app.py
    ```
2.  The app will open in your default web browser.

## CLI Usage

```bash
tf-explorer \
  --genome hg38 \
  --gene PWWP2A \
  --tf-list "E2F1,YY1,MYC,CREB1" \
  --jaspar-ids "MA0024.1,MA0095.1" \
  --promoter-up 2000 \
  --promoter-down 500 \
  --out results/
```

### Arguments

- `--gene`: Target gene symbol (e.g., "TP53").
- `--tf-list`: Comma-separated list of TFs to search in ENCODE (e.g., "E2F1,MYC").
- `--jaspar-ids`: Comma-separated list of JASPAR Matrix IDs for PWM scanning (e.g., "MA0139.1").
- `--genome`: Genome assembly (default: "hg38").
- `--promoter-up`: Base pairs upstream of TSS (default: 2000).
- `--promoter-down`: Base pairs downstream of TSS (default: 500).
- `--threshold`: PWM score threshold (default: 8.0).
- `--out`: Output directory.

## Output

The tool generates:
- `[GENE]_encode_hits.csv`: ENCODE peaks overlapping the promoter.
- `[GENE]_motif_predictions.csv`: Predicted binding sites from JASPAR.
- `[GENE]_combined_summary.csv`: Summary of findings.
- `[GENE]_tf_binding_plot.png`: Visualization.
