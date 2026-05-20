# 🧬 TF-Explorer v3.0 – Systems-Level Epigenetic Suite

A comprehensive Python suite and interactive **Streamlit Web Application** designed for high-resolution transcription factor (TF) ChIP-seq peak scanning, base-by-base evolutionary constraint profiling, and systems-level macromolecular interactome correlation for **ANY HUMAN GENE**.

---

## 🚀 Key Upgrades & Features in v3.0

### 1. Unified Multi-Track Promoter Track & Synergy Visualizer
* **Evolutionary Overlay:** View base-by-base `phastCons` / `phyloP` conservation plots overlaying ChIP-seq peaks, motif predictions, and regulatory spans.
* **Motif Synergy Scan:** Automatically highlights clusters where multiple different TFs bind in close proximity (within a 100 bp window) as shaded synergy hotspots.

### 2. UCSC Evolutionary Conservation Profiling
* **Constraint Metrics:** Computes average and maximum constraint values across the targeted promoter coordinates.
* **High-Confidence Conserved Peaks:** Automatically extracts and tabulates peaks overlapping genomic regions with $phastCons > 0.8$.
* **Nucleotide-Level Grid:** Explores base-by-base conservation values with full search, sorting, and CSV export.

### 3. GTEx Baseline Cross-Tissue Expression Correlation
* **TPM Expression Profiles:** Fetches and displays target gene and TF baseline mRNA expression levels across 54 non-diseased human tissues.
* **Log-Scale Grouped Visuals:** Plots comparative bar charts side-by-side using log-scaling to compare high-expression target genes and low-expression TFs seamlessly.
* **Subprocess Caching:** Leverages custom `st.session_state` caching for high-speed sub-millisecond tab switching.

### 4. STRING Epigenetic Interactome Visualization
* **Macromolecular Networks:** Renders a high-resolution macromolecular protein-protein network mapping physical and functional partners.
* **Interactions Table:** Maps physical database, neighborhood, co-occurrence, and textmining confidence values in a searchable grid.

### 5. Advanced Thermodynamics ChIP Primer Design
* **Thermodynamic Safety Checks:** Calculates primer hairpins, homodimers, and heterodimers using `primer3-py` (`calcHairpin`, `calcHomodimer`, `calcHeterodimer`).
* **Interactive Safety Badges:** Displays interactive color-coded status badges ("SAFE" vs "WARNING") for primer pairs.

---

## 💻 Installation & Setup

Install the required packages directly using pip:

```bash
pip install -r requirements.txt
pip install .
```

### 🧬 Science Skills Prerequisite
The GTEx, STRING, and UCSC Evolutionary databases leverage deep integration with locally pre-installed **Science Skills** under your agent's config folder. Make sure your shell has access to Python 3.12+ and `primer3-py`.

---

## 🏁 How to Run

### 🖥️ Streamlit Web Interface (Recommended)
To run the full Systems-level Epigenetic Suite:

```bash
py -m streamlit run tf_explorer/app.py
```
Or use the shortcut script:
```cmd
run_gui.bat
```

### ⌨️ Command Line Interface (CLI)
For quick general peak-scanning:

```bash
tf-explorer \
  --genome hg38 \
  --gene PWWP2A \
  --tf-list "E2F1,YY1,MYC" \
  --jaspar-ids "MA0024.1,MA0095.1" \
  --promoter-up 2000 \
  --promoter-down 500 \
  --out results/
```

### 🧪 Run the Verification Pipeline
To test the systems database connectivity, Windows subprocess execution, and session-state caching:
```bash
py C:\Users\rashi\.gemini\antigravity\brain\5d40b0db-f7c1-4512-baf3-529cf5ec2848/scratch/verify_upgrades.py
```

---

## 🔬 CLI Arguments & Options

* `--gene`: Target gene symbol (e.g., `"TP53"`, `"GAPDH"`).
* `--tf-list`: Comma-separated list of TFs to search in ENCODE (e.g., `"E2F1,YY1"`).
* `--jaspar-ids`: Comma-separated list of JASPAR Matrix IDs for PWM scanning (e.g., `"MA0139.1"`).
* `--genome`: Genome assembly (default: `"hg38"`).
* `--promoter-up`: Base pairs upstream of Transcription Start Site (default: `2000`).
* `--promoter-down`: Base pairs downstream of TSS (default: `500`).
* `--threshold`: JASPAR PWM score threshold (default: `8.0`).
* `--out`: Output directory.

---

## 📂 Output File Structure

Upon running an analysis, the suite generates:
* `[GENE]_encode_hits.csv`: ENCODE ChIP-seq peaks overlapping the promoter coordinates.
* `[GENE]_motif_predictions.csv`: JASPAR predicted binding motifs inside the region.
* `[GENE]_conservation.csv`: Base-by-base phyloP/phastCons UCSC constraint values.
* `[GENE]_synergy_hotspots.csv`: Transcription factor binding sites overlapping within 100 bp.
* `[GENE]_tcga_summary.csv`: Patient correlation scores across selected cancers.
