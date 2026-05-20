# Changelog

## [3.0.0] - 2026-05-20

### Added
- **Unified Multi-Track Promoter View:** Visualizes base-by-base UCSC Evolutionary Conservation (`phastCons`/`phyloP`) alongside ENCODE ChIP-seq peaks, JASPAR motifs, and shaded transcription factor synergy hotspots.
- **UCSC Evolutionary Conservation Profiling:** Integrates nucleotide-level constraint metrics, dynamic average/maximum score cards, high-confidence conserved peak filtering ($phastCons > 0.8$), and interactive database grids.
- **GTEx Baseline Tissue Expression:** Displays comparative baseline tissue mRNA profiles (TPM) across 54 non-diseased human tissues using side-by-side log-scaled bar charts for target genes and low-expression TFs.
- **STRING Epigenetic Interactome:** Maps macromolecular protein-protein networks and details interaction scores (neighborhood, database, co-occurrence, textmining) in interactive search grids.
- **Thermodynamic Primer Design:** Calculates hairpins, homodimers, and heterodimers using `primer3-py` thermodynamics (`calcHairpin`, `calcHomodimer`, `calcHeterodimer`), rendering interactive color-coded safety badges ("SAFE" vs "WARNING").
- **Subprocess Session-State Caching:** Implemented full query caching under Streamlit's `st.session_state` (`gtex_{gene}` and `string_{gene}`) for sub-millisecond tab switching.

## [1.3.0] - 2026-05-15

### Added
- **Transcript Selection:** High-resolution selection of genomic transcripts for precise promoter window targeting.
- **Primer Design Plot:** Stacks tiling strategies and region constraints visually.
- **Persistent Tabs:** Stores configuration and analysis inputs persistently in the UI session.

## [1.2.0] - 2025-11-30

### Added
- **Multi-Cell Line Comparison:**
    - Users can now select **multiple cell lines** (2 or more) to compare binding patterns simultaneously.
    - Added a **Jaccard Similarity Heatmap** to visualize the overlap similarity between different cell lines.
    - Added a **Multi-Track Plot** that stacks signal density and peak markers for all selected cell lines.
    - Added a **Unique Binding Sites** bar chart to show the number of base pairs unique to each cell line.
- **Multi-TF Comparison:**
    - New **"TF Comparison" Tab** allows users to compare binding patterns of **different Transcription Factors** on the same gene.
    - Supports selecting multiple TFs (e.g., YY1 vs CREB) to see how their binding sites overlap.
    - Includes the same advanced visualizations (Heatmap, Multi-Track Plot, Unique Counts) adapted for TFs.
- **Visualizations:**
    - **Bar Charts** for "Unique Binding Sites" replaced raw JSON output for better readability.
    - **Overview Chart** ("Binding Rates by Cell Line") now correctly uses detailed experiment statistics.

### Fixed
- **Overview Plot Error:** Resolved "Insufficient data for overview" error by correctly loading `_experiment_stats.csv`.
- **KeyError:** Fixed a crash caused by a column name mismatch (`num_overlapping_peaks_strict` vs `num_strict_peaks`).
- **Indentation:** Corrected indentation issues in `app.py` that caused syntax errors.

### Changed
- **UI Improvements:**
    - Replaced pairwise selection dropdowns with `st.multiselect` for greater flexibility.
    - Improved tab organization in the results section.
