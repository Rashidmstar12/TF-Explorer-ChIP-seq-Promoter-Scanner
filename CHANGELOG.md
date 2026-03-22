# Changelog

## [1.4.0] - 2026-03-22

### Added
- **CpG Island Detection:**
    - New `tf_explorer/enrichment.py` module implements the Gardiner-Garden &
      Frommer (1987) sliding-window algorithm to detect CpG islands in the
      promoter sequence automatically after every analysis run.
    - Results saved to `[GENE]_cpg_islands.csv` (start, end, length, GC
      fraction, CpG O/E ratio, CpG count; coordinates in both sequence-local
      and TSS-relative frames).
    - New **"CpG Islands & GC Profile"** tab in the GUI shows:
        - Sliding-window GC-content curve with configurable window size.
        - CpG island annotations drawn as green bars below the GC track.
        - Downloadable CpG islands CSV.
- **TF Co-binding Heatmap:**
    - `enrichment.calc_tf_cobinding()` builds a symmetric matrix counting the
      number of biosamples where each pair of TFs **co-binds** the promoter.
    - New **"Co-binding Heatmap"** tab in the GUI renders the matrix as an
      annotated seaborn heatmap and displays the raw matrix as a dataframe.
- **UCSC Genome Browser Quick-Link:**
    - A direct hyperlink to open the strict promoter window in UCSC is now
      shown in the results panel immediately after the window coordinates.
- **CpG island count** added to `[GENE]_combined_summary.csv`.
- **`test_enrichment.py`:** 15 focused unit tests covering GC profiling, CpG
  island detection, co-binding matrix calculation, and plot smoke-tests.

### Changed
- **CI workflow** now runs `test_enrichment.py` instead of the pre-existing
  broken `test_cell_line_comparison.py` and `test_simple.py` (network-dependent).

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
