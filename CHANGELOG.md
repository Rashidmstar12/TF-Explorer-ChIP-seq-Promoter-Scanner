# Changelog

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
