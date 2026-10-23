import pytest
import pandas as pd
import tempfile
import os
import matplotlib
matplotlib.use('Agg')
from tf_explorer.analysis import plot_cell_line_comparison, generate_bed_output


def test_plot_cell_line_comparison_schema_compatibility():
    """
    Verifies that plot_cell_line_comparison works cleanly with both:
    1. 'num_strict_peaks' (raw schema from experiment_stats.csv)
    2. 'num_overlapping_peaks_strict' (standardized schema from app.py)
    """
    # Schema 1: Standardized
    df1 = pd.DataFrame({
        'experiment_id': ['E1', 'E2'],
        'biosample': ['K562', 'HeLa'],
        'file_accession': ['F1', 'F2'],
        'num_overlapping_peaks_strict': [2, 0]
    })
    fig1 = plot_cell_line_comparison(df1)
    assert fig1 is not None, "Must generate figure with standardized schema"
    
    # Schema 2: Raw
    df2 = pd.DataFrame({
        'experiment_id': ['E1', 'E2'],
        'biosample': ['K562', 'HeLa'],
        'file_accession': ['F1', 'F2'],
        'num_strict_peaks': [2, 0]
    })
    fig2 = plot_cell_line_comparison(df2)
    assert fig2 is not None, "Must generate figure with raw schema"

def test_generate_bed_output_writes_both_peaks_and_motifs():
    """
    Verifies that generate_bed_output formats and exports:
    1. ENCODE strict peaks (red, 255,0,0)
    2. Predicted motif occurrences (blue, 0,0,255)
    """
    df_encode = pd.DataFrame([{
        'tf': 'CTCF',
        'experiment': 'EXP01',
        'biosample': 'K562',
        'file': 'F01',
        'peak_chrom': '1',
        'peak_start': 1000,
        'peak_end': 1200,
        'distance_to_tss': 100,
        'overlap': True
    }])
    
    df_motifs = pd.DataFrame([{
        'tf_name': 'CTCF',
        'jaspar_id': 'MA0139.1',
        'strand': '+',
        'start': 50,
        'end': 69,
        'score': 12.5,
        'sequence': 'CCGCGAGG'
    }])
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        generate_bed_output(
            df_encode=df_encode,
            df_motifs=df_motifs,
            chrom="1",
            promoter_start=900,
            gene_name="TEST_GENE",
            output_dir=tmp_dir
        )
        
        bed_file = os.path.join(tmp_dir, "TEST_GENE_overlaps.bed")
        assert os.path.exists(bed_file), "BED file must be created"
        
        with open(bed_file, "r") as f:
            lines = [l.strip() for l in f if l.strip()]
            
        # Verify track header is present
        assert any(l.startswith("track") for l in lines), "BED track header must be present"
        
        feature_lines = [l for l in lines if not l.startswith("track")]
        # Must have 2 feature lines: 1 peak line and 1 motif line
        assert len(feature_lines) == 2, f"Expected 2 feature lines in BED, got {len(feature_lines)}"
        
        # Check peak line format (9 tab-separated columns)
        peak_cols = feature_lines[0].split("\t")
        assert len(peak_cols) == 9, f"Expected 9 BED columns, got {len(peak_cols)}"
        assert peak_cols[0] == "chr1"
        assert peak_cols[1] == "1000"
        assert peak_cols[2] == "1200"
        assert peak_cols[8] == "255,0,0" # Red for overlapping peak
        
        # Check motif line format (9 tab-separated columns)
        motif_cols = feature_lines[1].split("\t")
        assert len(motif_cols) == 9, f"Expected 9 BED columns, got {len(motif_cols)}"
        assert motif_cols[0] == "chr1"
        # promoter_start (900) + start (50) = 950
        assert motif_cols[1] == "950"
        # promoter_start (900) + end (69) = 969
        assert motif_cols[2] == "969"
        assert motif_cols[5] == "+"
        assert motif_cols[8] == "0,0,255" # Blue for motif

def test_generate_bed_output_negative_strand():
    """
    Verifies that for negative strand promoters, BED export accurately translates
    sequence offsets from promoter_end in reverse direction.
    """
    df_motifs = pd.DataFrame([{
        'tf_name': 'TP53',
        'jaspar_id': 'MA0106.3',
        'strand': '+',
        'start': 100,
        'end': 120,
        'score': 15.0,
        'sequence': 'RRRCWWGYYY'
    }])
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        generate_bed_output(
            df_encode=pd.DataFrame(),
            df_motifs=df_motifs,
            chrom="17",
            promoter_start=24500,
            promoter_end=27000,
            gene_name="TP53",
            output_dir=tmp_dir,
            strand="-"
        )
        
        bed_file = os.path.join(tmp_dir, "TP53_overlaps.bed")
        with open(bed_file, "r") as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith("track")]
            
        assert len(lines) == 1
        cols = lines[0].split("\t")
        # 27000 - 120 = 26880 (chromStart)
        # 27000 - 100 = 26900 (chromEnd)
        assert cols[1] == "26880"
        assert cols[2] == "26900"
        # Since transcript is '-', '+' motif maps to '-' on chromosome
        assert cols[5] == "-"

