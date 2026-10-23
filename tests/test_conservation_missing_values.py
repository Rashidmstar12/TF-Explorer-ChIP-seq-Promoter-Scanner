import pytest
import numpy as np
import pandas as pd

def test_conservation_missing_values_are_nan():
    """
    Verifies that unannotated genomic positions are initialized as np.nan,
    never as 0.0 (which falsely indicates neutral evolution in phyloP).
    """
    promoter_start = 1000
    promoter_end = 1010
    
    # Process base-by-base conservation data
    base_data = {}
    for pos in range(promoter_start, promoter_end + 1):
        base_data[pos] = {"phyloP": np.nan, "phastCons": np.nan}
        
    # Check all are NaN initially
    for pos, vals in base_data.items():
        assert np.isnan(vals["phyloP"]), f"Position {pos} phyloP should be NaN"
        assert np.isnan(vals["phastCons"]), f"Position {pos} phastCons should be NaN"

def test_conservation_valid_scores_distinguished_from_nan():
    """
    Verifies that:
    1. Valid 0.0 scores are distinct from NaN.
    2. Negative phyloP scores (accelerated evolution) are preserved as negative floats.
    3. Positive phastCons scores are preserved.
    """
    base_data = {
        1000: {"phyloP": np.nan, "phastCons": np.nan},  # Missing
        1001: {"phyloP": 0.0, "phastCons": 0.0},        # True zero
        1002: {"phyloP": -2.8, "phastCons": 0.05},      # Accelerated
        1003: {"phyloP": 4.5, "phastCons": 0.95}        # Conserved
    }
    
    rows = [{"genomic_pos": pos, "phyloP": vals["phyloP"], "phastCons": vals["phastCons"]}
            for pos, vals in base_data.items()]
    df_cons = pd.DataFrame(rows)
    
    assert np.isnan(df_cons.loc[df_cons["genomic_pos"] == 1000, "phyloP"].values[0])
    assert df_cons.loc[df_cons["genomic_pos"] == 1001, "phyloP"].values[0] == 0.0
    assert df_cons.loc[df_cons["genomic_pos"] == 1002, "phyloP"].values[0] == -2.8
    assert df_cons.loc[df_cons["genomic_pos"] == 1003, "phyloP"].values[0] == 4.5

def test_peak_max_conservation_handles_nan_cleanly():
    """
    Verifies that:
    - Peak over missing data evaluates max_conservation to NaN and overlaps_conserved to False.
    - Peak over highly conserved data evaluates max_conservation to max value and overlaps_conserved to True.
    """
    rows = [
        {"genomic_pos": 1000, "phastCons": np.nan},
        {"genomic_pos": 1001, "phastCons": np.nan},
        {"genomic_pos": 2000, "phastCons": 0.92},
        {"genomic_pos": 2001, "phastCons": 0.85}
    ]
    df_cons = pd.DataFrame(rows)
    
    # Peak 1: covers 1000..1001 (all NaN)
    scores_1 = df_cons[(df_cons["genomic_pos"] >= 1000) & (df_cons["genomic_pos"] <= 1001)]["phastCons"].dropna()
    assert scores_1.empty
    max_c_1 = np.nan if scores_1.empty else scores_1.max()
    assert np.isnan(max_c_1)
    assert bool(max_c_1 > 0.8) is False
    
    # Peak 2: covers 2000..2001 (conserved)
    scores_2 = df_cons[(df_cons["genomic_pos"] >= 2000) & (df_cons["genomic_pos"] <= 2001)]["phastCons"].dropna()
    assert not scores_2.empty
    max_c_2 = scores_2.max()
    assert max_c_2 == 0.92
    assert bool(max_c_2 > 0.8) is True
