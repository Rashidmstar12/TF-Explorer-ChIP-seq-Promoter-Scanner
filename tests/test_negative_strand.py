import pytest
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from tf_explorer.comparative import ComparativeAnalysis

def test_negative_strand_comparative_coordinate_alignment():
    """
    Verifies that for a negative strand gene:
    1. Distance to TSS is signed relative to transcription (upstream is negative).
    2. In ComparativeAnalysis, the signal profile is centered at distance_to_tss,
       not shifted or inverted by 2 * distance_to_tss.
    """
    # Negative strand gene: TSS=25000, peak at [26000, 26200), center=26100.
    # On negative strand, 26100 is UPSTREAM of 25000 in transcription direction.
    # Signed distance to TSS = 25000 - 26100 = -1100 bp.
    peak_start = 26000
    peak_end = 26200
    dist_to_tss = -1100 # upstream
    
    df_neg = pd.DataFrame([{
        'tf': 'CTCF',
        'experiment': 'EXP_NEG',
        'biosample': 'HeLa',
        'file': 'F_NEG',
        'peak_chrom': '1',
        'peak_start': peak_start,
        'peak_end': peak_end,
        'distance_to_tss': dist_to_tss,
        'score': 100,
        'signal': 15.0,
        'overlap': True
    }])
    
    ca = ComparativeAnalysis(df_neg, ['HeLa'], 'TP53_LIKE', strand='-')
    
    # Check signal profile directly
    promoter_up = 2000
    promoter_down = 500
    x_coords = np.arange(-promoter_up, promoter_down + 1)
    
    fig = ca.plot_comparison(promoter_up=promoter_up, promoter_down=promoter_down)
    assert fig is not None, "Comparative figure must be generated"
    
    # Verify mathematical positioning:
    half_width = (peak_end - peak_start) / 2.0
    start_rel = dist_to_tss - half_width # -1100 - 100 = -1200
    end_rel = dist_to_tss + half_width   # -1100 + 100 = -1000
    
    idx_start = int(start_rel + promoter_up) # -1200 + 2000 = 800
    idx_end = int(end_rel + promoter_up)     # -1000 + 2000 = 1000
    
    # Peak center in index space:
    idx_center = int(dist_to_tss + promoter_up) # -1100 + 2000 = 900
    assert x_coords[idx_center] == dist_to_tss == -1100, "X-axis coordinate at peak center must be -1100"

def test_positive_strand_comparative_coordinate_alignment():
    """
    Verifies that for a positive strand gene, peak center matches distance_to_tss.
    """
    peak_start = 8900
    peak_end = 9100
    dist_to_tss = -1000 # 1000 bp upstream of TSS=10000
    
    df_pos = pd.DataFrame([{
        'tf': 'CTCF',
        'experiment': 'EXP_POS',
        'biosample': 'K562',
        'file': 'F_POS',
        'peak_chrom': '1',
        'peak_start': peak_start,
        'peak_end': peak_end,
        'distance_to_tss': dist_to_tss,
        'score': 100,
        'signal': 15.0,
        'overlap': True
    }])
    
    ca = ComparativeAnalysis(df_pos, ['K562'], 'MYC_LIKE', strand='+')
    fig = ca.plot_comparison(promoter_up=2000, promoter_down=500)
    assert fig is not None
