import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg') # Non-interactive backend
from tf_explorer import comparative
import os

def test_comparative_analysis():
    print("Testing ComparativeAnalysis...")
    
    # Mock Data
    data = {
        'biosample': ['CL1', 'CL1', 'CL2', 'CL2', 'CL3', 'CL3'],
        'tf':        ['TF_A', 'TF_B', 'TF_A', 'TF_B', 'TF_A', 'TF_B'],
        'experiment': ['E1', 'E2', 'E3', 'E4', 'E5', 'E6'],
        'peak_start': [100, 200, 150, 250, 100, 300],
        'peak_end':   [150, 250, 200, 300, 150, 350],
        'distance_to_tss': [-50, 50, 0, 100, -50, 150],
        'overlap': [True, True, True, True, True, True],
        'signal': [10, 20, 15, 25, 10, 30]
    }
    df = pd.DataFrame(data)
    
    # Test 1: Cell Line Comparison (group_by='biosample')
    print("\nTest 1: Cell Line Comparison (group_by='biosample')")
    cell_lines = ['CL1', 'CL2', 'CL3']
    comp_cl = comparative.ComparativeAnalysis(df, cell_lines, 'TEST_GENE', group_by='biosample')
    metrics_cl = comp_cl.calculate_metrics()
    print("Jaccard Matrix (Cell Lines):")
    print(metrics_cl['jaccard_matrix'])
    assert metrics_cl['jaccard_matrix'].shape == (3, 3)
    
    # Test 2: TF Comparison (group_by='tf')
    print("\nTest 2: TF Comparison (group_by='tf')")
    tfs = ['TF_A', 'TF_B']
    comp_tf = comparative.ComparativeAnalysis(df, tfs, 'TEST_GENE', group_by='tf')
    metrics_tf = comp_tf.calculate_metrics()
    print("Jaccard Matrix (TFs):")
    print(metrics_tf['jaccard_matrix'])
    
    # TF_A peaks: 
    # CL1: 100-150 (-50)
    # CL2: 150-200 (0)
    # CL3: 100-150 (-50)
    # Union bases: 100-200
    
    # TF_B peaks:
    # CL1: 200-250 (50)
    # CL2: 250-300 (100)
    # CL3: 300-350 (150)
    # Union bases: 200-350
    
    # Intersection between TF_A and TF_B: None (ranges don't overlap)
    # Jaccard should be 0
    assert metrics_tf['jaccard_matrix'].loc['TF_A', 'TF_B'] == 0.0
    
    # Test Interpretation
    text = comp_tf.generate_interpretation(metrics_tf)
    print("\nInterpretation (TF):")
    print(text)
    assert "Transcription Factors" in text
    
    # Test Plots
    fig_tracks = comp_tf.plot_comparison(2000, 500)
    fig_tracks.savefig("test_tracks_tf.png")
    
    print("\nTest Passed! Plots saved.")

if __name__ == "__main__":
    test_comparative_analysis()
