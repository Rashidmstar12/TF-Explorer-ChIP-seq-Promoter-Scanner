import pandas as pd
import sys
import os

# Add path to import app
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from tf_explorer.app import summarize_experiments

def test_summarize_experiments():
    print("Testing summarize_experiments (File-Based)...")
    
    # Case 1: Single experiment, single file, strict hit
    data1 = {
        'experiment_id': ['EXP1'],
        'num_strict_peaks': [1],
        'num_loose_peaks': [1],
        'total_peaks_in_file': [100]
    }
    df1 = pd.DataFrame(data1)
    n_strict, n_loose, t_strict, t_loose, t_exp, _ = summarize_experiments(df1)
    assert n_strict == 1
    assert n_loose == 1
    assert t_strict == 1
    assert t_exp == 1
    print("Case 1 Passed")
    
    # Case 2: Single experiment, multiple files, mixed hits
    # File 1: No hit, File 2: Strict hit
    # Should count as 2 files, 1 strict hit
    data2 = {
        'experiment_id': ['EXP2', 'EXP2'],
        'num_strict_peaks': [0, 1],
        'num_loose_peaks': [0, 1],
        'total_peaks_in_file': [50, 50]
    }
    df2 = pd.DataFrame(data2)
    n_strict, n_loose, t_strict, t_loose, t_exp, _ = summarize_experiments(df2)
    assert n_strict == 1 
    assert n_loose == 1
    assert t_strict == 1 
    assert t_exp == 2 # 2 files
    print("Case 2 Passed")
    
    # Case 3: Multiple experiments
    # EXP3: Strict hit
    # EXP4: Loose hit only
    # EXP5: No hit
    data3 = {
        'experiment_id': ['EXP3', 'EXP4', 'EXP5'],
        'num_strict_peaks': [1, 0, 0],
        'num_loose_peaks': [1, 1, 0],
        'total_peaks_in_file': [100, 100, 100]
    }
    df3 = pd.DataFrame(data3)
    n_strict, n_loose, t_strict, t_loose, t_exp, _ = summarize_experiments(df3)
    assert n_strict == 1 
    assert n_loose == 2 
    assert t_strict == 1
    assert t_loose == 2
    assert t_exp == 3
    print("Case 3 Passed")
    
    print("All tests passed!")

if __name__ == "__main__":
    test_summarize_experiments()
