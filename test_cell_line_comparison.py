import pandas as pd
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add path to import app
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from tf_explorer.app import plot_cell_line_comparison

def test_plot_cell_line_comparison():
    print("Testing plot_cell_line_comparison...")
    
    # Mock Data
    data = {
        'experiment_id': ['EXP1', 'EXP2', 'EXP3', 'EXP4', 'EXP5'],
        'biosample': ['HEK293', 'HEK293', 'MCF-7', 'MCF-7', 'K562'],
        'file_accession': ['F1', 'F2', 'F3', 'F4', 'F5'],
        'num_overlapping_peaks_strict': [1, 0, 1, 1, 0]
    }
    df = pd.DataFrame(data)
    
    # Call function
    fig = plot_cell_line_comparison(df)
    
    # Assertions
    assert fig is not None, "Figure should not be None"
    assert isinstance(fig, plt.Figure), "Result should be a matplotlib Figure"
    
    # Check if axes are correct
    ax = fig.axes[0]
    title = ax.get_title()
    assert "TF Binding Rate by Cell Line" in title, f"Unexpected title: {title}"
    
    print("Test Passed!")

if __name__ == "__main__":
    test_plot_cell_line_comparison()
