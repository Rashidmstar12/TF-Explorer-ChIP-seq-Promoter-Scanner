import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline
import os

def create_wave_plot(data_df, promoter_up=2000, promoter_down=500):
    """
    Creates a smooth wave plot of ChIP-seq signal density.
    """
    # Sort by distance
    df = data_df.sort_values('distance_to_tss')
    
    x = df['distance_to_tss'].values
    y = df['signal'].values
    
    # Generate smooth curve using CubicSpline
    # We need a dense x-axis for smoothness
    x_smooth = np.linspace(x.min(), x.max(), 500)
    cs = CubicSpline(x, y)
    y_smooth = cs(x_smooth)
    
    # Ensure no negative signal if that's physically impossible (usually yes for density)
    y_smooth = np.maximum(y_smooth, 0)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Green background for promoter window
    ax.axvspan(-promoter_up, promoter_down, color='green', alpha=0.1, label='Promoter Window')
    
    # Filled area under curve
    ax.fill_between(x_smooth, y_smooth, color='blue', alpha=0.3)
    ax.plot(x_smooth, y_smooth, color='blue', linewidth=2, label='Signal Density')
    
    # TSS Marker
    ax.axvline(x=0, color='black', linestyle='--', label='TSS')
    
    # Formatting
    ax.set_xlabel(f"Distance to TSS (bp)")
    ax.set_ylabel("Signal Intensity")
    ax.set_title("ChIP-seq Signal Density")
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Set x-limits to show context but focus on window
    # User asked for -2000 to +500 axis? Or just background?
    # "X-axis: Distance to TSS (-2000 to +500 bp)" implies limits.
    ax.set_xlim(-promoter_up, promoter_down)
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    # Input data
    data = {
        'distance_to_tss': [-165, -97, -74, 75, 87, 88],
        'signal': [27.8, 59.7, 36.6, 83.6, 199.98, 202.71]
    }
    df = pd.DataFrame(data)
    
    print("Generating wave plot...")
    fig = create_wave_plot(df)
    
    output_path = "signal_wave_plot.png"
    fig.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")
