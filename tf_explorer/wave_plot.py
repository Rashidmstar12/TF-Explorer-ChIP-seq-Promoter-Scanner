import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
import os

def create_wave_plot(data_df, promoter_up=2000, promoter_down=500):
    """
    Creates a smooth wave plot of ChIP-seq signal density.
    Robustly handles duplicate coordinates, small sample sizes, and empty inputs (DEF-016).
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axvspan(-promoter_up, promoter_down, color='green', alpha=0.1, label='Promoter Window')
    ax.axvline(x=0, color='black', linestyle='--', label='TSS')
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_ylabel("Signal Intensity")
    ax.set_title("ChIP-seq Signal Density")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-promoter_up, promoter_down)

    if data_df is None or data_df.empty or 'distance_to_tss' not in data_df.columns or 'signal' not in data_df.columns:
        ax.text(0.5, 0.5, "No ChIP-seq signal data available", transform=ax.transAxes,
                ha='center', va='center', fontsize=12, color='gray')
        ax.legend()
        plt.tight_layout()
        return fig

    # Clean missing values and sort
    df_clean = data_df.dropna(subset=['distance_to_tss', 'signal']).copy()
    if df_clean.empty:
        ax.text(0.5, 0.5, "No valid ChIP-seq signal data points", transform=ax.transAxes,
                ha='center', va='center', fontsize=12, color='gray')
        ax.legend()
        plt.tight_layout()
        return fig

    # DEF-016: Deduplicate identical TSS coordinates by averaging signal
    df_agg = df_clean.groupby('distance_to_tss', as_index=False)['signal'].mean().sort_values('distance_to_tss')
    x = df_agg['distance_to_tss'].values.astype(float)
    y = df_agg['signal'].values.astype(float)

    if len(x) == 1:
        # Single peak: draw discrete spike/scatter
        ax.scatter(x, y, color='blue', s=80, label='Peak Signal', zorder=5)
        ax.vlines(x=x[0], ymin=0, ymax=y[0], color='blue', linewidth=2)
        ax.set_ylim(bottom=0, top=max(y[0] * 1.25, 1.0))
    elif len(x) < 4:
        # 2 or 3 points: linear interpolation across smooth range
        x_min = max(x.min(), -promoter_up)
        x_max = min(x.max(), promoter_down)
        if x_min == x_max:
            x_min -= 10
            x_max += 10
        x_smooth = np.linspace(x_min, x_max, 300)
        y_smooth = np.interp(x_smooth, x, y)
        y_smooth = np.maximum(y_smooth, 0)
        ax.fill_between(x_smooth, y_smooth, color='blue', alpha=0.3)
        ax.plot(x_smooth, y_smooth, color='blue', linewidth=2, label='Signal Density (Linear)')
        ax.scatter(x, y, color='navy', s=30, zorder=5)
    else:
        # >= 4 points: use shape-preserving monotone PchipInterpolator to avoid Runge oscillations
        x_smooth = np.linspace(x.min(), x.max(), 500)
        pchip = PchipInterpolator(x, y)
        y_smooth = pchip(x_smooth)
        y_smooth = np.maximum(y_smooth, 0)
        ax.fill_between(x_smooth, y_smooth, color='blue', alpha=0.3)
        ax.plot(x_smooth, y_smooth, color='blue', linewidth=2, label='Signal Density (PCHIP)')
        ax.scatter(x, y, color='navy', s=25, alpha=0.6, zorder=5)

    ax.legend()
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
