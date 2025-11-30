import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
from typing import Dict, List, Tuple, Optional
import itertools

class ComparativeAnalysis:
    """
    Performs comparative analysis of TF binding between multiple groups (Cell Lines or TFs).
    """
    def __init__(self, df_encode: pd.DataFrame, selected_items: List[str], gene_name: str, group_by: str = 'biosample'):
        self.df = df_encode
        self.selected_items = selected_items
        self.gene_name = gene_name
        self.group_by = group_by
        
        # Filter data for each item
        self.data_map = {}
        self.peaks_map = {}
        
        for item in self.selected_items:
            # Filter by the group_by column
            if self.group_by in self.df.columns:
                item_df = self.df[self.df[self.group_by] == item].copy()
                self.data_map[item] = item_df
                # Filter for strict overlaps only for metrics
                self.peaks_map[item] = item_df[item_df['overlap'] == True]
            else:
                # Handle missing column gracefully? Or raise error?
                # For now, empty
                self.data_map[item] = pd.DataFrame()
                self.peaks_map[item] = pd.DataFrame()

    def calculate_metrics(self) -> Dict:
        """
        Calculates comparative metrics: 
        - Pairwise Jaccard indices
        - Shared by all count
        - Unique to each count
        """
        metrics = {}
        
        # 1. Pairwise Jaccard Matrix
        jaccard_matrix = pd.DataFrame(index=self.selected_items, columns=self.selected_items, dtype=float)
        
        def get_covered_bases(df):
            bases = set()
            if df.empty: return bases
            for _, row in df.iterrows():
                p_start = int(row['peak_start'])
                p_end = int(row['peak_end'])
                bases.update(range(p_start, p_end))
            return bases

        bases_map = {item: get_covered_bases(self.peaks_map[item]) for item in self.selected_items}
        
        for item1, item2 in itertools.combinations_with_replacement(self.selected_items, 2):
            b1 = bases_map[item1]
            b2 = bases_map[item2]
            
            intersection = len(b1.intersection(b2))
            union = len(b1.union(b2))
            jaccard = intersection / union if union > 0 else 0.0
            
            jaccard_matrix.loc[item1, item2] = jaccard
            jaccard_matrix.loc[item2, item1] = jaccard
            
        metrics['jaccard_matrix'] = jaccard_matrix
        
        # 2. Shared by All
        if len(self.selected_items) > 1:
            shared_bases = set.intersection(*bases_map.values())
            metrics['bases_shared_by_all'] = len(shared_bases)
        else:
            metrics['bases_shared_by_all'] = len(bases_map[self.selected_items[0]])
            
        # 3. Unique to Each (bases unique to that item)
        unique_counts = {}
        for item in self.selected_items:
            other_bases = set()
            for other_item in self.selected_items:
                if other_item != item:
                    other_bases.update(bases_map[other_item])
            
            unique_bases = bases_map[item] - other_bases
            unique_counts[item] = len(unique_bases)
            
        metrics['unique_bases_counts'] = unique_counts
        
        # 4. Basic Counts
        metrics['total_peaks'] = {item: len(self.peaks_map[item]) for item in self.selected_items}
        metrics['files_with_peaks'] = {item: self.peaks_map[item]['experiment'].nunique() if not self.peaks_map[item].empty else 0 for item in self.selected_items}
        
        return metrics

    def generate_interpretation(self, metrics: Dict) -> str:
        """Generates a biological interpretation of the results."""
        group_type = "Cell Lines" if self.group_by == 'biosample' else "Transcription Factors"
        
        text = f"Comparative analysis of **{self.gene_name}** across {len(self.selected_items)} {group_type}: **{', '.join(self.selected_items)}**.\n\n"
        
        # Check for empty results
        total_peaks = sum(metrics['total_peaks'].values())
        if total_peaks == 0:
            text += f"No binding observed in any of the selected {group_type}."
            return text
            
        # Jaccard Analysis
        matrix = metrics['jaccard_matrix']
        # Find highest similarity pair (excluding self)
        max_j = -1.0
        best_pair = None
        
        for item1, item2 in itertools.combinations(self.selected_items, 2):
            j = matrix.loc[item1, item2]
            if j > max_j:
                max_j = j
                best_pair = (item1, item2)
                
        if best_pair:
            text += f"**Highest Similarity:** {best_pair[0]} vs {best_pair[1]} (Jaccard: {max_j:.2f}). "
            if max_j > 0.5:
                text += f"These {group_type} share a significant portion of their binding sites, suggesting conserved regulation or co-binding.\n\n"
            elif max_j > 0.1:
                text += "They share some common binding sites but also have distinct regions.\n\n"
            else:
                text += "Their binding profiles are quite distinct.\n\n"
                
        # Specificity
        unique_counts = metrics['unique_bases_counts']
        most_unique_item = max(unique_counts, key=unique_counts.get)
        most_unique_count = unique_counts[most_unique_item]
        
        if most_unique_count > 0:
            text += f"**Specificity:** **{most_unique_item}** has the most unique binding footprint ({most_unique_count} bp unique), indicating potential specific regulatory functions."
            
        return text

    def plot_jaccard_heatmap(self, metrics: Dict) -> plt.Figure:
        """Plots the pairwise Jaccard similarity heatmap."""
        matrix = metrics['jaccard_matrix']
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(matrix, annot=True, cmap="YlGnBu", vmin=0, vmax=1, ax=ax, fmt=".2f")
        ax.set_title(f"Jaccard Similarity: {self.gene_name}")
        plt.tight_layout()
        return fig

    def plot_comparison(self, promoter_up: int, promoter_down: int) -> plt.Figure:
        """
        Generates a comparative track plot for N items.
        Top: Signal Density (Overlaid Lines)
        Bottom: Peak Regions (Stacked Tracks)
        """
        n_lines = len(self.selected_items)
        
        peak_panel_height = max(1, n_lines * 0.5)
        fig, (ax_sig, ax_peaks) = plt.subplots(2, 1, figsize=(12, 6 + peak_panel_height), sharex=True, 
                                             gridspec_kw={'height_ratios': [3, peak_panel_height]})
        
        # --- Signal Density Plot (Top) ---
        x_coords = np.arange(-promoter_up, promoter_down + 1)
        
        # Color palette
        colors = sns.color_palette("husl", n_lines)
        color_map = dict(zip(self.selected_items, colors))
        
        def get_signal_profile(peaks_df):
            profile = np.zeros_like(x_coords, dtype=float)
            if peaks_df.empty:
                return profile
            for _, row in peaks_df.iterrows():
                peak_center = (row['peak_start'] + row['peak_end']) // 2
                tss_inferred = peak_center - row['distance_to_tss']
                start_rel = row['peak_start'] - tss_inferred
                end_rel = row['peak_end'] - tss_inferred
                
                idx_start = int(start_rel + promoter_up)
                idx_end = int(end_rel + promoter_up)
                idx_start = max(0, idx_start)
                idx_end = min(len(profile), idx_end)
                
                if idx_start < idx_end:
                    val = row.get('signal', 0)
                    if pd.isna(val) or val == 0: val = row.get('score', 0)
                    if pd.isna(val) or val == 0: val = 1.0
                    profile[idx_start:idx_end] += val
            return profile

        for item in self.selected_items:
            profile = get_signal_profile(self.peaks_map[item])
            c = color_map[item]
            # Use line plot with slight fill
            ax_sig.plot(x_coords, profile, color=c, linewidth=1.5, label=item)
            ax_sig.fill_between(x_coords, profile, color=c, alpha=0.1)
            
        ax_sig.set_ylabel("Signal Density")
        ax_sig.legend(loc='upper right')
        ax_sig.set_title(f"Comparative Binding Profile: {self.gene_name}")
        ax_sig.grid(True, alpha=0.2)

        # --- Peak Regions Plot (Bottom) ---
        ax_peaks.axvline(x=0, color='black', linestyle='--', label='TSS')
        ax_peaks.set_xlim(-promoter_up, promoter_down)
        ax_peaks.set_xlabel("Distance to TSS (bp)")
        
        # Set y-ticks to item names
        ax_peaks.set_yticks(range(n_lines))
        ax_peaks.set_yticklabels(self.selected_items)
        ax_peaks.set_ylim(-0.5, n_lines - 0.5)
        
        for i, item in enumerate(self.selected_items):
            peaks = self.peaks_map[item]
            c = color_map[item]
            if not peaks.empty:
                for _, row in peaks.iterrows():
                    dist = row['distance_to_tss']
                    # Triangle marker
                    ax_peaks.plot(dist, i, marker='v', color=c, markersize=10, linestyle='None')
        
        ax_peaks.grid(True, axis='x', alpha=0.2)
        
        plt.tight_layout()
        return fig
