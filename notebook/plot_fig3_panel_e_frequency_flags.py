"""
Script: plot_fig3_panel_e_frequency_flags.py
Location: /notebook/
Description: 
    Generates Figure 3 Panel E.
    Plots the rank-abundance curve (UMI ratio) on a log scale and highlights 
    the total barcode richness with vertical drop-line flags.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
from matplotlib.ticker import LogLocator, PercentFormatter

# ==========================================
# Global Parameters & Formatting
# ==========================================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 1.2

CUSTOM_COLORS = {
    'HH + HdVRz': '#1F77B4',
    'HH + Twister': '#D62728',
    'RzB + Twister': '#FF7F0E'
}

# ==========================================
# 1. Calculation Module
# ==========================================
def get_rank_abundance_data(file_path):
    """
    Extracts UMI counts, sorts them, and calculates percentage frequency.
    """
    counts = []
    if not os.path.exists(file_path): 
        return np.array([]), np.array([])
        
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        first_line = f.readline().strip().split()
        if not (len(first_line) > 0 and first_line[0].lower() == 'count'):
            f.seek(0)
            
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5: continue
            try:
                bc = parts[3]
                if not set(bc.upper()).issubset({'A', 'T', 'C', 'G', 'N'}): continue
                umi = float(parts[4])
                if umi > 0:
                    counts.append(umi)
            except Exception:
                continue

    counts = np.array(counts)
    counts = np.sort(counts)[::-1] # Descending order
    total_umi = np.sum(counts)
    
    if total_umi == 0:
        return np.array([]), np.array([])
        
    umi_ratio_pct = (counts / total_umi) * 100
    ranks = np.arange(1, len(counts) + 1)
    
    return ranks, umi_ratio_pct

# ==========================================
# 2. Main Execution Flow
# ==========================================
def main():
    print("🚀 Starting Rank-Abundance Flag Analysis (Panel E)...")

    try:
        ROOT = Path(__file__).resolve().parents[1]
    except NameError:
        ROOT = Path.cwd().parent
        
    data_dir = ROOT / "data" / "processed" / "invitro_qc"
    out_dir = ROOT / "data" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    meta_path = data_dir / "Normalized_Metadata.csv"
    
    if not meta_path.exists():
        print(f"❌ Error: Metadata file not found at {meta_path}")
        return
        
    meta_df = pd.read_csv(meta_path)
    
    # Select only EnvA samples (Replicate 1 is sufficient for curve visualization)
    target_samples = meta_df[(meta_df['Stage'] == 'EnvA') & (meta_df['Replicate'] == 'Rep1')]

    # ==========================================
    # 3. Visualization
    # ==========================================
    fig, ax = plt.subplots(figsize=(8, 6))

    for _, row in target_samples.iterrows():
        group_name = row['Group']
        file_path = data_dir / row['Normalized_File']
        
        ranks, ratios = get_rank_abundance_data(file_path)
        if len(ranks) == 0: continue
            
        color = CUSTOM_COLORS.get(group_name, 'black')
        
        # Plot the main distribution line
        ax.plot(ranks, ratios, color=color, linewidth=2.5, label=group_name)

        # Render the Richness Flag (Vertical drop-line at the maximum rank)
        try:
            max_rank = ranks[-1]
            flag_y_positions = {'HH + HdVRz': 2.0, 'HH + Twister': 2.5, 'RzB + Twister': 1.5}
            flag_y = flag_y_positions.get(group_name, 2.0)
            
            # Draw vertical dashed line
            ax.plot([max_rank, max_rank], [0, flag_y], color=color, linestyle='--', linewidth=1.5, alpha=0.7)
            # Draw horizontal flag top
            ax.plot([max_rank * 0.8, max_rank * 1.2], [flag_y, flag_y], color=color, linewidth=1.5, alpha=0.7)
            # Add text annotation
            ax.text(max_rank, flag_y + 0.1, f'{max_rank:,} BCs', 
                    color=color, ha='center', va='bottom', 
                    fontsize=10, fontweight='bold', zorder=5)
                    
        except IndexError:
            continue

    # Axis Formatting
    ax.set_xscale('log')
    locmin = LogLocator(base=10.0, numticks=10)
    ax.xaxis.set_major_locator(locmin)
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=1))
    
    ax.set_xlabel('Barcode Rank (log10 scale)', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_ylabel('UMI Ratio per Barcode (% of Total)', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_title('Rank-Abundance Curve with Richness Flags', fontsize=16, fontweight='bold', pad=15)

    ax.grid(True, which="major", linestyle='-', alpha=0.3, color='grey')
    ax.grid(True, which="minor", linestyle=':', alpha=0.1, color='grey')
    
    # Extend Y-axis upper limit to prevent flag text clipping
    current_ymax = ax.get_ylim()[1]
    ax.set_ylim(-0.1, max(current_ymax, 3.2)) 
    
    ax.legend(title="Viral Library", frameon=False, fontsize=11, loc='upper right')
    sns.despine(trim=False)

    plt.tight_layout()
    save_path = out_dir / "Fig3_Panel_E_Frequency_Flags.svg"
    plt.savefig(save_path, format='svg', dpi=300)
    print(f"🎉 Plot saved to: {save_path}")
    plt.show()

if __name__ == "__main__":
    main()