"""
Script: plot_fig3_panel_f_rank_abundance.py
Location: /notebook/
Description: 
    Exploratory script for Figure 3.
    Generates a Rank-Abundance (Cumulative Distribution) curve to evaluate the 
    skewness of the viral libraries. A steeper curve indicates a severe jackpot 
    effect, while a flatter curve suggests a more uniform barcode distribution.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
from matplotlib.ticker import LogLocator, FuncFormatter

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
def get_cumulative_data(file_path):
    """
    Extracts UMI counts from normalized files and calculates the cumulative frequency.
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
                # Keep only valid barcodes
                if not set(bc.upper()).issubset({'A', 'T', 'C', 'G', 'N'}): continue
                umi = float(parts[4])
                if umi > 0:
                    counts.append(umi)
            except Exception:
                continue

    # Sort counts in descending order
    counts = np.sort(counts)[::-1]
    total_umi = np.sum(counts)
    
    if total_umi == 0:
        return np.array([]), np.array([])
        
    frequencies = (counts / total_umi) * 100
    cumulative_pct = np.cumsum(frequencies)
    ranks = np.arange(1, len(counts) + 1)
    
    return ranks, cumulative_pct

# ==========================================
# 2. Main Execution Flow
# ==========================================
def main():
    print("🚀 Starting Rank-Abundance Cumulative Analysis (Panel F)...")

    # Dynamic path resolution (GitHub safe)
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
    
    # Select EnvA samples, focusing on Replicate 1 for clear curve visualization
    target_samples = meta_df[(meta_df['Stage'] == 'EnvA') & (meta_df['Replicate'] == 'Rep1')]

    # ==========================================
    # 3. Visualization
    # ==========================================
    fig, ax = plt.subplots(figsize=(8, 6))

    # Store per-group 50%-crossing results for summary
    dom_clone_summary = {}

    for _, row in target_samples.iterrows():
        group_name = row['Group']
        file_path = data_dir / row['Normalized_File']
        
        print(f"  📊 Processing: {group_name}")
        ranks, cumulative = get_cumulative_data(file_path)
        
        if len(ranks) == 0: continue
            
        color = CUSTOM_COLORS.get(group_name, 'black')
        
        # Plot the cumulative curve
        ax.plot(ranks, cumulative, color=color, linewidth=3, label=group_name)

        # ---- Quantify dominant clones: rank at which 50% of reads is reached ----
        idx_50 = np.searchsorted(cumulative, 50)
        if idx_50 < len(ranks):
            n_dominant = ranks[idx_50]
            dom_clone_summary[group_name] = n_dominant
            print(f"      → {n_dominant:,} barcodes account for 50% of total reads")

            # Vertical drop-line from the curve down to the 50% threshold
            ax.plot([n_dominant, n_dominant], [0, 50], color=color,
                    linestyle=':', linewidth=1.2, alpha=0.8)
            # Anchor dot at the crossing point
            ax.scatter(n_dominant, 50, color=color, s=50, zorder=5)
            # Text label: shift right by ~40% in log₁₀ space
            txt_x = 10 ** (np.log10(n_dominant) + 0.15)
            ax.text(txt_x, 48, f'{n_dominant:,} BCs',
                    color=color, fontsize=9, fontweight='bold',
                    va='top', ha='left')
        else:
            print(f"      → 50% threshold not reached within {len(ranks):,} barcodes")

    # --- Axis Formatting ---
    ax.set_xscale('log')
    ax.set_xlim(1, 10**5.5) 
    ax.set_ylim(0, 105)
    
    # Configure Log ticks for X-axis
    locmin = LogLocator(base=10.0, numticks=10)
    ax.xaxis.set_major_locator(locmin)
    
    # Format Y-axis as percentage
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{int(y)}%'))
    
    # Add reference grids for readability
    ax.grid(True, which="major", ls="-", alpha=0.3, color='grey')
    ax.grid(True, which="minor", ls=":", alpha=0.1, color='grey')

    # Add 50% threshold reference line
    ax.axhline(50, color='gray', linestyle='--', alpha=0.7)
    ax.text(10, 52, '50% Reads Cumulative Threshold', fontsize=11, color='gray')

    # Print summary table
    if dom_clone_summary:
        print("\n  ═══════════════════════════════════════════")
        print("  Barcodes accounting for 50% of total reads:")
        for grp, n in dom_clone_summary.items():
            print(f"    {grp:<20s}  {n:>10,} BCs")
        print("  ═══════════════════════════════════════════\n")

    ax.set_xlabel('Barcode Rank (log scale)', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_ylabel('Cumulative Read Fraction', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_title('Rank-Abundance Cumulative Curve (EnvA)', fontsize=16, fontweight='bold', pad=15)
    
    ax.legend(title='Library Design', title_fontsize=12, fontsize=11, frameon=False, loc='lower right')
    sns.despine()

    plt.tight_layout()
    save_path = out_dir / "EDA_Rank_Abundance_Cumulative.svg"
    plt.savefig(save_path, format='svg', dpi=300)
    print(f"🎉 Plot saved to: {save_path}")
    plt.show()

if __name__ == "__main__":
    main()