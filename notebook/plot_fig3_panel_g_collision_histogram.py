"""
Script: plot_fig3_panel_g_collision_histogram.py
Location: /notebook/
Description: 
    Exploratory script for Figure 3.
    Performs Monte Carlo sampling on UMI-weighted barcode pools to calculate 
    pairwise Levenshtein edit distances. This visualizes the physical sequence 
    space distribution and the probability of barcode collisions.
"""

import pandas as pd
import numpy as np
import itertools
from Levenshtein import distance as lev_dist
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

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
# 1. Calculation Module (Monte Carlo Sampling)
# ==========================================
def get_weighted_distances(file_path, sample_size=2000):
    """
    Extracts UMI-weighted barcodes and calculates pairwise Levenshtein distances
    via random sampling to avoid memory overload.
    """
    print(f"  ⏳ Reading and sampling: {os.path.basename(file_path)}")
    barcodes, counts = [], []
    
    if not os.path.exists(file_path):
        print(f"  ❌ Error: File not found.")
        return []

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
                    barcodes.append(bc)
                    counts.append(umi)
            except Exception:
                continue

    barcodes = np.array(barcodes)
    counts = np.array(counts)
    total_umi = np.sum(counts)
    
    if total_umi == 0 or len(barcodes) < 2:
        return []

    # Calculate sampling probabilities based on UMI abundance
    probabilities = counts / total_umi
    
    # Fix random seed for reproducible sampling
    np.random.seed(42)
    
    # Sample barcodes considering their true abundance
    sampled_bcs = np.random.choice(barcodes, size=sample_size, p=probabilities, replace=True)
    
    print(f"  🔄 Calculating edit distances for {len(sampled_bcs)} sampled molecules...")
    
    # Calculate Levenshtein distance for all unique pairs in the sample
    distances = []
    for bc1, bc2 in itertools.combinations(sampled_bcs, 2):
        distances.append(lev_dist(bc1, bc2))
        
    return distances

# ==========================================
# 2. Main Execution Flow
# ==========================================
def main():
    print("🚀 Starting Sequence Collision Analysis (Panel G)...")

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
    
    # Focus analysis on the final EnvA pseudotyped library
    target_samples = meta_df[(meta_df['Stage'] == 'EnvA') & (meta_df['Replicate'] == 'Rep1')]

    # ==========================================
    # 3. Visualization
    # ==========================================
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for _, row in target_samples.iterrows():
        group_name = row['Group']
        file_path = data_dir / row['Normalized_File']
        
        distances = get_weighted_distances(file_path, sample_size=2000)
        if not distances: continue
            
        color = CUSTOM_COLORS.get(group_name, 'black')
        
        # Calculate key statistics
        mean_dist = np.mean(distances)
        collision_rate = (distances.count(0) / len(distances)) * 100
        print(f"    -> Mean Distance: {mean_dist:.2f} | Collision Rate: {collision_rate:.4f}%")
        
        # Use discrete step histogram for accurate Levenshtein representation
        sns.histplot(distances, bins=np.arange(-0.5, max(distances)+1.5, 1),
                     element="step", fill=True, stat="density", common_norm=False,
                     color=color, alpha=0.2, linewidth=2, 
                     label=f"{group_name}", ax=ax)

    # --- Highlight near-zero collisions ---
    ax.annotate('Founder Collisions', 
                xy=(0, 0.005), 
                xytext=(2.5, 0.008), 
                arrowprops=dict(facecolor='#1F77B4', shrink=0.05, width=1.5, headwidth=8),
                fontsize=11, fontweight='bold', color='#1F77B4',
                ha='left', va='center')

    # --- Axis Formatting ---
    ax.set_xlabel('Pairwise Levenshtein Distance', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_ylabel('Probability Density', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_title('Physical Collisions in Sequence Space (EnvA)', fontsize=16, fontweight='bold', pad=15)
    
    ax.set_xlim(-1, 22) 
    ax.xaxis.set_major_locator(plt.MultipleLocator(2))
    
    # Anchor the bottom spine strictly to Y=0
    ax.spines['bottom'].set_position('zero')
    
    ax.grid(True, linestyle='--', alpha=0.3)
    sns.despine(offset=5, trim=True)
    
    ax.legend(title='Library Design', title_fontsize=12, fontsize=11, loc='upper right', frameon=False)
    
    plt.tight_layout()
    save_path = out_dir / "EDA_Collision_Histogram.svg"
    plt.savefig(save_path, format='svg', dpi=300)
    print(f"🎉 Plot saved to: {save_path}")
    plt.show()

if __name__ == "__main__":
    main()