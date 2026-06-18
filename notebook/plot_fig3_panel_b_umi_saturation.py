"""
Script: plot_fig3_panel_b_umi_saturation.py
Location: /notebook/
Description: 
    Generates Figure 3 Panel B.
    Calculates UMI Rarefaction/Saturation curves using a high-performance 
    mathematical expectation algorithm (1 - (1-p)^d).
    This visualizes the untapped diversity of the Twister-optimized libraries.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import matplotlib.ticker as ticker

# ==========================================
# Global Plotting Styles (Nature Standards)
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
# 1. Math Module: Fast Rarefaction Calculation
# ==========================================
def calculate_rarefaction_curve(counts, n_points=50):
    """
    Computes rarefaction points using mathematical expectation.
    Significantly faster than stochastic downsampling for millions of reads.
    """
    total_reads = sum(counts)
    depths = np.linspace(0, total_reads, n_points).astype(int)
    unique_counts = []
    
    # Probability array
    freqs = counts / total_reads
    
    for d in depths:
        if d == 0:
            unique_counts.append(0)
            continue
        # Formula: Expected Unique = Sum(1 - (1 - frequency)^depth)
        expected_unique = np.sum(1 - (1 - freqs)**d)
        unique_counts.append(expected_unique)
        
    return depths, unique_counts

# ==========================================
# 2. Main Execution Block
# ==========================================
def main():
    print("🚀 Starting UMI Saturation Analysis (Panel B)...")
    
    # --- Dynamic Path Resolution ---
    # Since this script is in /notebook/, the project root is one level up.
    try:
        ROOT = Path(__file__).resolve().parents[1]
    except NameError:
        # Fallback for interactive Jupyter Notebook environments
        ROOT = Path.cwd().parent
        
    data_dir = ROOT / "data" / "processed" / "invitro_qc"
    out_dir = ROOT / "data" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    # File mapping for representative samples
    target_files = {
        'HH + HdVRz': 'HH + HdVRz -EnvA.target_normalized.txt',
        'HH + Twister': 'HH + Twister -EnvA.target_normalized.txt',
        'RzB + Twister': 'RzB + Twister -EnvA.target_normalized.txt'
    }

    fig, ax = plt.subplots(figsize=(7, 6))

    file_found = False
    for group_name, file_name in target_files.items():
        file_path = data_dir / file_name

        if file_path.exists():
            file_found = True
            print(f"  📊 Processing: {group_name}")
            df = pd.read_csv(file_path, sep='\t')
            
            # Robust column extraction
            if 'UMI_Count' in df.columns:
                counts = df['UMI_Count'].values
            else:
                counts = df.iloc[:, -1].values
                
            counts = counts[counts > 0].astype(int)

            # Perform calculation
            depths, species = calculate_rarefaction_curve(counts)
            
            # Actual sequencing depth and unique barcode count
            total_reads = int(depths[-1])
            total_unique = int(species[-1])
            print(f"      Total Reads: {total_reads:,}  ({total_reads/1e6:.2f} M)")
            print(f"      Unique Barcodes: {total_unique:,}")

            # Render plot
            color = CUSTOM_COLORS.get(group_name, 'black')
            ax.plot(depths, species, color=color, linewidth=3, label=group_name)
            # Add visual anchor point at the empirical maximum depth
            ax.scatter(depths[-1], species[-1], color=color, s=60, zorder=5)
            # Annotate endpoint with depth and unique count
            ax.annotate(
                f'{total_reads/1e6:.1f}M\n{total_unique:,}',
                xy=(depths[-1], species[-1]),
                xytext=(10, 10), textcoords='offset points',
                fontsize=8, color=color, fontweight='bold',
                va='bottom', ha='left'
            )
        else:
            print(f"⚠️ Warning: Data file {file_name} missing.")

    if not file_found:
        print(f"❌ Error: No valid data found in {data_dir}")
        return

    # ==========================================
    # 3. Aesthetics & Publication Formatting
    # ==========================================
    # X-axis: Format reads to Millions
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x/1e6:g}'))
    ax.set_xlabel('Sequencing Depth (Millions of Reads)', fontsize=14, fontweight='bold', labelpad=10)
    
    # Y-axis: Add comma separators for large integers
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
    ax.set_ylabel('Number of Unique Barcodes', fontsize=14, fontweight='bold', labelpad=10)
    
    ax.set_title('Library Complexity Saturation (Panel B)', fontsize=16, fontweight='bold', pad=15)
    ax.legend(frameon=False, loc='lower right', fontsize=11)
    
    sns.despine()
    plt.tight_layout()
    
    save_path = out_dir / "Fig3_Panel_B_UMI_Saturation.svg"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"🎉 Success! High-quality SVG saved to: {save_path}")
    plt.show()

if __name__ == "__main__":
    main()