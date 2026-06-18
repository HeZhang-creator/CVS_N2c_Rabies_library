"""
Figure Generator: Rarefaction Analysis (V1 vs LGN)
Description: Generates publication-ready rarefaction curves from pre-calculated simulation coordinates.
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import matplotlib as mpl

# ==========================================
# 1. Publication-quality plotting configuration
# ==========================================
mpl.rcParams['svg.fonttype'] = 'none' 

sns.set_theme(style="ticks", context="paper")
# Define the group color palette
CUSTOM_COLORS = {
    'HH + HdVRz': '#1F77B4',   # Muted Blue
    'HH + Twister': '#D62728', # Brick Red
    'RzB + Twister': '#FF7F0E' # Safety Orange
}

# Map replicate naming to line style. Adjust if your labels differ.
def determine_linestyle(Replicate):
    # Convert to string and lowercase for robust matching
    rep_lower = str(Replicate).lower()

    if '1' in rep_lower or 'rep1' in rep_lower:
        return '-'
    elif '2' in rep_lower or 'rep2' in rep_lower:
        return '--'
    elif '3' in rep_lower or 'rep3' in rep_lower:
        return ':'
    else:
        return '-'

# ==========================================
# 2. Load data from project root
# ==========================================
ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_PATH = ROOT_DIR / "data" / "processed" / "invivo_metrics" / "Rarefaction_Coordinates.csv"

if not DATA_PATH.exists():
    raise FileNotFoundError(
        f"Cannot find coordinates at {DATA_PATH}. Please run src/module2_invivo_tracing/03_rarefaction.py first."
    )

df_coords = pd.read_csv(DATA_PATH)

# ==========================================
# 3. Facet plot rendering
# ==========================================
# Expect the Brain_Region column to contain entries like V1_Cortex and LGN
regions = df_coords['Brain_Region'].unique()

fig, axes = plt.subplots(1, len(regions), figsize=(14, 6), sharey=False)

# Make sure axes is iterable when there is only one region
if len(regions) == 1:
    axes = [axes]

print(f"\n{'='*70}")
print(f"{'Region':<12} {'Group':<18} {'Replicate':<12} {'Total Pairs':>14} {'Unique Barcodes':>16}")
print(f"{'-'*70}")

for ax, region in zip(axes, regions):
    region_data = df_coords[df_coords['Brain_Region'] == region]
    
    # Collect all unique samples for this region
    unique_samples = region_data[['Group', 'Replicate']].drop_duplicates()
    
    for _, row in unique_samples.iterrows():
        group = row['Group']
        replicate = row['Replicate']
        
        # Extract the sample curve for this group and replicate
        sample_curve = region_data[(region_data['Group'] == group) & (region_data['Replicate'] == replicate)]
        
        # Endpoint stats
        total_pairs = int(sample_curve['Sampled_Pairs'].max())
        unique_bcs = int(sample_curve['Unique_Barcodes'].max())
        print(f"{region:<12} {group:<18} {str(replicate):<12} {total_pairs:>14,} {unique_bcs:>16,}")
        
        # Determine line style and color
        color = CUSTOM_COLORS.get(group, '#333333')
        linestyle = determine_linestyle(replicate)
        
        ax.plot(
            sample_curve['Sampled_Pairs'], 
            sample_curve['Unique_Barcodes'], 
            color=color, 
            linestyle=linestyle,
            linewidth=2, 
            alpha=0.85, 
            label=f"{group} ({replicate})"
        )
        # Endpoint scatter anchor
        ax.scatter(total_pairs, unique_bcs, color=color, s=40, zorder=5, alpha=0.9)
        
    ax.set_title(f'{region}', fontsize=14, fontweight='bold', pad=10)
    ax.set_xlabel('Sampled Cell-Barcode Observations', fontsize=12)
    ax.set_ylabel('Discovered Unique RV-Barcodes', fontsize=12)
    
    ax.grid(True, linestyle='--', alpha=0.4)
    sns.despine(ax=ax)

# ==========================================
# 4. Legend management and export
# ==========================================
# Collect unique legend handles and labels
handles, labels = axes[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))

fig.legend(
    by_label.values(), 
    by_label.keys(), 
    loc='center left', 
    bbox_to_anchor=(1.02, 0.5), 
    title="Experimental Group & Replicate",
    frameon=False
)

plt.tight_layout()

OUTPUT_IMAGE = ROOT_DIR / "data" / "processed" / "Fig4_Rarefaction_Analysis.svg"
plt.savefig(OUTPUT_IMAGE, format='svg', bbox_inches='tight')
print(f"✅ Publication-ready SVG saved to: {OUTPUT_IMAGE}")

plt.show()