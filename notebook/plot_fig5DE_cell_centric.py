"""
Figure Generator: Figure 5D-E (Cell-Centric Collision Analysis with GridSpec Layout)
Description: 
    Quantifies cell-centric collision status by plotting Collision-free vs. Collision-affected 
    cell percentages across brain regions (V1 and LGN) using a highly precise 2.3:1 global aspect ratio grid.
    Features robust error bars (SEM), Welch's t-test significance bars, and automated folder routing.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.ticker import PercentFormatter
from scipy.stats import ttest_ind

# ==========================================
# 1. Publication Aesthetics Configuration
# ==========================================
mpl.rcParams['svg.fonttype'] = 'none'
sns.set_theme(style="ticks", font_scale=1.1)

GROUP_ORDER = ['HH + HdVRz', 'HH + Twister', 'RzB + Twister']
REGIONS = ['V1_Cortex', 'LGN']

CUSTOM_COLORS = {
    'HH + HdVRz': '#1F77B4', 
    'HH + Twister': '#D62728', 
    'RzB + Twister': '#FF7F0E'
}

# ==========================================
# 2. Statistical Computation Engine (Welch's t-test)
# ==========================================
def draw_stats(ax, data, group1, group2, metric, x1, x2, y_pos, h_line=1.0):
    g1_data = data[data['Group'] == group1][metric]
    g2_data = data[data['Group'] == group2][metric]
    
    if len(g1_data) > 1 and len(g2_data) > 1:
        # Welch's t-test
        _, p_val = ttest_ind(g1_data, g2_data, equal_var=False)
        
        if p_val < 0.001: ast = '***'
        elif p_val < 0.01: ast = '**'
        elif p_val < 0.05: ast = '*'
        elif p_val <= 0.10: ast = '#'
        else: ast = 'ns'
        
        if ast != 'ns':
            ax.plot([x1, x1, x2, x2], [y_pos, y_pos + h_line, y_pos + h_line, y_pos], lw=1.2, c='black')
            ax.text((x1 + x2) * 0.5, y_pos + h_line, ast, ha='center', va='bottom', color='black', fontsize=12, fontweight='bold')
            return y_pos + h_line + (h_line * 1.5)
            
    return y_pos

# ==========================================
# 3. Pixel-Level Layout & Rendering Module
# ==========================================
def plot_pixel_perfect_panels(plot_data, output_prefix):
    print("🎨 Rendering Figure 5D-E with exact aspect ratios...")
    
    # Target 1: Global canvas aspect ratio forced to 2.3 : 1
    total_width_in = 20.0 
    total_height_in = total_width_in / 2.3 

    fig = plt.figure(figsize=(total_width_in, total_height_in))

    # Precise margins allocation via GridSpec
    left_margin = 0.06
    right_margin = 0.98
    bottom_margin = 0.15
    top_margin = 0.85
    wspace = 0.25 

    gs = fig.add_gridspec(1, 4, wspace=wspace, left=left_margin, right=right_margin, bottom=bottom_margin, top=top_margin)
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]

    plot_configs = []
    for region in REGIONS:
        # Collision-free cell: all barcodes globally unique (minority, low %)
        plot_configs.append({'region': region, 'metric': 'Collision_Free_Pct', 'title': 'Collision-free cell'})
        # Collision-affected cell: barcode shared with other cells (majority, high %)
        plot_configs.append({'region': region, 'metric': 'Collision_Affected_Pct', 'title': 'Collision-affected cell'})

    for idx, config in enumerate(plot_configs):
        ax = axes[idx]
        
        # Target 2: Force each individual subplot area to a 1.7 : 1 aspect ratio
        ax.set_box_aspect(1.7 / 1) 
        
        region = config['region']
        metric = config['metric']
        is_collision_affected = (metric == 'Collision_Affected_Pct')
        
        region_data = plot_data[plot_data['Brain_Region'] == region].copy()
        
        # Render Barplot with strictly isolated SEM tracking
        sns.barplot(
            data=region_data, x='Group', y=metric, ax=ax, 
            hue='Group', palette=CUSTOM_COLORS, order=GROUP_ORDER, dodge=False,
            capsize=0.15, errorbar='se', edgecolor='black', linewidth=1.5, alpha=0.85
        )
        
        if ax.get_legend() is not None: 
            ax.get_legend().remove()
        
        # Superimpose individual animal data points
        sns.stripplot(
            data=region_data, x='Group', y=metric, ax=ax,
            order=GROUP_ORDER, color='#222222', size=7, jitter=0.12, 
            edgecolor='white', linewidth=1.2, alpha=0.9
        )
        
        # 🌟 Stat bar height calculations
        base_offset = 1.5 if is_collision_affected else 3.0
        line_height = 1.0 if is_collision_affected else 2.0
        max_y_in_region = region_data[metric].max() + base_offset
        
        current_h = draw_stats(ax, region_data, GROUP_ORDER[0], GROUP_ORDER[1], metric, 0, 1, max_y_in_region, h_line=line_height)
        current_h = draw_stats(ax, region_data, GROUP_ORDER[0], GROUP_ORDER[2], metric, 0, 2, current_h, h_line=line_height)

        # Labels and floating structural boxes
        ax.set_title(config['title'], fontsize=14, fontweight='bold', pad=25)
        
        bg_color = '#e0f0ff' if 'V1' in region else '#fff0e0'
        ax.text(0.5, 0.98, f"{region.replace('_', ' ')}", transform=ax.transAxes, 
                fontsize=15, fontweight='bold', ha='center', va='top', 
                bbox=dict(facecolor=bg_color, edgecolor='black', boxstyle='round,pad=0.3', alpha=0.9))

        # X-axis configuration
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels([g.replace(' + ', ' +\n') for g in GROUP_ORDER], rotation=0, fontsize=12)
        for j, tick in enumerate(ax.get_xticklabels()):
            tick.set_color(CUSTOM_COLORS[GROUP_ORDER[j]])
            tick.set_fontweight('bold')
            
        ax.set_xlabel("")
        
        # Only preserve Y-axis title on the absolute leftmost panel
        if idx == 0:
            ax.set_ylabel("Proportion of Total Cells", fontsize=14, fontweight='bold')
        else:
            ax.set_ylabel("")
            
        ax.yaxis.set_major_formatter(PercentFormatter(100))
        
        # Y-axis domain: Collision-free cells are rare (~1-25%), Collision-affected dominate (~75-100%)
        if is_collision_affected:
            # Collision-affected cells are the majority → high Y-axis domain
            ax.set_yticks([70, 75, 80, 85, 90, 95, 100])
            ax.set_ylim(70, max(100, current_h + 2))
        else:
            # Collision-free cells are the minority → low Y-axis domain
            ax.set_yticks([0, 5, 10, 15, 20, 25, 30])
            ax.set_ylim(0, max(32, current_h + 2))
            
        ax.grid(axis='y', linestyle='--', alpha=0.4, color='gray')
        ax.set_axisbelow(True)

    sns.despine()
    for ax in axes:
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)

    # Export vector files and data tables to the consolidated results directory
    out_pdf = f"{output_prefix}.pdf"
    out_svg = f"{output_prefix}.svg"
    plt.savefig(out_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(out_svg, format='svg', dpi=300, bbox_inches='tight')
    print(f"✅ Success! Pixel-perfect figures saved to: {out_pdf}")
    plt.show()

# ==========================================
# 4. Runtime Main Execution
# ==========================================
if __name__ == "__main__":
    ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if ROOT_DIR not in sys.path:
        sys.path.insert(0, ROOT_DIR)

    CSV_PATH = os.path.join(ROOT_DIR, "data", "processed", "invivo_metrics", "Master_Corrected_Data.csv")
    
    # Establish uniform publication results routing
    RESULTS_DIR = os.path.join(ROOT_DIR, "data", "results")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    OUTPUT_PREFIX = os.path.join(RESULTS_DIR, "Fig5DE_Cell_Centric_Fidelity_PerfectGrid")
    
    try:
        # Read pre-computed metrics from 04_collision_resolution.py
        CELL_CSV = os.path.join(
            ROOT_DIR, "data", "processed", "invivo_metrics",
            "Fig5_Resolution_Metrics_Cell_Centric.csv"
        )
        if not os.path.exists(CELL_CSV):
            raise FileNotFoundError(
                f"Cell-centric metrics not found at {CELL_CSV}. "
                f"Run 04_collision_resolution.py first."
            )
        data = pd.read_csv(CELL_CSV)
        
        print("\n📊 === Loaded Plot Data (Source Data) ===")
        print(data.to_string(index=False))
        print("===============================================\n")

        # 2. Fire the custom layout plot engine
        plot_pixel_perfect_panels(data, OUTPUT_PREFIX)
        
    except FileNotFoundError:
        print(f"⚠️ Dataset missing at: {CSV_PATH}. Please verify project path mapping.")