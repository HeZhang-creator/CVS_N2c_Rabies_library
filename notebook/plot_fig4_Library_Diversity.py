"""
Figure Generator: Fig4_Library_Diversity
Description: 
    4-panel bar plots — V1 (dark) left, LGN (light) right.
    Bar positions are set manually with ax.bar() to guarantee ordering.
    Stripplot dots use ax.scatter() with tight jitter aligned to the same positions.
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from pathlib import Path
import os
import sys

mpl.rcParams['svg.fonttype'] = 'none'

# ==========================================
# 1. Colour palette & helpers
# ==========================================
def get_light_color(hex_color, white_mix=0.55):
    c = np.array(mcolors.to_rgb(hex_color))
    return tuple(c + (np.array([1., 1., 1.]) - c) * white_mix)

BASE_COLORS = ['#1F77B4', '#D62728', '#FF7F0E']
GROUPS     = ['HH + HdVRz', 'HH + Twister', 'RzB + Twister']
REGIONS    = ['V1', 'LGN']
N_GROUPS   = len(GROUPS)

DARK_COLORS  = [mcolors.to_rgba(c, alpha=0.85) for c in BASE_COLORS]
LIGHT_COLORS = [mcolors.to_rgba(get_light_color(c), alpha=0.85) for c in BASE_COLORS]

# ── Bar geometry (fixed, V1 always -offset, LGN always +offset) ──
X_CENTERS = np.arange(N_GROUPS)
BAR_WIDTH = 0.35
DODGE     = BAR_WIDTH * 0.50          # half the gap between paired bars
OFFSETS   = [-DODGE, +DODGE]          # V1 left, LGN right

def bar_x(gi, ri):
    """Return the bar centre x for group index gi, region index ri (0=V1, 1=LGN)."""
    return X_CENTERS[gi] + OFFSETS[ri]

def draw_sig_bracket(ax, x1, x2, y_base, label, gap_ratio=0.06, h_ratio=0.02):
    ylim = ax.get_ylim()
    yr = ylim[1] - ylim[0]
    y0 = y_base + yr * gap_ratio
    y1 = y0 + yr * h_ratio
    ax.plot([x1, x1, x2, x2], [y0, y1, y1, y0],
            lw=1.0, c='black', clip_on=False)
    ax.text((x1 + x2) / 2, y1, label, ha='center', va='bottom', fontsize=11)

# ==========================================
# 2. Main rendering
# ==========================================
def render_figure(metrics_csv, output_prefix):
    input_file = Path(metrics_csv)
    if not input_file.exists():
        raise FileNotFoundError(f"❌ Metrics table not found: {input_file}")

    audit_df = pd.read_csv(input_file)
    audit_df['Brain_Region'] = audit_df['Brain_Region'].replace('V1_Cortex', 'V1')

    print("🎨 Rendering publication-ready panels...")

    sns.set_theme(style="ticks", font_scale=1.1)
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.5))

    metrics = [
        ('Absolute_Diversity', 'Absolute Barcode Diversity', 'Total Unique Barcodes'),
        ('Effective_Capacity', 'Effective Barcode Capacity', 'Effective Capacity (from entropy)'),
        ('Information_Density', 'Sequencing Efficiency', 'Information Density\n(Barcodes per 10,000 UMIs)'),
        ('Total_UMI', 'Sequencing Depth Baseline', 'Total Valid UMIs')
    ]

    panel_brackets = {
        'Absolute_Diversity': [
            (0, 0, 1, 0, '*'),   # V1: canon vs HH+Twister
            (0, 0, 2, 0, '*'),   # V1: canon vs RzB+Twister
            (0, 1, 1, 1, '*'),   # LGN: canon vs HH+Twister
            (0, 1, 2, 1, '*'),   # LGN: canon vs RzB+Twister
        ],
        'Information_Density': [
            (0, 1, 1, 1, '*'),
            (0, 1, 2, 1, '*'),
        ],
        'Total_UMI': [
            (0, 0, 2, 0, '*'),   # V1: canon vs RzB+Twister (p=0.024)
        ],
    }

    for ax, (y_col, title, ylabel) in zip(axes, metrics):
        # ── Pre-compute means, SEMs, and raw values ──
        means   = np.zeros((N_GROUPS, 2))
        sems    = np.zeros((N_GROUPS, 2))
        raw_vals = [[None, None] for _ in range(N_GROUPS)]
        for gi, g in enumerate(GROUPS):
            for ri, r in enumerate(REGIONS):
                vals = audit_df[(audit_df['Group'] == g) & (audit_df['Brain_Region'] == r)][y_col]
                means[gi, ri] = vals.mean()
                sems[gi, ri]  = vals.sem()
                raw_vals[gi][ri] = vals.values

        # ── Bars ──
        for gi in range(N_GROUPS):
            for ri in range(2):
                color = DARK_COLORS[gi] if ri == 0 else LIGHT_COLORS[gi]
                bx = bar_x(gi, ri)
                ax.bar(bx, means[gi, ri], BAR_WIDTH,
                       color=color, edgecolor='black', linewidth=1.2,
                       yerr=sems[gi, ri], capsize=4,
                       error_kw={'linewidth': 1.5, 'color': 'black'},
                       zorder=2)

        # ── Stripplot dots (manual scatter, same bar_x positions) ──
        for gi in range(N_GROUPS):
            for ri in range(2):
                vals = raw_vals[gi][ri]
                if len(vals) == 0:
                    continue
                bx = bar_x(gi, ri)
                rng = np.random.default_rng(gi * 2 + ri + 42)
                x_jittered = np.full(len(vals), bx) + rng.uniform(-0.04, 0.04, len(vals))
                ax.scatter(x_jittered, vals, s=22, c='#333333',
                           edgecolors='white', linewidths=0.6, alpha=0.9, zorder=5)

        # ── Region labels below each bar ──
        for gi in range(N_GROUPS):
            for ri in range(2):
                bx = bar_x(gi, ri)
                ax.text(bx, -0.012, REGIONS[ri],
                        transform=ax.get_xaxis_transform(),
                        ha='center', va='top', fontsize=8.5,
                        fontweight='bold', color='#333333')

        # ── Axes styling ──
        ax.set_title(title, fontsize=15, fontweight='bold', pad=15)
        ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
        ax.set_xlabel('')
        ax.set_xticks(X_CENTERS)
        ax.set_xticklabels([g.replace(' + ', ' +\n') for g in GROUPS], fontweight='bold')

        if y_col == 'Total_UMI':
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: "{:,}".format(int(v))))

        sns.despine(ax=ax)
        ax.set_xlim(X_CENTERS[0] - 0.7, X_CENTERS[-1] + 0.7)

        # ── Significance brackets ──
        bracks = panel_brackets.get(y_col, [])
        if bracks:
            for gi1, ri1, gi2, ri2, label in bracks:
                x1 = bar_x(gi1, ri1)
                x2 = bar_x(gi2, ri2)
                y_vis = max(
                    means[gi1, ri1] + sems[gi1, ri1],
                    means[gi2, ri2] + sems[gi2, ri2],
                    raw_vals[gi1][ri1].max() if len(raw_vals[gi1][ri1]) else 0,
                    raw_vals[gi2][ri2].max() if len(raw_vals[gi2][ri2]) else 0,
                )
                draw_sig_bracket(ax, x1, x2, y_vis, label)
            ylim = ax.get_ylim()
            ax.set_ylim(ylim[0], ylim[1] * 1.12)

    # ── Panel letters ──
    for ax, letter in zip(axes, ['B', 'C', 'F', 'G']):
        ax.text(-0.10, 1.05, letter, transform=ax.transAxes, fontsize=20,
                fontweight='bold', va='top')

    plt.tight_layout(w_pad=3.0)
    plt.subplots_adjust(bottom=0.12)

    out_svg = Path(f"{output_prefix}.svg")
    out_pdf = Path(f"{output_prefix}.pdf")
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_svg, format='svg', dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, format='pdf', dpi=300, bbox_inches='tight')
    print(f"✅ Panels saved to {out_pdf}")
    plt.show()

if __name__ == "__main__":
    ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if ROOT_DIR not in sys.path:
        sys.path.insert(0, ROOT_DIR)
    CSV_PATH = os.path.join(ROOT_DIR, "data", "processed", "invivo_metrics",
                            "Diversity_Metrics_Summary.csv")
    print("📁 Project root directory:", ROOT_DIR)
    print("📄 Data path:", CSV_PATH)
    print("✅ File exists:", os.path.exists(CSV_PATH))
    OUTPUT_PREFIX = os.path.join(ROOT_DIR, "data", "processed",
                                 "Fig4_Library_Diversity_V1_LGN")
    render_figure(CSV_PATH, OUTPUT_PREFIX)
