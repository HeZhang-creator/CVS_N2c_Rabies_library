"""
Figure Generator: Fig6_Topology_BF
Description:
    Generates a 5-panel publication-ready figure (Panels B–F) for Figure 6,
    visualizing cell-centric network topology metrics for V1→LGN trans-synaptic
    tracing across three viral library groups.

    Panel B: Pure Network Edges — total trans-synaptic yield from collision-free V1 cells.
    Panel C: LGN Convergence (In-Degree) — number of V1 barcodes converging per LGN cell.
    Panel D: Out-Degree (All cells) — collision-inflated V1 divergence.
    Panel E: Out-Degree (Collision-free cells) — true V1 divergence on collision-free cells.
    Panel F: Distribution overlay — all cells (solid) vs collision-free cells (dashed).

    Input:
    - Fig6_Topology_Metrics_Summary.csv
    - Fig6_Topology_Metrics_OutDegree_Dist.csv
    - Fig6_Topology_Metrics_SingletOutDegree_Dist.csv

    Output:
    - Fig6_Topology_BF.pdf / .svg / .png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from matplotlib.lines import Line2D
import os

# ==========================================
# 1. Configuration
# ==========================================
GROUP_ORDER = ['HH + HdVRz', 'HH + Twister', 'RzB + Twister']
CUSTOM_COLORS = {
    'HH + HdVRz': '#1F77B4',
    'HH + Twister': '#D62728',
    'RzB + Twister': '#FF7F0E'
}

matplotlib.rcParams.update({
    'svg.fonttype': 'none',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.labelsize': 10,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'figure.dpi': 150,
})

# ==========================================
# 2. Load Data
# ==========================================
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

SUMMARY_CSV = os.path.join(
    ROOT, "data", "processed", "invivo_metrics", "Fig6_Topology_Metrics_Summary.csv"
)
OUTDEGREE_CSV = os.path.join(
    ROOT, "data", "processed", "invivo_metrics", "Fig6_Topology_Metrics_OutDegree_Dist.csv"
)
SINGLET_CSV = os.path.join(
    ROOT, "data", "processed", "invivo_metrics", "Fig6_Topology_Metrics_SingletOutDegree_Dist.csv"
)

df_summary = pd.read_csv(SUMMARY_CSV)
df_outdist = pd.read_csv(OUTDEGREE_CSV)
df_singlet = pd.read_csv(SINGLET_CSV)

# ==========================================
# 3. Helper Functions
# ==========================================
def draw_stats_bar(ax, data, metric, group1, group2, x1, x2, y_pos, h_line=1.0):
    """Draw significance brackets between two bar groups (Welch's t-test)."""
    g1 = data[data['Group'] == group1][metric].dropna().values
    g2 = data[data['Group'] == group2][metric].dropna().values
    if len(g1) > 1 and len(g2) > 1:
        _, p = ttest_ind(g1, g2, equal_var=False)
        if p < 0.001:   ast = '***'
        elif p < 0.01:  ast = '**'
        elif p < 0.05:  ast = '*'
        elif p <= 0.10: ast = '#'
        else: ast = 'ns'
        if ast != 'ns':
            ax.plot([x1, x1, x2, x2],
                    [y_pos, y_pos + h_line, y_pos + h_line, y_pos],
                    lw=1.0, c='black', clip_on=False)
            ax.text((x1 + x2) * 0.5, y_pos + h_line + h_line * 0.15, ast,
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
            return y_pos + h_line + h_line * 1.6
    return y_pos


def bar_with_points(ax, data, metric, groups, colors, y_label, y_top, title,
                    do_stats=True):
    """Bar chart with individual replicate points (pure matplotlib, no seaborn)."""
    n = len(groups)
    spacing = 0.70   # tight group spacing
    width = 0.52     # bar width
    x_pos = np.arange(n) * spacing

    means, sems = [], []
    for g in groups:
        vals = data[data['Group'] == g][metric].dropna().values
        means.append(vals.mean())
        sems.append(vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)

    ax.bar(x_pos, means, width, color=[colors[g] for g in groups],
           edgecolor='black', linewidth=1.1, capstyle='round')
    ax.errorbar(x_pos, means, yerr=sems, fmt='none', ecolor='black',
                elinewidth=1.3, capsize=5, capthick=1.3)

    for i, g in enumerate(groups):
        vals = data[data['Group'] == g][metric].dropna().values
        jitter = np.random.RandomState(42).uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), x_pos[i]) + jitter, vals,
                   s=30, color='#333333', edgecolors='white', linewidth=0.7,
                   zorder=10, alpha=0.9)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([g.replace(' + ', ' +\n') for g in groups],
                       fontweight='bold', fontsize=9)
    for i, tick in enumerate(ax.get_xticklabels()):
        tick.set_color(colors[groups[i]])
    ax.set_xlim(x_pos[0] - 0.50, x_pos[-1] + 0.50)
    ax.set_ylabel(y_label, fontsize=11, fontweight='bold')
    ax.set_ylim(0, y_top)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)
    ax.grid(axis='y', linestyle='--', alpha=0.35, color='gray')
    ax.set_axisbelow(True)

    if do_stats:
        h = draw_stats_bar(ax, data, metric, groups[0], groups[1],
                          x_pos[0], x_pos[1], max(means) * 1.05)
        _ = draw_stats_bar(ax, data, metric, groups[0], groups[2],
                          x_pos[0], x_pos[2], h)

    return ax


# ==========================================
# 4. Main Figure: 5 Panels (B–F)
# ==========================================
def plot_fig6_topology_BF(output_prefix):
    fig, axes = plt.subplots(1, 5, figsize=(24, 5.0))
    panel_letters = ['B', 'C', 'D', 'E', 'F']

    # ── Panel B: Pure Network Edges ──
    metric = 'Pure_Network_Edges'
    max_y = df_summary[metric].max()
    bar_with_points(axes[0], df_summary, metric, GROUP_ORDER, CUSTOM_COLORS,
                    'Total Trans-Synaptic Edges', max_y * 1.30,
                    'Pure Network Edges')

    # ── Panel C: LGN Convergence (In-Degree) ──
    metric = 'Mean_In_Degree'
    max_y = df_summary[metric].max()
    bar_with_points(axes[1], df_summary, metric, GROUP_ORDER, CUSTOM_COLORS,
                    'Mean In-Degree\n(V1 barcodes per LGN cell)', max_y * 1.15,
                    'LGN Convergence (In-Degree)')
    axes[1].set_ylim(1.0, max_y * 1.15)

    # ── Panel D: Out-Degree — All cells ──
    metric = 'Mean_Out_Degree'
    max_y = df_summary[metric].max()
    bar_with_points(axes[2], df_summary, metric, GROUP_ORDER, CUSTOM_COLORS,
                    'Mean Out-Degree\n(LGN targets per V1 cell)', max_y * 1.35,
                    'Out-Degree (All cells)')

    # ── Panel E: Out-Degree — Collision-free cells ──
    metric = 'Mean_Singlet_Out_Degree'
    max_y_s = df_summary[metric].max()
    bar_with_points(axes[3], df_summary, metric, GROUP_ORDER, CUSTOM_COLORS,
                    'Mean Out-Degree\n(LGN targets per V1 cell)\n[collision-free cells]',
                    max_y_s * 1.40, 'Out-Degree (Collision-free cells)')

    # ── Panel F: Distribution Overlay ──
    ax_f = axes[4]
    ax_f.set_box_aspect(1 / 1.45)
    max_bin = 40
    bins = np.arange(0, max_bin + 1, 1)
    for g in GROUP_ORDER:
        # All cells — solid
        vals_all = df_outdist[df_outdist['Group'] == g]['Out_Degree'].values
        vals_all = vals_all[vals_all <= max_bin]
        if len(vals_all) > 0:
            weights_all = np.ones(len(vals_all)) / len(
                df_outdist[df_outdist['Group'] == g])
            ax_f.hist(vals_all, bins=bins, weights=weights_all,
                      histtype='step', color=CUSTOM_COLORS[g],
                      linewidth=2.0, linestyle='-')
        # Collision-free cells — dashed
        vals_s = df_singlet[df_singlet['Group'] == g][
            'Singlet_Out_Degree'].values
        vals_s = vals_s[vals_s <= max_bin]
        if len(vals_s) > 0:
            weights_s = np.ones(len(vals_s)) / len(
                df_singlet[df_singlet['Group'] == g])
            ax_f.hist(vals_s, bins=bins, weights=weights_s,
                      histtype='step', color=CUSTOM_COLORS[g],
                      linewidth=1.6, linestyle='--')

    legend_elements = (
        [Line2D([0], [0], color=CUSTOM_COLORS[g], lw=2.0, linestyle='-',
                label=g) for g in GROUP_ORDER]
        + [Line2D([0], [0], color='black', lw=2.0, linestyle='-',
                  label='All cells')]
        + [Line2D([0], [0], color='black', lw=1.6, linestyle='--',
                  label='Collision-free cells')]
    )
    ax_f.legend(handles=legend_elements, loc='upper right', frameon=True,
                fontsize=7, handlelength=1.3)
    ax_f.set_xlim(0, max_bin)
    ax_f.set_title('Out-Degree Distribution\n(All cells vs Collision-free cells)',
                   fontsize=13, fontweight='bold', pad=12)
    ax_f.set_xlabel('Out-Degree (Targets / V1 Neuron)',
                    fontsize=11, fontweight='bold')
    ax_f.set_ylabel('Proportion of Neurons', fontsize=11, fontweight='bold')
    ax_f.spines['bottom'].set_linewidth(1.2)
    ax_f.spines['left'].set_linewidth(1.2)
    ax_f.grid(axis='both', linestyle='--', alpha=0.3)
    ax_f.set_axisbelow(True)

    # ── Panel letters & spine styling ──
    for ax, letter in zip(axes, panel_letters):
        ax.spines['bottom'].set_linewidth(1.2)
        ax.spines['left'].set_linewidth(1.2)
        ax.text(-0.08, 1.04, letter, transform=ax.transAxes,
                fontsize=18, fontweight='bold', va='top')

    plt.tight_layout(w_pad=2.0, rect=[0, 0, 1, 1])

    # ── Export ──
    for fmt in ['pdf', 'svg', 'png']:
        out_path = f"{output_prefix}.{fmt}"
        fig.savefig(out_path, format=fmt, dpi=300, bbox_inches='tight')
        print(f"  ✓ {out_path}")

    plt.close('all')
    return fig


# ==========================================
# 5. Runtime Entry Point
# ==========================================
if __name__ == "__main__":
    OUTPUT_PREFIX = os.path.join(
        ROOT, "data", "processed", "Fig6_Topology_BF"
    )
    print(f"📁 Project root: {ROOT}")
    print(f"📄 Summary: {SUMMARY_CSV} — {'✓' if os.path.exists(SUMMARY_CSV) else '✗'}")
    print(f"📄 OutDegree: {OUTDEGREE_CSV} — {'✓' if os.path.exists(OUTDEGREE_CSV) else '✗'}")
    print(f"📄 Singlet: {SINGLET_CSV} — {'✓' if os.path.exists(SINGLET_CSV) else '✗'}")
    print("🎨 Rendering Fig 6 B–F ...")
    plot_fig6_topology_BF(OUTPUT_PREFIX)
    print("🎉 Done!")
