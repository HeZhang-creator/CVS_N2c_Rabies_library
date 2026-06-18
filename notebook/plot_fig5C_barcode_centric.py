"""
Figure Generator: Figure 5C (Barcode-Centric Collision Profiles by Replicate)
Description: 
    Generates a stacked bar chart illustrating the structural composition of the viral libraries.
    Displays EACH biological replicate independently to demonstrate inter-animal reproducibility.
    Categorizes RV-barcode species into four collision levels based on host cell infection count.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import os
import sys

# ==========================================
# 1. Publication Aesthetics Configuration
# ==========================================
# Ensure text is editable in Adobe Illustrator
mpl.rcParams['svg.fonttype'] = 'none'
sns.set_theme(style="ticks", font_scale=1.1)

# Exact group ordering from the manuscript
GROUP_ORDER = ['HH + HdVRz', 'HH + Twister', 'RzB + Twister']
BRAIN_REGIONS = ['V1_Cortex', 'LGN']

# 🎯 定制颜色字典 (来自您的需求)
CUSTOM_COLORS = {
    'HH + HdVRz': '#1F77B4', 
    'HH + Twister': '#D62728', 
    'RzB + Twister': '#FF7F0E'
}

# 🎯 碰撞等级颜色与名称 (来自您的需求)
COLLISION_COLORS = ["#3C8DBC", "#00A087", "#F39C12", "#E64B35"] 
LEVEL_NAMES = ['1 (Unique)', '2 (Low)', '3 (Med)', '>3 (High)']

# ==========================================
# 2. Plotting Module
# ==========================================
def plot_stacked_collision(plot_data, output_prefix):
    print("🎨 Rendering Figure 5C (Replicate-resolved)...")
    
    # 增加图片宽度以适应更多的柱子 (e.g., 9 bars per region)
    fig, axes = plt.subplots(1, len(BRAIN_REGIONS), figsize=(16, 6), sharey=True)
    if len(BRAIN_REGIONS) == 1: 
        axes = [axes]

    for i, region in enumerate(BRAIN_REGIONS):
        region_data = plot_data[plot_data['Brain_Region'] == region].copy()
        
        # 强制规定 Group 的排序，并按 Group 和 Rep 排序
        region_data['Group'] = pd.Categorical(region_data['Group'], categories=GROUP_ORDER, ordered=True)
        region_data = region_data.sort_values(['Group', 'Replicate'])
        
        # 🌟 构建独立的 X 轴标签： "组名 \n Rep X"
        region_data['X_Label'] = region_data['Group'].astype(str) + "\nRep " + region_data['Replicate'].astype(str)
        
        # Pivot data for stacked bar plotting
        pivot_df = region_data.pivot(index='X_Label', columns='Collision_Level', values='Proportion')
        
        # 提取排序好的 X 轴列表，安全地 Reindex 并填充缺失值为 0
        ordered_x_labels = region_data['X_Label'].drop_duplicates().tolist()
        pivot_df = pivot_df.reindex(index=ordered_x_labels, columns=LEVEL_NAMES).fillna(0)
        
        # 绘图
        ax = axes[i]
        pivot_df.plot(
            kind='bar', stacked=True, ax=ax, color=COLLISION_COLORS,
            edgecolor='black', linewidth=0.8, width=0.8, legend=False
        )
        
        ax.set_title(f"Brain Region: {region}", fontsize=15, fontweight='bold', pad=10)
        # 将文字倾斜 45 度，防止拥挤
        ax.set_xticklabels(ordered_x_labels, rotation=45, ha='right', fontweight='bold', fontsize=11)
        ax.set_xlabel("")
        
        if i == 0:
            ax.set_ylabel("Proportion of Barcode Species (%)", fontsize=14, fontweight='bold')
        
        sns.despine(ax=ax, top=True, right=True)

        # 🌟 绝佳美学修饰：将底部的文字按组别上色，呼应您的 custom_colors
        for tick_label in ax.get_xticklabels():
            # 提取原组名（截掉 "\nRep X"）
            group_name = tick_label.get_text().split('\nRep')[0]
            if group_name in CUSTOM_COLORS:
                tick_label.set_color(CUSTOM_COLORS[group_name])

    # 统一图例
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], title='Collision Severity', 
               bbox_to_anchor=(1.02, 0.85), loc='upper left', frameon=False, fontsize=11, title_fontsize=12)

    plt.tight_layout()
    
    out_pdf = f"{output_prefix}.pdf"
    out_svg = f"{output_prefix}.svg"
    plt.savefig(out_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(out_svg, format='svg', dpi=300, bbox_inches='tight')
    print(f"✅ Saved to {out_pdf}")

if __name__ == "__main__":
    # ===================== Project Root Configuration =====================
    ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if ROOT_DIR not in sys.path:
        sys.path.insert(0, ROOT_DIR)

    CSV_PATH = os.path.join(ROOT_DIR, "data", "processed", "invivo_metrics", "Master_Corrected_Data.csv")
    
    RESULTS_DIR = os.path.join(ROOT_DIR, "data", "results")
    os.makedirs(RESULTS_DIR, exist_ok=True)

    OUTPUT_PREFIX = os.path.join(RESULTS_DIR, "Fig5C_Barcode_Centric_Collision_ByRep")
    
    try:
        # Read pre-computed metrics from 04_collision_resolution.py
        BARCODE_CSV = os.path.join(
            ROOT_DIR, "data", "processed", "invivo_metrics",
            "Fig5_Resolution_Metrics_Barcode_Centric.csv"
        )
        if not os.path.exists(BARCODE_CSV):
            raise FileNotFoundError(
                f"Barcode-centric metrics not found at {BARCODE_CSV}. "
                f"Run 04_collision_resolution.py first."
            )
        data = pd.read_csv(BARCODE_CSV)
        
   
        print("\n📊 === Loaded Plot Data (Source Data) ===")
        print(data.to_string(index=False))
        print("===============================================\n")

        plot_stacked_collision(data, OUTPUT_PREFIX)
        
    except FileNotFoundError:
        print(f"⚠️ Required files not found. Please run 04_collision_resolution.py first.")