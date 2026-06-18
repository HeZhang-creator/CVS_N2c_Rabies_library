# 导入绘图依赖
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ==========================================
# 全局格式设置 (Nature Methods 标准)
# ==========================================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 1.2

CUSTOM_COLORS = {
    'HH + HdVRz': '#1F77B4',
    'HH + Twister': '#D62728',
    'RzB + Twister': '#FF7F0E'
}

LINESTYLES = {
    'Plasmid': '-',
    'EnvA': '--',
    'B19G': ':',
    'P1': '-.'
}

# ==========================================
# 读取 03 脚本预计算好的数据
# ==========================================
metrics_dir = Path("../data/processed/invitro_metrics")
plot_dir = Path("../plots/supplementary")
plot_dir.mkdir(parents=True, exist_ok=True)

df_dist = pd.read_csv(metrics_dir / "levenshtein_distances.csv")

# ==========================================
# 绘图 1：每个 Group 内部的跨 Stage 对比 (Individual Plots)
# ==========================================
groups = df_dist['Group'].unique()

for group in groups:
    plt.figure(figsize=(7, 5))
    df_group = df_dist[df_dist['Group'] == group]
    color = CUSTOM_COLORS.get(group, '#333333')
    
    for stage in df_group['Stage'].unique():
        df_stage = df_group[df_group['Stage'] == stage]['Levenshtein_Distance']
        mean_dist = df_stage.mean()
        col_rate = (df_stage == 0).sum() / len(df_stage)
        ls = LINESTYLES.get(stage, '-')
        
        # 使用 seaborn 的 histplot 绘制阶梯密度图 (比纯 matplotlib 更丝滑)
        sns.histplot(
            df_stage, bins=np.arange(-0.5, df_stage.max() + 1.5, 1), 
            element="step", fill=False, stat="density", common_norm=False,
            color=color, linestyle=ls, linewidth=2.5, alpha=0.9,
            label=f"{stage} (Mean: {mean_dist:.1f}, Col: {col_rate:.1%})"
        )

    plt.xlabel("Pairwise Levenshtein Distance", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.title(f"Sequence Space Analysis: {group}", fontsize=14, fontweight='bold', pad=10)
    plt.legend(frameon=False, loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig(plot_dir / f"Levenshtein_Individual_{group.replace(' ', '_')}.svg", dpi=300)
    plt.show()

# ==========================================
# 绘图 2：所有 Group 和 Stage 的终极汇总大图 (Combined Plot)
# ==========================================
plt.figure(figsize=(10, 6))

for group in groups:
    df_group = df_dist[df_dist['Group'] == group]
    color = CUSTOM_COLORS.get(group, '#333333')
    
    for stage in df_group['Stage'].unique():
        df_stage = df_group[df_group['Stage'] == stage]['Levenshtein_Distance']
        mean_dist = df_stage.mean()
        ls = LINESTYLES.get(stage, '-')
        
        sns.histplot(
            df_stage, bins=np.arange(-0.5, df_stage.max() + 1.5, 1), 
            element="step", fill=False, stat="density", common_norm=False,
            color=color, linestyle=ls, linewidth=2, alpha=0.8,
            label=f"{group} - {stage} (Mean: {mean_dist:.1f})"
        )

plt.xlabel("Pairwise Levenshtein Distance", fontsize=12)
plt.ylabel("Probability Density", fontsize=12)
plt.title("Combined Sequence Space Collision Analysis", fontsize=14, fontweight='bold')

# 把图例放在外侧
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
sns.despine()
plt.tight_layout()
plt.savefig(plot_dir / "Levenshtein_Combined_All.svg", dpi=300)
plt.show()

print("🎉 编辑距离渲染完成！所有高精度 SVG 已保存在 plots 文件夹下。")