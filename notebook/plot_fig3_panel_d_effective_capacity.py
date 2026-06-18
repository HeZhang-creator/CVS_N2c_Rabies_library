"""
Script: plot_fig3_panel_d_effective_capacity.py
Location: /notebook/
Description: 
    Generates Figure 3 Panel D.
    Calculates Shannon entropy and effective barcode capacity (E = e^H) 
    to visualize the library bottleneck during viral packaging.
"""

import pandas as pd
import numpy as np
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
# 1. Calculation Module
# ==========================================
def calculate_shannon_entropy(file_path):
    """
    Calculates Shannon Entropy (H) from the UMI distribution.
    """
    counts = []
    if not os.path.exists(file_path): 
        return None
        
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5 or parts[0].lower() == 'count': 
                continue
            try:
                umi = float(parts[4])
                if umi > 0: 
                    counts.append(umi)
            except Exception: 
                continue
                
    if not counts: 
        return 0.0
        
    counts = np.array(counts)
    probabilities = counts / np.sum(counts)
    probabilities = probabilities[probabilities > 0]
    
    return -np.sum(probabilities * np.log2(probabilities))

# ==========================================
# 2. Main Execution Flow
# ==========================================
def main():
    print("Starting Effective Capacity Analysis (Panel D)...")

    try:
        ROOT = Path(__file__).resolve().parents[1]
    except NameError:
        ROOT = Path.cwd().parent
        
    data_dir = ROOT /  "data" / "processed" / "invitro_qc"
    out_dir = ROOT / "data" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    meta_path = data_dir / "Normalized_Metadata.csv"
    
    if not meta_path.exists():
        print("❌ Error: Metadata file not found at {meta_path}")
        return
        
    meta_df = pd.read_csv(meta_path)
    
    # Unified stage mapping
    stage_map = {
        'plasmid': 'Plasmid',
        'enva': 'EnvA',
        'b19g': 'P1',
        'p1': 'P1'
    }
    
    plot_data = []
    for index, row in meta_df.iterrows():
        raw_stage = str(row['Stage']).lower()
        if raw_stage not in stage_map: 
            continue
            
        stage_label = stage_map[raw_stage]
        file_name = row['Normalized_File']
        file_path = data_dir / file_name
        
        entropy = calculate_shannon_entropy(file_path)
        
        if entropy is not None:
            # Calculate Effective Capacity: E = e^H
            # Note: converting bits (log2) to nats (ln) requires division by log2(e)
            effective_bcs = np.exp(entropy * np.log(2)) 
            
            plot_data.append({
                'Group': row['Group'],
                'Stage': stage_label,
                'Effective_BCs': effective_bcs,
                'Replicate': row.get('Replicate', f'Rep{index}')
            })

    if not plot_data:
        print("❌ No valid data extracted.")
        return

    df = pd.DataFrame(plot_data)

    # ==========================================
    # 3. Visualization
    # ==========================================
    fig, ax = plt.subplots(figsize=(7, 6))
    stage_order = ['Plasmid', 'P1', 'EnvA']

    sns.pointplot(data=df, x='Stage', y='Effective_BCs', hue='Group',
                  palette=CUSTOM_COLORS, order=stage_order, 
                  markers=['o', 's', '^'], markersize=10, 
                  err_kws={'linewidth': 2}, capsize=0.1, ax=ax)

    ax.set_yscale('log')
    
    # Add text annotations for data points
    mean_div = df.groupby(['Group', 'Stage'])['Effective_BCs'].mean().reset_index()

    for _, row in mean_div.iterrows():
        group = row['Group']
        stage = row['Stage']
        val = row['Effective_BCs']
        if stage not in stage_order: continue
        x_pos = stage_order.index(stage)
        
        if stage == 'Plasmid':
            y_offset = val * 1.3 if group == 'HH + HdVRz' else val * 0.65
        else:
            if group == 'HH + Twister': y_offset = val * 1.3
            elif group == 'RzB + Twister': y_offset = val * 0.65
            else: y_offset = val * 0.5
                
        color = CUSTOM_COLORS[group]
        ax.text(x_pos, y_offset, f"{int(val):,}", 
                ha='center', va='center', fontsize=10, color=color, fontweight='bold')

    # Formatting
    sns.despine(trim=False)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='Library Design', frameon=False, loc='lower left')

    ax.set_xlabel('Viral Production Stage', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_ylabel('Effective Library Capacity (Effective BCs)', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_title('Viral Library Bottleneck During Packaging', fontsize=16, fontweight='bold', pad=15)

    plt.tight_layout()
    save_path = out_dir /"Fig3_Panel_D_Effective_Capacity.svg"
    plt.savefig(save_path, format='svg', dpi=300)
    print("🎉 Plot saved to: {save_path}")
    plt.show()
if __name__ == "__main__":
    main()
