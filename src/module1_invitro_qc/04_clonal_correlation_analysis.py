"""
Script: 04_clonal_correlation_analysis.py
Description: 
    Data processing engine for clonal fidelity analysis (Plasmid vs EnvA).
    Computes Spearman rank correlations and prepares aligned frequency data.
    *Note: Plotting is delegated to Jupyter Notebooks for visual fine-tuning.*
"""

import pandas as pd
import os
from pathlib import Path
from scipy.stats import spearmanr

# ==========================================
# 1. Core Computation Module: Data Extraction and Alignment
# ==========================================
def get_barcode_data(file_path):
    """Read normalized files from script 01, extract Barcode and percentage Frequency"""
    if not os.path.exists(file_path):
        return pd.DataFrame(columns=['Barcode', 'Frequency'])
    
    df = pd.read_csv(file_path, sep='\t')
    total_umi = df['UMI_Count'].sum()
    if total_umi == 0: return pd.DataFrame(columns=['Barcode', 'Frequency'])
    
    df['Frequency'] = (df['UMI_Count'] / total_umi) * 100
    return df[['Barcode', 'Frequency']]

def main():
    print("🚀 Starting Clonal Fidelity data calculation...")
    
    # Path settings (adjust according to your ROOT)
   # Path settings
    ROOT = Path(__file__).resolve().parents[2]
    data_dir = ROOT / "data" / "processed" / "invitro_qc"
    meta_path = data_dir / "Normalized_Metadata.csv"
    out_dir = ROOT / "data" / "processed" / "invitro_metrics"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not meta_path.exists():
        print(f"❌ Metadata file not found: {meta_path}")
        return

    meta_df = pd.read_csv(meta_path)
    # Normalize stage text to avoid casing mismatches between metadata and code expectations
    meta_df['Stage_norm'] = meta_df['Stage'].astype(str).str.strip().str.lower()
    groups = meta_df['Group'].unique()
    
    all_scatter_data = []
    stats_records = []

    for group in groups:
        print(f"  ⏳ Aligning and calculating group: {group}")
        try:
            # Extract corresponding files based on normalized stage values
            p_row = meta_df[(meta_df['Group'] == group) & (meta_df['Stage_norm'] == 'plasmid')]
            e_row = meta_df[(meta_df['Group'] == group) & (meta_df['Stage_norm'] == 'enva')]

            p_file = data_dir / p_row['Normalized_File'].values[0]
            e_file = data_dir / e_row['Normalized_File'].values[0]
            
            df_p = get_barcode_data(p_file)
            df_e = get_barcode_data(e_file)
            
            # Execute Inner Join to get Shared BCs
            merged = pd.merge(df_p, df_e, on='Barcode', suffixes=('_Plasmid', '_EnvA'), how='inner')
            
            if not merged.empty:
                rho, _ = spearmanr(merged['Frequency_Plasmid'], merged['Frequency_EnvA'])
                
                # Record statistical results
                stats_records.append({
                    'Group': group,
                    'Spearman_Rho': rho,
                    'Shared_BCs': len(merged)
                })
                
                # Append scatter data
                merged['Group'] = group
                all_scatter_data.append(merged)
                
        except IndexError:
            print(f"  ⚠️ Warning: {group} missing corresponding files, skipping.")

    # ==========================================
    # Save calculation results for Notebook reading
    # ==========================================
    if all_scatter_data:
        df_scatter = pd.concat(all_scatter_data, ignore_index=True)
        df_scatter.to_csv(out_dir / "fidelity_scatter_data.csv", index=False)
        
        df_stats = pd.DataFrame(stats_records)
        df_stats.to_csv(out_dir / "fidelity_stats.csv", index=False)
        
        print(f"\n🎉 Calculation complete! Data exported to: {out_dir}")
        print("👉 Next step: Open the plotting script in notebooks/ for visualization.")

if __name__ == "__main__":
    main()