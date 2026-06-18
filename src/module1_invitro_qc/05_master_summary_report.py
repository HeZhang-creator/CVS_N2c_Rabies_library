"""
Script: 06_master_summary_report.py
Description: 
    Aggregates metrics from individual pipeline outputs to create a comprehensive 
    master summary table containing diversity, capacity, and collision statistics.
    Perfect for Supplementary Tables in publications.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Biological order mapping for elegant table sorting
STAGE_ORDER_MAP = {'plasmid': 1, 'p1': 2, 'enva': 3, 'b19g': 4}

def main():
    print("\n📊 Compiling Master Summary Report...")
    
    # Dynamic path resolution
    try:
        ROOT = Path(__file__).resolve().parents[2]
    except NameError:
        ROOT = Path.cwd().parents[1]
        
    qc_dir = ROOT / "data" / "processed" / "invitro_qc"
    metrics_dir = ROOT / "data" / "processed" / "invitro_metrics"
    results_dir = ROOT / "data" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    meta_path = qc_dir / "Normalized_Metadata.csv"
    lev_path = metrics_dir / "fig3_levenshtein.csv"
    
    if not meta_path.exists():
        print("❌ Error: Normalized_Metadata.csv not found.")
        return
        
    df_meta = pd.read_csv(meta_path)
    
    # Try to load Levenshtein data if it exists
    df_lev = pd.read_csv(lev_path) if lev_path.exists() else None
        
    summary_records = []
    
    for idx, row in df_meta.iterrows():
        group = row['Group']
        stage = row['Stage']
        replicate = row.get('Replicate', f'Rep{idx}')
        
        file_path = qc_dir / row['Normalized_File']
        if not file_path.exists():
            continue
            
        # 1. Read normalized counts
        df_counts = pd.read_csv(file_path, sep='\t')
        if 'UMI_Count' in df_counts.columns:
            counts = df_counts['UMI_Count'].values
        else:
            counts = df_counts.iloc[:, -1].values
            
        counts = counts[counts > 0]
        if len(counts) == 0:
            continue
            
        raw_unique = len(counts)
        total_umi = np.sum(counts)
        
        # 2. Calculate Shannon & Effective BCs (E = 2^H, bits)
        probs = counts / total_umi
        shannon = -np.sum(probs * np.log2(probs))  # log2 → bits
        effective = 2 ** shannon 
        
        # 3. Extract Levenshtein & Collision Rate
        mean_lev = np.nan
        col_rate = np.nan
        
        if df_lev is not None:
            sub_lev = df_lev[(df_lev['Group'] == group) & (df_lev['Stage'] == stage)]
            if not sub_lev.empty:
                dists = sub_lev['Levenshtein_Distance'].values
                mean_lev = np.mean(dists)
                col_rate = (np.sum(dists == 0) / len(dists)) * 100
                
        # Append to records
        summary_records.append({
            'Group': group,
            'Stage': stage,
            'Replicate': replicate,
            'Raw_Unique_BCs': raw_unique,
            'Total_UMI_Weight': int(total_umi),
            'Shannon_Index': round(shannon, 3),
            'Effective_BCs': round(effective, 1),
            'Mean_Levenshtein': round(mean_lev, 2) if not np.isnan(mean_lev) else "N/A",
            'Collision_Rate(%)': round(col_rate, 3) if not np.isnan(col_rate) else "N/A"
        })
        
    if summary_records:
        df_summary = pd.DataFrame(summary_records)
        
        # Sort biologically and logically
        df_summary['Stage_Order'] = df_summary['Stage'].str.lower().map(STAGE_ORDER_MAP).fillna(99)
        df_summary = df_summary.sort_values(['Group', 'Stage_Order', 'Replicate']).drop(columns=['Stage_Order'])
        
        out_file = results_dir / "Master_Library_Summary_Table.csv"
        df_summary.to_csv(out_file, index=False)
        print(f"🎉 Master summary table generated successfully: {out_file}")
    else:
        print("⚠️ No data was processed for the summary table.")

if __name__ == "__main__":
    main()