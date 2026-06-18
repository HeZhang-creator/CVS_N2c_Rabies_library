"""
Script: 04_collision_resolution.py
Description: 
    Executes the Dual-Resolution Collision Analysis.
    - Resolution 1 (Barcode-Centric): Evaluates structural purity prior to amplification.
    - Resolution 2 (Cell-Centric): Evaluates biological footprint and host cell collision status.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time

def compute_dual_resolution_collision(input_csv, output_prefix):
    start_time = time.time()
    print(f"🚀 Loading dataset: {input_csv}")
    
    df = pd.read_csv(input_csv)
    barcode_col = 'RV_Barcode_Corrected' if 'RV_Barcode_Corrected' in df.columns else 'RV_Barcode'
    cell_col = 'Cell_Barcode'
    
    groups = ['HH + HdVRz', 'HH + Twister', 'RzB + Twister']
    regions = ['V1_Cortex', 'LGN']
    
    barcode_metrics = []
    cell_metrics = []
    
    print("🧠 Crunching Barcode-centric and Cell-centric collision profiles...")
    
    for group in groups:
        group_df = df[df['Group'] == group]
        for rep in group_df['Replicate'].unique():
            rep_df = group_df[group_df['Replicate'] == rep]
            
            for region in regions:
                region_df = rep_df[rep_df['Brain_Region'] == region]
                if region_df.empty: continue
                
                # ==========================================
                # 1. Barcode-level statistics (count how many cells each RV-Barcode infected)
                # ==========================================
                bc_infection_counts = region_df.groupby(barcode_col)[cell_col].nunique()
                
                # Categorize: 1 (Unique), 2 (Low), 3 (Med), >3 (High)
                bc_cat = pd.cut(
                    bc_infection_counts, 
                    bins=[0, 1, 2, 3, np.inf], 
                    labels=['1 (Unique)', '2 (Low)', '3 (Med)', '>3 (High)']
                )
                bc_distribution = bc_cat.value_counts(normalize=True) * 100
                
                for cat, pct in bc_distribution.items():
                    barcode_metrics.append({
                        'Group': group, 'Replicate': rep, 'Brain_Region': region,
                        'Collision_Level': cat, 'Proportion': pct
                    })
                
                # ==========================================
                # 2. Cell-level statistics (assess host cell collision status)
                # ==========================================
                # Map each RV-Barcode's collision count back to each row
                region_df = region_df.copy()
                region_df['BC_Infection_Count'] = region_df[barcode_col].map(bc_infection_counts)
                
                # Evaluate each host cell: if any RV-Barcode in it infected other cells (Count>1),
                # the cell is Collision-affected. It is only Collision-free when all RV-Barcodes
                # in that cell are globally unique to it (Count==1).
                cell_max_conflict = region_df.groupby(cell_col)['BC_Infection_Count'].max()
                
                total_cells = len(cell_max_conflict)
                collision_free_cells = (cell_max_conflict == 1).sum()
                collision_affected_cells = (cell_max_conflict > 1).sum()
                
                collision_free_pct = (collision_free_cells / total_cells) * 100 if total_cells > 0 else 0
                collision_affected_pct = (collision_affected_cells / total_cells) * 100 if total_cells > 0 else 0
                
                cell_metrics.append({
                    'Group': group, 'Replicate': rep, 'Brain_Region': region,
                    'Total_Cells': total_cells,
                    'Collision_Free_Pct': collision_free_pct,
                    'Collision_Affected_Pct': collision_affected_pct
                })

    df_barcode = pd.DataFrame(barcode_metrics)
    df_cell = pd.DataFrame(cell_metrics)
    
    out_dir = Path(output_prefix).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    
    df_barcode.to_csv(f"{output_prefix}_Barcode_Centric.csv", index=False)
    df_cell.to_csv(f"{output_prefix}_Cell_Centric.csv", index=False)
    
    print(f"🎉 Analysis complete in {time.time()-start_time:.2f}s!")
    return df_barcode, df_cell

if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[2]
    INPUT = ROOT / "data" / "processed" / "invivo_metrics" / "Master_Corrected_Data.csv"
    OUTPUT_PREFIX = ROOT / "data" / "processed" / "invivo_metrics" / "Fig5_Resolution_Metrics"
    compute_dual_resolution_collision(INPUT, OUTPUT_PREFIX)