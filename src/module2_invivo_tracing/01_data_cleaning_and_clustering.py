"""
Script: 01_data_cleaning_and_clustering.py
Description: 
    Executes a conservative "anchor-first" preprocessing pipeline for in vivo 
    trans-synaptic scRNA-seq data to prevent artificial noise daisy-chaining.
    
    Phase 1: Dual-Threshold Noise Filtering.
             - Absolute threshold: UMI_Kinds >= 3 (Mitigates template switching).
             - Relative threshold: Intra-cellular fraction >= 5% (Eliminates ambient RNA).
    Phase 2: Clonal Consolidation (1bp Levenshtein Clustering).
             - Performs top-down greedy merging exclusively on the surviving 
               high-confidence "anchor" sequences to resolve discrete point mutations.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time
from tqdm import tqdm

# ==========================================
# 🧬 Fast 1bp Levenshtein Neighbor Generator
# ==========================================
def generate_1bp_neighbors(seq):
    """Generates all theoretical sequences within 1bp Levenshtein distance in O(N) time."""
    bases = ['A', 'C', 'G', 'T', 'N']
    neighbors = set()
    # Substitutions
    for i in range(len(seq)):
        for b in bases:
            if b != seq[i]: neighbors.add(seq[:i] + b + seq[i+1:])
    # Deletions
    for i in range(len(seq)):
        neighbors.add(seq[:i] + seq[i+1:])
    # Insertions
    for i in range(len(seq) + 1):
        for b in bases:
            neighbors.add(seq[:i] + b + seq[i:])
    return neighbors

def perform_anchor_clustering(sub_df, barcode_col='RV_Barcode', count_col='UMI_Kinds'):
    """Performs greedy clustering exclusively on high-confidence anchor barcodes."""
    bc_counts = sub_df.groupby(barcode_col)[count_col].sum().sort_values(ascending=False)
    sorted_barcodes = bc_counts.index.tolist()
    
    mapping = {}
    processed = set()
    
    for core_bc in sorted_barcodes:
        if core_bc in processed: 
            continue
        # Assign core barcode to itself
        mapping[core_bc] = core_bc
        processed.add(core_bc)
        
        # Search for and absorb 1bp variants that survived the noise filters
        for neighbor in generate_1bp_neighbors(core_bc):
            if neighbor in bc_counts and neighbor not in processed:
                mapping[neighbor] = core_bc
                processed.add(neighbor)
                
    return mapping

# ==========================================
# 🚀 Main Pipeline Execution
# ==========================================
def run_anchor_first_pipeline(input_csv, output_csv, umi_thresh=3, frac_thresh=0.05):
    start_time = time.time()
    print(f"📥 [Step 1] Loading raw dataset: {input_csv}")
    
    input_file = Path(input_csv)
    if not input_file.exists():
        raise FileNotFoundError(f"❌ Raw dataset not found: {input_file}")

    df_raw = pd.read_csv(input_file)
    total_initial = len(df_raw)
    print(f"📊 Initial raw connections: {total_initial:,}")

    # ---------------------------------------------------------
    # Phase 1: Stringent Noise Filtering (Anchor Selection)
    # ---------------------------------------------------------
    print(f"\n🛡️ [Phase 1] Selecting High-Confidence Anchors...")
    
    # Filter 1: UMI >= 3
    df_pass_1 = df_raw[df_raw['UMI_Kinds'] >= umi_thresh].copy()
    drop_1 = total_initial - len(df_pass_1)
    print(f"   - Threshold 1 (UMI >= {umi_thresh}): Removed {drop_1:,} singleton/doubleton artifacts.")
    
    # Filter 2: Fraction >= 5%
    cell_group_cols = ['Group', 'Brain_Region', 'Replicate', 'Cell_Barcode']
    
    # Safely compute total UMIs per host cell
    cell_totals = df_pass_1.groupby(cell_group_cols)['UMI_Kinds'].transform('sum')
    df_pass_1['Cell_Internal_Fraction'] = df_pass_1['UMI_Kinds'] / cell_totals
    
    # Retain only the solid anchors
    df_anchors = df_pass_1[df_pass_1['Cell_Internal_Fraction'] >= frac_thresh].copy()
    drop_2 = len(df_pass_1) - len(df_anchors)
    print(f"   - Threshold 2 (Fraction >= {frac_thresh*100}%): Removed {drop_2:,} ambient noise events.")
    print(f"   ✅ Solid anchors established: {len(df_anchors):,} events.")

    # ---------------------------------------------------------
    # Phase 2: Clonal Consolidation (1bp Clustering)
    # ---------------------------------------------------------
    print("\n🧬 [Phase 2] Executing 1bp Clonal Consolidation on Anchors...")
    
    # Isolate clustering boundaries to prevent artificial merging across subjects
    biol_boundaries = ['Group', 'Brain_Region', 'Replicate'] 
    
    corrected_dfs = []
    for (group, region, rep), sub in tqdm(df_anchors.groupby(biol_boundaries), desc="Clustering Replicates"):
        mapping_dict = perform_anchor_clustering(sub)
        
        sub = sub.copy()
        # Map to corrected barcode, defaulting to original if unmapped
        sub['RV_Barcode_Corrected'] = sub['RV_Barcode'].map(mapping_dict).fillna(sub['RV_Barcode'])
        corrected_dfs.append(sub)

    df_corrected = pd.concat(corrected_dfs, ignore_index=True)
    
    # Consolidate UMI counts for variants merged within the same host cell
    merge_cols = ['Group', 'Brain_Region', 'Replicate', 'Cell_Barcode', 'RV_Barcode_Corrected']
    df_final = df_corrected.groupby(merge_cols, as_index=False)['UMI_Kinds'].sum()
    
    # ---------------------------------------------------------
    # Export
    # ---------------------------------------------------------
    print("\n✅ Pipeline Completed Successfully!")
    print(f"🏆 Final consolidated biological connections: {len(df_final):,} ({(len(df_final)/total_initial)*100:.2f}% retained)")
    
    output_file = Path(output_csv)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df_final.to_csv(output_file, index=False)
    
    print(f"💾 Master Corrected Data saved to: {output_file}")
    print(f"⏱️ Total execution time: {time.time()-start_time:.2f} seconds.")

if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[2]
    INPUT_PATH = ROOT / "data" / "raw_in_vivo" / "Master_Raw_Data.csv"
    OUTPUT_PATH = ROOT / "data" / "processed" / "invivo_metrics" / "Master_Corrected_Data.csv"
    
    # Execute the Anchor-First pipeline
    run_anchor_first_pipeline(INPUT_PATH, OUTPUT_PATH, umi_thresh=3, frac_thresh=0.05)