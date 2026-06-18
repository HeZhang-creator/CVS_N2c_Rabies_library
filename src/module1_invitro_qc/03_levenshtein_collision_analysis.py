"""
Script: 03_levenshtein_collision_analysis.py
Description: 
    Data processing engine for Sequence Space Analysis (Collision/Jackpot effect).
    Calculates abundance-weighted pairwise Levenshtein distances using Monte Carlo sampling.
    *Note: Plotting is delegated to Jupyter Notebooks for visual fine-tuning.*
"""

import pandas as pd
import numpy as np
import itertools
from Levenshtein import distance as lev_dist
import os
from pathlib import Path

# ================= 1. Core: Weighted sampling and distance calculation =================
def calculate_weighted_levenshtein(file_path, sample_size=1500, random_seed=42):
    """Read normalized data, sample based on UMI abundance, and calculate pairwise edit distances"""
    if not os.path.exists(file_path):
        return []

    # Fix random seed to ensure absolute reproducibility
    np.random.seed(random_seed) 
    
    # Read normalized, clean data from script 01
    df = pd.read_csv(file_path, sep='\t')
    
    if df.empty or 'UMI_Count' not in df.columns:
        return []
    
    # Filter valid Barcodes (pure ATCGN)
    df = df[df['Barcode'].str.match(r'^[ATCGN]+$', case=False, na=False)]
    
    barcodes = df['Barcode'].values
    counts = df['UMI_Count'].values
    total_umis = np.sum(counts)
    
    if total_umis == 0: return []

    # 1. Calculate abundance weight probabilities
    probabilities = counts / total_umis
    
    # 2. Monte Carlo sampling (simulate real physical collision probabilities)
    sampled_bcs = np.random.choice(barcodes, size=sample_size, p=probabilities)
    
    # 3. Calculate pairwise Levenshtein distances
    distances = [lev_dist(bc1, bc2) for bc1, bc2 in itertools.combinations(sampled_bcs, 2)]
    
    return distances

# ================= 2. Main Execution Flow =================
def main():
    print("🚀 Starting Levenshtein Distance collision analysis engine...")
    
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
    
    # Clean potential whitespaces
    meta_df['Group'] = meta_df['Group'].astype(str).str.strip()
    meta_df['Stage'] = meta_df['Stage'].astype(str).str.strip()
    
    all_distances = []

    for _, row in meta_df.iterrows():
        group = row['Group']
        stage = row['Stage']
        file_name = str(row.get('Normalized_File', ''))
        
        file_path = data_dir / file_name
        
        if not file_path.exists() or not file_name:
            continue
            
        print(f"  ⏳ Sampling and calculating: [{group}] - [{stage}]")
        
        # To ensure statistical power and speed, sample 1500 molecules by default (generates ~1.12 million distance pairs)
        distances = calculate_weighted_levenshtein(file_path, sample_size=1500)
        
        if distances:
            # Flatten distances and append to list
            for dist in distances:
                all_distances.append({
                    'Group': group,
                    'Stage': stage,
                    'Levenshtein_Distance': dist
                })

    # ================= 3. Export final table for Notebook rendering =================
    if all_distances:
        print("\n💾 Writing millions of distance records to CSV (this may take a few seconds)...")
        df_out = pd.DataFrame(all_distances)
        df_out.to_csv(out_dir / "levenshtein_distances.csv", index=False)
        print(f"🎉 Calculation complete! Data exported to: {out_dir / 'levenshtein_distances.csv'}")
        print("👉 Next step: Open the plotting script in notebooks/ for elegant visualization.")

if __name__ == "__main__":
    main()