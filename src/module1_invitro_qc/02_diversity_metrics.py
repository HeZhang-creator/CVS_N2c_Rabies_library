"""
Script: 02_diversity_metrics.py
Description:
    Reads normalized in vitro QC outputs from 01_invitro_qc_pipeline.py and computes
    library diversity metrics for Figure 3.

    Outputs:
    - fig3_freq_dist.csv: rank-abundance and cumulative frequency data.
    - fig3_effective_capacity.csv: Shannon entropy / effective barcode capacity.
    - fig3_levenshtein.csv: sampled pairwise Levenshtein distances for stage/group combinations.
"""

import pandas as pd
import numpy as np
import random
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# ==========================================
# Core utility: Levenshtein edit distance
# ==========================================
def levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

# ==========================================
# Core utility: Effective capacity
# ==========================================
def calculate_effective_capacity(counts):
    """
    Computes Shannon entropy (bits) and effective barcode capacity (E = 2^H).
    Returns (effective_capacity, shannon_bits).
    """
    total = np.sum(counts)
    if total <= 0:
        return 0, 0
    probs = counts / total
    probs = probs[probs > 0]
    shannon_bits = -np.sum(probs * np.log2(probs))  # log2 → bits
    effective_capacity = 2 ** shannon_bits
    return effective_capacity, shannon_bits

# ==========================================
# Core utility: Sampled Levenshtein computation
# ==========================================
def sample_levenshtein_distances(barcodes, counts, sample_pairs=1000, sample_size=5000, seed=42):
    np.random.seed(seed)
    total_counts = np.sum(counts)
    if total_counts <= 0 or len(barcodes) < 2:
        return []
        
    probs = counts / total_counts
    actual_sample_size = min(sample_size, len(barcodes))
    
    # Sample based on abundance weight to reflect jackpot clones
    sampled_bcs = np.random.choice(barcodes, size=actual_sample_size, p=probs, replace=True)
    
    distances = []
    for _ in range(sample_pairs):
        bc1, bc2 = np.random.choice(sampled_bcs, 2, replace=False)
        distances.append(levenshtein_distance(bc1, bc2))
        
    return distances

def main():
    print("🚀 Starting 02_diversity_metrics.py: Library diversity calculation...")
    
    input_dir = ROOT / "data" / "processed" / "invitro_qc"
    out_dir = ROOT / "data" / "processed" / "invitro_metrics"
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_path = input_dir / "Normalized_Metadata.csv"

    if not meta_path.exists():
        print(f"❌ Metadata file not found: {meta_path}. Please run step 01 first.")
        return

    meta_df = pd.read_csv(meta_path)
    
    effective_capacity_rows = []
    freq_rows = []
    distance_rows = []

    # Sampling parameters
    sample_pairs = 1000
    sample_limit = 5000

    for idx, row in meta_df.iterrows():
        stage = row['Stage']
        group = row['Group']
        replicate = row['Replicate']
        file_name = row['Normalized_File']
        
        file_path = input_dir / file_name
        if not file_path.exists():
            continue
            
        print(f"  ⏳ Analyzing: {group} | {stage} | {replicate}")
        
        df = pd.read_csv(file_path, sep='\t')
        if df.empty or 'UMI_Count' not in df.columns:
            continue
            
        barcode_counter = dict(zip(df['Barcode'], df['UMI_Count']))
        counts = np.array(list(barcode_counter.values()))
        
        eff_cap, shannon_bits = calculate_effective_capacity(counts)
        effective_capacity_rows.append({
            'Stage': stage, 
            'Group': group, 
            'Replicate': replicate,
            'Shannon_Index_bits': round(shannon_bits, 3),
            'Effective_BCs': round(eff_cap, 1)
        })

        counts_sorted = np.sort(counts)[::-1]
        total_counts = np.sum(counts_sorted)
        frequencies = (counts_sorted / total_counts) * 100
        cumulative = np.cumsum(frequencies)
        
        # Limit to top 10000 ranks to reduce file size
        for rank, (freq, cum) in enumerate(zip(frequencies[:10000], cumulative[:10000]), start=1):
            freq_rows.append({
                'Stage': stage,
                'Group': group,
                'Rank': rank,
                'Frequency': freq,
                'Cumulative_Freq': cum
            })

        distances = sample_levenshtein_distances(
            np.array(list(barcode_counter.keys()), dtype=object),
            np.array(list(barcode_counter.values()), dtype=float),
            sample_pairs=sample_pairs,
            sample_size=sample_limit,
            seed=42
        )
        for dist in distances:
            distance_rows.append({
                'Stage': stage,
                'Group': group,
                'Levenshtein_Distance': dist
            })

    pd.DataFrame(effective_capacity_rows).to_csv(out_dir / "fig3_effective_capacity.csv", index=False)
    pd.DataFrame(freq_rows).to_csv(out_dir / "fig3_freq_dist.csv", index=False)
    pd.DataFrame(distance_rows).to_csv(out_dir / "fig3_levenshtein.csv", index=False)

    print(f"\n🎉 Perfect! Calculation complete. Three core data tables have been generated in {out_dir}.")
    print("👉 Next step: You can directly run the plotting scripts to render these tables into charts!")

if __name__ == "__main__":
    main()