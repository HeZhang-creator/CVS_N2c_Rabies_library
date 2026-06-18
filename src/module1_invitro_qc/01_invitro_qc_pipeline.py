"""
Script: 01_invitro_qc_pipeline.py
Description: 
    A unified quality control pipeline for in vitro bulk sequencing data (e.g., Plasmid, EnvA).
    
    Step 1: Greedy Directional Adjacency Network Clustering (1bp Error Correction).
            Collapses PCR and sequencing errors by absorbing 1bp variants into dominant core clones.
    Step 2: Stage-specific Multinomial Downsampling.
            Normalizes sequencing depth across all replicates within the same biological stage 
            to eliminate sequencing 'floor effects' when comparing library diversity.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time

# 🌟 设置全局随机种子，确保结果可复现，且消除局部调用的“0和5”聚类假象
np.random.seed(42)

ROOT = Path(__file__).resolve().parents[2]

# ==========================================
# 🧬 Core algorithm: 1bp neighbor generator
# ==========================================
def generate_1bp_neighbors(seq):
    bases = ['A', 'C', 'G', 'T', 'N']
    neighbors = []
    # Substitutions
    for i in range(len(seq)):
        for b in bases:
            if b != seq[i]: neighbors.append(seq[:i] + b + seq[i+1:])
    # Deletions
    for i in range(len(seq)):
        neighbors.append(seq[:i] + seq[i+1:])
    # Insertions
    for i in range(len(seq) + 1):
        for b in bases:
            neighbors.append(seq[:i] + b + seq[i:])
    return neighbors

# ==========================================
# 🧬 Core algorithm: greedy 1bp network correction
# ==========================================
def error_correction_1bp(barcode_counts):
    """
    Input: dict of {barcode: count}
    Output: dict of {corrected_barcode: aggregated_count}
    Logic: Sort by abundance descending. Core clones absorb 1bp mutated neighbors.
    """
    sorted_bcs = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)
    corrected_dict = {}
    processed = set()
    
    for bc, count in sorted_bcs:
        if bc in processed:
            continue
            
        core_count = count
        processed.add(bc)
        neighbors = generate_1bp_neighbors(bc)
        
        for neighbor in neighbors:
            if neighbor in barcode_counts and neighbor not in processed:
                core_count += barcode_counts[neighbor]
                processed.add(neighbor)
                
        corrected_dict[bc] = core_count
        
    return corrected_dict

# ==========================================
# 🧬 Core algorithm: Multinomial Downsampling
# ==========================================
def downsample_reads(counts, target_depth):
    """
    使用多项式分布将 read counts 下采样到目标深度。
    利用了全局 np.random.seed，保证科学天然性与结果可复现。
    """
    probabilities = counts / np.sum(counts)
    downsampled_counts = np.random.multinomial(target_depth, probabilities)
    return downsampled_counts

# ==========================================
# 🚀 Main Pipeline Execution
# ==========================================
def main():
    print("🚀 Starting In Vitro sequencing data cleaning and error correction pipeline...")
    start_time = time.time()
    
    raw_path_dir = ROOT / "data" / "raw_in_vitro"
    out_path_dir = ROOT / "data" / "processed" / "invitro_qc"
    meta_path = ROOT / "data" / "metadata" / "RV_library_metadata.csv"
    
    out_path_dir.mkdir(parents=True, exist_ok=True)
    
    if not meta_path.exists():
        print(f"❌ Metadata file not found: {meta_path}")
        return
        
    meta_df = pd.read_csv(meta_path)
    
    # Clean whitespace in metadata
    meta_df['Group'] = meta_df['Group'].astype(str).str.strip()
    meta_df['Stage'] = meta_df['Stage'].astype(str).str.strip()

    # Data structure to store intermediate results
    processed_data = []

    print("\n[Phase 1/2] Executing 1bp network correction by abundance...")
    for index, row in meta_df.iterrows():
        
        # --- SMART PATH RESOLUTION START ---
        file_name = str(row.get('File_Name', '')).strip()
        
        if 'Full_Path' in row and pd.notna(row['Full_Path']) and str(row['Full_Path']).strip() != "":
            raw_path_str = str(row['Full_Path']).strip()
            # If the CSV has an absolute path (C:\...), use it. Otherwise, append to ROOT
            if Path(raw_path_str).is_absolute():
                file_path = Path(raw_path_str)
            else:
                file_path = ROOT / raw_path_str
        else:
            file_path = raw_path_dir / file_name
        # --- SMART PATH RESOLUTION END ---
        
        group = row.get('Group', 'Unknown')
        stage = row.get('Stage', 'Unknown')
        
        # Safely handle 'Replicate' if it's missing from your CSV
        replicate = str(row.get('Replicate', 'Rep1')).strip()
        
        if not file_path.exists():
            print(f"⚠️ ERROR: Cannot find raw data file at -> {file_path}")
            continue
        if not file_name:
            continue
            
        # Parse Raw Count File
        bc_counts = {}
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 5 or parts[0].lower() == 'count': continue
                try:
                    bc = parts[3]
                    umi = float(parts[4])
                    if umi > 0 and set(bc.upper()).issubset({'A', 'T', 'C', 'G', 'N'}): 
                        bc_counts[bc] = bc_counts.get(bc, 0) + umi
                except Exception:
                    continue
                    
        total_umi_before = sum(bc_counts.values())
        print(f"  👉 {len(bc_counts)} Unique BCs, preparing for 1bp correction... [{file_name}]")
        
        # Apply Correction
        corrected_counts = error_correction_1bp(bc_counts)
        
        processed_data.append({
            'metadata_row': row,
            'corrected_dict': corrected_counts,
            'total_umi': total_umi_before,
            'stage_key': f"{group}_{stage}",
            'replicate': replicate
        })

    print("\n[Phase 2/2] Executing Stage-specific multinomial downsampling...")
    # Calculate minimum UMI depth for each stage
    stage_min_umi = {}
    for data in processed_data:
        stage = data['metadata_row']['Stage']
        current_total = sum(data['corrected_dict'].values())
        if stage not in stage_min_umi:
            stage_min_umi[stage] = current_total
        else:
            stage_min_umi[stage] = min(stage_min_umi[stage], current_total)
            
    new_metadata_records = []

    for data in processed_data:
        row = data['metadata_row']
        bc_dict = data['corrected_dict']
        stage_label = row['Stage']
        stage_key = data['stage_key']
        group = row['Group']
        replicate = data['replicate']
        
        if stage_label not in stage_min_umi:
            print(f"  ⚠️ Not enough data in stage {stage_label} for downsampling.")
            continue
            
        target_depth = int(stage_min_umi[stage_label])
        
        bcs = list(bc_dict.keys())
        counts = np.array(list(bc_dict.values()))
        
        # 🌟 调用外层的下采样函数 (直接传入 counts 和 target_depth)
        downsampled_counts = downsample_reads(counts, target_depth)
        
        final_bcs = 0
        output_file = out_path_dir / f"{Path(row.get('File_Name', 'unknown')).stem}_normalized.txt"
            
        with open(output_file, 'w', encoding='utf-8') as fout:
            fout.write("Count\tPercent\tLength\tBarcode\tUMI_Count\n")
            for i, bc in enumerate(bcs):
                if downsampled_counts[i] > 0:
                    final_bcs += 1
                    percent = downsampled_counts[i] / target_depth * 100
                    fout.write(f"{downsampled_counts[i]}\t{percent:.4f}\t{len(bc)}\t{bc}\t{downsampled_counts[i]}\n")
                    
        print(f"  ✅ [{group}] {row.get('File_Name', 'unknown')} -> Retained {final_bcs:,} unique BCs. (Target: {target_depth:,})")
            
        new_row = row.copy()
        new_row['Original_File'] = str(row.get('File_Name', '')).strip()
        new_row['Normalized_File'] = output_file.name
        new_row['Stage_Key'] = stage_key
        new_row['Stage'] = stage_label
        new_row['Replicate'] = replicate
        new_row['Total_UMI_before'] = int(data['total_umi'])
        new_row['Downsampled_UMI'] = int(downsampled_counts.sum())
        new_row['Target_Depth'] = int(target_depth)
        new_metadata_records.append(new_row)
            
    # Save updated Metadata
    if new_metadata_records:
        pd.DataFrame(new_metadata_records).to_csv(out_path_dir / "Normalized_Metadata.csv", index=False)
    print(f"\n🎉 Invitro QC completed in {time.time()-start_time:.2f} seconds!")

if __name__ == "__main__":
    main()