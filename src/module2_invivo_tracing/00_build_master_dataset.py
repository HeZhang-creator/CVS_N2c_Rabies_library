"""
Script: 00_build_master_dataset.py
Description: 
    Consolidates raw UMI count matrices from individual sequencing runs 
    into a single master dataset, driven by a metadata manifest file.
    This ensures full data provenance and reproducibility across operating systems.
"""

import pandas as pd
from pathlib import Path
import time

def consolidate_raw_data(metadata_csv_path, output_csv_path):
    start_time = time.time()
    print(f"🚀 [Step 0] Loading metadata manifest from: {metadata_csv_path}")
    
    metadata_file = Path(metadata_csv_path)
    if not metadata_file.exists():
        raise FileNotFoundError(f"❌ Manifest not found: {metadata_file}. Please check directory structure.")

    meta_df = pd.read_csv(metadata_file, encoding='utf-8-sig', sep=None, engine='python')
    meta_df.columns = meta_df.columns.str.strip()
    print(f"📊 DEBUG: 成功读取的表头是: {meta_df.columns.tolist()}")
    
    if 'Sample_ID' not in meta_df.columns:
        raise KeyError(f"❌ 关键列 'Sample_ID' 未找到！请检查 CSV 分隔符。当前识别出的表头为: {meta_df.columns.tolist()}")
    master_results_list = []
    
    print(f"📦 Total biological samples to process: {len(meta_df)}")
    
    # Dynamically locate the 'raw' data directory relative to the metadata file
    # (Assuming metadata is in data/metadata/ and raw data is in data/raw_in_vivo/)
    raw_data_dir = metadata_file.parent.parent / 'raw_in_vivo'
    
    for index, row in meta_df.iterrows():
        # 🌟 修改 1: 完美对接您最新的 Metadata 表头 (Sample_ID)
        sample_name = str(row['Sample_ID']).strip()
        group = str(row['Group']).strip()
        region = str(row['Brain_Region']).strip()
        replicate = str(row['Replicate']).strip()
        
        # 🌟 修改 2: 直接构建扁平化的独立文件路径，摒弃嵌套文件夹
        # 假设您的文件后缀为 .txt。如果是 .csv，请在此处修改。
        raw_file_path = raw_data_dir / f"{sample_name}.txt"
        
        print(f"  ▶️ Processing [{index+1:02d}/{len(meta_df)}]: {sample_name} ({group} | {region} | Rep{replicate})")
        
        if not raw_file_path.exists():
            print(f"  ❌ WARNING: File missing at {raw_file_path}, skipping this sample.")
            continue
            
        try:
            # Load the raw text file (assuming tab-separated based on the .txt extension)
            # Note: Change sep='\\t' to sep=',' if your txt file is actually comma-separated.
            temp_df = pd.read_csv(raw_file_path, sep='\t')
            
            # Tag every row with its biological provenance (crucial for downstream groupby operations)
            temp_df['Sample_ID'] = sample_name
            temp_df['Group'] = group
            temp_df['Brain_Region'] = region
            temp_df['Replicate'] = int(replicate)
            
            master_results_list.append(temp_df)
        except Exception as e:
            print(f"  ❌ ERROR reading {sample_name}: {e}")

    # Consolidate and Export
    if master_results_list:
        final_master_df = pd.concat(master_results_list, ignore_index=True)
        print(f"\n🎉 Consolidation complete! Total trans-synaptic connections recorded: {len(final_master_df):,}")
        
        output_file = Path(output_csv_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        final_master_df.to_csv(output_file, index=False)
        print(f"💾 Master Raw Dataset saved to: {output_file}")
        print(f"⏱️ Execution time: {time.time()-start_time:.2f} seconds.")
        return final_master_df
    else:
        print("\n⚠️ No valid data consolidated. Please check your raw files.")
        return None

if __name__ == "__main__":
    # Define paths relative to this script's location
    # Assuming script is run from src/module2_invivo_tracing/
    
    MANIFEST_PATH = "../../data/metadata/Master_Sample_Metadata.csv"
    OUTPUT_PATH = "../../data/raw_in_vivo/Master_Raw_Data.csv"
    
    consolidate_raw_data(MANIFEST_PATH, OUTPUT_PATH)