"""
Script: 02_diversity_metrics.py
Description: 
    The centralized calculation engine for all library complexity metrics 
    from in vivo trans-synaptic data.
    
    Metrics Computed:
    - Absolute Diversity: Total unique RV-Barcodes.
    - Total UMI: Overall sequencing depth baseline.
    - Shannon Index & Effective Capacity (e^H): True information carrying capacity.
    - Information Density: Unique barcodes per 10,000 UMIs (Sequencing Efficiency).
"""

import pandas as pd
import numpy as np
from scipy.stats import entropy
from pathlib import Path
import time

def calculate_diversity_metrics(cleaned_csv_path, output_metrics_path):
    print(f"🚀 [Step 2] Loading cleaned data: {cleaned_csv_path}")
    
    input_file = Path(cleaned_csv_path)
    if not input_file.exists():
        raise FileNotFoundError(f"❌ Dataset not found: {input_file}")

    df_clean = pd.read_csv(input_file)
    # Support both possible barcode column names
    barcode_col = 'RV_Barcode_Corrected' if 'RV_Barcode_Corrected' in df_clean.columns else 'RV_Barcode'
    
    print("\n🧬 Calculating comprehensive diversity metrics (Entropy, Capacity, Density)...")
    metrics_list = []
    
    # Assume metadata includes Replicate and Brain_Region columns
    # Group using the most detailed available keys
    grouped = df_clean.groupby(['Group', 'Replicate', 'Brain_Region'])
    
    for (group, rep, region), sample_df in grouped:
        barcode_abundance = sample_df.groupby(barcode_col)['UMI_Kinds'].sum()
        total_umi = barcode_abundance.sum()
        
        if total_umi > 0:
            # 1. Absolute diversity count
            unique_bcs = len(barcode_abundance)
            
            # 2. Shannon entropy and effective capacity (base e)
            probabilities = barcode_abundance / total_umi
            shannon_index = entropy(probabilities, base=np.e)
            effective_capacity = np.exp(shannon_index)
            
            # 3. Information density (unique barcodes per 10,000 UMIs)
            info_density = (unique_bcs / total_umi) * 10000
            
            metrics_list.append({
                'Group': group,
                'Replicate': rep,
                'Brain_Region': region,
                'Total_UMI': total_umi,
                'Absolute_Diversity': unique_bcs,
                'Shannon_Index': round(shannon_index, 4),
                'Effective_Capacity': round(effective_capacity, 2),
                'Information_Density': round(info_density, 2)
            })

    df_metrics = pd.DataFrame(metrics_list)
    df_metrics = df_metrics.sort_values(by=['Group', 'Brain_Region', 'Replicate'])
    
    output_file = Path(output_metrics_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df_metrics.to_csv(output_file, index=False)
    
    print(f"🎉 Unified diversity metrics saved to: {output_file}")
    return df_metrics

if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[2]
    INPUT_PATH = ROOT / "data" / "processed" / "invivo_metrics" / "Master_Corrected_Data.csv"
    OUTPUT_PATH = ROOT / "data" / "processed" / "invivo_metrics" / "Diversity_Metrics_Summary.csv"
    calculate_diversity_metrics(INPUT_PATH, OUTPUT_PATH)