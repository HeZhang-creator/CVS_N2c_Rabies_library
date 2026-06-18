"""
Script: 03_rarefaction.py
Description: 
    Performs rarefaction analysis on cleaned in vivo trans-synaptic data.
    
    This script evaluates library saturation by simulating random subsampling 
    (without replacement) of valid trans-synaptic connections. It calculates 
    the expected number of unique RV-Barcodes discovered at increasing sampling depths.
    Steeper, non-asymptotic curves indicate high library complexity and severe 
    under-sequencing.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time

def calculate_rarefaction_curve(data, target_col='RV_Barcode', step_size=500, random_seed=42):
    """
    Simulates random sampling to generate rarefaction coordinates.
    """
    # Randomize the data order to simulate random sampling
    shuffled_data = data.sample(frac=1, random_state=random_seed).reset_index(drop=True)
    total_records = len(shuffled_data)
    
    if total_records == 0:
        return [], []

    # Generate sampling points (e.g., 0, 500, 1000...)
    sample_indices = np.arange(0, total_records, step_size)
    if len(sample_indices) == 0 or sample_indices[-1] != total_records - 1:
        sample_indices = np.append(sample_indices, total_records - 1)
        
    unique_counts = []
    seen_barcodes = set()
    barcode_values = shuffled_data[target_col].values
    current_pos = 0
    
    # Cumulatively count the number of unique barcodes discovered at each sampling point
    for idx in sample_indices:
        seen_barcodes.update(barcode_values[current_pos : idx + 1])
        unique_counts.append(len(seen_barcodes))
        current_pos = idx + 1
        
    return sample_indices, unique_counts

def generate_rarefaction_data(cleaned_csv_path, output_csv_path):
    print(f"🚀 [Step 4] Loading cleaned trans-synaptic data: {cleaned_csv_path}")
    
    input_file = Path(cleaned_csv_path)
    if not input_file.exists():
        raise FileNotFoundError(f"❌ Cleaned dataset not found: {input_file}.")

    df_clean = pd.read_csv(input_file)
    print("\n📈 Simulating rarefaction curves (this may take a moment)...")
    
    start_time = time.time()
    rarefaction_results = []
    
    # Ensure required columns are present
    barcode_col = 'RV_Barcode_Corrected' if 'RV_Barcode_Corrected' in df_clean.columns else 'RV_Barcode'
    required_cols = ['Group', 'Brain_Region', 'Replicate', barcode_col]
    for col in required_cols:
        if col not in df_clean.columns:
            raise KeyError(f"❌ Missing critical column: {col}")

    # Group by biological replicate for simulation
    grouped = df_clean.groupby(['Group', 'Brain_Region', 'Replicate'])
    
    for (group, brain_region, Replicate), sample_df in grouped:
        print(f"  ▶️ Processing: [{group} | {brain_region}] {Replicate} (n={len(sample_df)})")
        
        x_indices, y_uniques = calculate_rarefaction_curve(sample_df, target_col=barcode_col)
        
        # Serialize the calculated coordinates for export
        for step, unique_count in zip(x_indices, y_uniques):
            rarefaction_results.append({
                'Group': group,
                'Brain_Region': brain_region,
                'Replicate': Replicate,
                'Sampled_Pairs': int(step + 1), # Actual sampled count is index + 1
                'Unique_Barcodes': int(unique_count)
            })

    # Export the results
    df_results = pd.DataFrame(rarefaction_results)
    
    output_file = Path(output_csv_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(output_file, index=False)
    
    print(f"\n🎉 Rarefaction simulation completed in {time.time()-start_time:.2f} seconds!")
    print(f"💾 Coordinates saved to: {output_file}")
    
    return df_results

if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[2]
    INPUT_PATH = ROOT / "data" / "processed" / "invivo_metrics" / "Master_Corrected_Data.csv"
    OUTPUT_PATH = ROOT / "data" / "processed" / "invivo_metrics" / "Rarefaction_Coordinates.csv"
    
    generate_rarefaction_data(INPUT_PATH, OUTPUT_PATH)