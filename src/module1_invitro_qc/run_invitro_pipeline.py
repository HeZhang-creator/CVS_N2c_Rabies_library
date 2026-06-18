"""
Pipeline runner for in vitro QC and diversity metrics.
This script executes the complete heavy-lifting workflow:
  1) 01_invitro_qc_pipeline: 1bp Error correction & Multinomial downsampling
  2) 02_diversity_metrics: Shannon entropy & Rank abundance calculation
  3) 03_levenshtein_collision_analysis: Monte Carlo sequence space collision
  4) 04_clonal_correlation_analysis: Cross-stage fidelity and Spearman correlation
  
Note: All plotting is handled separately in Jupyter Notebooks to ensure UI flexibility.
"""

import importlib.util
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
MODULE_DIR = ROOT / "src" / "module1_invitro_qc"

def load_module_from_path(module_name, path):
    """Dynamically load modules to avoid relative import errors"""
    if not path.exists():
        raise FileNotFoundError(f"❌ Critical module file not found: {path}")
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def main():
    print("="*65)
    print("🧬 CVS-N2c Library Pipeline: Starting In Vitro Master Workflow")
    print("="*65)
    
    start_time = time.time()

    # ---------------------------------------------------------
    # Step 1: Data Cleaning and Normalized Downsampling
    # ---------------------------------------------------------
    print("\n▶️ [STEP 1/4] Running QC and Normalization Pipeline...")
    qc_module = load_module_from_path("mod_01", MODULE_DIR / "01_invitro_qc_pipeline.py")
    qc_module.main()

    # ---------------------------------------------------------
    # Step 2: Diversity and Effective Capacity Calculation
    # ---------------------------------------------------------
    print("\n▶️ [STEP 2/4] Calculating Diversity Metrics (Shannon & Effective BCs)...")
    metrics_module = load_module_from_path("mod_02", MODULE_DIR / "02_diversity_metrics.py")
    metrics_module.main()

    # ---------------------------------------------------------
    # Step 3: Sequence Space Collision Analysis (Levenshtein)
    # ---------------------------------------------------------
    print("\n▶️ [STEP 3/4] Executing Monte Carlo Sequence Collision Analysis...")
    lev_module = load_module_from_path("mod_03", MODULE_DIR / "03_levenshtein_collision_analysis.py")
    lev_module.main()

    # ---------------------------------------------------------
    # Step 4: Clonal Fidelity and Correlation Analysis
    # ---------------------------------------------------------
    print("\n▶️ [STEP 4/4] Performing Clonal Fidelity & Correlation Analysis...")
    corr_module = load_module_from_path("mod_04", MODULE_DIR / "04_clonal_correlation_analysis.py")
    corr_module.main()
    

    # --------------------------------------------------------
    # Step 5: Master Summary Table Compilation
    # --------------------------------------------------------- 
    print("\n▶️ [STEP 5/5] Compiling Master Summary Table...")
    summary_module = load_module_from_path("mod_05", MODULE_DIR / "05_master_summary_report.py")
    summary_module.main()
    # 👆 ---------------------------------- 👆

    # ---------------------------------------------------------
    # Summary and Next Steps
    # ---------------------------------------------------------
    elapsed_time = time.time() - start_time
    print("\n" + "="*65)
    # ... 后面的保持不变 ...
    # ---------------------------------------------------------
    # Summary and Next Steps
    # ---------------------------------------------------------
    elapsed_time = time.time() - start_time
    print("\n" + "="*65)
    print(f"✅ Master Pipeline Completed Successfully in {elapsed_time:.2f} seconds!")
    print("💾 All hard calculations are done. Metric tables are saved in:")
    print("   -> data/processed/invitro_metrics/")
    print("🎨 NEXT STEP: Open the Jupyter Notebooks in 'notebooks/' to render publication figures.")
    print("="*65)

if __name__ == "__main__":
    main()