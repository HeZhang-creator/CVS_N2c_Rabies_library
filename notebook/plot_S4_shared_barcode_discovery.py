"""
shared_barcode_discovery.py
============================
Joint V1–LGN downsampling to evaluate shared-barcode discovery rate
as a function of sequencing depth.

Input:   data/processed/invivo_metrics/Master_Corrected_Data.csv
Output:  data/processed/invivo_metrics/Downsampling_Robustness.csv
         data/processed/invivo_metrics/Figure_SharedBC_Discovery.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import time
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 0. 全局配置
# ============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'axes.unicode_minus': False,
    'axes.linewidth': 1.5,
})

CUSTOM_COLORS = {
    'HH + HdVRz':  '#1F77B4',
    'HH + Twister': '#D62728',
    'RzB + Twister': '#FF7F0E',
}
GROUPS = list(CUSTOM_COLORS.keys())

FRACTIONS = np.arange(0.1, 1.05, 0.1)   # 10% → 100%
N_TRIALS  = 10
RANDOM_SEED = 42

# ============================================================================
# 1. 路径
# ============================================================================
# Resolve repo root: when placed at CVS_N2c_Rabies_library/notebook/ this
# script, parents[1] gives the repo root.  Adjust as needed.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent if SCRIPT_DIR.name == 'notebook' else Path.cwd()

INPUT_CSV      = REPO_ROOT / "data" / "processed" / "invivo_metrics" / "Master_Corrected_Data.csv"
OUTPUT_DIR     = REPO_ROOT / "data" / "processed" / "invivo_metrics"
DOWNSAMPLE_CSV = OUTPUT_DIR / "Downsampling_Robustness.csv"
OUTPUT_FIG_PDF = OUTPUT_DIR / "Figure_SharedBC_Discovery.pdf"


# ============================================================================
# 2. 加载
# ============================================================================
def load_data(path):
    if not path.exists():
        raise FileNotFoundError(f"❌ {path}")
    df = pd.read_csv(path)
    bc_col = 'RV_Barcode_Corrected' if 'RV_Barcode_Corrected' in df.columns else 'RV_Barcode'
    required = ['Group', 'Brain_Region', 'Replicate', bc_col]
    for c in required:
        if c not in df.columns:
            raise KeyError(f"Missing column: {c}")
    df['Region_Short'] = df['Brain_Region'].apply(
        lambda x: 'V1' if 'V1' in str(x) else ('LGN' if 'LGN' in str(x) else str(x))
    )
    return df, bc_col


# ============================================================================
# 3. Joint V1+LGN downsampling
# ============================================================================
def joint_downsampling(df, bc_col):
    """
    For each Group × Replicate, pool V1+LGN, subsample at increasing fractions,
    and count shared barcodes (V1 ∩ LGN).  10 independent trials per fraction.
    """
    print("🔬 Joint V1+LGN downsampling …")
    results = []
    rng = np.random.default_rng(RANDOM_SEED)

    groups_reps = df[['Group', 'Replicate']].drop_duplicates().values
    n_total = len(groups_reps) * len(FRACTIONS) * N_TRIALS
    n_done = 0

    for grp, rep in groups_reps:
        rep_data = df[(df['Group'] == grp) & (df['Replicate'] == rep)]
        v1  = rep_data[rep_data['Region_Short'] == 'V1'].reset_index(drop=True)
        lgn = rep_data[rep_data['Region_Short'] == 'LGN'].reset_index(drop=True)

        if v1.empty or lgn.empty:
            print(f"  ⚠️ [{grp}] Rep {rep} — missing region, skipping")
            continue

        pooled = pd.concat([v1, lgn], ignore_index=True)
        n_pooled = len(pooled)

        for frac in FRACTIONS:
            n_sample = max(1, int(n_pooled * frac))
            for trial in range(N_TRIALS):
                idxs = rng.choice(n_pooled, size=n_sample, replace=False)
                sub = pooled.iloc[idxs]

                v1_bcs  = set(sub.loc[sub['Region_Short'] == 'V1', bc_col])
                lgn_bcs = set(sub.loc[sub['Region_Short'] == 'LGN', bc_col])
                shared = len(v1_bcs & lgn_bcs)

                results.append({
                    'Group': grp,
                    'Replicate': rep,
                    'Fraction': round(frac, 2),
                    'Trial': trial,
                    'Sampled_Cells': n_sample,
                    'Shared_Barcodes': shared,
                })
                n_done += 1
        print(f"  ✅ [{grp}] Rep {rep} — done")

    print(f"  🎉 {n_done} trials completed")
    return pd.DataFrame(results)


# ============================================================================
# 4. 绘图
# ============================================================================
def plot_shared_bc_discovery(df_down, out_pdf):
    agg = df_down.groupby(['Group', 'Fraction']).agg(
        Shared_Mean=('Shared_Barcodes', 'mean'),
        Shared_SEM=('Shared_Barcodes', 'sem'),
    ).reset_index()

    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)

    for grp in GROUPS:
        gdata = agg[agg['Group'] == grp]
        ax.errorbar(gdata['Fraction'] * 100,
                    gdata['Shared_Mean'],
                    yerr=gdata['Shared_SEM'],
                    marker='o', markersize=6, linewidth=2,
                    color=CUSTOM_COLORS[grp], capsize=3.5, label=grp,
                    zorder=3)

    ax.set_xlabel('Subsampling Fraction of UMIs (%)', fontsize=12)
    ax.set_ylabel('Shared Barcodes (V1 ∩ LGN)', fontsize=12)
    ax.set_title('Trans-synaptic Barcode Discovery Rate', fontsize=13, fontweight='bold', pad=12)
    ax.set_xlim(0, 105)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(frameon=False, fontsize=10)
    sns.despine()

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches='tight', transparent=True)
    print(f"💾 Figure → {out_pdf}")
    plt.show()


# ============================================================================
# 5. Main
# ============================================================================
if __name__ == "__main__":
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df, bc_col = load_data(INPUT_CSV)
    print(f"✅ {len(df)} records  |  barcode column = {bc_col}")

    df_down = joint_downsampling(df, bc_col)
    df_down.to_csv(DOWNSAMPLE_CSV, index=False)
    print(f"💾 CSV → {DOWNSAMPLE_CSV}  ({len(df_down)} rows)")

    plot_shared_bc_discovery(df_down, OUTPUT_FIG_PDF)

    print(f"\n⏱️  {time.time() - t0:.1f} s")
