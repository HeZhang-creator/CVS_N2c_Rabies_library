"""
Script: 05_network_topology_cell_centric.py
Description: 
    Reconstructs the trans-synaptic neural network topology (V1 -> LGN) using
    a STRICTLY CELL-CENTRIC framework for Figure 6.
    
    Key principle: each V1 cell (not each barcode) casts exactly one vote 
    in the out-degree calculation. Collision barcodes that infect N V1 cells 
    contribute their inflated out-degree N times, systematically amplifying 
    the divergence of collision-heavy libraries (HH + HdVRz).
    
    Metrics Computed:
    - Pure Network Edges (Yield): Total trans-synaptic edges from singlet V1 cells.
    - In-degree (LGN): Number of distinct V1 barcodes per LGN cell.
    - Out-degree (V1): Number of distinct LGN targets per V1 source cell (cell-centric).
    - Singlet Out-degree (V1): Out-degree using only singlet barcodes (collision-free).
    
    Output:
    - Per-replicate table (stdout)
    - Group-level mean ± SD summary (stdout)
    - Fig6_Topology_Metrics_Summary.csv
    - Fig6_Topology_Metrics_OutDegree_Dist.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time


def compute_network_topology_cell_centric(input_csv, output_prefix):
    start_time = time.time()
    print(f"🚀 [Fig 6] Loading connectomics dataset: {input_csv}")
    
    df = pd.read_csv(input_csv)
    barcode_col = 'RV_Barcode_Corrected' if 'RV_Barcode_Corrected' in df.columns else 'RV_Barcode'
    cell_col = 'Cell_Barcode'
    
    groups = ['HH + HdVRz', 'HH + Twister', 'RzB + Twister']
    
    summary_metrics = []
    outdegree_distributions = []
    singlet_outdegree_distributions = []
    
    print("🕸️  Reconstructing neural network topology (Cell-Centric)...\n")
    
    # ──────────────────────────────────────────
    # Per-replicate table
    # ──────────────────────────────────────────
    HEADER = (f"{'Group':<16s} {'Rep':>4s} {'Pure_Edges':>12s} "
              f"{'In_Deg':>8s} {'Out_Deg':>8s} {'Sngl_Out':>8s} "
              f"{'V1_Cells':>10s} {'LGN_Cells':>10s}")
    print(HEADER)
    print("-" * 76)
    
    for group in groups:
        group_df = df[df['Group'] == group]
        if group_df.empty:
            continue
        for rep in sorted(group_df['Replicate'].unique()):
            rep_df = group_df[group_df['Replicate'] == rep]
            
            # Separate brain region nodes
            v1_df = rep_df[rep_df['Brain_Region'].str.contains('V1', case=False, na=False)]
            lgn_df = rep_df[rep_df['Brain_Region'].str.contains('LGN', case=False, na=False)]
            
            if v1_df.empty or lgn_df.empty:
                continue
            
            # ──────────────────────────────────────
            # Prep: deduplicate (V1 cell, barcode) and (LGN cell, barcode) pairs
            # ──────────────────────────────────────
            v1_pairs = v1_df[[cell_col, barcode_col]].drop_duplicates()
            lgn_pairs = lgn_df[[cell_col, barcode_col]].drop_duplicates()
            
            # ──────────────────────────────────────
            # 1. Identify globally unique (singlet) barcodes in V1
            # ──────────────────────────────────────
            v1_bc_counts = v1_pairs.groupby(barcode_col)[cell_col].nunique()
            singlet_bcs = v1_bc_counts[v1_bc_counts == 1].index
            
            # ──────────────────────────────────────
            # 2. Pure Network Edges — CELL-CENTRIC
            # ──────────────────────────────────────
            v1_singlet = v1_pairs[v1_pairs[barcode_col].isin(singlet_bcs)]
            pure_edges_df = v1_singlet.merge(
                lgn_pairs, on=barcode_col, suffixes=('_V1', '_LGN')
            )
            pure_edges = len(pure_edges_df.drop_duplicates(
                subset=[f'{cell_col}_V1', f'{cell_col}_LGN']
            ))
            
            # ──────────────────────────────────────
            # 3. In-degree Convergence — CELL-CENTRIC
            # ──────────────────────────────────────
            in_degrees = lgn_pairs.groupby(cell_col)[barcode_col].nunique()
            mean_indegree = in_degrees.mean()
            
            # ──────────────────────────────────────
            # 4. Out-degree Divergence — CELL-CENTRIC
            # ──────────────────────────────────────
            v1_to_lgn = v1_pairs.merge(
                lgn_pairs, on=barcode_col, suffixes=('_V1', '_LGN')
            )
            v1_cell_outdegree = v1_to_lgn.groupby(f'{cell_col}_V1')[f'{cell_col}_LGN'].nunique()
            mean_outdegree = v1_cell_outdegree.mean()
            
            n_v1 = v1_pairs[cell_col].nunique()
            n_lgn = lgn_pairs[cell_col].nunique()

            # ──────────────────────────────────────
            # 5. Singlet-Only Out-degree — COLLISION-FREE
            # ──────────────────────────────────────
            singlet_outdegree = None
            if len(singlet_bcs) > 0:
                v1_singlet_only = v1_singlet.copy()
                singlet_v1_to_lgn = v1_singlet_only.merge(
                    lgn_pairs, on=barcode_col, suffixes=('_V1', '_LGN')
                )
                if not singlet_v1_to_lgn.empty:
                    singlet_out = singlet_v1_to_lgn.groupby(
                        f'{cell_col}_V1'
                    )[f'{cell_col}_LGN'].nunique()
                    singlet_outdegree = singlet_out.mean()
                    # Save per-cell singlet-only out-degree for distribution
                    for v1_cell, degree in singlet_out.items():
                        singlet_outdegree_distributions.append({
                            'Group': group, 'Replicate': rep,
                            'V1_Cell_Barcode': v1_cell,
                            'Singlet_Out_Degree': degree
                        })

            # Per-replicate print
            s_out_str = f"{singlet_outdegree:>8.2f}" if singlet_outdegree is not None else f"{'N/A':>8s}"
            print(f"{group:<16s} {rep:>4d} {pure_edges:>12d} "
                  f"{mean_indegree:>8.4f} {mean_outdegree:>8.2f} "
                  f"{s_out_str} "
                  f"{n_v1:>10d} {n_lgn:>10d}")
            
            # Save per-cell out-degree for distribution plotting
            for v1_cell, degree in v1_cell_outdegree.items():
                outdegree_distributions.append({
                    'Group': group, 'Replicate': rep,
                    'V1_Cell_Barcode': v1_cell, 'Out_Degree': degree
                })
            
            # Save summary metrics
            summary_metrics.append({
                'Group': group,
                'Replicate': rep,
                'Pure_Network_Edges': pure_edges,
                'Mean_In_Degree': round(mean_indegree, 4),
                'Mean_Out_Degree': round(mean_outdegree, 4),
                'Mean_Singlet_Out_Degree': round(singlet_outdegree, 4) if singlet_outdegree is not None else None
            })
    
    # ──────────────────────────────────────────
    # Group-level summary
    # ──────────────────────────────────────────
    df_summary = pd.DataFrame(summary_metrics)
    
    print()
    print("Group-Level Summary (mean ± SD, n=3):")
    for grp in groups:
        grp_df = df_summary[df_summary['Group'] == grp]
        if len(grp_df) == 0:
            continue
        for col in ['Pure_Network_Edges', 'Mean_In_Degree', 'Mean_Out_Degree', 'Mean_Singlet_Out_Degree']:
            vals = grp_df[col]
            print(f"  {grp:15s} | {col:20s}: {vals.mean():.2f} ± {vals.std(ddof=1):.2f}")
    
    # ──────────────────────────────────────────
    # Save outputs
    # ──────────────────────────────────────────
    df_outdist = pd.DataFrame(outdegree_distributions)
    
    out_dir = Path(output_prefix).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    
    summary_path = f"{output_prefix}_Summary.csv"
    dist_path = f"{output_prefix}_OutDegree_Dist.csv"
    singlet_dist_path = f"{output_prefix}_SingletOutDegree_Dist.csv"
    df_summary.to_csv(summary_path, index=False)
    df_outdist.to_csv(dist_path, index=False)

    # Save singlet-only out-degree distribution
    if singlet_outdegree_distributions:
        df_singlet_outdist = pd.DataFrame(singlet_outdegree_distributions)
        df_singlet_outdist.to_csv(singlet_dist_path, index=False)
    else:
        df_singlet_outdist = pd.DataFrame()

    print(f"\n🎉  Cell-centric network reconstruction complete in {time.time()-start_time:.2f}s!")
    print(f"💾  Summary: {summary_path}")
    print(f"💾  Distributions: {dist_path}")
    if not df_singlet_outdist.empty:
        print(f"💾  Singlet Out-Degree Distributions: {singlet_dist_path}")
    
    return df_summary, df_outdist


if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[2]
    INPUT = ROOT / "data" / "processed" / "invivo_metrics" / "Master_Corrected_Data.csv"
    OUTPUT_PREFIX = ROOT / "data" / "processed" / "invivo_metrics" / "Fig6_Topology_Metrics"
    compute_network_topology_cell_centric(INPUT, OUTPUT_PREFIX)
