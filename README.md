# CVS-N2c Rabies Virus Barcode Library — Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-%E2%89%A53.9-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

Bioinformatics analysis pipeline for **ribozyme‑optimized, hyper‑diverse barcoded CVS‑N2cΔG rabies virus** single‑neuron projection connectivity mapping.

This repository accompanies:

_He Zhang1,3,*, Daoyin Liu1,2,3, Xueping Gao1, Xiangyu Ge1, Yahui Li1, Xia Zhang1, Lei Jin1,*_
1 Lingang Laboratory, Shanghai 200031, China.
2 School of Life Science and Technology, ShanghaiTech University, 201210, Shanghai, China.
3 These authors contributed equally: He Zhang, Daoyin Liu
*Correspondence: zhanghe@lglab.ac.cn	; jinlei@lglab.ac.cn

---

## Overview

We developed a systematic framework combining **Twister‑ribozyme‑engineered hyper‑diverse barcode libraries** with **collision‑free, cell‑centric topology reconstruction** to achieve high‑fidelity monosynaptic projection mapping at single‑neuron resolution.

This repository contains all custom analysis scripts for:

- In vitro barcode quality control, error correction, and diversity profiling
- In vivo single‑cell barcode preprocessing (anchor‑first pipeline)
- Dual‑resolution collision analysis (barcode‑centric + cell‑centric)
- Cell‑centric network topology reconstruction
- All manuscript figure generation

---

## Pipeline Architecture

```
                               ┌──────────────────────────┐
                               │   Raw sequencing reads    │
                               │   (Plasmid / P1 / EnvA)  │
                               └────────────┬─────────────┘
                                            │
                               ┌────────────▼─────────────┐
                               │  01_invitro_qc_pipeline  │
                               │  • Greedy 1‑bp error     │
                               │    correction            │
                               │  • Multinomial depth     │
                               │    normalization         │
                               └────────────┬─────────────┘
                                            │
              ┌─────────────────────────────┼─────────────────────────────┐
              │                             │                             │
   ┌──────────▼──────────┐       ┌─────────▼─────────┐       ┌───────────▼──────────┐
   │  02_diversity_      │       │  03_levenshtein_   │       │  04_clonal_          │
   │  metrics.py         │       │  collision_        │       │  correlation_        │
   │  • Shannon entropy  │       │  analysis.py       │       │  analysis.py         │
   │  • Effective BCs    │       │  • Edit distance   │       │  • Plasmid ↔ EnvA    │
   │  • Rank‑abundance   │       │  • Collision rate  │       │  • Spearman ρ        │
   └──────────┬──────────┘       └─────────┬─────────┘       └───────────┬──────────┘
              │                             │                             │
              └─────────────────────────────┼─────────────────────────────┘
                                            │
                              ┌─────────────▼─────────────┐
                              │      In Vivo Pipeline      │
                              └─────────────┬─────────────┘
                                            │
                              ┌─────────────▼─────────────┐
                              │  00_build_master_dataset  │
                              │  01_data_cleaning_        │
                              │     clustering.py         │
                              │  • Anchor‑first QC        │
                              │  • UMI‑Kinds ≥ 3          │
                              │  • 5% fraction filter     │
                              │  • 1‑bp clonal merge      │
                              └─────────────┬─────────────┘
                                            │
              ┌─────────────────────────────┼─────────────────────────────┐
              │                             │                             │
   ┌──────────▼──────────┐       ┌─────────▼─────────┐       ┌───────────▼──────────┐
   │  02_diversity_      │       │  03_rarefaction.py │       │  Dual‑Resolution      │
   │  metrics.py         │       │  • Without         │       │  Collision Analysis   │
   │  • Absolute BCs     │       │    replacement     │       │  plot_fig5C/DE        │
   │  • Shannon H        │       │  • Step = 500      │       │  • BC‑centric         │
   │  • Info density     │       │  • Seed = 42       │       │  • Cell‑centric       │
   └──────────┬──────────┘       └─────────┬─────────┘       └───────────┬──────────┘
              │                             │                             │
              └─────────────────────────────┼─────────────────────────────┘
                                            │
                              ┌─────────────▼─────────────┐
                              │  Cell‑Centric Network      │
                              │  Topology Reconstruction   │
                              │  • Pure network edges      │
                              │  • V1 in‑degree            │
                              │  • LGN out‑degree          │
                              └───────────────────────────┘
```

---

## Repository Structure

```
CVS_N2c_Rabies_library/
│
├── in_vitro/                          # In vitro library characterization
│   ├── 01_invitro_qc_pipeline.py      # Error correction + depth normalization
│   ├── 02_diversity_metrics.py        # Shannon entropy, effective capacity
│   ├── 03_levenshtein_collision_analysis.py   # Edit-distance collision analysis
│   ├── 04_clonal_correlation_analysis.py      # Plasmid–EnvA clonal fidelity
│   │
│   ├── plot_fig3_panel_b_umi_saturation.py    # UMI rarefaction curves (Fig. 3B)
│   ├── plot_fig3_panel_e_frequency_flags.py   # Rank‑abundance frequency (Fig. 3E)
│   ├── plot_fig3_panel_f_rank_abundance.py    # Cumulative dominance (Fig. 3F)
│   ├── plot_fig3_panel_g_collision_histogram.py # Levenshtein histogram (Fig. 3G)
│   └── plot_S3_Clonal_Correlation.ipynb       # Clonal correlation scatter (Fig. S2)
│
├── in_vivo/                           # In vivo single‑cell barcode analysis
│   ├── 00_build_master_dataset.py     # Metadata‑driven data consolidation
│   ├── 01_data_cleaning_and_clustering.py  # Anchor‑first QC pipeline
│   ├── 02_diversity_metrics.py        # In vivo diversity (Fig. 4B, C, F, G)
│   ├── 03_rarefaction.py              # Cell‑barcode rarefaction (Fig. 4D, E)
│   │
│   ├── plot_fig4_Library_Diversity.py          # Grouped diversity bar plots
│   ├── plot_fig4_Rarefaction_Analysis.py       # Rarefaction curves
│   ├── plot_fig5C_barcode_centric.py           # BC‑centric collision profiles
│   ├── plot_fig5DE_cell_centric.py             # Cell‑centric collision profiles
│   └── plot_fig6_network_topology.py           # Network topology reconstruction
│
├── metadata/
│   ├── Normalized_Metadata.csv         # In vitro sample metadata
│   └── MasterCorrectedData.csv         # In vivo master corrected dataset
│
├── data/                               # Example / test data (see below)
│
├── requirements.txt                    # Python dependencies
├── README.md
└── LICENSE
```

---

## Key Methods

### In Vitro Pipeline

#### Error Correction (`01_invitro_qc_pipeline.py`)
- **Greedy 1‑bp adjacency network** error correction
- Exhaustive enumeration of single‑nucleotide substitutions, deletions, insertions
- Barcodes sorted by descending UMI abundance; core clones absorb 1‑bp neighbors
- Stage‑specific multinomial downsampling to minimum UMI depth (seed = 42)
- Target depths: Plasmid 736,197 | P1‑B19G 2,945,431 | EnvA 1,775,959

#### Diversity Metrics (`02_diversity_metrics.py`)
- Shannon entropy: \( H = -\Sigma\ p_i \ln(p_i) \)
- Effective library capacity: \( E = e^{H} \)
- Computed on stage‑normalized, error‑corrected data

#### Collision Analysis (`03_levenshtein_collision_analysis.py`)
- UMI‑weighted Monte Carlo sampling (n = 2,000 molecules, seed = 42)
- Pairwise Levenshtein edit distances; collision rate = fraction of zero‑distance pairs

### In Vivo Pipeline

#### Anchor‑First Preprocessing (`01_data_cleaning_and_clustering.py`)
- **Phase 1** — Stringent quality filters:
  - `UMI‑Kinds ≥ 3` (removes singleton/doubleton artifacts)
  - Intra‑cellular fraction ≥ 5% (removes ambient RNA contamination)
- **Phase 2** — 1‑bp Levenshtein greedy clonal consolidation, applied within biological boundaries (Group × BrainRegion × Replicate)

#### Dual‑Resolution Collision Analysis
- **Barcode‑centric**: classify barcode species by number of infected host cells (Unique / Low / Medium / High)
- **Cell‑centric**: quantify proportion of collision‑free vs. collision‑affected host neurons

#### Cell‑Centric Network Topology
- **Pure network edges**: unique (V1, LGN) pairs mediated exclusively by globally singlet barcodes
- **V1 in‑degree**: number of distinct LGN inputs per V1 starter neuron
- **LGN out‑degree**: number of distinct V1 sources per LGN target cell

---

## Installation

```bash
git clone https://github.com/HeZhang‑creator/CVS_N2c_Rabies_library.git
cd CVS_N2c_Rabies_library
pip install -r requirements.txt
```

### Dependencies

```
python ≥ 3.9
numpy
pandas
scipy
matplotlib
seaborn
python-Levenshtein
jupyter          # for .ipynb notebooks
```

---

## Usage

### In Vitro Analysis

```bash
# 1. Error correction and depth normalization
python in_vitro/01_invitro_qc_pipeline.py

# 2. Diversity metrics
python in_vitro/02_diversity_metrics.py

# 3. Collision analysis
python in_vitro/03_levenshtein_collision_analysis.py

# 4. Clonal correlation
python in_vitro/04_clonal_correlation_analysis.py

# Generate figures
python in_vitro/plot_fig3_panel_b_umi_saturation.py
python in_vitro/plot_fig3_panel_e_frequency_flags.py
python in_vitro/plot_fig3_panel_f_rank_abundance.py
python in_vitro/plot_fig3_panel_g_collision_histogram.py
```

### In Vivo Analysis

```bash
# 1. Build master dataset
python in_vivo/00_build_master_dataset.py

# 2. Anchor‑first QC and clonal consolidation
python in_vivo/01_data_cleaning_and_clustering.py

# 3. Diversity and rarefaction
python in_vivo/02_diversity_metrics.py
python in_vivo/03_rarefaction.py

# Generate figures
python in_vivo/plot_fig4_Library_Diversity.py
python in_vivo/plot_fig4_Rarefaction_Analysis.py
python in_vivo/plot_fig5C_barcode_centric.py
python in_vivo/plot_fig5DE_cell_centric.py
python in_vivo/plot_fig6_network_topology.py
```

---

## Reproducibility

All random sampling operations use a **fixed seed (42)**. Analysis parameters:

| Parameter | Value |
|-----------|-------|
| Random seed | 42 |
| Levenshtein error‑correction radius | 1 bp |
| In vitro sampling (EnVA) | 2,000 molecules |
| In vivo rarefaction step | 500 observations |
| Anchor QC: UMI‑Kinds threshold | ≥ 3 |
| Anchor QC: intra‑cellular fraction | ≥ 5% |

---

## Data Availability

Raw and processed data are available upon reasonable request. The plasmid and viral barcode libraries described in this study are available from the corresponding authors under standard material transfer agreements.

---

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@article{Zhang2026CVSN2c,
  title     = {Ribozyme‑optimized CVS‑N2c{Δ}G Rabies Virus with Hyper‑diverse Barcodes
               for High‑fidelity Single‑Neuron Projection Connectivity},
  author    = {Zhang, He and Liu, Daoyin and Gao, Xueping and Ge, Xiangyu and
               Li, Yahui and Zhang, Xia and Jin, Lei},
  journal   = {Lingang Laboratory},
  year      = {2026},
}
```

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

- **He Zhang** — [zhanghe@lglab.ac.cn](mailto:zhanghe@lglab.ac.cn)
- **Lei Jin** — [jinlei@lglab.ac.cn](mailto:jinlei@lglab.ac.cn)

Lingang Laboratory, Shanghai 200031, China
