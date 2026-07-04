# Cross-Model Structural Consensus for Evaluating De Novo Protein Binders

**MSc Bioinformatics Project 1 - University of Copenhagen**
Neha Jay · Supervisors: Albert Jelke Kooistra, Melika Keshavarz (Dept. of Drug Design & Pharmacology), Thomas Wim Hamelryck (Dept. of Biology, SCARB)

## Overview

Experimental validation of de novo-designed protein binders is slow and expensive, and individual AI structure-prediction models don't always transfer their confidence metrics reliably to novel (non-natural) sequences. This project asks: **can agreement between multiple independent structure-prediction models act as a proxy for experimental binding, without running a wet-lab assay?**

Binder sequences and predicted complexes were pooled from **both rounds of AdaptyivBio's EGFR Binder Design Competition** (579 de novo binders total - 61 experimentally confirmed binders, 518 non-binders) and run through **9 structure-prediction models**: AlphaFold2, AlphaFold3, OpenFold2, OpenFold3, Boltz-2, Chai-1, HelixFold3, Protenix, and SeedFold.

For each binder, three structural consensus metrics were computed across all rank-1, pairwise model combinations:

- **Jaccard index** - overlap of predicted target-binder interface contacts between model pairs
- **Pocket-aligned RMSD** - binder positional consistency relative to the EGFR binding pocket
- **Binder-aligned Cα RMSD** - binder fold consistency between models

These were z-score normalized and combined into a **composite consensus score**, then validated against the known experimental outcomes with Mann-Whitney U tests (rank-biserial effect sizes, FDR-corrected).

## Key findings

- The composite consensus score separated true binders from non-binders better than any single metric alone (rank-biserial r = 0.34, p = 1.5x10^-5).
- Interface contact agreement (Jaccard index) was a significant, positive predictor of binding (r = 0.24, p = 0.006).
- Binder-fold divergence (binder-aligned RMSD) was, surprisingly, *higher* among true binders (r = 0.29, p = 0.001) - the opposite of the initial hypothesis.
- Pocket-aligned RMSD showed no significant association (p = 0.63), likely because all binders in this dataset target the same EGFR epitope.
- Effect sizes throughout were modest - this is an exploratory, single-target, class-imbalanced study, not a validated predictor.

See `report/EGFR_binder_consensus_report.pdf` for the full write-up, figures, and discussion.

## Repository structure

```
├ report/
│   └ EGFR_binder_consensus_report.pdf   # Full write-up
├ models/                       # Step 0: per-model structure prediction
│   ├ 01_alphafold2/
│   ├ 02_alphafold3/
│   ├ 03_openfold2/             # + Amber relaxation (OpenFold2 doesn't relax by default)
│   ├ 04_openfold3/
│   ├ 05_boltz2/
│   ├ 06_chai1/
│   ├ 07_helixfold3/
│   ├ 08_protenix/
│   ├ 09_seedfold/
│   ├ run_all_models.sh
│   └ README.md                 # per-model notes, envs, upstream repos
├ scripts/                      # Steps 2-5: analysis pipeline
│   ├ 01_build_records.py       #assembles the raw per-model prediction outputs into the `prediction_records.parquet` the analysis pipeline starts from.
|   ├ 02_extract_metrics.py     # Pull native confidence metrics (pLDDT, iPTM, PTM) per model
│   ├ 03_rmsd_rank1.py          # Rank-1 filtering + pairwise binder/pocket RMSD
│   ├ 04_contacts.py            # Interface contact extraction (5A cutoff) -> Jaccard input
│   └ 05_plots.py               # Summary figures
├ notebooks/
│   └ consensus_analysis.ipynb  # Composite score, statistics, publication figures
├ run_all_scripts.sh
└ requirements.txt
```

## Pipeline

| Step | Location | Description |
|------|----------|-------------|
| 0 | `models/` | Run each of the 9 structure-prediction models on the merged binder set |
| 1 | `scripts/01_build_records.py` | Compile all model outputs + metadata into `prediction_records.parquet` |
| 2 | `scripts/02_extract_metrics.py` | Extract native confidence metrics per method |
| 3 | `scripts/03_rmsd_rank1.py` | Rank-1 filtering, binder CA RMSD, pocket-aligned RMSD |
| 4 | `scripts/04_contacts.py` | Rank-1 interface contact extraction |
| 5 | `scripts/05_plots.py` | Summary figures |
| - | `notebooks/consensus_analysis.ipynb` | Composite score + statistics + publication figures |

Model inference (`models/`) needs a separate environment per model - see `models/README.md`. Once step 1 is added, run the analysis pipeline with:

```bash
bash run_all_scripts.sh
```

## Setup

```bash
pip install -r requirements.txt
```

`scripts/amber_relaxation.py` additionally requires an AlphaFold2 installation (`alphafold`, `absl-py`) and is run separately.

## Data

Raw prediction outputs and structure files are not tracked in this repo (they're large - 579 binders x 9 models). `data/`, `plots/`, `*.parquet`, and structure files are git-ignored by default; adjust `.gitignore` if you want to track a subset (e.g. the final summary parquets).

## Data source

- Adaptyv Bio, *AdaptyivBio EGFR Binder Design Competition*, Rounds 1 & 2, via [ProteinBase](https://proteinbase.com).
