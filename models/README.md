# Model Inference Scripts

Scripts used to run each of the 9 structure-prediction models on the merged EGFR binder set (both AdaptyivBio competition rounds) on the shared GPU workstation. These produce the raw structure predictions consumed by `scripts/01_build_records.py` downstream.

| Folder | Model | Notes |
|---|---|---|
| `01_alphafold2/` | AlphaFold2 | Multimer preset, reduced DBs, GPU relax on |
| `02_alphafold3/` | AlphaFold3 | Uses precomputed MSAs (`--norun_data_pipeline`) |
| `03_openfold2/` | OpenFold2 | `reorganize_alignments.sh` restructures precomputed alignments into OpenFold2's expected layout before running; `amber_relaxation.py` (adapted from [deepmind/alphafold#721](https://github.com/deepmind/alphafold/issues/721)) is run afterward to Amber-relax OpenFold2's unrelaxed outputs, since OpenFold2 doesn't relax by default |
| `04_openfold3/` | OpenFold3 | Two variants (`_gpu0` / `_gpu1`) for splitting the query set across two GPUs |
| `05_boltz2/` | Boltz-2 | `run_boltz_fix_cpu.py` is a thin wrapper patching `torch.linalg.svd` to drop the unsupported `driver` kwarg on this GPU setup, then calls Boltz's own CLI |
| `06_chai1/` | Chai-1 | Reduced trunk recycles / diffusion steps to fit GPU memory |
| `07_helixfold3/` | HelixFold3 | `setup_helix_inputs.sh` assembles the MSA folder structure HelixFold3 expects; `run_helixfold3.sh` then runs inference with skip/resume logic per input |
| `08_protenix/` | Protenix | |
| `09_seedfold/` | SeedFold | Predictions were generated via SeedFold's web dashboard (https://seedfold.io/dashboard/proteinPrediction) - closed, commercial, no local script; results downloaded as `.tar.gz` archives per run. `extract_seedfold_metrics.py` parses the returned `confidence_*.json` files into a metrics table |

## Before running

Set up your machine config once:

```bash
cp models/env.example.sh models/env.sh
# then edit models/env.sh with your own paths (project root, conda envs, shared databases, etc.)
```

Each script sources `env.sh` and uses those variables instead of hardcoded paths. `env.sh` itself is git-ignored, so your local paths never get committed - only `env.example.sh` (the template) is tracked.

Each model also has its own installation and database requirements (e.g. AlphaFold2/3 need their respective sequence databases; OpenFold2 needs precomputed alignments; HelixFold3 needs its own conda env). See each model's own repository for setup instructions:

- AlphaFold2: https://github.com/google-deepmind/alphafold
- AlphaFold3: https://github.com/google-deepmind/alphafold3
- OpenFold2: https://github.com/aqlaboratory/openfold
- OpenFold3: https://github.com/aqlaboratory/openfold-3
- Boltz-2: https://github.com/jwohlwend/boltz
- Chai-1: https://github.com/chaidiscovery/chai-lab
- HelixFold3: https://github.com/PaddlePaddle/PaddleHelix
- Protenix: https://github.com/bytedance/Protenix
- SeedFold: https://seedfold.io (web dashboard, closed/commercial - no public repo)

## Order

Run `run_all_models.sh` for a documented walkthrough of the order these were executed in (it doesn't auto-run everything end-to-end, since each model needs its own environment activated first - see comments inside).
