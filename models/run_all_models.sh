#!/bin/bash
# Model inference - run order
# Each model needs its own conda/venv environment, so this script
# is a documented walkthrough rather than a single one-shot pipeline.
# Uncomment the block(s) you want to run, after activating the
# matching environment and editing the hardcoded paths inside
# each script to match your own setup.
set -euo pipefail
cd "$(dirname "$0")"

if [[ ! -f env.sh ]]; then
  echo "Missing models/env.sh - copy env.example.sh to env.sh and fill in your paths first."
  exit 1
fi
source env.sh

echo "=== 1/9: AlphaFold2 ==="
# bash 01_alphafold2/run_af2.sh

echo "=== 2/9: AlphaFold3 ==="
# bash 02_alphafold3/run_af3.sh

echo "=== 3/9: OpenFold2 ==="
# bash 03_openfold2/reorganize_alignments.sh   # restructure precomputed alignments first
# bash 03_openfold2/run_openfold2.sh
# python 03_openfold2/amber_relaxation.py <unrelaxed_pdb_dir> -output_dir <relaxed_dir>

echo "=== 4/9: OpenFold3 ==="
# bash 04_openfold3/run_of3_gpu0.sh
# bash 04_openfold3/run_of3_gpu1.sh

echo "=== 5/9: Boltz-2 ==="
# bash 05_boltz2/run_boltz2.sh   # internally calls run_boltz_fix_cpu.py

echo "=== 6/9: Chai-1 ==="
# bash 06_chai1/run_chai1.sh

echo "=== 7/9: HelixFold3 ==="
# bash 07_helixfold3/setup_helix_inputs.sh
# bash 07_helixfold3/run_helixfold3.sh

echo "=== 8/9: Protenix ==="
# bash 08_protenix/run_protenix.sh

echo "=== 9/9: SeedFold ==="
# Submit binder x EGFR predictions via https://seedfold.io/dashboard/proteinPrediction,
# download the resulting .tar.gz archives, unpack them, then:
# python 09_seedfold/extract_seedfold_metrics.py --results_dir <unpacked_results_dir>

echo "All model steps documented above - uncomment as needed."
