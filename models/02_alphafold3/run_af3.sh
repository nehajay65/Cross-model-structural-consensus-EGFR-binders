#!/bin/bash
source "$(dirname "$0")/../env.sh"

export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

HMMER_BIN="$ALPHAFOLD3_CODE_DIR/local_hmmer/bin"
export PATH="$HMMER_BIN:$PATH"

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENVS_ROOT/af3_blackwell"

AF3_DIR="$ALPHAFOLD3_CODE_DIR"
AF3_DB="$SHARED_DB_ROOT/AF3/alphafold3/public_databases"
AF3_MODEL="$ALPHAFOLD3_CODE_DIR/AF3_models"

# paths
INPUT_DIR="$PROJECT_ROOT/af3_inputs"
OUTPUT="$PROJECT_ROOT/af3_outputs_27_5"

CUDA_VISIBLE_DEVICES=0 \
JAX_TRACEBACK_FILTERING=off \
python "${AF3_DIR}/run_alphafold.py" \
  --model_dir="${AF3_MODEL}" \
  --db_dir="${AF3_DB}" \
  --output_dir="${OUTPUT}" \
  --input_dir="${INPUT_DIR}" \
  --jackhmmer_binary_path="$HMMER_BIN/jackhmmer" \
  --hmmsearch_binary_path="$HMMER_BIN/hmmsearch" \
  --hmmbuild_binary_path="$HMMER_BIN/hmmbuild" \
  --jax_compilation_cache_dir="$PROJECT_ROOT/jax_cache" \
  --norun_data_pipeline   # skips MSA generation, uses precomputed MSAs
