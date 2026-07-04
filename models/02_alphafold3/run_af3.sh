#!/bin/bash
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

HMMER_BIN="/home/postyr/alphafold3/local_hmmer/bin"
export PATH="$HMMER_BIN:$PATH"

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate /home/postyr/af3_blackwell

AF3_DIR="/home/postyr/alphafold3"
AF3_DB="/mnt/dsdd_share/AF3/alphafold3/public_databases"
AF3_MODEL="/home/postyr/alphafold3/AF3_models"

# paths
INPUT_DIR="/home/nehajay/project_new/af3_inputs"
OUTPUT="/home/nehajay/project_new/af3_outputs_27_5"

CUDA_VISIBLE_DEVICES=0 \
JAX_TRACEBACK_FILTERING=off \
python ${AF3_DIR}/run_alphafold.py \
  --model_dir=${AF3_MODEL} \
  --db_dir=${AF3_DB} \
  --output_dir=${OUTPUT} \
  --input_dir=${INPUT_DIR} \
  --jackhmmer_binary_path="$HMMER_BIN/jackhmmer" \
  --hmmsearch_binary_path="$HMMER_BIN/hmmsearch" \
  --hmmbuild_binary_path="$HMMER_BIN/hmmbuild" \
  --jax_compilation_cache_dir="/home/nehajay/project_new/jax_cache" \
  --norun_data_pipeline   # ✅ skips MSA generation, uses your precomputed MSAs!
