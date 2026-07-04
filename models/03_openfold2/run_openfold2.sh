#!/bin/bash


set -euo pipefail

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

FASTA_DIR="$HOME/project_new/openfold2_fastas"
OUTPUT_DIR="$HOME/project_new/openfold2_outputs_27_5"
PRECOMPUTED_ALIGNMENT_DIR="$HOME/project_new/openfold2_alignments"
MMCIF_DIR="/mnt/dsdd_share/AF3/alphafold3/public_databases/mmcif_files"
PARAMS_DIR="$HOME/openfold/openfold/resources/params"
OPENFOLD_SCRIPT="$HOME/openfold/run_pretrained_openfold.py"

MODEL="model_3_multimer_v3"
DEVICE="cuda:1"
SEED=42


for dir in "$FASTA_DIR" "$MMCIF_DIR" "$PRECOMPUTED_ALIGNMENT_DIR" "$PARAMS_DIR"; do
  if [[ ! -d "$dir" ]]; then
    echo "[ERROR] Directory not found: $dir" >&2
    exit 1
  fi
done

PARAM_FILE="${PARAMS_DIR}/params_${MODEL}.npz"
[[ ! -f "$PARAM_FILE" ]] && { echo "[ERROR] Not found: $PARAM_FILE" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"
LOG_FILE="${OUTPUT_DIR}/run.log"


echo "============================================"
echo " Model  : $MODEL"
echo " Seed   : $SEED"
echo " FASTAs : $FASTA_DIR"
echo " Aligns : $PRECOMPUTED_ALIGNMENT_DIR"
echo " Output : $OUTPUT_DIR"
echo " Started: $(date)"
echo "============================================"

START=$(date +%s)

python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$MMCIF_DIR" \
  --output_dir                 "$OUTPUT_DIR" \
  --config_preset              "$MODEL" \
  --model_device               "$DEVICE" \
  --data_random_seed           "$SEED" \
  --jax_param_path             "$PARAM_FILE" \
  --use_precomputed_alignments "$PRECOMPUTED_ALIGNMENT_DIR" \
  --multimer_ri_gap            200 \
  --skip_relaxation \
  2>&1 | tee "$LOG_FILE"


echo " Log: $LOG_FILE"
