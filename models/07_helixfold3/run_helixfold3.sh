#!/bin/bash
source "$(dirname "$0")/../env.sh"
set -euo pipefail
shopt -s nullglob

# Workstation config
HF3_HOME="$HELIXFOLD3_HOME"
INPUT_GLOB="$PROJECT_ROOT/hf3_inputs_2/*.json"
OUTPUT="$PROJECT_ROOT/outputs/helix_fold_outputs_27_5"
RUN_INFER="${HF3_HOME}/run_infer_no_dbs_test.sh"

# GPU and template behavior
export CUDA_VISIBLE_DEVICES=0
export NO_PDB_TEMPLATES=TRUE

# Match this with --infer_times in run_infer_no_dbs_test.sh
N_PREDS=1

cd "$HF3_HOME"

# Check essentials
if [[ ! -x "$RUN_INFER" ]]; then
  echo "ERROR: run script not executable or missing:"
  echo "  $RUN_INFER"
  echo "Fix with: chmod +x $RUN_INFER"
  exit 1
fi

if [[ ! -d "$OUTPUT" ]]; then
  mkdir -p "$OUTPUT"
fi

mapfile -t inputs < <(ls -1 ${INPUT_GLOB} 2>/dev/null | sort)
n=${#inputs[@]}

if (( n == 0 )); then
  echo "ERROR: No inputs found matching:"
  echo "  ${INPUT_GLOB}"
  exit 1
fi

echo "HF3_HOME: $HF3_HOME"
echo "INPUTS:   $INPUT_GLOB"
echo "OUTPUT:   $OUTPUT"
echo "RUNNER:   $RUN_INFER"
echo "GPU:      CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
echo "NO_PDB_TEMPLATES=${NO_PDB_TEMPLATES}"
echo "Total inputs: $n"
echo

# Completion check
is_complete() {
  local name="$1"
  local base="${OUTPUT}/${name}"

  for k in $(seq 1 "$N_PREDS"); do
    local d="${base}/${name}-pred-${k}-1"
    [[ -d "$d" ]] || return 1

    find "$d" -maxdepth 1 -type f -name "*.cif" -print -quit 2>/dev/null | grep -q . || return 1
  done

  return 0
}

clean_partial() {
  local name="$1"
  local base="${OUTPUT}/${name}"

  echo "CLEANING partial outputs: $name"
  rm -rf "${base}/${name}-pred-"[1-5]"-1" \
         "${base}/${name}-rank"[1-5] 2>/dev/null || true
}

# Run sequentially
for input_json in "${inputs[@]}"; do
  name="$(basename "$input_json" .json)"

  if is_complete "$name"; then
    echo "SKIP complete: $name"
    continue
  fi

  echo "RUN: $name"

  if [[ -d "${OUTPUT}/${name}" ]]; then
    clean_partial "$name"
  fi

  bash "$RUN_INFER" "$input_json" "$OUTPUT"

  echo "DONE: $name"
  echo
done

echo "All jobs finished."
