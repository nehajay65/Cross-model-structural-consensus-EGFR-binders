#!/bin/bash
source "$(dirname "$0")/../env.sh"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENVS_ROOT/protenix"

export CUDA_VISIBLE_DEVICES=0

mkdir -p "$PROJECT_ROOT/protenix_outputs_27_5"

for json_file in "$PROJECT_ROOT"/protenix_inputs/*.json; do
    name=$(basename "$json_file" .json)
    echo "Running: $name"
    protenix pred \
        -i "$json_file" \
        -o "$PROJECT_ROOT/protenix_outputs_27_5/" \
        -n protenix_base_default_v1.0.0 \
        --seed 42 && echo "done: $name" || echo "FAILED: $name"
done
