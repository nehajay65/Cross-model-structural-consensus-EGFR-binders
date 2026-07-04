#!/bin/bash
source "$(dirname "$0")/../env.sh"

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENVS_ROOT/boltz-clean"
export CUDA_VISIBLE_DEVICES=0
mkdir -p "$PROJECT_ROOT/boltz_outputs"

python "$(dirname "$0")/run_boltz_fix_cpu.py" predict "$PROJECT_ROOT/boltz_yaml_inputs/" \
    --out_dir "$PROJECT_ROOT/boltz_outputs_27_5/" \
    --seed 42
