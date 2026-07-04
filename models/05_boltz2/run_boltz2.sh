#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /home/postyr/.conda/envs/boltz-clean
export CUDA_VISIBLE_DEVICES=0
mkdir -p /home/nehajay/project_new/boltz_outputs

python /home/nehajay/project_new/run_boltz_fix_cpu.py predict /home/nehajay/project_new/boltz_yaml_inputs/ \
    --out_dir /home/nehajay/project_new/boltz_outputs_27_5/ \
    --seed 42