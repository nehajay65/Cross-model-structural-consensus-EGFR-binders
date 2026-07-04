#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /home/postyr/.conda/envs/protenix

export CUDA_VISIBLE_DEVICES=0

mkdir -p /home/nehajay/project_new/protenix_outputs_27_5

for json_file in /home/nehajay/project_new/protenix_inputs/*.json; do
    name=$(basename "$json_file" .json)
    echo "Running: $name"
    protenix pred \
        -i "$json_file" \
        -o /home/nehajay/project_new/protenix_outputs_27_5/ \
        -n protenix_base_default_v1.0.0 \
        --seed 42 && echo "done: $name" || echo "FAILED: $name"
done
