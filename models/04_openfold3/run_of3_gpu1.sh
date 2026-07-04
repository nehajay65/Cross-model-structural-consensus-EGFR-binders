#!/bin/bash
#eval "$(conda shell.bash hook)"
#conda activate openfold3

export CUDA_HOME=/usr/local/cuda-13.2/
export PATH="$CUDA_HOME/bin:$PATH"
export LIBRARY_PATH="$CUDA_HOME/targets/x86_64-linux/lib:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$CUDA_HOME/targets/x86_64-linux/lib:${LD_LIBRARY_PATH:-}"

export CUDA_VISIBLE_DEVICES=1

export TORCH_EXTENSIONS_DIR="/home/postyr/.cache/torch_extensions_gpu1"
export XDG_CACHE_HOME="/home/postyr/.cache_gpu1"

export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_CLIENT_MEM_FRACTION=0.95

export CUTLASS_PATH=$(python - << 'PY'
import cutlass_library, pathlib
print(pathlib.Path(cutlass_library.__file__).resolve().parent.joinpath("source"))
PY
)

mkdir -p /home/nehajay/project_new/of3_outputs

run_openfold predict \
    --query_json /home/nehajay/project_new/of3_queries_fixed.json \
    --use_msa_server=False \
    --output_dir /home/nehajay/project_new/of3_outputs/ \
    --runner_yaml /home/postyr/openfold-3/inference_precomputed.yml \
    2>&1 | tee /home/nehajay/project_new/of3_run.log