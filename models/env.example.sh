#!/bin/bash
# Machine-specific configuration for model inference scripts.
#
# Copy this file to `env.sh` (in the same models/ folder) and
# fill in the paths for your own workstation. `env.sh` is
# git-ignored, so your local paths never get committed.
#
#   cp models/env.example.sh models/env.sh
#   # then edit models/env.sh

# Root folder where all project inputs/outputs live
# (fasta/json inputs, MSAs, per-model outputs, logs)
export PROJECT_ROOT="${PROJECT_ROOT:-$HOME/project_new}"

# Root folder containing your conda environments for each model
# e.g. if AlphaFold2's env lives at /some/path/af2_env,
# set CONDA_ENVS_ROOT=/some/path
export CONDA_ENVS_ROOT="${CONDA_ENVS_ROOT:-$HOME/.conda/envs}"

# Root of shared/mounted reference databases (UniRef90, MGnify, BFD,
# mmCIF templates, etc.) used by AlphaFold2/3 and OpenFold2
export SHARED_DB_ROOT="${SHARED_DB_ROOT:-/mnt/shared_dbs}"

# Local clones of each model's code repository
export ALPHAFOLD2_CODE_DIR="${ALPHAFOLD2_CODE_DIR:-$HOME/alphafold}"
export ALPHAFOLD3_CODE_DIR="${ALPHAFOLD3_CODE_DIR:-$HOME/alphafold3}"
export OPENFOLD2_CODE_DIR="${OPENFOLD2_CODE_DIR:-$HOME/openfold}"
export OPENFOLD3_CODE_DIR="${OPENFOLD3_CODE_DIR:-$HOME/openfold-3}"
export HELIXFOLD3_HOME="${HELIXFOLD3_HOME:-$HOME/helixfold3}"

# CUDA toolkit location (used by OpenFold3)
export CUDA_TOOLKIT_HOME="${CUDA_TOOLKIT_HOME:-/usr/local/cuda}"
