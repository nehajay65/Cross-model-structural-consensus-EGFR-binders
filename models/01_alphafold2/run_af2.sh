#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /home/postyr/af2_blackwell

export CUDA_VISIBLE_DEVICES=1

DB=/mnt/dsdd_share/AF3/alphafold3/public_databases
HF=/mnt/dsdd_share/weights_databases/helixfold
SRC=/mnt/dsdd_share/src

mkdir -p /home/nehajay/project_new/af2_outputs

for fasta_file in /home/nehajay/project_new/af2_inputs/*/*.fasta; do
    name=$(basename "$fasta_file" .fasta)
    echo "Running: $name"
    python /home/nehajay/alphafold/run_alphafold.py \
        --fasta_paths="$fasta_file" \
        --output_dir=/home/nehajay/project_new/af2_outputs/ \
        --use_precomputed_msas=True \
        --model_preset=multimer \
        --db_preset=reduced_dbs \
        --max_template_date=2021-09-30 \
        --obsolete_pdbs_path=/home/nehajay/project_new/obsolete.dat \
        --use_gpu_relax=True \
        --data_dir=$DB \
        --uniref90_database_path=$DB/uniref90_2022_05.fa \
        --mgnify_database_path=$DB/mgy_clusters_2022_05.fa \
        --small_bfd_database_path=$HF/small_bfd/bfd-first_non_consensus_sequences.fasta \
        --template_mmcif_dir=$DB/mmcif_files \
        --pdb_seqres_database_path=$DB/pdb_seqres_2022_09_28.fasta \
        --uniprot_database_path=$HF/uniprot/uniprot.fasta \
        --jackhmmer_binary_path=$SRC/hmmer-3.4/src/jackhmmer \
        --hmmbuild_binary_path=$SRC/hmmer-3.4/src/hmmbuild \
        --hmmsearch_binary_path=$SRC/hmmer-3.4/src/hmmsearch \
        --hhblits_binary_path=$SRC/hh-suite/bin/hhblits \
        --hhsearch_binary_path=$SRC/hh-suite/bin/hhsearch \
        --kalign_binary_path=$SRC/kalign-3.3.2/src/kalign \
        && echo "yes $name" || echo "no $name"
done