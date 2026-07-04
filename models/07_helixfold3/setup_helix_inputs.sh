#!/bin/bash
source "$(dirname "$0")/../env.sh"
set -euo pipefail

BINDER_MSAS="$PROJECT_ROOT/a3m_files/alignments_openfold_new"
EGFR_MSAS="$PROJECT_ROOT/egfr_msa_converted"
INPUT_DIR="$PROJECT_ROOT/af3_inputs"
OUT="$PROJECT_ROOT/outputs/helix_fold_outputs"

STO_FILES="mgnify_hits.sto pdb_hits.sto reduced_bfd_hits.sto uniprot_hits.sto uniref90_hits.sto"

count=0; skipped=0; missing=0

for json in "$INPUT_DIR"/*.json; do
    name=$(basename "$json" .json)

    binder_src="$BINDER_MSAS/$name"
    if [[ ! -d "$binder_src" ]]; then
        echo "MISSING binder MSA: $name"
        (( missing++ )) || true
        continue
    fi

    # Create folder structure
    binder_dst="$OUT/$name/msas/protein_A-1/A-1"
    egfr_dst="$OUT/$name/msas/protein_B-1/B-1"
    mkdir -p "$binder_dst" "$egfr_dst"

    # Copy binder .sto files
    for f in $STO_FILES; do
        src="$binder_src/$f"
        [[ -f "$src" ]] && cp "$src" "$binder_dst/$f" || touch "$binder_dst/$f"
    done
    touch "$binder_dst/no_pdb"

    # Copy EGFR .sto files
    for f in $STO_FILES small_bfd_hits.sto no_pdb; do
        [[ -f "$EGFR_MSAS/$f" ]] && cp "$EGFR_MSAS/$f" "$egfr_dst/$f"
    done

    (( count++ )) || true
    echo "OK: $name"
done

echo ""
echo "Done - set up: $count | missing MSAs: $missing"
