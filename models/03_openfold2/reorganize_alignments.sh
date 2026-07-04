#!/bin/bash
source "$(dirname "$0")/../env.sh"

ALIGN_DIR="$PROJECT_ROOT/a3m_files/alignments_openfold_new"
NEW_ALIGN_DIR="$PROJECT_ROOT/a3m_files/alignments_openfold_restructured"

mkdir -p "$NEW_ALIGN_DIR"

for protein_dir in "$ALIGN_DIR"/*/; do
  protein=$(basename "$protein_dir")
  chain_id="${protein%%-*}"   # extracts just the first token, e.g. "swift" from "swift-crane-cedar"

  echo "Fixing: $protein -> subfolder: $chain_id"
  mkdir -p "${NEW_ALIGN_DIR}/${protein}/${chain_id}"
  cp "${protein_dir}"/* "${NEW_ALIGN_DIR}/${protein}/${chain_id}/" 2>/dev/null
done

echo "Done!"
