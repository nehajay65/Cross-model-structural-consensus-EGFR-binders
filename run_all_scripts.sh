#!/bin/bash
set -e  # stop on any error

cd "$(dirname "$0")"

echo "=== Step 1: build_records.py ==="
python scripts/01_build_records.py \
    --root /home/nehajay/project_new/a_predictions/organized_predictions \
    --out-dir ./data \
    --target-length 621

echo "=== Step 2: extract_metrics.py ==="
python scripts/02_extract_metrics.py \
    --records ./data/prediction_records.parquet \
    --out-dir ./data

echo "=== Step 3: rmsd_rank1.py ==="
python scripts/03_rmsd_rank1.py \
    --records ./data/prediction_records_with_metrics.parquet \
    --out-dir ./data

echo "=== Step 4: contacts.py ==="
python scripts/04_contacts.py

echo "=== Step 5: plots.py ==="
python scripts/05_plots.py \
    --data-dir ./data \
    --out-dir ./plots

echo "=== ALL DONE ==="
