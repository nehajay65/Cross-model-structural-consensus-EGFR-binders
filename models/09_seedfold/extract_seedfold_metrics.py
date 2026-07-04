 #first unzip all the tar,gz files using 
# for f in *.tar.gz; do tar -xzf "$f"; done
 
 
 

"""
SeedFold Metrics Extractor
Extracts pLDDT, ptm, iptm, ranking_score, and per-chain metrics
from confidence_*.json files across all unpacked SeedFold results.

Usage:
    python extract_seedfold_metrics.py --results_dir ./seedfold_results/batch_1
    python extract_seedfold_metrics.py --results_dir ./seedfold_results/batch_1 --out metrics.csv
"""

import json
import argparse
import glob
import os
import numpy as np
import pandas as pd
from pathlib import Path

# Helpers
def load_confidence(json_path: str) -> dict:
    with open(json_path) as f:
        return json.load(f)

def extract_metrics(data: dict) -> dict:
    """Pull every useful scalar and array-derived metric from one confidence JSON."""
    metrics = {}

    #  Top-level scalars (present in most SeedFold / AF3-style outputs) 
    for key in ("ptm", "iptm", "ranking_score", "fraction_disordered",
                "has_clash", "num_recycles"):
        if key in data:
            metrics[key] = data[key]

    #  pLDDT  (stored as per-atom or per-residue array) 
    plddt_key = next(
        (k for k in ("atom_plddts", "plddt", "per_residue_plddt") if k in data),
        None
    )
    if plddt_key:
        arr = np.array(data[plddt_key])
        metrics["mean_plddt"]   = float(arr.mean())
        metrics["median_plddt"] = float(np.median(arr))
        metrics["min_plddt"]    = float(arr.min())
        metrics["max_plddt"]    = float(arr.max())
        # % high-confidence residues
        metrics["pct_plddt_gt70"] = float((arr > 70).mean() * 100)
        metrics["pct_plddt_gt90"] = float((arr > 90).mean() * 100)

    #  Per-chain pLDDT (chain_plddts dict) 
    if "chain_plddts" in data:
        for chain, vals in data["chain_plddts"].items():
            arr = np.array(vals)
            metrics[f"chain_{chain}_mean_plddt"] = float(arr.mean())

    #  PAE-derived metrics 
    pae_key = next(
        (k for k in ("predicted_aligned_error", "pae") if k in data),
        None
    )
    if pae_key:
        pae = np.array(data[pae_key])
        metrics["mean_pae"]   = float(pae.mean())
        metrics["median_pae"] = float(np.median(pae))
        metrics["max_pae"]    = float(pae.max())

    #  Chain-pair PAE (interface quality) 
    if "chain_pair_pae_min" in data:
        for pair, val in data["chain_pair_pae_min"].items():
            metrics[f"pair_pae_min_{pair}"] = float(val)

    if "chain_pair_iptm" in data:
        for pair, val in data["chain_pair_iptm"].items():
            metrics[f"pair_iptm_{pair}"] = float(val)

    return metrics

def parse_name_and_model(json_path: str):
    """
    confidence_amber-crane-onyx_model_0.json
    → run_name = 'amber-crane-onyx', model_idx = 0
    """
    stem = Path(json_path).stem                        # confidence_..._model_N
    stem = stem.removeprefix("confidence_")            # ..._model_N
    parts = stem.rsplit("_model_", 1)
    run_name  = parts[0] if len(parts) == 2 else stem
    model_idx = int(parts[1]) if len(parts) == 2 else -1
    return run_name, model_idx

# Main
def collect_metrics(results_dir: str) -> pd.DataFrame:
    pattern = os.path.join(results_dir, "**", "confidence_*.json")
    json_files = sorted(glob.glob(pattern, recursive=True))

    if not json_files:
        raise FileNotFoundError(
            f"No confidence_*.json files found under: {results_dir}\n"
            "Make sure you've unpacked the .tar.gz archives first."
        )

    print(f"Found {len(json_files)} confidence JSON files.\n")
    rows = []

    for jf in json_files:
        run_name, model_idx = parse_name_and_model(jf)
        try:
            data    = load_confidence(jf)
            metrics = extract_metrics(data)
            row = {"run_name": run_name, "model": model_idx, "json_path": jf}
            row.update(metrics)
            rows.append(row)
        except Exception as e:
            print(f"  [!] Skipping {jf}: {e}")

    df = pd.DataFrame(rows).sort_values(["run_name", "model"]).reset_index(drop=True)
    return df

def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Best model per run (by ranking_score, falling back to mean_plddt)."""
    score_col = "ranking_score" if "ranking_score" in df.columns else "mean_plddt"
    best = (
        df.sort_values(score_col, ascending=False)
          .groupby("run_name", sort=False)
          .first()
          .reset_index()
    )
    best.insert(1, "best_model_by", score_col)
    return best

# CLI
def main():
    parser = argparse.ArgumentParser(description="Extract SeedFold confidence metrics.")
    parser.add_argument("--results_dir", default=".", help="Root dir of unpacked results")
    parser.add_argument("--out",         default="seedfold_metrics.csv",
                        help="Output CSV path  (default: seedfold_metrics.csv)")
    parser.add_argument("--summary",     default="seedfold_best_models.csv",
                        help="Best-per-run summary CSV")
    args = parser.parse_args()

    df = collect_metrics(args.results_dir)

    #  Print a quick preview 
    preview_cols = ["run_name", "model", "ranking_score", "mean_plddt",
                    "ptm", "iptm", "mean_pae"]
    show_cols = [c for c in preview_cols if c in df.columns]
    print(" All models ")
    print(df[show_cols].to_string(index=False))

    summary = summarise(df)
    print("\n Best model per run ")
    print(summary[show_cols].to_string(index=False))

    #  Save 
    df.to_csv(args.out, index=False)
    summary.to_csv(args.summary, index=False)
    print(f"\n[OK] Full table  → {args.out}")
    print(f"[OK] Best models → {args.summary}")

if __name__ == "__main__":
    main()