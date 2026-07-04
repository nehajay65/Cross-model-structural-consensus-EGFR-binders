"""
STEP 2 - Extract native confidence metrics
Run this after step1. It reads prediction_records.parquet, pulls per-method
confidence numbers (pLDDT, iPTM, PTM, ranking scores) from the JSON/pkl/npz
files next to each structure, and saves an enriched parquet.

Usage:
    python step2_extract_metrics.py \
        --records ./data/prediction_records.parquet \
        --out-dir ./data

Output:
    data/prediction_records_with_metrics.parquet
"""

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

class NativeMetricExtractor:

    def _safe_float(self, value):
        try:
            if isinstance(value, (list, tuple, dict)):
                return None
            if isinstance(value, np.ndarray):
                return float(value.item()) if value.size == 1 else None
            return float(value)
        except Exception:
            return None

    def _load_json(self, p):
        with open(p) as f:
            return json.load(f)

    def _load_pkl(self, p):
        with open(p, "rb") as f:
            return pickle.load(f)

    def _load_npz(self, p):
        return dict(np.load(p, allow_pickle=True))

    def _flatten(self, d, parent=""):
        out = {}
        if isinstance(d, dict):
            for k, v in d.items():
                key = f"{parent}.{k}" if parent else str(k)
                if isinstance(v, dict):
                    out.update(self._flatten(v, key))
                else:
                    out[key] = v
        return out

    def _pick_numeric(self, flat, keys, prefix):
        out = {}
        for k in keys:
            if k in flat:
                v = self._safe_float(flat[k])
                if v is not None:
                    out[f"{prefix}_{k}"] = v
        return out

    def _all_numeric(self, flat, prefix, exclude=()):
        out = {}
        for k, v in flat.items():
            if any(x in k.lower() for x in exclude):
                continue
            val = self._safe_float(v)
            if val is not None:
                out[f"{prefix}_{k.replace('.', '__')}"] = val
        return out

    #  per-method extractors 

    def _af2(self, row):
        out = {}
        pkl_path = row.get("pkl_path")
        if pd.notna(pkl_path) and Path(pkl_path).exists():
            data = self._load_pkl(Path(pkl_path))
            if "plddt" in data:
                plddt = np.asarray(data["plddt"])
                out["af2_plddt_mean"]   = float(np.mean(plddt))
                out["af2_plddt_median"] = float(np.median(plddt))
                out["af2_plddt_min"]    = float(np.min(plddt))
                out["af2_plddt_max"]    = float(np.max(plddt))
            for k in ["ptm", "iptm"]:
                v = self._safe_float(data.get(k))
                if v is not None:
                    out[f"af2_{k}"] = v
        return out

    def _af3(self, row):
        out = {}
        jp = row.get("json_path")
        if pd.notna(jp) and Path(jp).exists():
            flat = self._flatten(self._load_json(Path(jp)))
            out.update(self._pick_numeric(flat,
                ["ranking_score","ptm","iptm","fraction_disordered","has_clash"], "af3"))
            for k, v in self._all_numeric(flat, "af3",
                    exclude=["seed","sample","token","atom","residue","chain"]).items():
                out.setdefault(k, v)
        return out

    def _helixfold3(self, row):
        out = {}
        jp = row.get("json_path")
        if pd.notna(jp) and Path(jp).exists():
            flat = self._flatten(self._load_json(Path(jp)))
            out.update(self._pick_numeric(flat, ["ranking_score","ptm","iptm","plddt","mean_plddt"], "helixfold3"))
        return out

    def _chai1(self, row):
        out = {}
        npz = row.get("npz_path")
        if pd.notna(npz) and Path(npz).exists():
            data = self._load_npz(Path(npz))
            for k, v in data.items():
                val = self._safe_float(v)
                if val is not None:
                    out[f"chai1_{k}"] = val
        return out

    def _boltz2(self, row):
        out = {}
        candidates = []
        jp = row.get("json_path")
        if pd.notna(jp) and Path(jp).exists():
            candidates.append(Path(jp))
        path = row.get("path")
        if pd.notna(path):
            candidates.extend(sorted(Path(path).parent.glob("*.json")))
        seen = set()
        for p in candidates:
            sp = str(p)
            if sp in seen:
                continue
            seen.add(sp)
            try:
                flat = self._flatten(self._load_json(p))
                numeric = self._all_numeric(flat, "boltz2",
                    exclude=["seed","sample","token","atom","residue","chain"])
                if numeric:
                    out.update(numeric)
                    break
            except Exception:
                continue
        return out

    def _seedfold(self, row):
        jp = row.get("json_path")
        if pd.notna(jp) and Path(jp).exists():
            return self._all_numeric(self._flatten(self._load_json(Path(jp))), "seedfold")
        return {}

    def _protenix(self, row):
        jp = row.get("json_path")
        if pd.notna(jp) and Path(jp).exists():
            return self._all_numeric(self._flatten(self._load_json(Path(jp))), "protenix")
        return {}

    def _of3(self, row):
        jp = row.get("json_path")
        if pd.notna(jp) and Path(jp).exists():
            return self._all_numeric(self._flatten(self._load_json(Path(jp))), "of3")
        return {}

    #  dispatcher 

    def extract_row(self, row):
        dispatch = {
            "AF2": self._af2,
            "AF3": self._af3,
            "HelixFold3": self._helixfold3,
            "Chai-1": self._chai1,
            "Boltz-2": self._boltz2,
            "SeedFold": self._seedfold,
            "Protenix": self._protenix,
            "OF3": self._of3,
            "OF2": lambda _: {},
        }
        return dispatch.get(row["method"], lambda _: {})(row)

# CLI
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", required=True,
                        help="Path to prediction_records.parquet from step 1")
    parser.add_argument("--out-dir", default="./data")
    args = parser.parse_args()

    out_path = Path(args.out_dir) / "prediction_records_with_metrics.parquet"
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    print("Loading records...")
    df = pd.read_parquet(args.records)
    print(f"  {len(df)} rows, {df['method'].nunique()} methods")

    print("Extracting native metrics (this may take a few minutes)...")
    extractor = NativeMetricExtractor()
    metrics = df.apply(extractor.extract_row, axis=1, result_type="expand")

    df_out = pd.concat([df, metrics], axis=1)
    df_out.to_parquet(out_path, index=False)

    print(f"\nDone. Saved → {out_path}")
    print(f"Added {len(metrics.columns)} metric columns.")

    # Quick coverage summary
    method_cols = [c for c in df_out.columns if any(
        c.startswith(p) for p in
        ["af2_","af3_","helixfold3_","boltz2_","chai1_","seedfold_","protenix_","of3_"]
    )]
    print("\nMetric coverage per method (non-null count):")
    for method in sorted(df_out["method"].unique()):
        sub = df_out[df_out["method"] == method]
        filled = {c: sub[c].notna().sum() for c in method_cols if sub[c].notna().sum() > 0}
        if filled:
            print(f"  {method}: {filled}")

#execut in terminal:
""" python extract_metrics_2.py \
    --records ./data/prediction_records.parquet \
    --out-dir ./data """