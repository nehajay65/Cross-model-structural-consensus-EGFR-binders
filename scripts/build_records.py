"""
STEP 1-  Building prediction_records.parquet
It crawls your outputs/ folder, parses every structure file,
extracts binder/target chain info and confidence scores, and saves a parquet.

Usage:
    python build_records.py \
        --root ./outputs \
        --out-dir ./data \
        --target-length 621

Output:
    data/prediction_records.parquet
"""

import argparse
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Any
import re
import json
import pickle
import tempfile
import os

import numpy as np
import pandas as pd
import gemmi
from Bio.PDB import PDBParser, MMCIFParser
from Bio.Data.IUPACData import protein_letters_3to1
from joblib import Parallel, delayed


# data model

@dataclass
class StructureRecord:
    path: Path
    file_type: str
    method: str
    complex_name: str | None = None
    prediction_id: str | None = None
    rank: int | None = None
    target_chain_id: str | None = None
    binder_chain_id: str | None = None
    target_length: int | None = None
    binder_length: int | None = None
    binder_sequence: str | None = None
    target_sequence: str | None = None
    score: float | None = None
    score_name: str | None = None
    json_path: Path | None = None
    pkl_path: Path | None = None
    npz_path: Path | None = None
    parse_success: bool = True
    error_message: str | None = None


#dataset class

class PredictionDataset:

    KNOWN_METHODS = {
        "AF2", "AF3", "Boltz-2", "Chai-1",
        "HelixFold3", "OF2", "OF3", "Protenix", "SeedFold",
    }

    def __init__(self, root_dir, target_length=621):
        self.root_dir = Path(root_dir)
        self.target_length = target_length
        self.records = []

    # method detection 

    def detect_method(self, path):
        parts = Path(path).relative_to(self.root_dir).parts
        return parts[0] if parts[0] in self.KNOWN_METHODS else "unknown"

    #file filtering

    def should_keep(self, path):
        path = Path(path)
        if not path.is_file():
            return False
        if path.suffix.lower() not in {".pdb", ".cif"}:
            return False

        method = self.detect_method(path)
        stem = path.stem.lower()
        parts_lower = [p.lower() for p in path.parts]

        if method == "AF2":
            return stem.startswith("ranked_")

        if method == "HelixFold3":
            return any("rank" in p for p in parts_lower)

        if method == "AF3":
            return path.parent.name.lower().startswith("seed-")

        if method == "OF2":
            return "_relaxed" in stem and "unrelaxed" not in stem
        return True
    #can put the same for OF2

    def find_structure_files(self):
        all_files = list(self.root_dir.rglob("*.pdb")) + list(self.root_dir.rglob("*.cif"))
        return sorted(f for f in all_files if self.should_keep(f))

    # structure loading 

    def _cif_to_tmp_pdb(self, cif_path):
        doc = gemmi.cif.read_file(str(cif_path))
        structure = gemmi.make_structure_from_block(doc.sole_block())
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        structure.write_pdb(tmp.name)
        tmp.close()
        return tmp.name

    def load_structure(self, path):
        path = Path(path)
        method = self.detect_method(path)
        if path.suffix.lower() == ".pdb":
            return PDBParser(QUIET=True).get_structure(path.stem, str(path))
        elif path.suffix.lower() == ".cif":
            if method == "OF3":
                tmp = self._cif_to_tmp_pdb(path)
                try:
                    return PDBParser(QUIET=True).get_structure(path.stem, tmp)
                finally:
                    os.unlink(tmp)
            return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
        raise ValueError(f"Unsupported file: {path}")

    # chain identification

    def identify_chains(self, structure):
        target_chain = binder_chain = None
        for chain in structure.get_chains():
            residues = [r for r in chain.get_residues() if r.id[0] == " "]
            if len(residues) == self.target_length:
                target_chain = chain
            else:
                binder_chain = chain
        return target_chain, binder_chain

    def extract_sequence(self, chain):
        seq = []
        for res in chain.get_residues():
            if res.id[0] != " ":
                continue
            seq.append(protein_letters_3to1.get(res.resname.title(), "X"))
        return "".join(seq)

    #score extraction (per model specifically)

    def _safe_float(self, v):
        try:
            return float(v)
        except Exception:
            return None

    def _load_json(self, p):
        with open(p) as f:
            return json.load(f)

    def extract_ranking_info(self, path, method):
        path = Path(path)
        stem = path.stem

        if method == "AF2":
            ranking_json = path.parent / "ranking_debug.json"
            score = pred_id = pkl_path = None
            m = re.match(r"ranked_(\d+)$", stem.lower())
            if m and ranking_json.exists():
                data = self._load_json(ranking_json)
                idx = int(m.group(1))
                model_name = data.get("order", [None] * (idx + 1))[idx]
                if model_name:
                    pred_id = model_name
                    score = self._safe_float((data.get("iptm+ptm") or {}).get(model_name))
                    candidate = path.parent / f"result_{model_name}.pkl"
                    if candidate.exists():
                        pkl_path = candidate
            return {"prediction_id": pred_id or stem, "score": score,
                    "score_name": "iptm+ptm" if score is not None else None,
                    "json_path": ranking_json if ranking_json.exists() else None,
                    "pkl_path": pkl_path, "npz_path": None}

        if method == "AF3":
            summary = next((j for j in sorted(path.parent.glob("*.json"))
                            if "summary" in j.name.lower()), None)
            score = self._safe_float(self._load_json(summary).get("ranking_score")) if summary else None
            return {"prediction_id": path.parent.name, "score": score,
                    "score_name": "ranking_score" if score is not None else None,
                    "json_path": summary, "pkl_path": None, "npz_path": None}

        if method == "Boltz-2":
            jsons = sorted(path.parent.glob("*.json"))
            jp = jsons[0] if jsons else None
            score = self._safe_float(self._load_json(jp).get("confidence_score")) if jp else None
            return {"prediction_id": stem, "score": score,
                    "score_name": "confidence_score" if score is not None else None,
                    "json_path": jp, "pkl_path": None, "npz_path": None}

        if method == "Chai-1":
            suffix = stem.replace("pred.", "")
            npz = path.with_name(f"scores.{suffix}.npz")
            score = None
            if npz.exists():
                data = np.load(npz, allow_pickle=True)
                if "aggregate_score" in data:
                    score = self._safe_float(np.array(data["aggregate_score"]).item())
            return {"prediction_id": stem, "score": score,
                    "score_name": "aggregate_score" if score is not None else None,
                    "json_path": None, "pkl_path": None,
                    "npz_path": npz if npz.exists() else None}

        if method == "HelixFold3":
            all_results = path.with_name("all_results.json")
            score = None
            if all_results.exists():
                data = self._load_json(all_results)
                score = self._safe_float(data.get("ranking_confidence"))
            return {"prediction_id": path.parent.name, "score": score,
                "score_name": "ranking_confidence" if score is not None else None,
                "json_path": all_results if all_results.exists() else None,
                "pkl_path": None, "npz_path": None}

        if method == "OF3":
            jp = None
            if stem.endswith("_model"):
                jp = path.with_name(f"{stem[:-6]}_confidences_aggregated.json")
            score = self._safe_float(self._load_json(jp).get("sample_ranking_score")) if (jp and jp.exists()) else None
            return {"prediction_id": stem, "score": score,
                    "score_name": "sample_ranking_score" if score is not None else None,
                    "json_path": jp if (jp and jp.exists()) else None,
                    "pkl_path": None, "npz_path": None}

        if method == "Protenix":
            m = re.match(r"(.+)_sample_(\d+)$", stem)
            jp = None
            if m:
                prefix, sid = m.groups()
                jp = path.with_name(f"{prefix}_summary_confidence_sample_{sid}.json")
            score = self._safe_float(self._load_json(jp).get("ranking_score")) if (jp and jp.exists()) else None
            return {"prediction_id": stem, "score": score,
                    "score_name": "ranking_score" if score is not None else None,
                    "json_path": jp if (jp and jp.exists()) else None,
                    "pkl_path": None, "npz_path": None}

        if method == "SeedFold":
            jp = path.with_name(f"confidence_{stem}.json")
            score = self._safe_float(self._load_json(jp).get("confidence_score")) if jp.exists() else None
            return {"prediction_id": stem, "score": score,
                    "score_name": "confidence_score" if score is not None else None,
                    "json_path": jp if jp.exists() else None,
                    "pkl_path": None, "npz_path": None}

        if method == "OF2":
            return {"prediction_id": stem, "score": None, "score_name": None, "json_path": None, "pkl_path": None, "npz_path": None}

        # OF2 + unknown
        return {"prediction_id": stem, "score": None, "score_name": None,
                "json_path": None, "pkl_path": None, "npz_path": None}

    

    def extract_complex_name(self, path, method):
        if method == "OF2":
             # amber_crane_onyx-EGFR_model_3... → amber-crane-onyx
            stem = path.stem  # amber_crane_onyx-EGFR_model_3_multimer_v3_relaxed
            binder_part = stem.split("-EGFR")[0]  # amber_crane_onyx
            return binder_part.replace("_", "-")  # amber-crane-onyx

        if method == "Boltz-2":
            return Path(path).parent.name  # amber-crane-onyx
        parts = Path(path).relative_to(self.root_dir).parts
        return parts[1] if len(parts) >= 2 else Path(path).parent.name

    # record builder

    def extract_record(self, path):
        path = Path(path)
        method = "unknown"
        try:
            method = self.detect_method(path)
            structure = self.load_structure(path)
            target_chain, binder_chain = self.identify_chains(structure)

            target_seq = binder_seq = None
            target_id = binder_id = None
            target_len = binder_len = None

            if target_chain:
                target_id = target_chain.id
                target_seq = self.extract_sequence(target_chain)
                target_len = len(target_seq)
            if binder_chain:
                binder_id = binder_chain.id
                binder_seq = self.extract_sequence(binder_chain)
                binder_len = len(binder_seq)

            ranking = self.extract_ranking_info(path, method)

            return StructureRecord(
                path=path, file_type=path.suffix.lower(), method=method,
                complex_name=self.extract_complex_name(path, method),
                prediction_id=ranking["prediction_id"],
                score=ranking["score"], score_name=ranking["score_name"],
                json_path=ranking["json_path"], pkl_path=ranking["pkl_path"],
                npz_path=ranking["npz_path"],
                target_chain_id=target_id, binder_chain_id=binder_id,
                target_length=target_len, binder_length=binder_len,
                target_sequence=target_seq, binder_sequence=binder_seq,
                parse_success=True, error_message=None,
            )
        except Exception as e:
            return StructureRecord(
                path=path, file_type=path.suffix.lower(), method=method,
                parse_success=False, error_message=str(e),
            )

    # rank assignment

    def assign_ranks(self):
        groups = {}
        for rec in self.records:
            groups.setdefault((rec.method, rec.complex_name), []).append(rec)

        for (method, _), group in groups.items():
            if method == "AF2":
                def af2_key(r):
                    m = re.match(r"ranked_(\d+)$", Path(r.path).stem.lower())
                    return int(m.group(1)) if m else 999
                group.sort(key=af2_key)
                for i, r in enumerate(group, 1):
                    r.rank = i
                continue
            if method == "OF2":          # add this
                for r in group:
                   r.rank = 1
                continue

            scored = sorted([r for r in group if r.score is not None],
                    key=lambda r: r.score, reverse=True)
            if len(scored) == 1:
                scored[0].rank = 1
            elif len(scored) > 1:
                for i, r in enumerate(scored, 1):
                    r.rank = i
            if not scored and len(group) == 1:
                group[0].rank = 1

    # run + save

    def run(self, out_path, n_jobs=-2):
        files = self.find_structure_files()
        print(f"Found {len(files)} structure files across all methods.")

        self.records = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(self.extract_record)(f) for f in files
        )

        self.assign_ranks()

        df = pd.DataFrame([asdict(r) for r in self.records])
        for col in df.columns:
            df[col] = df[col].apply(lambda x: str(x) if isinstance(x, Path) else x)

        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(out_path, index=False)

        ok = df["parse_success"].sum()
        fail = (~df["parse_success"]).sum()
        print(f"\nDone. {ok} records parsed OK, {fail} failed.")
        print(f"Saved → {out_path}")
        print("\nRecords per method:")
        print(df.groupby("method").size().to_string())
        return df


# CLI

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, help="Root folder containing method subfolders (e.g. ./outputs)")
    parser.add_argument("--out-dir", default="./data", help="Where to save outputs (default: ./data)")
    parser.add_argument("--target-length", type=int, default=621, help="Target protein residue count (default: 621)")
    parser.add_argument("--n-jobs", type=int, default=-2, help="Parallel jobs (default: -2 = all but one CPU)")
    args = parser.parse_args()

    out_path = Path(args.out_dir) / "prediction_records.parquet"

    dataset = PredictionDataset(
        root_dir=args.root,
        target_length=args.target_length,
    )
    dataset.run(out_path=out_path, n_jobs=args.n_jobs)





#run it by using the command line, for example:
""" python build_records.py \
    --root  /home/nehajay/project_new/a_predictions \ #path to your directory where your predictions lies
    --out-dir ./data \
    --target-length 621 (length of your target protein)"""