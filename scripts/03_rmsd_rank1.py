"""
STEP 3 — Rank-1 filtering + pairwise binder RMSD
==================================================
Run this after step2. It:
  1. Filters to rank-1 (best) model per binder - method
  2. Builds all method pairs for the same binder
  3. Computes binder CA RMSD (binder-on-binder alignment) for every pair
  4. Computes binder CA RMSD aligned on the target binding pocket (5 Å interface)

Usage:
    python step3_rank1_and_rmsd.py \
        --records ./data/prediction_records_with_metrics.parquet \
        --out-dir ./data

Outputs:
    data/rank1_structures.parquet
    data/pairwise_binder_ca_rmsd.parquet
    data/pairwise_pocket_align_rmsd.parquet
"""

import argparse
import os
import tempfile
from itertools import combinations, permutations
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, Superimposer, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from joblib import Parallel, delayed


# ─── Structure loading ────────────────────────────────────────────────────────

def _cif_to_tmp_pdb(cif_path):
    doc = gemmi.cif.read_file(str(cif_path))
    st = gemmi.make_structure_from_block(doc.sole_block())
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    st.write_pdb(tmp.name)
    tmp.close()
    return tmp.name

def load_structure(path, file_type=None, sid=None):
    path = Path(path)
    sid = sid or path.stem
    ft = (file_type or path.suffix).lower().lstrip(".")
    if ft == "pdb":
        return PDBParser(QUIET=True).get_structure(sid, str(path))
    if ft in {"cif", "mmcif"}:
        try:
            return MMCIFParser(QUIET=True).get_structure(sid, str(path))
        except Exception:
            tmp = _cif_to_tmp_pdb(path)
            try:
                return PDBParser(QUIET=True).get_structure(sid, tmp)
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)
    raise ValueError(f"Unsupported file type: {ft}")

def get_chain(structure, chain_id):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain
    available = [c.id for m in structure for c in m]
    raise KeyError(f"Chain '{chain_id}' not found. Available: {sorted(set(available))}")

def get_protein_residues(chain):
    return [r for r in chain if r.id[0] == " " and is_aa(r, standard=False)]

def get_ca_atoms(chain):
    return [r["CA"] for r in get_protein_residues(chain) if "CA" in r]


# ─── Rank-1 table builder ─────────────────────────────────────────────────────

def build_rank1_df(df):
    required = ["binder_sequence", "method", "path", "rank",
                "file_type", "binder_chain_id", "target_chain_id"]
    work = df[df["rank"] == 1].copy()
    work = work.dropna(subset=required)

    dupes = work.duplicated(subset=["binder_sequence", "method"], keep=False)
    if dupes.any():
        print(f"  WARNING: {dupes.sum()} duplicate rank-1 rows — keeping first.")
        work = work.drop_duplicates(subset=["binder_sequence", "method"], keep="first")

    return work.sort_values(["binder_sequence", "method"]).reset_index(drop=True)


# ─── Pair builder ─────────────────────────────────────────────────────────────

def build_pairs(rank1_df, directed=False):
    rows = []
    for binder_seq, group in rank1_df.groupby("binder_sequence", sort=False):
        if len(group) < 2:
            continue
        records = group.to_dict("records")
        pair_iter = permutations(records, 2) if directed else combinations(records, 2)
        for r1, r2 in pair_iter:
            m1_col = "reference_method" if directed else "method_1"
            m2_col = "query_method"     if directed else "method_2"
            rows.append({
                "binder_sequence": binder_seq,
                m1_col: r1["method"], m2_col: r2["method"],
                "path_1": r1["path"], "path_2": r2["path"],
                "file_type_1": r1["file_type"], "file_type_2": r2["file_type"],
                "binder_chain_1": r1["binder_chain_id"],
                "binder_chain_2": r2["binder_chain_id"],
                "target_chain_1": r1["target_chain_id"],
                "target_chain_2": r2["target_chain_id"],
            })
    sort_cols = (["binder_sequence", "reference_method", "query_method"] if directed
                 else ["binder_sequence", "method_1", "method_2"])
    return pd.DataFrame(rows).sort_values(sort_cols).reset_index(drop=True)


# ─── Binder CA RMSD (aligned on binder) ──────────────────────────────────────

def binder_ca_rmsd_for_row(row):
    try:
        s1 = load_structure(row["path_1"], row.get("file_type_1"), "ref")
        s2 = load_structure(row["path_2"], row.get("file_type_2"), "mov")
        ca1 = get_ca_atoms(get_chain(s1, row["binder_chain_1"]))
        ca2 = get_ca_atoms(get_chain(s2, row["binder_chain_2"]))
        if not ca1 or not ca2:
            return {"binder_ca_rmsd": np.nan, "n_ca": 0, "status": "failed",
                    "error": "No CA atoms found"}
        if len(ca1) != len(ca2):
            return {"binder_ca_rmsd": np.nan, "n_ca": min(len(ca1), len(ca2)),
                    "status": "failed",
                    "error": f"CA count mismatch: {len(ca1)} vs {len(ca2)}"}
        sup = Superimposer()
        sup.set_atoms(ca1, ca2)
        return {"binder_ca_rmsd": float(sup.rms), "n_ca": len(ca1),
                "status": "ok", "error": None}
    except Exception as e:
        return {"binder_ca_rmsd": np.nan, "n_ca": 0, "status": "failed", "error": str(e)}


# ─── Pocket-aligned RMSD ─────────────────────────────────────────────────────

def residue_key(res):
    return (int(res.id[1]), res.id[2].strip())

def get_interface_keys(target_chain, binder_chain, cutoff):
    binder_atoms = [a for r in get_protein_residues(binder_chain) for a in r.get_atoms()]
    if not binder_atoms:
        return []
    ns = NeighborSearch(binder_atoms)
    keys = []
    for res in get_protein_residues(target_chain):
        for atom in res.get_atoms():
            if ns.search(atom.coord, cutoff):
                keys.append(residue_key(res))
                break
    return keys

def get_matched_ca(chain1, chain2, keys):
    d1 = {residue_key(r): r for r in get_protein_residues(chain1)}
    d2 = {residue_key(r): r for r in get_protein_residues(chain2)}
    a1, a2 = [], []
    for k in keys:
        r1, r2 = d1.get(k), d2.get(k)
        if r1 and r2 and "CA" in r1 and "CA" in r2:
            a1.append(r1["CA"]); a2.append(r2["CA"])
    return a1, a2

def calc_rmsd_with_transform(atoms1, atoms2, rot, tran):
    c1 = np.array([a.coord for a in atoms1])
    c2 = np.array([a.coord for a in atoms2])
    c2t = c2 @ rot + tran
    return float(np.sqrt(np.mean(np.sum((c1 - c2t) ** 2, axis=1))))

def pocket_rmsd_for_row(row, cutoff=5.0):
    try:
        s1 = load_structure(row["path_1"], row.get("file_type_1"), "ref")
        s2 = load_structure(row["path_2"], row.get("file_type_2"), "mov")
        b1 = get_chain(s1, row["binder_chain_1"])
        b2 = get_chain(s2, row["binder_chain_2"])
        t1 = get_chain(s1, row["target_chain_1"])
        t2 = get_chain(s2, row["target_chain_2"])

        keys = get_interface_keys(t1, b1, cutoff)
        if len(keys) < 3:
            return {"pocket_binder_ca_rmsd": np.nan, "n_pocket_residues": len(keys),
                    "status_pocket": "failed",
                    "error_pocket": f"Too few interface residues: {len(keys)}"}

        shell_ca1, shell_ca2 = get_matched_ca(t1, t2, keys)
        if len(shell_ca1) < 3:
            return {"pocket_binder_ca_rmsd": np.nan, "n_pocket_residues": len(keys),
                    "status_pocket": "failed",
                    "error_pocket": f"Too few matched shell CA: {len(shell_ca1)}"}

        sup = Superimposer()
        sup.set_atoms(shell_ca1, shell_ca2)
        rot, tran = sup.rotran

        bca1, bca2 = get_ca_atoms(b1), get_ca_atoms(b2)
        if len(bca1) != len(bca2):
            return {"pocket_binder_ca_rmsd": np.nan, "n_pocket_residues": len(keys),
                    "status_pocket": "failed",
                    "error_pocket": f"Binder CA mismatch: {len(bca1)} vs {len(bca2)}"}

        rmsd = calc_rmsd_with_transform(bca1, bca2, rot, tran)
        return {"pocket_binder_ca_rmsd": rmsd, "n_pocket_residues": len(keys),
                "n_shell_ca_aligned": len(shell_ca1),
                "status_pocket": "ok", "error_pocket": None}
    except Exception as e:
        return {"pocket_binder_ca_rmsd": np.nan, "n_pocket_residues": np.nan,
                "status_pocket": "failed", "error_pocket": str(e)}


# ─── CLI ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", required=True,
                        help="prediction_records_with_metrics.parquet from step 2")
    parser.add_argument("--out-dir", default="./data")
    parser.add_argument("--cutoff", type=float, default=5.0,
                        help="Interface cutoff in Angstrom for pocket alignment (default: 5.0)")
    parser.add_argument("--n-jobs", type=int, default=-2)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading records...")
    df = pd.read_parquet(args.records)

    # ── Step 3a: rank-1 table ─────────────────────────────────────────────────
    print("\nBuilding rank-1 table...")
    rank1 = build_rank1_df(df)
    rank1_path = out_dir / "rank1_structures.parquet"
    rank1.to_parquet(rank1_path, index=False)
    print(f"  {len(rank1)} rank-1 rows ({rank1['binder_sequence'].nunique()} unique binders)")
    print(f"  Saved → {rank1_path}")

    # ── Step 3b: triangular pairs + binder CA RMSD ───────────────────────────
    print("\nBuilding triangular method pairs...")
    pairs = build_pairs(rank1, directed=False)
    print(f"  {len(pairs)} pairs")

    print("Computing binder CA RMSD (binder-on-binder alignment)...")
    results = Parallel(n_jobs=args.n_jobs, verbose=5)(
        delayed(binder_ca_rmsd_for_row)(row) for _, row in pairs.iterrows()
    )
    rmsd_df = pd.concat([pairs.reset_index(drop=True), pd.DataFrame(results)], axis=1)
    rmsd_path = out_dir / "pairwise_binder_ca_rmsd.parquet"
    rmsd_df.to_parquet(rmsd_path, index=False)
    print(f"  Status: {rmsd_df['status'].value_counts().to_dict()}")
    print(f"  Saved → {rmsd_path}")

    # ── Step 3c: directed pairs + pocket-aligned RMSD ────────────────────────
    print("\nBuilding directed pairs for pocket-aligned RMSD...")
    dir_pairs = build_pairs(rank1, directed=True)
    print(f"  {len(dir_pairs)} directed pairs")

    print(f"Computing pocket-aligned binder RMSD (cutoff={args.cutoff} Å)...")
    pocket_results = Parallel(n_jobs=args.n_jobs, verbose=5)(
        delayed(pocket_rmsd_for_row)(row, args.cutoff)
        for _, row in dir_pairs.iterrows()
    )
    pocket_df = pd.concat([dir_pairs.reset_index(drop=True),
                           pd.DataFrame(pocket_results)], axis=1)
    pocket_path = out_dir / "pairwise_pocket_align_rmsd.parquet"
    pocket_df.to_parquet(pocket_path, index=False)
    print(f"  Status: {pocket_df['status_pocket'].value_counts().to_dict()}")
    print(f"  Saved → {pocket_path}")

    print("\n✓ Step 3 complete. Next: run step4_plots.py")



""" Execute this by using this command in terminal:
python rmsd_rank1_3.py \
    --records ./data/prediction_records_with_metrics.parquet \
    --out-dir ./data """