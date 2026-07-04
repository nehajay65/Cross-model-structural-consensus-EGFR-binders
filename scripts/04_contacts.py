"""
Contact interface script — rank-1 binder/target contacts
Run after step2 (prediction_records_with_metrics.parquet)
"""

from pathlib import Path
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import os
import tempfile

import gemmi
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa


# ─── Structure loading ────────────────────────────────────────────────────────

def _cif_to_tmp_pdb(cif_path):
    doc = gemmi.cif.read_file(str(cif_path))
    st = gemmi.make_structure_from_block(doc.sole_block())
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    st.write_pdb(tmp.name)
    tmp.close()
    return tmp.name

def load_structure_biopython(path, file_type, structure_id="s"):
    path = Path(path)
    ft = file_type.lower().lstrip(".")
    if ft == "pdb":
        return PDBParser(QUIET=True).get_structure(structure_id, str(path))
    if ft in {"cif", "mmcif"}:
        try:
            return MMCIFParser(QUIET=True).get_structure(structure_id, str(path))
        except Exception:
            tmp = _cif_to_tmp_pdb(path)
            try:
                return PDBParser(QUIET=True).get_structure(structure_id, tmp)
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)
    raise ValueError(f"Unsupported file type: {ft}")

def get_chain_by_id(structure, chain_id):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain
    raise KeyError(f"Chain '{chain_id}' not found.")

def get_protein_residues(chain):
    return [r for r in chain if r.id[0] == " " and is_aa(r, standard=False)]

def residue_key(res):
    return (int(res.id[1]), res.id[2].strip())


# ─── Rank-1 filter ────────────────────────────────────────────────────────────

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


# ─── Contact extraction ───────────────────────────────────────────────────────

def make_residue_record(residue, prefix):
    hetflag, resseq, icode = residue.id
    return {
        f"{prefix}_chain_id": residue.get_parent().id,
        f"{prefix}_resname": residue.get_resname(),
        f"{prefix}_resseq": int(resseq),
        f"{prefix}_icode": "" if icode is None or icode == " " else str(icode),
        f"{prefix}_hetflag": hetflag,
    }

def extract_interface_contacts(structure, binder_chain_id, target_chain_id, cutoff=5.0):
    binder_chain = get_chain_by_id(structure, binder_chain_id)
    target_chain = get_chain_by_id(structure, target_chain_id)
    binder_residues = get_protein_residues(binder_chain)
    target_residues = get_protein_residues(target_chain)

    target_atoms = [atom for res in target_residues for atom in res]
    if not target_atoms:
        raise ValueError("No target atoms found.")

    ns = NeighborSearch(target_atoms)
    target_contacts, binder_contacts, contact_pairs = {}, {}, {}

    for binder_res in binder_residues:
        if binder_res.id[0] != " " or not is_aa(binder_res, standard=False):
            continue
        binder_key = residue_key(binder_res)
        for atom in binder_res:
            for target_res in ns.search(atom.coord, cutoff, level="R"):
                if target_res.get_parent().id != target_chain_id:
                    continue
                if target_res.id[0] != " " or not is_aa(target_res, standard=False):
                    continue
                target_key = residue_key(target_res)
                target_contacts[target_key] = target_res
                binder_contacts[binder_key] = binder_res
                contact_pairs[(target_key, binder_key)] = (target_res, binder_res)

    return target_contacts, binder_contacts, contact_pairs

def compute_contacts_for_row(row, cutoff=5.0):
    base = {
        "binder_sequence": row["binder_sequence"],
        "method": row["method"],
        "path": row["path"],
        "file_type": row["file_type"],
        "binder_chain_id": row["binder_chain_id"],
        "target_chain_id": row["target_chain_id"],
        "score": row.get("score"),
        "score_name": row.get("score_name"),
        "contact_cutoff": cutoff,
    }
    try:
        structure = load_structure_biopython(row["path"], row["file_type"], "contact_struct")
        target_contacts, binder_contacts, contact_pairs = extract_interface_contacts(
            structure, row["binder_chain_id"], row["target_chain_id"], cutoff
        )
        summary = {**base,
            "n_target_contact_residues": len(target_contacts),
            "n_binder_contact_residues": len(binder_contacts),
            "n_residue_contact_pairs": len(contact_pairs),
            "status_contacts": "ok", "error_contacts": None}
        target_rows = [{**base, **make_residue_record(r, "target")}
                       for _, r in sorted(target_contacts.items())]
        binder_rows = [{**base, **make_residue_record(r, "binder")}
                       for _, r in sorted(binder_contacts.items())]
        pair_rows   = [{**base, **make_residue_record(tr, "target"), **make_residue_record(br, "binder")}
                       for _, (tr, br) in sorted(contact_pairs.items())]
        return summary, target_rows, binder_rows, pair_rows
    except Exception as e:
        summary = {**base,
            "n_target_contact_residues": 0, "n_binder_contact_residues": 0,
            "n_residue_contact_pairs": 0,
            "status_contacts": "failed", "error_contacts": str(e)}
        return summary, [], [], []


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── CHANGE THESE PATHS ────────────────────────────────────────────────────
    INPUT_PATH  = Path("./data/prediction_records_with_metrics.parquet")
    OUTPUT_DIR  = Path("./data/contacts")
    CUTOFF      = 5.0
    # ─────────────────────────────────────────────────────────────────────────

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading records...")
    df = pd.read_parquet(INPUT_PATH)

    print("Filtering to rank-1...")
    rank1_df = build_rank1_df(df)
    print(f"  {len(rank1_df)} rank-1 rows | {rank1_df['binder_sequence'].nunique()} unique binders")
    print(f"  Methods: {sorted(rank1_df['method'].unique())}")

    rank1_df.to_parquet(OUTPUT_DIR / "rank1_structures.parquet", index=False)

    print(f"\nComputing contacts at {CUTOFF}Å cutoff...")
    results = Parallel(n_jobs=-2, verbose=10)(
        delayed(compute_contacts_for_row)(row, cutoff=CUTOFF)
        for _, row in rank1_df.iterrows()
    )

    summary_df     = pd.DataFrame([x[0] for x in results])
    target_long_df = pd.DataFrame([r for x in results for r in x[1]])
    binder_long_df = pd.DataFrame([r for x in results for r in x[2]])
    pair_long_df   = pd.DataFrame([r for x in results for r in x[3]])

    summary_df.to_parquet(    OUTPUT_DIR / "rank1_contact_summary_5A.parquet",      index=False)
    target_long_df.to_parquet(OUTPUT_DIR / "rank1_target_contacts_long_5A.parquet", index=False)
    binder_long_df.to_parquet(OUTPUT_DIR / "rank1_binder_contacts_long_5A.parquet", index=False)
    pair_long_df.to_parquet(  OUTPUT_DIR / "rank1_contact_pairs_long_5A.parquet",   index=False)

    print("\nContact status:")
    print(summary_df["status_contacts"].value_counts(dropna=False))
    print(f"\nSaved all parquets to: {OUTPUT_DIR}")