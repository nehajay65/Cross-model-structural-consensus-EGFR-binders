"""
STEP 4 — Generate analysis plots

Usage:
    python plots.py \
        --data- ./data \
        --out- ./plots

Outputs (all saved to ./plots/):
    01_records_per_method.png         — how many structures per method
    02_score_distributions.png        — confidence score distributions per method
    03_iptm_per_method.png            — iPTM boxplot across methods (rank-1)
    04_plddt_per_method.png           — pLDDT boxplot across methods (rank-1)
    05_binder_length_dist.png         — binder length distribution
    06_binder_ca_rmsd_heatmap.png     — median binder CA RMSD between method pairs
    07_pocket_rmsd_heatmap.png        — median pocket-aligned RMSD between method pairs
    08_rmsd_by_method_pair.png        — RMSD boxplot per method pair
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)

METHOD_ORDER = ["AF2", "AF3", "Boltz-2", "Chai-1",
                "HelixFold3", "OF2", "OF3", "Protenix", "SeedFold"]

IPTM_COLS = {
    "AF2":       "af2_iptm",
    "AF3":       "af3_iptm",
    "Boltz-2":   "boltz2_iptm",
    "Chai-1":    "chai1_iptm",
    "HelixFold3":"helixfold3_iptm",
    "OF3":       "of3_iptm",
    "Protenix":  "protenix_iptm",
    "SeedFold":  "seedfold_iptm",
}

PLDDT_COLS = {
    "AF2":       "af2_plddt_mean",
    "AF3":       "af3_plddt_mean",
    "Boltz-2":   "boltz2_complex_plddt",
    "HelixFold3":"helixfold3_mean_plddt",
    "OF3":       "of3_avg_plddt",
    "Protenix":  "protenix_plddt",
    "SeedFold":  "seedfold_complex_plddt",
}


def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")


# ── Plot 1: records per method ────────────────────────────────────────────────

def plot_records_per_method(df, out_dir):
    counts = df.groupby("method").size().reindex(METHOD_ORDER).dropna()
    fig, ax = plt.subplots(figsize=(9, 4))
    counts.plot(kind="bar", ax=ax, color="#6e7de8", edgecolor="white", width=0.65)
    ax.set_xlabel("Method"); ax.set_ylabel("Number of structure files")
    ax.set_title("Total structure files per method")
    ax.tick_params(axis="x", rotation=35)
    sns.despine()
    save(fig, out_dir / "01_records_per_method.png")


# ── Plot 2: score distributions ───────────────────────────────────────────────

def plot_score_distributions(df, out_dir):
    ok = df[df["score"].notna() & df["parse_success"]]
    methods_with_scores = [m for m in METHOD_ORDER if m in ok["method"].values]
    if not methods_with_scores:
        print("  Skipping score distributions — no scores found.")
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(data=ok[ok["method"].isin(methods_with_scores)],
                x="method", y="score", order=methods_with_scores,
                color="#d9e6f2", width=0.55, showfliers=False, ax=ax)
    sns.stripplot(data=ok[ok["method"].isin(methods_with_scores)],
                  x="method", y="score", order=methods_with_scores,
                  color="#1f2937", alpha=0.3, size=2.5, jitter=0.2, ax=ax)
    ax.set_xlabel("Method"); ax.set_ylabel("Confidence score (method-native scale)")
    ax.set_title("Confidence score distributions per method (all ranks)")
    ax.tick_params(axis="x", rotation=35)
    sns.despine()
    save(fig, out_dir / "02_score_distributions.png")


# ── Plot 3: iPTM (rank-1) ─────────────────────────────────────────────────────

def plot_iptm(rank1, out_dir):
    rows = []
    for method, col in IPTM_COLS.items():
        sub = rank1[rank1["method"] == method]
        if col not in sub.columns:
            continue
        vals = sub[col].dropna()
        for v in vals:
            rows.append({"method": method, "iptm": v})
    if not rows:
        print("  Skipping iPTM plot — no data found.")
        return
    plot_df = pd.DataFrame(rows)
    order = [m for m in METHOD_ORDER if m in plot_df["method"].values]
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(data=plot_df, x="method", y="iptm", order=order,
                color="#d9e6f2", width=0.55, showfliers=False, ax=ax)
    sns.stripplot(data=plot_df, x="method", y="iptm", order=order,
                  color="#1f2937", alpha=0.35, size=3, jitter=0.2, ax=ax)
    ax.set_xlabel("Method"); ax.set_ylabel("iPTM")
    ax.set_title("iPTM distribution per method (rank-1 predictions)")
    ax.tick_params(axis="x", rotation=35)
    sns.despine()
    save(fig, out_dir / "03_iptm_per_method.png")


# ── Plot 4: pLDDT (rank-1) ───────────────────────────────────────────────────

def plot_plddt(rank1, out_dir):
    rows = []
    for method, col in PLDDT_COLS.items():
        sub = rank1[rank1["method"] == method]
        if col not in sub.columns:
            continue
        vals = sub[col].dropna()
        if method == "Boltz-2":
            vals = vals * 100  # ← add this to convert Boltz-2 pLDDT to 0-100 scale for better comparison   
        for v in vals:
            rows.append({"method": method, "plddt": v})
    if not rows:
        print("  Skipping pLDDT plot — no data found.")
        return
    plot_df = pd.DataFrame(rows)
    order = [m for m in METHOD_ORDER if m in plot_df["method"].values]
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(data=plot_df, x="method", y="plddt", order=order,
                color="#d9eed9", width=0.55, showfliers=False, ax=ax)
    sns.stripplot(data=plot_df, x="method", y="plddt", order=order,
                  color="#1f2937", alpha=0.35, size=3, jitter=0.2, ax=ax)
    ax.set_xlabel("Method"); ax.set_ylabel("pLDDT (mean)")
    ax.set_title("pLDDT distribution per method (rank-1 predictions)")
    ax.tick_params(axis="x", rotation=35)
    sns.despine()
    save(fig, out_dir / "04_plddt_per_method.png")


# ── Plot 5: binder length distribution ───────────────────────────────────────

def plot_binder_lengths(rank1, out_dir):
    binders = rank1.drop_duplicates("binder_sequence")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(binders["binder_length"].dropna(), bins=30, color="#6e7de8",
            edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Binder length (residues)")
    ax.set_ylabel("Count")
    ax.set_title(f"Binder length distribution ({len(binders)} unique binders)")
    sns.despine()
    save(fig, out_dir / "05_binder_length_dist.png")


# ── Plot 6: binder CA RMSD heatmap ───────────────────────────────────────────

def _rmsd_heatmap(rmsd_df, value_col, method1_col, method2_col, title, path):
    ok = rmsd_df[rmsd_df[value_col].notna()]
    methods = sorted(set(ok[method1_col].unique()) | set(ok[method2_col].unique()))
    mat = pd.DataFrame(np.nan, index=methods, columns=methods)
    for _, row in ok.iterrows():
        m1, m2, v = row[method1_col], row[method2_col], row[value_col]
        mat.loc[m1, m2] = v
        mat.loc[m2, m1] = v
    mat_arr = mat.values.copy()
    np.fill_diagonal(mat_arr, 0)
    mat = pd.DataFrame(mat_arr, index=mat.index, columns=mat.columns)
    fig, ax = plt.subplots(figsize=(9, 7))
    sns.heatmap(mat.astype(float), annot=True, fmt=".2f", cmap="YlOrRd",
                linewidths=0.5, ax=ax, vmin=0)
    ax.set_title(title)
    plt.tight_layout()
    save(fig, path)

def plot_rmsd_heatmaps(rmsd_df, pocket_df, out_dir):
    if rmsd_df is not None and "binder_ca_rmsd" in rmsd_df.columns:
        med = rmsd_df[rmsd_df["status"] == "ok"].groupby(
            ["method_1", "method_2"])["binder_ca_rmsd"].median().reset_index()
        _rmsd_heatmap(med, "binder_ca_rmsd", "method_1", "method_2",
                      "Median binder CA RMSD — binder-on-binder alignment (Å)",
                      out_dir / "06_binder_ca_rmsd_heatmap.png")

    if pocket_df is not None and "pocket_binder_ca_rmsd" in pocket_df.columns:
        m1_col = "reference_method" if "reference_method" in pocket_df.columns else "method_1"
        m2_col = "query_method"     if "query_method" in pocket_df.columns else "method_2"
        med = pocket_df[pocket_df["status_pocket"] == "ok"].groupby(
            [m1_col, m2_col])["pocket_binder_ca_rmsd"].median().reset_index()
        _rmsd_heatmap(med, "pocket_binder_ca_rmsd", m1_col, m2_col,
                      "Median binder CA RMSD — pocket-aligned (Å)",
                      out_dir / "07_pocket_rmsd_heatmap.png")


# ── Plot 7: RMSD by method pair boxplot ──────────────────────────────────────

def plot_rmsd_by_pair(rmsd_df, out_dir):
    if rmsd_df is None or "binder_ca_rmsd" not in rmsd_df.columns:
        return
    ok = rmsd_df[rmsd_df["status"] == "ok"].copy()
    ok["pair"] = ok["method_1"] + " vs " + ok["method_2"]
    order = ok.groupby("pair")["binder_ca_rmsd"].median().sort_values().index
    fig, ax = plt.subplots(figsize=(max(10, len(order) * 0.6), 5))
    sns.boxplot(data=ok, x="pair", y="binder_ca_rmsd", order=order,
                color="#d9e6f2", width=0.55, showfliers=False, ax=ax)
    ax.set_xlabel("Method pair"); ax.set_ylabel("Binder CA RMSD (Å)")
    ax.set_title("Binder CA RMSD per method pair (rank-1, binder-on-binder alignment)")
    ax.tick_params(axis="x", rotation=55)
    sns.despine()
    save(fig, out_dir / "08_rmsd_by_method_pair.png")


# ─── CLI ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="./data",
                        help="Folder containing parquets from steps 1-3")
    parser.add_argument("--out-dir", default="./plots")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    out_dir  = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load whatever is available
    def try_load(name):
        p = data_dir / name
        if p.exists():
            print(f"  Loading {name}...")
            return pd.read_parquet(p)
        print(f"  Skipping {name} (not found).")
        return None

    df      = try_load("prediction_records_with_metrics.parquet")
    rank1   = try_load("rank1_structures.parquet")
    rmsd_df = try_load("pairwise_binder_ca_rmsd.parquet")
    pocket  = try_load("pairwise_pocket_align_rmsd.parquet")

    if df is None:
        df = try_load("prediction_records.parquet")  # fallback to step1 output

    print("\nGenerating plots...")

    if df is not None:
        plot_records_per_method(df, out_dir)
        plot_score_distributions(df, out_dir)

    if rank1 is not None:
        plot_iptm(rank1, out_dir)
        plot_plddt(rank1, out_dir)
        plot_binder_lengths(rank1, out_dir)
    elif df is not None:
        # If no rank1 parquet yet, build on the fly
        rank1_tmp = df[df["rank"] == 1].drop_duplicates(
            subset=["binder_sequence", "method"])
        plot_iptm(rank1_tmp, out_dir)
        plot_plddt(rank1_tmp, out_dir)
        plot_binder_lengths(rank1_tmp, out_dir)

    plot_rmsd_heatmaps(rmsd_df, pocket, out_dir)
    plot_rmsd_by_pair(rmsd_df, out_dir)

    print(f"\n✓ All plots saved to {out_dir}/")
