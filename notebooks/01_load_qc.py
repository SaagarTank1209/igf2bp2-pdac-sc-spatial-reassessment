# %% [markdown]
# # 01 Load and QC
# Reads the GSE111672 tables, applies the cleaning rules locked in docs/input_manifest.md,
# writes one h5ad per array and QC summary tables.
#
# Run: python notebooks/01_load_qc.py --counts_dir DIR --raw_dir DIR --outdir DIR

# %% imports
import argparse
import gzip
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# %% locked inputs (see docs/input_manifest.md)
SPATIAL = {
    "A1": ("A", "GSM3036911_PDAC-A-ST1-filtered.txt.gz"),
    "A2": ("A", "GSM4100721_PDAC-A-st2.tsv.gz"),
    "A3": ("A", "GSM4100722_PDAC-A-st3.tsv.gz"),
    "B1": ("B", "GSM3405534_PDAC-B-ST1-filtered.txt.gz"),
    "B2": ("B", "GSM4100723_PDAC-B-st2.tsv.gz"),
    "B3": ("B", "GSM4100724_PDAC-B-st3.tsv.gz"),
    "D": ("D", "GSM4100725_PDAC-D-st1.tsv.gz"),
    "E": ("E", "GSM4100726_PDAC-E-st1.tsv.gz"),
    "G": ("G", "GSM4100728_PDAC-G-st1.tsv.gz"),
}
SC = {
    "A": "GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz",
    "B": "GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz",
}
KEY_GENES = ["IGF2BP2", "SGMS2", "CD274", "KRT19", "EPCAM", "PTPRC", "COL1A1", "CD68"]


# %% readers
def read_matrix(path):
    """Read a genes x columns count table. Returns (counts int32 array genes x cols, genes, labels)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        first_row = fh.readline().rstrip("\n").split("\t")
    if len(header) == len(first_row):      # header carries a label for the gene column
        labels = header[1:]
    elif len(header) == len(first_row) - 1:  # header has no gene-column label
        labels = header
    else:
        raise SystemExit(f"{path.name}: header has {len(header)} fields, first row {len(first_row)}")
    df = pd.read_csv(path, sep="\t", header=None, skiprows=1, index_col=0)
    if df.shape[1] != len(labels):
        raise SystemExit(f"{path.name}: {df.shape[1]} data columns but {len(labels)} labels")
    vals = df.to_numpy()
    if np.isnan(vals).any() or (vals < 0).any() or not np.array_equal(vals, np.round(vals)):
        raise SystemExit(f"{path.name}: values are not non-negative integers")
    return vals.astype(np.int32), df.index.astype(str).to_numpy(), np.asarray(labels, dtype=str)


def drop_duplicate_labels(X, labels, name):
    """Keep the first column of any duplicated label. Stop if the duplicates are not identical."""
    seen, keep, n_dup = {}, [], 0
    for j, lab in enumerate(labels):
        if lab in seen:
            n_dup += 1
            if not np.array_equal(X[:, j], X[:, seen[lab]]):
                raise SystemExit(f"{name}: duplicated label {lab} with different counts")
            continue
        seen[lab] = j
        keep.append(j)
    return X[:, keep], labels[keep], n_dup


def parse_xy(labels, name):
    try:
        xy = np.array([[int(v) for v in lab.split("x")] for lab in labels])
    except ValueError:
        raise SystemExit(f"{name}: spot labels are not of the form XxY, e.g. {labels[:3]}")
    if xy.shape[1] != 2:
        raise SystemExit(f"{name}: spot labels do not have exactly two coordinates")
    return xy


# %% QC helpers
def spot_geometry(xy):
    coords = {tuple(c) for c in xy}
    has_nb = sum(
        any((x + dx, y + dy) in coords for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)))
        for x, y in coords
    )
    xr, yr = xy[:, 0].max() - xy[:, 0].min() + 1, xy[:, 1].max() - xy[:, 1].min() + 1
    return {
        "x_min": xy[:, 0].min(), "x_max": xy[:, 0].max(),
        "y_min": xy[:, 1].min(), "y_max": xy[:, 1].max(),
        "bbox_fill": round(len(coords) / (xr * yr), 3),
        "frac_with_grid_neighbor": round(has_nb / len(coords), 3),
    }


def key_gene_table(adata, label):
    rows = []
    genes = list(adata.var_names)
    for g in KEY_GENES:
        if g in genes:
            col = adata[:, g].X
            col = col.toarray().ravel() if sparse.issparse(col) else np.asarray(col).ravel()
            rows.append({"dataset": label, "gene": g, "present": True,
                         "frac_detected": round(float((col > 0).mean()), 4),
                         "total_counts": int(col.sum())})
        else:
            rows.append({"dataset": label, "gene": g, "present": False,
                         "frac_detected": np.nan, "total_counts": np.nan})
    return rows


def q(a, p):
    return int(np.percentile(a, p))


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts_dir", required=True)
    ap.add_argument("--raw_dir", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    counts_dir, raw_dir, out = Path(a.counts_dir), Path(a.raw_dir), Path(a.outdir)
    (out / "h5ad").mkdir(parents=True, exist_ok=True)

    sc_rows, sp_rows, gene_rows, ref_genes, notes = [], [], [], None, []

    # single-cell references
    for tumor, fname in SC.items():
        X, genes, cells = read_matrix(raw_dir / fname)
        if len(set(genes)) != len(genes):
            raise SystemExit(f"{fname}: duplicated gene names")
        if ref_genes is None:
            ref_genes = genes
        adata = ad.AnnData(X=sparse.csr_matrix(X.T), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
        adata.obs["tumor"] = tumor
        umi = np.asarray(adata.X.sum(axis=1)).ravel()
        ngene = np.asarray((adata.X > 0).sum(axis=1)).ravel()
        mt = np.array([g.startswith("MT-") for g in genes])
        pct_mt = 100 * np.asarray(adata[:, mt].X.sum(axis=1)).ravel() / np.maximum(umi, 1)
        sc_rows.append({"tumor": tumor, "n_cells": adata.n_obs, "n_genes": adata.n_vars,
                        "umi_q05": q(umi, 5), "umi_median": q(umi, 50), "umi_q95": q(umi, 95),
                        "genes_median": q(ngene, 50), "n_mt_genes": int(sum(mt)),
                        "pct_mt_median": round(float(np.median(pct_mt)), 2),
                        "pct_mt_max": round(float(pct_mt.max()), 2)})
        gene_rows += key_gene_table(adata, f"sc_{tumor}")
        adata.write_h5ad(out / "h5ad" / f"sc_{tumor}.h5ad", compression="gzip")

    # spatial arrays
    coords_for_plot = {}
    for arr, (tumor, fname) in SPATIAL.items():
        X, genes, labels = read_matrix(counts_dir / fname)
        n_raw = X.shape[1]
        X, labels, n_dup = drop_duplicate_labels(X, labels, arr)
        xy = parse_xy(labels, arr)
        same_set = set(genes) == set(ref_genes)
        same_order = bool(same_set and np.array_equal(genes, ref_genes))
        adata = ad.AnnData(X=sparse.csr_matrix(X.T), obs=pd.DataFrame(index=labels), var=pd.DataFrame(index=genes))
        adata.obs["tumor"], adata.obs["array"] = tumor, arr
        adata.obsm["spatial"] = xy
        umi = np.asarray(adata.X.sum(axis=1)).ravel()
        ngene = np.asarray((adata.X > 0).sum(axis=1)).ravel()
        row = {"array": arr, "tumor": tumor, "n_spots_raw": n_raw, "n_duplicates_dropped": n_dup,
               "n_spots": adata.n_obs, "n_genes": adata.n_vars,
               "umi_q05": q(umi, 5), "umi_median": q(umi, 50), "umi_q95": q(umi, 95),
               "genes_median": q(ngene, 50), "genes_same_set_as_sc": same_set,
               "genes_same_order_as_sc": same_order}
        row.update(spot_geometry(xy))
        sp_rows.append(row)
        gene_rows += key_gene_table(adata, f"st_{arr}")
        coords_for_plot[arr] = (xy, umi)
        adata.write_h5ad(out / "h5ad" / f"st_{arr}.h5ad", compression="gzip")

    sc_df, sp_df, gene_df = pd.DataFrame(sc_rows), pd.DataFrame(sp_rows), pd.DataFrame(gene_rows)
    sc_df.to_csv(out / "sc_qc.tsv", sep="\t", index=False)
    sp_df.to_csv(out / "spatial_qc.tsv", sep="\t", index=False)
    gene_df.to_csv(out / "key_gene_detection.tsv", sep="\t", index=False)

    # spot maps
    fig, axes = plt.subplots(3, 3, figsize=(12, 12))
    for ax, (arr, (xy, umi)) in zip(axes.ravel(), coords_for_plot.items()):
        sc_ = ax.scatter(xy[:, 0], xy[:, 1], c=np.log10(umi + 1), s=14, cmap="viridis")
        ax.set_title(f"{arr}  ({len(umi)} spots)")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        fig.colorbar(sc_, ax=ax, label="log10 UMI")
    fig.tight_layout()
    fig.savefig(out / "spot_maps.png", dpi=120)

    pd.set_option("display.width", 200)
    print("== single-cell ==")
    print(sc_df.to_string(index=False))
    print("\n== spatial ==")
    print(sp_df.drop(columns=["tumor"]).to_string(index=False))
    print("\n== key gene detection (fraction of cells or spots with count > 0) ==")
    print(gene_df.pivot(index="gene", columns="dataset", values="frac_detected").to_string())
    print(f"\nWrote files to {out}")


if __name__ == "__main__":
    main()