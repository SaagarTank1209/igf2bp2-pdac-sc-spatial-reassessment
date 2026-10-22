# %% [markdown]
# # 02 Single-cell annotation (exploratory)
# Normalizes and clusters each single-cell reference, then writes cluster marker tables so that cell
# types and malignant clusters can be called from evidence. No labels are assigned in this script.
#
# Run: python notebooks/02_sc_annotation.py --h5ad_dir DIR --outdir DIR [--resolution 0.8]

# %% imports
import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# %% settings
# Spreadsheet-damaged or renamed gene labels (see docs/input_manifest.md): excluded from HVGs and marker ranking.
DAMAGED = re.compile(r"^[0-9]+-(Mar|Sep|Dec)|__[0-9]+$")

MARKERS = {
    "epithelial / ductal": ["KRT19", "KRT7", "KRT8", "KRT18", "EPCAM", "MUC1"],
    "cancer clusters (Moncada)": ["TM4SF1", "S100A4"],
    "hypoxia": ["CA9", "ERO1A"],
    "terminal duct": ["TFF1", "TFF2", "TFF3"],
    "centroacinar": ["CRISP3", "CFTR", "AQP3"],
    "antigen presenting duct": ["CD74", "HLA-DRA"],
    "acinar": ["PRSS1", "CPA1", "CELA3A", "CTRB1"],
    "endocrine": ["INS", "GCG"],
    "fibroblast": ["COL1A1", "DCN", "LUM"],
    "stellate / pericyte": ["RGS5", "ACTA2"],
    "endothelial": ["PECAM1", "VWF", "CDH5"],
    "macrophage": ["CD68", "CD163", "C1QA", "MS4A4A"],
    "T / NK": ["CD3E", "CD2", "NKG7"],
    "B / plasma": ["MS4A1", "MZB1", "JCHAIN"],
    "mast": ["TPSAB1", "KIT"],
    "target genes": ["IGF2BP2", "SGMS2"],
}
UMAP_GENES = ["KRT19", "EPCAM", "TM4SF1", "S100A4", "COL1A1", "PTPRC", "CD68", "IGF2BP2"]


# %% helpers (plain numpy / pandas, so they can be tested without scanpy)
def damaged_mask(genes):
    return np.array([bool(DAMAGED.search(g)) for g in genes])


def to_dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)


def marker_table(X, var_names, clusters, panel):
    """Mean expression (genes x clusters) for the panel genes present. Returns (table, missing genes)."""
    var_names = list(var_names)
    wanted = [(grp, g) for grp, genes in panel.items() for g in genes]
    present = [(grp, g) for grp, g in wanted if g in var_names]
    missing = [g for grp, g in wanted if g not in var_names]
    sub = to_dense(X[:, [var_names.index(g) for _, g in present]])
    clusters = np.asarray(clusters).astype(str)
    order = sorted(set(clusters), key=int)
    table = pd.DataFrame({c: sub[clusters == c].mean(axis=0) for c in order},
                         index=pd.MultiIndex.from_tuples(present, names=["group", "gene"]))
    return table, missing


def frac_positive(C, var_names, clusters, genes):
    """Fraction of cells per cluster with a count above zero (genes x clusters)."""
    var_names = list(var_names)
    genes = [g for g in genes if g in var_names]
    sub = to_dense(C[:, [var_names.index(g) for g in genes]])
    clusters = np.asarray(clusters).astype(str)
    order = sorted(set(clusters), key=int)
    return pd.DataFrame({c: (sub[clusters == c] > 0).mean(axis=0) for c in order}, index=genes)


def gene_vector(adata, gene):
    return to_dense(adata[:, gene].X).ravel()


# %% one tumor
def annotate_one(path, tumor, out, resolution):
    adata = ad.read_h5ad(path)
    adata.X = adata.X.astype(np.float32)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    bad = damaged_mask(adata.var_names)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata.var.loc[bad, "highly_variable"] = False
    hv = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=30, random_state=0)
    sc.pp.neighbors(hv, n_neighbors=15, n_pcs=30, random_state=0)
    sc.tl.umap(hv, random_state=0)
    sc.tl.leiden(hv, resolution=resolution, random_state=0, flavor="igraph", n_iterations=2, directed=False)
    adata.obs["leiden"] = hv.obs["leiden"].astype(str).to_numpy()
    adata.obsm["X_umap"] = hv.obsm["X_umap"]
    adata.obsm["X_pca"] = hv.obsm["X_pca"]

    clusters = adata.obs["leiden"].to_numpy()
    sizes = adata.obs["leiden"].value_counts()
    sizes = sizes.reindex(sorted(sizes.index, key=int))

    mean, missing = marker_table(adata.X, adata.var_names, clusters, MARKERS)
    frac = frac_positive(adata.layers["counts"], adata.var_names, clusters, ["IGF2BP2", "SGMS2"])

    ranked = adata[:, ~bad].copy()
    sc.tl.rank_genes_groups(ranked, "leiden", method="wilcoxon", n_genes=10, use_raw=False)
    top = pd.DataFrame(ranked.uns["rank_genes_groups"]["names"]).head(8)

    adata.obs[["leiden"]].to_csv(out / f"clusters_{tumor}.tsv", sep="\t")
    mean.round(3).to_csv(out / f"marker_means_{tumor}.tsv", sep="\t")
    frac.round(3).to_csv(out / f"target_gene_fraction_{tumor}.tsv", sep="\t")
    top.to_csv(out / f"top_genes_{tumor}.tsv", sep="\t", index=False)
    adata.write_h5ad(out / f"sc_{tumor}_clustered.h5ad", compression="gzip")

    # UMAP overview
    umap = adata.obsm["X_umap"]
    cl = adata.obs["leiden"].astype(int).to_numpy()
    fig, axes = plt.subplots(3, 3, figsize=(13, 12))
    ax = axes[0, 0]
    ax.scatter(umap[:, 0], umap[:, 1], c=cl, cmap="tab20", s=6)
    for c in np.unique(cl):
        x, y = np.median(umap[cl == c], axis=0)
        ax.text(x, y, str(c), fontsize=11, weight="bold")
    ax.set_title(f"{tumor}: leiden clusters (resolution {resolution})")
    for ax, g in zip(axes.ravel()[1:], UMAP_GENES):
        if g in adata.var_names:
            v = gene_vector(adata, g)
            s = ax.scatter(umap[:, 0], umap[:, 1], c=v, cmap="viridis", s=6)
            fig.colorbar(s, ax=ax)
        ax.set_title(g)
    fig.tight_layout()
    fig.savefig(out / f"umap_{tumor}.png", dpi=120)
    plt.close(fig)

    # compact printout
    print(f"\n######## tumor {tumor}: {adata.n_obs} cells, {len(sizes)} clusters ########")
    print("cluster sizes:", ", ".join(f"{c}:{n}" for c, n in sizes.items()))
    print("\nmean log-normalized expression (rows genes, columns clusters)")
    print(mean.droplevel("group").round(2).to_string())
    print("\nfraction of cells with IGF2BP2 / SGMS2 counts")
    print(frac.round(2).to_string())
    print("\ntop genes per cluster (damaged labels excluded)")
    for c in top.columns:
        print(f"  cluster {c}: " + ", ".join(top[c].tolist()))
    if missing:
        print("\npanel genes absent from the data:", ", ".join(missing))


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad_dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--resolution", type=float, default=0.8)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 220)
    for tumor in ("A", "B"):
        annotate_one(Path(a.h5ad_dir) / f"sc_{tumor}.h5ad", tumor, out, a.resolution)
    print(f"\nWrote files to {out}")


if __name__ == "__main__":
    main()