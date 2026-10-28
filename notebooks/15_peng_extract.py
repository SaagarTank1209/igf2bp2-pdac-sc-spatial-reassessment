# %% [markdown]
# # 15 Second external dataset (Peng et al. 2019; PRJCA001063 / CRA001160; Besca reprocessed copy on Zenodo)
# Extracts the cancer cells (the original authors' "Ductal cell type 2" label, from tumor samples), checks the label,
# applies the eligibility rule of the plan in docs/input_manifest.md, and writes a cache in the same layout as the
# GSE202051 cache, so that steps 11 and 12 run on it unchanged. Only the NUMBER of cells with IGF2BP2 counts is read.
#
# Run: python notebooks/15_peng_extract.py --h5ad FILE --outdir DIR --cache_dir DIR

# %% imports
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

# %% locked settings (second external dataset plan)
TARGET = "IGF2BP2"
CANCER_LABEL = "Ductal cell type 2"
MIN_DETECTED = 30
MIN_UNDETECTED = 30
MIN_ELIGIBLE = 8
MIN_TUMOR_SHARE = 0.90     # label check: share of the cancer label that comes from tumor samples


def to_dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cache_dir", required=True)
    args = ap.parse_args()
    out, cache_dir = Path(args.outdir), Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(args.h5ad)
    obs = a.obs
    print("loaded", a.n_obs, "cells x", a.n_vars, "genes")
    X = sparse.csr_matrix(a.X)
    nz = X.data[: 2_000_000]
    print("integer counts:", bool(np.allclose(nz, np.round(nz))), "| maximum", float(X.data.max()))

    print("\ncell type by sample group (CONDITION):")
    print(pd.crosstab(obs["Cell_type"], obs["CONDITION"]).to_string())
    mal = (obs["Cell_type"] == CANCER_LABEL).to_numpy()
    tumor = (obs["CONDITION"] == "T").to_numpy()
    share = float(tumor[mal].mean())
    print(f"\ncells labeled '{CANCER_LABEL}': {int(mal.sum())}; share from tumor samples: {share:.3f} (rule: at least {MIN_TUMOR_SHARE})")
    if share < MIN_TUMOR_SHARE:
        raise SystemExit("label check failed: the cancer label is not mostly from tumor samples; stop and report")

    keep = mal & tumor
    total = np.asarray(X.sum(axis=1)).ravel()
    ngenes = np.asarray((X > 0).sum(axis=1)).ravel()
    j = list(a.var_names).index(TARGET)
    detected = to_dense(X[:, j]).ravel() > 0

    pid = obs["Patient"].astype(str).to_numpy()
    rows = []
    for p in sorted(set(pid[keep]), key=lambda s: (len(s), s)):
        m = keep & (pid == p)
        n, k = int(m.sum()), int(detected[m].sum())
        rows.append({"pid": p, "n_cancer_cells": n, "n_detected": k, "n_undetected": n - k,
                     "eligible": bool(k >= MIN_DETECTED and (n - k) >= MIN_UNDETECTED)})
    elig = pd.DataFrame(rows)
    elig.to_csv(out / "eligibility.tsv", sep="\t", index=False)
    n_el = int(elig["eligible"].sum())
    pd.set_option("display.width", 200)
    print(f"\ntumors with cancer cells: {len(elig)}; eligible: {n_el} (threshold {MIN_ELIGIBLE})")
    print("second external check feasible:", n_el >= MIN_ELIGIBLE)
    print(elig.to_string(index=False))

    Xk = X[keep].astype(np.float32)
    norm = sparse.diags((1e4 / total[keep]).astype(np.float32)) @ Xk
    Xlog = norm.tocsr().log1p()
    new_obs = pd.DataFrame({"pid": pid[keep], "sampleid": pid[keep], "treatment_status": "Untreated",
                            "total_counts": total[keep], "n_genes": ngenes[keep],
                            "Cell_type": obs["Cell_type"].astype(str).to_numpy()[keep]},
                           index=obs["CELL"].astype(str).to_numpy()[keep])
    cache = ad.AnnData(X=Xlog, obs=new_obs, var=pd.DataFrame(index=a.var_names))
    cache.layers["counts"] = Xk
    cache.write_h5ad(cache_dir / "Peng2019_cancer_cells.h5ad", compression="gzip")
    print("\nwrote", cache_dir / "Peng2019_cancer_cells.h5ad", "with", int(keep.sum()), "cancer cells")


if __name__ == "__main__":
    main()
