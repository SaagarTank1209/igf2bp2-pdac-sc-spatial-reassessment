# %% [markdown]
# # 11 State extension, part 1: comparison genes and the patient-corrected map (GSE202051 cancer nuclei)
# Chooses IGF2BP2 plus 100 detection-matched comparison genes BEFORE any map exists, removes all 101 genes
# from the map inputs, then builds the map (3,000 variable genes, PCA 30, Harmony on patient, 30 neighbours,
# Leiden at 0.2, 0.3, 0.5) and writes the cluster labels. Rules: docs/input_manifest.md (state extension plan).
#
# Run: python notebooks/11_state_map.py --cache FILE --outdir DIR --population untreated --variant main

# %% imports
import argparse
import importlib.util
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def load_module(name, fname):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m8 = load_module("step08", "08_external_analysis.py")
m3 = m8.m3

TARGET = "IGF2BP2"
N_COMPARISON = 100
N_HVG, N_PCS, N_NEIGHBORS, SEED = 3000, 30, 30, 0
K3_GRID = [0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0]
K3_MIN_SHARE = 0.05
NUCLEAR = ["MALAT1", "NEAT1", "FTX"]
NONCODING = re.compile(r"^(LINC|MIR|SNHG)|-AS[0-9]+$|-IT[0-9]+$|^(AC|AL|AP|AF|AJ|BX|CT|CU|FP|Z)[0-9]{5,}\.[0-9]+$")


# %% data and populations
def load_cache(path):
    a = ad.read_h5ad(path)
    genes = np.array(a.var_names)
    X = sparse.csr_matrix(a.X)
    counts = sparse.csr_matrix(a.layers["counts"])
    j = list(genes).index(TARGET)
    counts_t = m8.to_dense(counts[:, j]).ravel() > 0
    return a, genes, X, counts, a.obs, counts_t


def population_pids(obs, counts_t, which):
    pid_all = obs["pid"].astype(str).to_numpy()
    treat = obs.assign(pid=pid_all).groupby("pid")["treatment_status"].first()
    n_det = pd.Series(counts_t).groupby(pid_all).sum()
    n_tot = pd.Series(counts_t).groupby(pid_all).size()
    eligible = [p for p in sorted(treat.index) if n_det[p] >= m8.MIN_N and (n_tot[p] - n_det[p]) >= m8.MIN_N]
    return [p for p in eligible if treat[p] == "Untreated"] if which == "untreated" else eligible


# %% comparison genes (chosen before any map is built)
def choose_comparison_genes(counts, X, genes, pid_all, pids, n_comp=N_COMPARISON):
    """100 genes closest to IGF2BP2 in per-patient detection fraction AND mean log expression among detected nuclei."""
    rows_by = {p: np.where(pid_all == p)[0] for p in pids}
    pooled = np.concatenate(list(rows_by.values()))
    frac_pool = np.asarray((counts[pooled] > 0).mean(axis=0)).ravel()
    ok = (frac_pool >= m8.MIN_GENE_FRAC) & ~m3.damaged_mask(genes)
    for i, g in enumerate(genes):
        if g == TARGET or g == "IGF2BP2-AS1" or g in NUCLEAR or NONCODING.search(g):
            ok[i] = False
    cand = np.where(ok)[0]
    j = list(genes).index(TARGET)
    frac, level, n, f0, l0 = [], [], [], [], []
    for p in pids:
        rows = rows_by[p]
        nd = np.asarray((counts[rows][:, cand] > 0).sum(axis=0)).ravel()
        sums = np.asarray(X[rows][:, cand].sum(axis=0)).ravel()
        frac.append(nd / len(rows))
        level.append(np.where(nd > 0, sums / np.maximum(nd, 1), np.nan))
        n.append(len(rows))
        det_t = (counts[rows][:, j].toarray().ravel() > 0)
        xt = X[rows][:, j].toarray().ravel()
        f0.append(det_t.mean())
        l0.append(xt[det_t].mean() if det_t.any() else np.nan)
    frac, level = np.vstack(frac), np.vstack(level)
    n, f0, l0 = np.array(n)[:, None], np.array(f0), np.array(l0)
    qualifies = ((frac * n >= m8.MIN_N) & ((1 - frac) * n >= m8.MIN_N)).all(axis=0)
    d1 = np.abs(frac - f0[:, None]).mean(axis=0)
    d2 = np.nanmean(np.abs(level - l0[:, None]), axis=0)
    good = qualifies & np.isfinite(d1) & np.isfinite(d2)
    dist = np.full(len(cand), np.inf)
    dist[good] = d1[good] / np.median(d1[good]) + d2[good] / np.median(d2[good])
    order = np.argsort(dist)[:n_comp]
    order = order[np.isfinite(dist[order])]
    return genes[cand[order]], dist[order]


# %% the map
def choose_k3(cluster_fn, grid=K3_GRID, min_share=K3_MIN_SHARE):
    """Smallest resolution in the grid that gives exactly three clusters, each with at least min_share of the nuclei."""
    for r in grid:
        lab = cluster_fn(r)
        share = lab.value_counts(normalize=True)
        if len(share) == 3 and share.min() >= min_share:
            return r, lab
    return None, None


def build_map(X, obs, genes, exclude, keep_rows, resolutions):
    import scanpy as sc          # imported here so that the helpers above can be used without scanpy
    import harmonypy

    keep_genes = np.where(~np.isin(genes, list(exclude)))[0]
    adm = ad.AnnData(X=X[keep_rows][:, keep_genes].copy(),
                     obs=pd.DataFrame({"pid": obs["pid"].astype(str).to_numpy()[keep_rows]}, index=obs.index[keep_rows]),
                     var=pd.DataFrame(index=genes[keep_genes]))
    sc.pp.highly_variable_genes(adm, n_top_genes=N_HVG, flavor="seurat", batch_key="pid")
    adm = adm[:, adm.var["highly_variable"]].copy()
    sc.pp.scale(adm, max_value=10)
    sc.tl.pca(adm, n_comps=N_PCS, random_state=SEED)
    ho = harmonypy.run_harmony(adm.obsm["X_pca"], adm.obs[["pid"]], "pid")
    Z = np.asarray(ho.Z_corr)
    adm.obsm["X_harmony"] = Z if Z.shape[0] == adm.n_obs else Z.T       # harmonypy 2.1.0: orientation guard
    sc.pp.neighbors(adm, n_neighbors=N_NEIGHBORS, use_rep="X_harmony", random_state=SEED)
    cache = {}

    def leiden_at(r):
        if r not in cache:
            key = f"tmp_{r}"
            sc.tl.leiden(adm, resolution=r, random_state=SEED, flavor="igraph", n_iterations=2, directed=False, key_added=key)
            cache[r] = adm.obs[key].astype(str)
        return cache[r]

    out = {f"r{r}": leiden_at(r) for r in resolutions}
    r3, lab3 = choose_k3(leiden_at)
    if lab3 is not None:
        out["k3"] = lab3
        print(f"k3: resolution {r3} gives exactly three clusters of at least {K3_MIN_SHARE:.0%} of the nuclei each")
    else:
        print("k3: no resolution in the grid gives exactly three clusters of at least 5 percent each; variant S4 not formed")
    return pd.DataFrame(out)


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--population", choices=["untreated", "all_eligible"], default="untreated")
    ap.add_argument("--variant", choices=["main", "S1"], default="main")
    ap.add_argument("--resolutions", type=float, nargs="+", default=[0.2, 0.3, 0.5])
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    a, genes, X, counts, obs, counts_t = load_cache(args.cache)
    pid_all = obs["pid"].astype(str).to_numpy()
    pids = population_pids(obs, counts_t, args.population)
    print(f"{args.population}: {len(pids)} patients")

    gfile = out / f"exposure_genes_{args.population}.txt"
    if gfile.exists():
        expo = [g.strip() for g in gfile.read_text().splitlines() if g.strip()]
        print(f"reusing {len(expo)} exposure genes from {gfile}")
    else:
        comp, dist = choose_comparison_genes(counts, X, genes, pid_all, pids)
        expo = [TARGET] + list(comp)
        gfile.write_text("\n".join(expo) + "\n")
        print(f"chose {len(comp)} comparison genes (match distance {dist.min():.2f} to {dist.max():.2f}); wrote {gfile}")

    in_pop = np.isin(pid_all, pids)
    keep = in_pop.copy()
    if args.variant == "S1":
        if "cnv_score" not in obs.columns:
            raise SystemExit("S1 needs the cnv_score column in the cache. Add it to KEEP_OBS in step 07 and rerun step 07.")
        cnv = obs["cnv_score"].astype(float).to_numpy()
        keep = in_pop & (cnv > np.median(cnv[in_pop]))
        print(f"S1: {int(keep.sum())} of {int(in_pop.sum())} nuclei above the median copy-number score")
    rows = np.where(keep)[0]
    labels = build_map(X, obs, genes, expo, rows, args.resolutions)
    labels.to_csv(out / f"clusters_{args.population}_{args.variant}.tsv", sep="\t")
    for col in labels.columns:
        print(col, "clusters:", labels[col].nunique(), "| sizes:", labels[col].value_counts().sort_index().to_dict())
    print("wrote", out / f"clusters_{args.population}_{args.variant}.tsv")


if __name__ == "__main__":
    main()