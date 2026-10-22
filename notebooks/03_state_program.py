# %% [markdown]
# # 03 IGF2BP2-associated program in candidate cancer cells
# For each tumor, finds genes that differ between IGF2BP2-detected and IGF2BP2-undetected candidate cancer
# cells (adjusting for depth and cell cycle), then tests whether the program found in one tumor reproduces
# in the other. The null shuffles IGF2BP2 detection among cells of similar depth and repeats the whole
# derivation, so selection effects and depth artifacts are part of the null.
#
# Run: python notebooks/03_state_program.py --h5ad_dir DIR --labels FILE --outdir DIR

# %% imports
import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.stats import hypergeom, rankdata, spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# %% settings
DAMAGED = re.compile(r"^[0-9]+-(Mar|Sep|Dec)|__[0-9]+$")
CYCLE_GENES = ["MKI67", "TOP2A", "CENPF", "CCNB1", "CCNB2", "UBE2C", "PCNA", "MCM2", "MCM4", "MCM6",
               "TYMS", "RRM2", "CDK1", "BIRC5", "AURKB", "PLK1", "CDC20", "HMGB2", "STMN1"]
TARGET = "IGF2BP2"
MIN_DETECTED = 0.10   # a gene must be detected in at least this fraction of candidate cancer cells, in both tumors
TOP_N = 100
N_PERM = 1000
N_STRATA = 10
SEED = 0


# %% helpers (plain numpy / pandas / scipy, tested without scanpy)
def damaged_mask(genes):
    return np.array([bool(DAMAGED.search(g)) for g in genes])


def to_dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)


def zscore_columns(Y):
    sd = Y.std(axis=0)
    sd[sd == 0] = np.nan
    return (Y - Y.mean(axis=0)) / sd


def cycle_score(expr, genes):
    idx = [list(genes).index(g) for g in CYCLE_GENES if g in set(genes)]
    return np.nan_to_num(zscore_columns(expr[:, idx])).mean(axis=1)


def ols_t(Y, design, col):
    """Least squares of each column of Y on the design. Returns (beta, t) for design column `col`."""
    xtx_inv = np.linalg.inv(design.T @ design)
    coef = xtx_inv @ design.T @ Y
    resid = Y - design @ coef
    dof = design.shape[0] - design.shape[1]
    s2 = (resid ** 2).sum(axis=0) / dof
    se = np.sqrt(s2 * xtx_inv[col, col])
    se[se == 0] = np.nan
    return coef[col], coef[col] / se


def residualize(y, covariates):
    d = np.column_stack([np.ones(len(y))] + list(covariates))
    b, *_ = np.linalg.lstsq(d, y, rcond=None)
    return y - d @ b


def auc(score, label):
    """Probability that a random positive cell scores above a random negative cell."""
    r = rankdata(score)
    n1 = int(label.sum())
    n0 = len(label) - n1
    return (r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def stratified_permutation(pos, strata, rng):
    """Shuffle detection labels among cells of similar depth."""
    new = pos.copy()
    for k in np.unique(strata):
        idx = np.where(strata == k)[0]
        new[idx] = rng.permutation(pos[idx])
    return new


def derive_t(expr_u, design, pos):
    d = design.copy()
    d[:, 1] = pos.astype(float)
    return ols_t(expr_u, d, 1)


def top_indices(t, n):
    return np.argsort(-np.nan_to_num(t, nan=-np.inf))[:n]


def program_auc(top_idx, z_dst, pos_dst, covs_dst):
    """Score dst cells with a program, remove depth and cycle, and return AUC for detected vs undetected."""
    return auc(residualize(z_dst[:, top_idx].mean(axis=1), covs_dst), pos_dst)


def statistics(tA, tB, posA, posB, zA, zB, covA, covB, n_top):
    topA, topB = top_indices(tA, n_top), top_indices(tB, n_top)
    return {"spearman_t": spearmanr(tA, tB, nan_policy="omit")[0],
            "top_overlap": len(set(topA.tolist()) & set(topB.tolist())),
            "auc_A_to_B": program_auc(topA, zB, posB, covB),
            "auc_B_to_A": program_auc(topB, zA, posA, covA)}


def reproducibility(data, universe, n_perm, rng, n_top=TOP_N):
    """Observed statistics and a null from depth-stratified shuffles of IGF2BP2 detection in both tumors."""
    A, B = data["A"], data["B"]
    zA = np.nan_to_num(zscore_columns(A["expr"][:, universe]))
    zB = np.nan_to_num(zscore_columns(B["expr"][:, universe]))
    covA, covB = [A["log_umi"], A["cycle"]], [B["log_umi"], B["cycle"]]
    sA = pd.qcut(A["log_umi"], N_STRATA, labels=False, duplicates="drop")
    sB = pd.qcut(B["log_umi"], N_STRATA, labels=False, duplicates="drop")
    A["beta"], A["t"] = derive_t(A["expr"][:, universe], A["design"], A["pos"])
    B["beta"], B["t"] = derive_t(B["expr"][:, universe], B["design"], B["pos"])
    obs = statistics(A["t"], B["t"], A["pos"], B["pos"], zA, zB, covA, covB, n_top)
    null = {k: [] for k in obs}
    for _ in range(n_perm):
        pA, pB = stratified_permutation(A["pos"], sA, rng), stratified_permutation(B["pos"], sB, rng)
        _, tA = derive_t(A["expr"][:, universe], A["design"], pA)
        _, tB = derive_t(B["expr"][:, universe], B["design"], pB)
        for k, v in statistics(tA, tB, pA, pB, zA, zB, covA, covB, n_top).items():
            null[k].append(v)
    rows = []
    for k, v in obs.items():
        arr = np.array(null[k])
        rows.append({"statistic": k, "observed": round(float(v), 3), "null_mean": round(float(arr.mean()), 3),
                     "null_sd": round(float(arr.std()), 3),
                     "p_upper": round((1 + (arr >= v).sum()) / (1 + n_perm), 4)})
    return pd.DataFrame(rows), len(universe)


# %% data loading
def load_tumor(h5ad_dir, labels, tumor):
    adata = ad.read_h5ad(Path(h5ad_dir) / f"sc_{tumor}_clustered.h5ad")
    cls = adata.obs["leiden"].astype(int).map(labels[tumor]).to_numpy()
    sub = adata[cls == "malignant"].copy()
    return {"genes": np.array(sub.var_names), "expr": to_dense(sub.X).astype(np.float64),
            "counts": to_dense(sub.layers["counts"])}


def prepare(d, j_target):
    counts = d["counts"]
    d["pos"] = counts[:, j_target] > 0
    d["log_umi"] = np.log(counts.sum(axis=1))
    d["cycle"] = cycle_score(d["expr"], d["genes"])
    d["design"] = np.column_stack([np.ones(len(d["pos"])), d["pos"].astype(float),
                                   d["log_umi"] - d["log_umi"].mean(), d["cycle"]])
    d["detected"] = (counts > 0).mean(axis=0)
    return d


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad_dir", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    with open(a.labels) as fh:
        labels = yaml.safe_load(fh)
    data = {t: load_tumor(a.h5ad_dir, labels, t) for t in ("A", "B")}
    genes = data["A"]["genes"]
    if not np.array_equal(genes, data["B"]["genes"]):
        raise SystemExit("gene lists differ between tumors A and B")
    j_target = list(genes).index(TARGET)
    for t in data:
        prepare(data[t], j_target)

    keep = ~damaged_mask(genes)
    keep[j_target] = False
    for t in data:
        keep &= data[t]["detected"] >= MIN_DETECTED
    universe = np.where(keep)[0]
    ugenes = genes[universe]

    repro, n_universe = reproducibility(data, universe, N_PERM, rng)
    repro.to_csv(out / "reproducibility.tsv", sep="\t", index=False)
    for t in data:
        data[t]["top_idx"] = top_indices(data[t]["t"], TOP_N)
    rho = float(repro.loc[repro["statistic"] == "spearman_t", "observed"].iloc[0])

    stats = pd.DataFrame({"gene": ugenes,
                          "t_A": data["A"]["t"], "beta_A": data["A"]["beta"],
                          "t_B": data["B"]["t"], "beta_B": data["B"]["beta"]})
    stats["top_A"] = False
    stats["top_B"] = False
    stats.loc[data["A"]["top_idx"], "top_A"] = True
    stats.loc[data["B"]["top_idx"], "top_B"] = True
    stats.round(4).to_csv(out / "state_stats.tsv", sep="\t", index=False)
    for t in ("A", "B"):
        top = stats.loc[stats[f"top_{t}"]].sort_values(f"t_{t}", ascending=False)["gene"]
        top.to_csv(out / f"program_top{TOP_N}_{t}.txt", index=False, header=False)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(stats["t_A"], stats["t_B"], s=4, c="lightgrey")
    both = stats["top_A"] & stats["top_B"]
    ax.scatter(stats.loc[stats["top_A"] & ~both, "t_A"], stats.loc[stats["top_A"] & ~both, "t_B"], s=10, c="tab:blue", label="top A")
    ax.scatter(stats.loc[stats["top_B"] & ~both, "t_A"], stats.loc[stats["top_B"] & ~both, "t_B"], s=10, c="tab:orange", label="top B")
    ax.scatter(stats.loc[both, "t_A"], stats.loc[both, "t_B"], s=16, c="tab:red", label="top in both")
    ax.axhline(0, c="k", lw=0.5)
    ax.axvline(0, c="k", lw=0.5)
    ax.set_xlabel(f"t statistic, IGF2BP2 detected vs not, tumor A")
    ax.set_ylabel(f"t statistic, tumor B")
    ax.set_title(f"Spearman rho = {rho:.2f} over {len(universe)} genes")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "program_scatter.png", dpi=120)

    pd.set_option("display.width", 220)
    print("candidate cancer cells:", {t: len(data[t]["pos"]) for t in data},
          "| IGF2BP2 detected:", {t: int(data[t]["pos"].sum()) for t in data},
          "| universe genes:", len(universe))
    print("\n== reproducibility (null: IGF2BP2 detection shuffled among cells of similar depth, both tumors) ==")
    print(repro.to_string(index=False))
    print(f"expected top-{TOP_N} overlap by chance without any shuffling: {TOP_N * TOP_N / len(universe):.1f}")
    for t in ("A", "B"):
        top = stats.loc[stats[f"top_{t}"]].sort_values(f"t_{t}", ascending=False)["gene"].tolist()
        print(f"\ntop 30 program genes, tumor {t}:", ", ".join(top[:30]))
    print(f"\nWrote files to {out}")


if __name__ == "__main__":
    main()