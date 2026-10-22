# %% [markdown]
# # 04 Spatial scores and tests
# Builds per-spot signature scores for each array and tests, per array and per tumor:
#   H1 lipid score and H2 macrophage score against IGF2BP2 detection (primary),
#   secondary scores, and a positive control (stress score against stromal score).
# Models adjust for log total UMI, an epithelial score and (for H1, H2) a stromal score.
# Significance uses a spatially aware null: the exposure is shuffled within spatial blocks, which keeps
# large-scale gradients and tests only local association. Random gene sets matched for expression are run
# through the same procedure as a negative control.
#
# Run: python notebooks/04_spatial_tests.py --h5ad_dir DIR --outdir DIR

# %% imports
import argparse
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import t as t_dist

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# %% settings (signature lists are locked in docs/input_manifest.md)
SIGS = {
    "lipid": ["FASN", "SCD", "ACACA", "ACLY", "ELOVL6", "DGAT1", "SREBF1", "SREBF2", "HMGCR", "HMGCS1",
              "INSIG1", "LDLR", "PLIN2", "FABP5", "CD36"],
    "macrophage": ["CD68", "CSF1R", "C1QA", "C1QB", "C1QC", "MS4A4A", "CD163", "MRC1", "AIF1"],
    "lam": ["TREM2", "APOE", "APOC1", "LPL", "CD9", "GPNMB", "SPP1"],
    "exhaustion": ["PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2", "CD8A", "FOXP3"],
    "stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "epithelial": ["KRT19", "KRT7", "KRT8", "KRT18", "MUC1", "EPCAM"],
    "stress": ["CEBPB", "DNAJA1", "FOS", "FOSB", "JUN", "JUNB", "KLF2", "KLF4", "MAP2K3"],
}
TARGET = "IGF2BP2"
ARRAYS = ["A1", "A2", "A3", "B1", "B2", "B3", "D", "E", "G"]
TUMOR_OF = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B", "B3": "B", "D": "D", "E": "E", "G": "G"}
BLOCK = 5          # block edge in array grid units
N_PERM = 1000
N_RANDOM = 200
MIN_EXPOSED = 8    # fewer exposed spots than this: the test is skipped
SEED = 0


# %% helpers (plain numpy / pandas / scipy)
def to_dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)


def zscore_columns(Y):
    sd = Y.std(axis=0)
    sd[sd == 0] = np.nan
    return (Y - Y.mean(axis=0)) / sd


def zs(v):
    sd = np.std(v)
    return (v - np.mean(v)) / sd if sd > 0 else np.zeros_like(v)


def score_from_z(Z, idx):
    """Mean z-score over the signature genes that vary. Returns None if none do."""
    if len(idx) == 0:
        return None
    sub = Z[:, idx]
    ok = ~np.isnan(sub).all(axis=0)
    return np.nanmean(sub[:, ok], axis=1) if ok.any() else None


def ols_t(Y, design, col):
    xtx_inv = np.linalg.inv(design.T @ design)
    coef = xtx_inv @ design.T @ Y
    resid = Y - design @ coef
    dof = design.shape[0] - design.shape[1]
    se = np.sqrt((resid ** 2).sum(axis=0) / dof * xtx_inv[col, col])
    se[se == 0] = np.nan
    return coef[col], coef[col] / se


def block_ids(xy, size):
    bx, by = xy[:, 0] // size, xy[:, 1] // size
    return pd.factorize(pd.Series(list(zip(bx, by))))[0]


def block_permute(values, blocks, rng):
    """Shuffle values among the spots of each spatial block."""
    new = values.copy()
    for b in np.unique(blocks):
        idx = np.where(blocks == b)[0]
        new[idx] = rng.permutation(values[idx])
    return new


def matched_random_sets(mean_expr, eligible, sig_idx, n_sets, rng, n_bins=10):
    """Random gene sets of the same size, each gene drawn from the same expression decile as a signature gene."""
    bins = pd.qcut(pd.Series(mean_expr).rank(method="first"), n_bins, labels=False).to_numpy()
    pools = {b: np.where(eligible & (bins == b))[0] for b in range(n_bins)}
    return [np.array([rng.choice(pools[bins[g]]) for g in sig_idx], dtype=int) for _ in range(n_sets)]


def run_group(exposure, covs, outcomes, random_sets, blocks, rng):
    """Observed and block-permutation statistics for several outcomes sharing one design.
    outcomes: dict name -> vector. random_sets: dict name -> (n_spots x n_sets) score matrix."""
    names = list(outcomes)
    cols = [zs(outcomes[n]) for n in names]
    owner = names.copy()
    for n, mat in random_sets.items():
        for j in range(mat.shape[1]):
            cols.append(zs(mat[:, j]))
            owner.append("random:" + n)
    Y = np.column_stack(cols)
    design = np.column_stack([np.ones(len(exposure)), exposure.astype(float)] + covs)
    beta, t_obs = ols_t(Y, design, 1)
    exceed = np.zeros(Y.shape[1])
    for _ in range(N_PERM):
        d = design.copy()
        d[:, 1] = block_permute(exposure.astype(float), blocks, rng)
        _, t_p = ols_t(Y, d, 1)
        exceed += np.abs(t_p) >= np.abs(t_obs)
    p = (1 + exceed) / (1 + N_PERM)
    res = pd.DataFrame({"outcome": owner, "beta": beta, "t": t_obs, "p_block": p})
    # second, stricter comparison: how often do matched random gene sets show an association as large?
    res["p_vs_random"] = np.nan
    for n in names:
        rt = np.abs(res.loc[res["outcome"] == "random:" + n, "t"].to_numpy())
        if len(rt):
            res.loc[res["outcome"] == n, "p_vs_random"] = (1 + (rt >= abs(res.loc[res["outcome"] == n, "t"].iloc[0])).sum()) / (1 + len(rt))
    return res


def analyze_array(arr, adata, rng):
    counts = to_dense(adata.X).astype(float)
    genes = list(adata.var_names)
    gi = {g: i for i, g in enumerate(genes)}
    total = counts.sum(axis=1)
    logx = np.log1p(counts / total[:, None] * 1e4)
    Z = zscore_columns(logx)
    xy = np.asarray(adata.obsm["spatial"])
    blocks = block_ids(xy, BLOCK)
    sig_idx = {n: [gi[g] for g in s if g in gi] for n, s in SIGS.items()}
    scores = {n: score_from_z(Z, idx) for n, idx in sig_idx.items()}
    log_umi = zs(np.log(total))
    exposed = counts[:, gi[TARGET]] > 0

    mean_expr = logx.mean(axis=0)
    eligible = ((counts > 0).mean(axis=0) >= 0.01) & ~np.isnan(Z).all(axis=0)
    exclude = set(i for idx in sig_idx.values() for i in idx) | {gi[TARGET]}
    for i in exclude:
        eligible[i] = False

    varies = ~np.isnan(Z).all(axis=0)

    def random_scores(name):
        match_idx = [i for i in sig_idx[name] if varies[i]]
        sets = matched_random_sets(mean_expr, eligible, match_idx, N_RANDOM, rng)
        return np.column_stack([np.nan_to_num(Z[:, s]).mean(axis=1) for s in sets])

    rows = []
    epi, stro = zs(scores["epithelial"]), zs(scores["stromal"])
    if exposed.sum() >= MIN_EXPOSED:
        outs = {"H1_lipid": scores["lipid"], "H2_macrophage": scores["macrophage"],
                "S1_lam": scores["lam"], "S2_exhaustion": scores["exhaustion"]}
        outs = {k: v for k, v in outs.items() if v is not None}
        rnd = {"H1_lipid": random_scores("lipid"), "H2_macrophage": random_scores("macrophage")}
        rnd = {k: v for k, v in rnd.items() if k in outs}
        res = run_group(exposed, [log_umi, epi, stro], outs, rnd, blocks, rng)
        res["test"] = "IGF2BP2 detected"
        rows.append(res)
    if scores["stress"] is not None:
        rnd = {"C1_stress": random_scores("stress")}
        res = run_group(stro, [log_umi, epi], {"C1_stress": scores["stress"]}, rnd, blocks, rng)
        res["test"] = "stromal score"
        rows.append(res)
    res = pd.concat(rows, ignore_index=True)
    res.insert(0, "array", arr)
    info = {"array": arr, "tumor": TUMOR_OF[arr], "n_spots": len(total), "n_exposed": int(exposed.sum()),
            "genes_used": {n: int(len(idx)) for n, idx in sig_idx.items()}}
    spot = pd.DataFrame({"array": arr, "spot": adata.obs_names, "x": xy[:, 0], "y": xy[:, 1],
                         "log_umi": log_umi, "IGF2BP2_detected": exposed})
    for n, v in scores.items():
        spot[n] = v if v is not None else np.nan
    return res, info, spot


def tumor_summary(real):
    """Tumor-level effect per outcome: mean over the tumor's arrays, then a t interval over the 5 tumors.
    meets_criterion also requires a median p versus matched random gene sets below 0.05 (see manifest)."""
    out = []
    for outcome, g in real.groupby("outcome"):
        per_tumor = g.groupby("tumor")["beta"].mean()
        n = len(per_tumor)
        mean, sd = per_tumor.mean(), per_tumor.std(ddof=1)
        half = t_dist.ppf(0.975, n - 1) * sd / np.sqrt(n) if n > 1 else np.nan
        n_pos = int((per_tumor > 0).sum())
        med_p = g["p_vs_random"].median()
        supports_pos = bool(n_pos >= 4 and mean - half > 0)
        supports_neg = bool((n - n_pos) >= 4 and mean + half < 0)
        out.append({"outcome": outcome, "n_tumors": n, "n_positive": n_pos, "mean_beta": round(mean, 3),
                    "ci_low": round(mean - half, 3), "ci_high": round(mean + half, 3),
                    "median_p_vs_random": round(med_p, 3) if med_p == med_p else np.nan,
                    "supports_positive": supports_pos, "supports_negative": supports_neg,
                    "meets_criterion": bool((supports_pos or supports_neg) and med_p == med_p and med_p < 0.05)})
    return pd.DataFrame(out)


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad_dir", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    results, infos, spots = [], [], []
    for arr in ARRAYS:
        adata = ad.read_h5ad(Path(a.h5ad_dir) / f"st_{arr}.h5ad")
        res, info, spot = analyze_array(arr, adata, rng)
        res["tumor"] = info["tumor"]
        results.append(res)
        infos.append(info)
        spots.append(spot)
        print(f"{arr}: {info['n_spots']} spots, {info['n_exposed']} with IGF2BP2 detected, genes used {info['genes_used']}")
    allres = pd.concat(results, ignore_index=True)
    real = allres[~allres["outcome"].str.startswith("random:")].copy()
    rand = allres[allres["outcome"].str.startswith("random:")]

    real.round(4).to_csv(out / "array_tests.tsv", sep="\t", index=False)
    pd.concat(spots, ignore_index=True).round(4).to_csv(out / "spot_scores.tsv", sep="\t", index=False)
    summ = tumor_summary(real)
    summ.to_csv(out / "tumor_summary.tsv", sep="\t", index=False)
    calib = (rand.groupby(["array", "outcome"])
             .agg(n_sets=("t", "size"), frac_abs_t_gt_2=("t", lambda s: float((s.abs() > 2).mean())),
                  frac_p_block_lt_05=("p_block", lambda s: float((s < 0.05).mean())))
             .round(3).reset_index())
    calib.to_csv(out / "negative_control_calibration.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7, 4 + 0.35 * real["outcome"].nunique()))
    for i, (outcome, g) in enumerate(real.groupby("outcome")):
        ax.scatter(g["beta"], np.full(len(g), i), s=22, c="tab:blue", alpha=0.7)
        ax.scatter([g.groupby("tumor")["beta"].mean().mean()], [i], marker="D", c="k", s=30)
        ax.text(ax.get_xlim()[1], i, "", va="center")
    ax.set_yticks(range(real["outcome"].nunique()))
    ax.set_yticklabels(sorted(real["outcome"].unique()))
    ax.axvline(0, c="k", lw=0.5)
    ax.set_xlabel("adjusted effect (SD of outcome per unit of exposure); dots = arrays, diamond = mean over tumors")
    fig.tight_layout()
    fig.savefig(out / "effects.png", dpi=120)

    pd.set_option("display.width", 220)
    print("\n== per array (block-permutation p) ==")
    piv = real.pivot_table(index="outcome", columns="array", values="beta").round(2)
    print(piv[[c for c in ARRAYS if c in piv.columns]].to_string())
    pv = real.pivot_table(index="outcome", columns="array", values="p_block").round(3)
    print(pv[[c for c in ARRAYS if c in pv.columns]].to_string())
    print("\n== per array: p versus matched random gene sets (blank if no random sets) ==")
    pr = real.pivot_table(index="outcome", columns="array", values="p_vs_random").round(3)
    print(pr[[c for c in ARRAYS if c in pr.columns]].to_string())
    print("\n== tumor-level summary (5 tumors) ==")
    print(summ.to_string(index=False))
    print("\n== negative control: share of matched random gene sets with block p < 0.05 (expect about 0.05) ==")
    print(calib.pivot_table(index="outcome", columns="array", values="frac_p_block_lt_05").round(2).to_string())
    print(f"\nWrote files to {out}")


if __name__ == "__main__":
    main()