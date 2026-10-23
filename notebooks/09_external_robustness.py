# %% [markdown]
# # 09 External check: robustness and specificity of E1 (post-hoc, untreated patients, primary exposure)
# E1 passed. This script asks whether the pass could come from a technical axis rather than from IGF2BP2:
#   V0 baseline; V1 adds log(detected genes per nucleus); V2 adds a nuclear-lncRNA score and removes
#   lncRNA-like genes from the programs; V3 does both; V4 repeats E1 with 100 prevalence-matched random genes
#   as the exposure. All of it is post-hoc and is reported as such.
#
# Run: python notebooks/09_external_robustness.py --cache FILE --outdir DIR [--n_perm 200]

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
m3, m4 = m8.m3, m8.m4

N_MATCH = 100
SEED = 1
NUCLEAR = ["MALAT1", "NEAT1", "FTX"]
NONCODING = re.compile(r"^(LINC|MIR|SNHG)|-AS[0-9]+$|-IT[0-9]+$|^(AC|AL|AP|AF|AJ|BX|CT|CU|FP|Z)[0-9]{5,}\.[0-9]+$")


# %% model with extra covariates
def t_x(P, e, extra):
    X = np.column_stack([np.ones(P.n), e.astype(float), P.ld, P.cycle] + extra).astype(np.float32)
    k = X.shape[1]
    inv = np.linalg.inv((X.T @ X).astype(np.float64))
    xty = (X.T @ P.Y[:, :P.U]).astype(np.float64)
    coef = inv @ xty
    rss = np.maximum(P.ssy - (coef * xty).sum(axis=0), 0)
    se = np.sqrt(rss / (P.n - k) * inv[1, 1])
    return np.where(se > 0, coef[1] / np.where(se > 0, se, 1.0), np.nan)


def loto_x(patients, Ts, exposures, extras, mask_cols):
    """Like step 08, with extra covariates removed from the held-out score and some genes barred from programs."""
    Ts = Ts.copy()
    if len(mask_cols):
        Ts[:, mask_cols] = np.nan
    ok = ~np.isnan(Ts)
    S, C = np.nansum(Ts, axis=0), ok.sum(axis=0)
    aucs = []
    for i, P in enumerate(patients):
        c = C - ok[i]
        others = np.where(c > 0, (S - np.nan_to_num(Ts[i])) / np.where(c > 0, c, 1), np.nan)
        top = np.argsort(-np.nan_to_num(others, nan=-np.inf))[:m8.TOP_N]
        score = ((P.Y[:, top] - P.mu[top]) * P.inv_sd[top]).mean(axis=1)
        D = np.column_stack([np.ones(P.n), P.ld, P.cycle] + extras[i])
        resid = score - D @ np.linalg.lstsq(D, score, rcond=None)[0]
        aucs.append(m3.auc(resid, exposures[i]))
    return np.array(aucs)


def run_variant(patients, extras, mask_cols, n_perm, rng):
    exp = [P.exposure["primary"] for P in patients]
    Ts = np.vstack([t_x(P, e, ex) for P, e, ex in zip(patients, exp, extras)])
    aucs = loto_x(patients, Ts, exp, extras, mask_cols)
    null = []
    for _ in range(n_perm):
        perm = [m3.stratified_permutation(e, P.strata, rng) for P, e in zip(patients, exp)]
        Tp = np.vstack([t_x(P, e, ex) for P, e, ex in zip(patients, perm, extras)])
        null.append(loto_x(patients, Tp, perm, extras, mask_cols).mean())
    null = np.array(null)
    mean, lo, hi = m8.t_interval(aucs)
    p = (1 + (null >= aucs.mean()).sum()) / (1 + n_perm)
    share = float((aucs > 0.5).mean())
    return {"mean_auc": mean, "ci_low": lo, "ci_high": hi, "share_auc_gt_half": share, "null_mean": null.mean(),
            "null_sd": null.std(), "p_upper": p, "passes": bool(mean > 0.5 and p <= 0.05 and share >= m8.E1_SHARE)}


# %% V4: prevalence-matched random exposures
def matched_exposures(patients, counts, rows_by, universe, genes, n_match):
    f0 = np.array([P.exposure["primary"].mean() for P in patients])
    frac = np.vstack([np.asarray((counts[rows_by[P.pid]][:, universe] > 0).mean(axis=0)).ravel() for P in patients])
    n = np.array([P.n for P in patients])[:, None]
    qualifies = ((frac * n >= m8.MIN_N) & ((1 - frac) * n >= m8.MIN_N)).all(axis=0)
    dist = np.abs(frac - f0[:, None]).mean(axis=0)
    dist[~qualifies] = np.inf
    for i, g in enumerate(genes[universe]):
        if g == "IGF2BP2-AS1":
            dist[i] = np.inf          # overlaps the IGF2BP2 locus
    chosen = np.argsort(dist)[:n_match]
    chosen = chosen[np.isfinite(dist[chosen])]
    det = [(counts[rows_by[P.pid]][:, universe[chosen]] > 0).toarray() for P in patients]
    rows = []
    for k, u in enumerate(chosen):
        exp = [d[:, k] for d in det]
        Ts = np.vstack([m8.t_for_exposure(P, e) for P, e in zip(patients, exp)])
        aucs = loto_x(patients, Ts, exp, [[] for _ in patients], [u])
        rows.append({"gene": genes[universe[u]], "distance": dist[u], "mean_auc": aucs.mean()})
    return pd.DataFrame(rows)


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n_perm", type=int, default=m8.N_PERM)
    ap.add_argument("--n_match", type=int, default=N_MATCH)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    a = ad.read_h5ad(args.cache)
    genes = np.array(a.var_names)
    X = sparse.csr_matrix(a.X)
    counts = sparse.csr_matrix(a.layers["counts"])
    j = list(genes).index(m8.TARGET)
    counts_t = m8.to_dense(counts[:, j]).ravel() > 0
    obs = a.obs
    total_counts = obs["total_counts"].to_numpy(dtype=float)
    n_genes = obs["n_genes"].to_numpy(dtype=float)
    pid_all = obs["pid"].astype(str).to_numpy()
    treat = obs.assign(pid=pid_all).groupby("pid")["treatment_status"].first()
    n_det = pd.Series(counts_t).groupby(pid_all).sum()
    n_tot = pd.Series(counts_t).groupby(pid_all).size()
    pids = [p for p in sorted(treat.index)
            if treat[p] == "Untreated" and n_det[p] >= m8.MIN_N and (n_tot[p] - n_det[p]) >= m8.MIN_N]
    patients, universe = m8.prepare_population("untreated", pids, obs, X, counts_t, total_counts, genes, m8.MIN_PATIENTS)
    patients = [P for P in patients if P.ok["primary"]]
    rows_by = {p: np.where(pid_all == p)[0] for p in pids}

    nuc_idx = [list(genes).index(g) for g in NUCLEAR if g in set(genes)]
    ug = genes[universe]
    mask_nc = np.array([i for i, g in enumerate(ug) if g in NUCLEAR or NONCODING.search(g)], dtype=int)
    print(f"nuclear lncRNA genes used as covariate: {[genes[i] for i in nuc_idx]}")
    print(f"genes barred from programs in V2 and V3: {len(mask_nc)}")
    log_ng, nuc = [], []
    for P in patients:
        ng = np.log(n_genes[rows_by[P.pid]])
        log_ng.append(ng - ng.mean())
        pos = [P.pos[i] for i in nuc_idx if i in P.pos]
        sc = m4.score_from_z(m4.zscore_columns(P.Y[:, pos]), list(range(len(pos)))) if pos else None
        nuc.append(np.nan_to_num(sc) if sc is not None else np.zeros(P.n))

    variants = {"V0 baseline": ([[] for _ in patients], []),
                "V1 + log detected genes": ([[x] for x in log_ng], []),
                "V2 + nuclear lncRNA score, lncRNA-like genes removed": ([[x] for x in nuc], mask_nc),
                "V3 both": ([[x, y] for x, y in zip(log_ng, nuc)], mask_nc)}
    rows = []
    for name, (ex, mask) in variants.items():
        rows.append({"variant": name, **run_variant(patients, ex, mask, args.n_perm, rng)})
        print("done:", name)
    var = pd.DataFrame(rows)
    var.round(4).to_csv(out / "e1_variants.tsv", sep="\t", index=False)

    rnd = matched_exposures(patients, counts, rows_by, universe, genes, args.n_match)
    rnd.round(4).to_csv(out / "e1_matched_random_exposures.tsv", sep="\t", index=False)
    obs_auc = float(var.loc[0, "mean_auc"])
    summary = {"igf2bp2_mean_auc": obs_auc, "n_random_genes": len(rnd), "random_mean_auc": rnd["mean_auc"].mean(),
               "random_sd": rnd["mean_auc"].std(), "random_max": rnd["mean_auc"].max(),
               "share_random_at_or_above_igf2bp2": float((rnd["mean_auc"] >= obs_auc).mean())}
    pd.Series(summary).round(4).to_csv(out / "e1_matched_random_summary.tsv", sep="\t", header=False)

    pd.set_option("display.width", 220)
    print("\n== E1 variants (untreated, primary exposure) ==")
    print(var.round(3).to_string(index=False))
    print("\n== E1 with prevalence-matched random genes as the exposure ==")
    for k, v in summary.items():
        print(f"  {k}: {v:.4g}" if isinstance(v, float) else f"  {k}: {v}")
    print("\nmatched genes with the highest AUC:")
    print(rnd.sort_values("mean_auc", ascending=False).head(8).round(3).to_string(index=False))
    print(f"\nWrote files to {out}")


if __name__ == "__main__":
    main()