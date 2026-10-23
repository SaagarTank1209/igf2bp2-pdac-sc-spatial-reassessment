# %% [markdown]
# # 08 External check on GSE202051 (E1 and E2, as locked in docs/input_manifest.md)
# E1: is an IGF2BP2-linked program reproducible across tumors in cancer nuclei (leave one tumor out)?
# E2: is IGF2BP2 linked to the lipid score within cancer nuclei, beyond matched random gene sets?
# Primary exposure: IGF2BP2 detected. Amendment A1: IGF2BP2 expression above the patient's median.
# Primary population: untreated patients. Secondary: all eligible patients.
#
# Run: python notebooks/08_external_analysis.py --cache FILE --outdir DIR [--n_perm 200]

# %% imports
import argparse
import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import t as t_dist


def load_module(name, fname):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m3 = load_module("step03", "03_state_program.py")   # cycle genes, auc, stratified_permutation
m4 = load_module("step04", "04_spatial_tests.py")   # signatures, matched random gene sets, ols_t, zs

# %% locked settings (see the external check plan in docs/input_manifest.md)
TARGET = "IGF2BP2"
MIN_N = 30               # nuclei with and without the exposure, per patient
MIN_PATIENTS = 8
MIN_GENE_FRAC = 0.05     # gene universe: detected in at least this fraction of pooled cancer nuclei
TOP_N = 100
N_PERM = 200
N_STRATA = 10
N_RANDOM = 200
E1_SHARE, E2_SHARE = 0.60, 0.70
SEED = 0


def to_dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)


def t_interval(values):
    v = np.asarray(values, dtype=float)
    n = len(v)
    half = t_dist.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    return v.mean(), v.mean() - half, v.mean() + half


# %% per-patient objects
class Patient:
    pass


def build_patient(pid, rows, X, detected, total_counts, universe, extras, genes_index):
    """Dense log-normalized matrix with the universe genes first, then the extra genes (signatures, cycle, target)."""
    cols = np.concatenate([universe, extras])
    P = Patient()
    P.pid = pid
    P.Y = X[rows][:, cols].toarray().astype(np.float32)
    P.U = len(universe)
    pos = {g: i for i, g in enumerate(cols)}
    P.pos = pos
    n = P.Y.shape[0]
    P.n = n
    depth = np.log(total_counts[rows])
    P.ld = depth - depth.mean()
    cyc_pos = [pos[g] for g in genes_index["cycle"] if g in pos]
    P.cycle = (np.nan_to_num(m4.zscore_columns(P.Y[:, cyc_pos])).mean(axis=1) if len(cyc_pos) >= 3 else np.zeros(n))
    P.strata = pd.qcut(P.ld, N_STRATA, labels=False, duplicates="drop")
    YU = P.Y[:, :P.U]
    P.ssy = (YU.astype(np.float64) ** 2).sum(axis=0)
    P.mu = YU.mean(axis=0, dtype=np.float64)
    sd = YU.std(axis=0, dtype=np.float64)
    P.inv_sd = np.where(sd > 0, 1.0 / np.where(sd > 0, sd, 1.0), 0.0)
    D3 = np.column_stack([np.ones(n), P.ld, P.cycle])
    P.D3, P.pinv = D3, np.linalg.pinv(D3)
    P.lipid_idx = genes_index["lipid"]
    expr_t = P.Y[:, pos[genes_index["target"]]]
    P.exposure = {"primary": detected[rows], "A1": expr_t > np.median(expr_t)}
    P.ok = {k: bool(e.sum() >= MIN_N and (~e).sum() >= MIN_N) for k, e in P.exposure.items()}
    return P


def t_for_exposure(P, e):
    """t statistic of the exposure for every universe gene: expression on exposure, log depth and cycle score."""
    X = np.column_stack([np.ones(P.n), e.astype(float), P.ld, P.cycle]).astype(np.float32)
    inv = np.linalg.inv((X.T @ X).astype(np.float64))
    xty = (X.T @ P.Y[:, :P.U]).astype(np.float64)
    coef = inv @ xty
    rss = np.maximum(P.ssy - (coef * xty).sum(axis=0), 0)
    se = np.sqrt(rss / (P.n - 4) * inv[1, 1])
    return np.where(se > 0, coef[1] / np.where(se > 0, se, 1.0), np.nan)


def loto_aucs(patients, Ts, exposures):
    """Leave-one-tumor-out: program from the other patients, scored in the held-out patient (depth and cycle removed)."""
    ok = ~np.isnan(Ts)
    S, C = np.nansum(Ts, axis=0), ok.sum(axis=0)
    aucs = []
    for i, P in enumerate(patients):
        c = C - ok[i]
        others = np.where(c > 0, (S - np.nan_to_num(Ts[i])) / np.where(c > 0, c, 1), np.nan)
        top = np.argsort(-np.nan_to_num(others, nan=-np.inf))[:TOP_N]
        score = ((P.Y[:, top] - P.mu[top]) * P.inv_sd[top]).mean(axis=1)
        resid = score - P.D3 @ (P.pinv @ score)
        aucs.append(m3.auc(resid, exposures[i]))
    return np.array(aucs)


# %% E1
def run_e1(patients, key, n_perm, rng):
    exp = [P.exposure[key] for P in patients]
    Ts = np.vstack([t_for_exposure(P, e) for P, e in zip(patients, exp)])
    aucs = loto_aucs(patients, Ts, exp)
    obs = aucs.mean()
    null = []
    for _ in range(n_perm):
        perm = [m3.stratified_permutation(e, P.strata, rng) for P, e in zip(patients, exp)]
        Tp = np.vstack([t_for_exposure(P, e) for P, e in zip(patients, perm)])
        null.append(loto_aucs(patients, Tp, perm).mean())
    null = np.array(null)
    mean, lo, hi = t_interval(aucs)
    p = (1 + (null >= obs).sum()) / (1 + n_perm)
    share = float((aucs > 0.5).mean())
    summary = {"n_patients": len(patients), "mean_auc": mean, "ci_low": lo, "ci_high": hi,
               "share_auc_gt_half": share, "null_mean": null.mean(), "null_sd": null.std(), "p_upper": p,
               "passes": bool(mean > 0.5 and p <= 0.05 and share >= E1_SHARE)}
    per = pd.DataFrame({"pid": [P.pid for P in patients], "auc": aucs})
    consensus = np.nanmean(Ts, axis=0)
    return summary, per, consensus


# %% E2
def run_e2(patients, key, rng, universe_n):
    rows = []
    for P in patients:
        sig_pos = [P.pos[g] for g in P.lipid_idx if g in P.pos]
        Z = m4.zscore_columns(P.Y)
        score = m4.score_from_z(Z, sig_pos)
        if score is None:
            continue
        varies = ~np.isnan(Z).all(axis=0)
        match = [i for i in sig_pos if varies[i]]
        eligible = np.zeros(P.Y.shape[1], dtype=bool)
        eligible[:universe_n] = True
        eligible[sig_pos] = False
        sets = m4.matched_random_sets(P.Y.mean(axis=0), eligible, match, N_RANDOM, rng)
        rand = np.column_stack([np.nan_to_num(Z[:, s]).mean(axis=1) for s in sets])
        Yo = np.column_stack([m4.zs(score)] + [m4.zs(rand[:, j]) for j in range(rand.shape[1])])
        design = np.column_stack([np.ones(P.n), P.exposure[key].astype(float), P.ld, P.cycle])
        beta, t = m4.ols_t(Yo, design, 1)
        p_rand = (1 + (np.abs(t[1:]) >= abs(t[0])).sum()) / (1 + len(t[1:]))
        rows.append({"pid": P.pid, "beta": beta[0], "t": t[0], "p_vs_random": p_rand,
                     "rand_mean_t": t[1:].mean(), "lipid_genes_used": len(match)})
    per = pd.DataFrame(rows)
    if per.empty:
        return {"n_patients": 0, "mean_effect": np.nan, "ci_low": np.nan, "ci_high": np.nan, "share_positive": np.nan,
                "median_p_vs_random": np.nan, "meets_criterion": False}, per
    mean, lo, hi = t_interval(per["beta"])
    share = float((per["beta"] > 0).mean())
    med_p = float(per["p_vs_random"].median())
    summary = {"n_patients": len(per), "mean_effect": mean, "ci_low": lo, "ci_high": hi,
               "share_positive": share, "median_p_vs_random": med_p,
               "meets_criterion": bool(share >= E2_SHARE and lo > 0 and med_p < 0.05)}
    return summary, per


# %% populations
def prepare_population(label, pids, obs, X, counts_t, total_counts, genes, n_min_patients):
    pid_all = obs["pid"].astype(str).to_numpy()
    j = list(genes).index(TARGET)
    rows_by = {p: np.where(pid_all == p)[0] for p in pids}
    pooled = np.concatenate(list(rows_by.values()))
    frac = np.asarray((X[pooled] > 0).sum(axis=0)).ravel() / len(pooled)
    keep = (frac >= MIN_GENE_FRAC) & ~m3.damaged_mask(genes)
    keep[j] = False
    universe = np.where(keep)[0]
    gset = list(genes)
    cyc = [gset.index(g) for g in m3.CYCLE_GENES if g in gset]
    lip = [gset.index(g) for g in m4.SIGS["lipid"] if g in gset]
    extras = np.array(sorted(set(cyc + lip + [j]) - set(universe.tolist())), dtype=int)
    # the cycle genes are looked up by gene index inside build_patient
    gi = {"cycle": cyc, "target": j, "lipid": lip}
    patients = [build_patient(p, rows_by[p], X, counts_t, total_counts, universe, extras, gi) for p in pids]
    print(f"{label}: {len(patients)} patients, {len(universe)} universe genes, "
          f"{sum(P.n for P in patients)} cancer nuclei")
    return patients, universe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n_perm", type=int, default=N_PERM)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    a = ad.read_h5ad(args.cache)
    genes = np.array(a.var_names)
    X = sparse.csr_matrix(a.X)
    counts_t = to_dense(a.layers["counts"][:, list(genes).index(TARGET)]).ravel() > 0
    obs = a.obs
    total_counts = obs["total_counts"].to_numpy(dtype=float)
    pid_all = obs["pid"].astype(str).to_numpy()
    treat = obs.assign(pid=pid_all).groupby("pid")["treatment_status"].first()
    pids_all = sorted(treat.index)
    n_det = pd.Series(counts_t).groupby(pid_all).sum()
    n_tot = pd.Series(counts_t).groupby(pid_all).size()
    eligible_all = [p for p in pids_all if n_det[p] >= MIN_N and (n_tot[p] - n_det[p]) >= MIN_N]
    populations = {"untreated": [p for p in eligible_all if treat[p] == "Untreated"], "all_eligible": eligible_all}

    e1_rows, e2_rows, e1_pat, e2_pat = [], [], [], []
    for label, pids in populations.items():
        if len(pids) < MIN_PATIENTS:
            print(f"{label}: only {len(pids)} eligible patients, skipped (threshold {MIN_PATIENTS})")
            continue
        patients, universe = prepare_population(label, pids, obs, X, counts_t, total_counts, genes, MIN_PATIENTS)
        for key in ("primary", "A1"):
            use = [P for P in patients if P.ok[key]]
            if len(use) < MIN_PATIENTS:
                print(f"{label} / {key}: only {len(use)} patients qualify for this exposure, skipped")
                continue
            s1, per1, cons = run_e1(use, key, args.n_perm, rng)
            e1_rows.append({"population": label, "exposure": key, **s1})
            e1_pat.append(per1.assign(population=label, exposure=key))
            s2, per2 = run_e2(use, key, rng, len(universe))
            e2_rows.append({"population": label, "exposure": key, **s2})
            e2_pat.append(per2.assign(population=label, exposure=key))
            if label == "untreated" and key == "primary":
                order = np.argsort(-np.nan_to_num(cons, nan=-np.inf))[:TOP_N]
                pd.Series(genes[universe[order]]).to_csv(out / "consensus_program_top100_untreated.txt", index=False, header=False)
            print(f"done: {label} / {key}")

    e1, e2 = pd.DataFrame(e1_rows), pd.DataFrame(e2_rows)
    e1.round(4).to_csv(out / "e1_summary.tsv", sep="\t", index=False)
    e2.round(4).to_csv(out / "e2_summary.tsv", sep="\t", index=False)
    pd.concat(e1_pat).round(4).to_csv(out / "e1_patients.tsv", sep="\t", index=False)
    pd.concat(e2_pat).round(4).to_csv(out / "e2_patients.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n== E1: leave-one-tumor-out reproducibility of the IGF2BP2-linked program ==")
    print(e1.round(3).to_string(index=False))
    print("\n== E2: lipid score versus IGF2BP2 within cancer nuclei ==")
    print(e2.round(3).to_string(index=False))
    print(f"\nWrote files to {out}")


if __name__ == "__main__":
    main()