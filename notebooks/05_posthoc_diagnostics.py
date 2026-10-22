# %% [markdown]
# # 05 Post-hoc diagnostics for the spatial analysis (not part of the locked primary analysis)
# 1. Spike-in: add a known shift to the lipid score of IGF2BP2-detected spots and find the smallest shift at
#    which the locked criterion would have been met on these data.
# 2. Split-half plumbing control: split a signature's genes into two random halves and test whether one half
#    predicts the other through the same models. Shows whether real co-variation survives the adjustments.
# 3. Signed random-set statistics: shows whether matched random gene sets have a systematic bias in each model.
# Everything here is post-hoc and is reported as such.
#
# Run: python notebooks/05_posthoc_diagnostics.py --h5ad_dir DIR --outdir DIR

# %% imports
import argparse
import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import t as t_dist

# reuse the helpers and the locked settings of step 04
_spec = importlib.util.spec_from_file_location("step04", Path(__file__).with_name("04_spatial_tests.py"))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75]   # shifts in SD of the lipid score
N_SPLITS = 30
SEED = 1


# %% per-array preparation (same formulas as step 04)
def prepare(adata, rng):
    counts = m.to_dense(adata.X).astype(float)
    genes = list(adata.var_names)
    gi = {g: i for i, g in enumerate(genes)}
    total = counts.sum(axis=1)
    logx = np.log1p(counts / total[:, None] * 1e4)
    Z = m.zscore_columns(logx)
    xy = np.asarray(adata.obsm["spatial"])
    sig_idx = {n: [gi[g] for g in s if g in gi] for n, s in m.SIGS.items()}
    scores = {n: m.score_from_z(Z, idx) for n, idx in sig_idx.items()}
    varies = ~np.isnan(Z).all(axis=0)
    eligible = ((counts > 0).mean(axis=0) >= 0.01) & varies
    for i in set(i for idx in sig_idx.values() for i in idx) | {gi[m.TARGET]}:
        eligible[i] = False
    mean_expr = logx.mean(axis=0)

    def random_scores(name):
        match = [i for i in sig_idx[name] if varies[i]]
        sets = m.matched_random_sets(mean_expr, eligible, match, m.N_RANDOM, rng)
        return np.column_stack([np.nan_to_num(Z[:, s]).mean(axis=1) for s in sets])

    return {"Z": Z, "scores": scores, "sig_idx": sig_idx, "varies": varies,
            "blocks": m.block_ids(xy, m.BLOCK), "log_umi": m.zs(np.log(total)),
            "epi": m.zs(scores["epithelial"]), "stro": m.zs(scores["stromal"]),
            "exposed": counts[:, gi[m.TARGET]] > 0,
            "random": {n: random_scores(n) for n in ("lipid", "macrophage", "stress")}}


def tumor_level(df, value):
    per = df.groupby("tumor")[value].mean()
    n = len(per)
    half = t_dist.ppf(0.975, n - 1) * per.std(ddof=1) / np.sqrt(n)
    return {"n_tumors": n, "n_positive": int((per > 0).sum()), "mean": round(per.mean(), 3),
            "ci_low": round(per.mean() - half, 3), "ci_high": round(per.mean() + half, 3)}


# %% 1. spike-in
def spike_in(a, rng):
    covs = [a["log_umi"], a["epi"], a["stro"]]
    base = m.zs(a["scores"]["lipid"])
    rows = []
    for d in DELTAS:
        res = m.run_group(a["exposed"], covs, {"H1_lipid": base + d * a["exposed"]},
                          {"H1_lipid": a["random"]["lipid"]}, a["blocks"], rng)
        r = res[res["outcome"] == "H1_lipid"].iloc[0]
        rows.append({"delta": d, "outcome": "H1_lipid", "beta": r["beta"], "p_vs_random": r["p_vs_random"]})
    return pd.DataFrame(rows)


# %% 2. split-half plumbing control
def split_half(a, name, covs, rng):
    idx = [i for i in a["sig_idx"][name] if a["varies"][i]]
    betas = []
    for _ in range(N_SPLITS):
        perm = rng.permutation(idx)
        h = len(perm) // 2
        x = np.nan_to_num(a["Z"][:, perm[:h]]).mean(axis=1)
        y = np.nan_to_num(a["Z"][:, perm[h:]]).mean(axis=1)
        design = np.column_stack([np.ones(len(x)), m.zs(x)] + covs)
        beta, _ = m.ols_t(m.zs(y)[:, None], design, 1)
        betas.append(beta[0])
    return float(np.median(betas))


# %% 3. signed random-set statistics
def random_bias(a, name, outcome, exposure, covs):
    mat = a["random"][name]
    Y = np.column_stack([m.zs(outcome)] + [m.zs(mat[:, j]) for j in range(mat.shape[1])])
    design = np.column_stack([np.ones(len(exposure)), exposure.astype(float)] + covs)
    _, t = m.ols_t(Y, design, 1)
    rt = t[1:]
    return {"t_real": t[0], "rand_mean_t": rt.mean(), "rand_sd_t": rt.std(),
            "z_vs_random": (t[0] - rt.mean()) / rt.std(), "frac_random_t_positive": float((rt > 0).mean())}


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad_dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    spike_rows, split_rows, bias_rows = [], [], []
    for arr in m.ARRAYS:
        a = prepare(ad.read_h5ad(Path(args.h5ad_dir) / f"st_{arr}.h5ad"), rng)
        tumor = m.TUMOR_OF[arr]
        s = spike_in(a, rng)
        s.insert(0, "array", arr)
        s["tumor"] = tumor
        spike_rows.append(s)
        umi, epi, stro = a["log_umi"], a["epi"], a["stro"]
        for name, covs in (("lipid", [umi, epi, stro]), ("macrophage", [umi, epi, stro]),
                           ("stromal", [umi, epi]), ("stress", [umi, epi])):
            split_rows.append({"array": arr, "tumor": tumor, "signature": name,
                               "median_beta": split_half(a, name, covs, rng)})
        for model, name, outcome, exposure, covs in (
                ("H1_lipid", "lipid", a["scores"]["lipid"], a["exposed"], [umi, epi, stro]),
                ("H2_macrophage", "macrophage", a["scores"]["macrophage"], a["exposed"], [umi, epi, stro]),
                ("C1_stress", "stress", a["scores"]["stress"], stro, [umi, epi])):
            b = random_bias(a, name, outcome, exposure, covs)
            bias_rows.append({"array": arr, "tumor": tumor, "model": model, **b})
        print(f"{arr} done")

    spike = pd.concat(spike_rows, ignore_index=True)
    spike.round(4).to_csv(out / "spike_in_arrays.tsv", sep="\t", index=False)
    srows = []
    for d, g in spike.groupby("delta"):
        summ = m.tumor_summary(g)
        srows.append(summ.assign(delta=d))
    spike_summary = pd.concat(srows, ignore_index=True)[
        ["delta", "n_positive", "mean_beta", "ci_low", "ci_high", "median_p_vs_random",
         "supports_positive", "meets_criterion"]]
    spike_summary.to_csv(out / "spike_in_summary.tsv", sep="\t", index=False)

    split = pd.DataFrame(split_rows)
    split.round(4).to_csv(out / "split_half_arrays.tsv", sep="\t", index=False)
    split_summary = pd.DataFrame([{"signature": n, **tumor_level(g, "median_beta")} for n, g in split.groupby("signature")])
    split_summary.to_csv(out / "split_half_summary.tsv", sep="\t", index=False)

    bias = pd.DataFrame(bias_rows)
    bias.round(4).to_csv(out / "random_set_bias_arrays.tsv", sep="\t", index=False)
    bias_summary = pd.DataFrame([{"model": n, "mean_rand_t": round(g["rand_mean_t"].mean(), 3),
                                  "mean_rand_sd_t": round(g["rand_sd_t"].mean(), 3),
                                  **{"z_" + k: v for k, v in tumor_level(g, "z_vs_random").items()}}
                                 for n, g in bias.groupby("model")])
    bias_summary.to_csv(out / "random_set_bias_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 220)
    base = spike_summary[spike_summary["delta"] == 0].iloc[0]
    se = (base["ci_high"] - base["ci_low"]) / 2 / t_dist.ppf(0.975, 4)
    mde = (t_dist.ppf(0.975, 4) + t_dist.ppf(0.8, 4)) * se
    print("\n== 1. spike-in: shift added to the lipid score of IGF2BP2-detected spots (SD of the score) ==")
    print(spike_summary.to_string(index=False))
    print(f"approximate minimum detectable tumor-level effect at 80% power (t test over 5 tumors, ignoring the "
          f"random-set condition): {mde:.2f} SD")
    print("\n== 2. split-half plumbing control: tumor-level effect of one random half of a signature on the other ==")
    print(split_summary.to_string(index=False))
    print("\n== 3. signed random-set statistics (t of matched random gene sets in each model) ==")
    print(bias_summary.to_string(index=False))
    print("\nper array, mean t of random sets:")
    print(bias.pivot_table(index="model", columns="array", values="rand_mean_t").round(2)[
        [c for c in m.ARRAYS]].to_string())
    print(f"\nWrote files to {out}")


if __name__ == "__main__":
    main()