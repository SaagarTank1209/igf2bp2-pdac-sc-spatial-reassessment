# %% [markdown]
# # 12 State extension, part 2: held-out test of the IGF2BP2-enriched cluster (GSE202051 cancer nuclei)
# Implements criteria C1 to C4 of the state extension plan in docs/input_manifest.md. Patients are split at random
# into discovery and validation halves (20 splits, fixed seeds). For IGF2BP2 and each of 100 comparison genes the
# enriched cluster and a 100-gene signature are found in discovery patients and tested in validation patients.
#
# Run: python notebooks/12_state_test.py --cache FILE --outdir DIR --population untreated --variant main --resolution r0.3

# %% imports
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


def load_module(name, fname):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m11 = load_module("step11", "11_state_map.py")
m8, m3, m4 = m11.m8, m11.m3, m11.m8.m4

# %% locked settings (state extension plan)
N_SPLITS = 20
MIN_IN = 30              # nuclei inside and outside a cluster (and with and without the gene), per patient
MIN_PATIENTS = 3
SHARE_MIN = 0.05
C1_SHARE = 0.60
RANK_SHARE = 0.95
SPLIT_SHARE = 0.70
SIG_N = 100


# %% patients
def prepare(obs, X, counts, counts_t, total_counts, genes, pids, keep, labels, expo_idx):
    pid_all = obs["pid"].astype(str).to_numpy()
    j = list(genes).index(m11.TARGET)
    rows_by = {p: np.where((pid_all == p) & keep)[0] for p in pids}
    pooled = np.concatenate(list(rows_by.values()))
    frac = np.asarray((X[pooled] > 0).sum(axis=0)).ravel() / len(pooled)
    ok = (frac >= m8.MIN_GENE_FRAC) & ~m3.damaged_mask(genes)
    ok[j] = False
    universe = np.where(ok)[0]
    gset = list(genes)
    cyc = [gset.index(g) for g in m3.CYCLE_GENES if g in gset]
    lip = [gset.index(g) for g in m4.SIGS["lipid"] if g in gset]
    extras = np.array(sorted(set(cyc + lip + [j]) - set(universe.tolist())), dtype=int)
    gi = {"cycle": cyc, "target": j, "lipid": lip}
    patients = []
    for p in pids:
        rows = rows_by[p]
        P = m8.build_patient(p, rows, X, counts_t, total_counts, universe, extras, gi)
        P.cl = labels.loc[obs.index[rows]].to_numpy()
        P.det = (counts[rows][:, expo_idx] > 0).toarray()
        patients.append(P)
    return patients, universe


# %% per patient and cluster summaries
def summarize(patients, clusters, n_genes):
    Np, K, U = len(patients), len(clusters), patients[0].U
    valid = np.zeros((Np, K), dtype=bool)
    share = np.zeros((Np, K))
    Dm = np.zeros((Np, K, U), dtype=np.float32)
    Fd = np.full((Np, K, n_genes), np.nan)
    for i, P in enumerate(patients):
        YU = P.Y[:, :P.U]
        tot = YU.sum(axis=0, dtype=np.float64)
        for k, c in enumerate(clusters):
            m = P.cl == c
            nin = int(m.sum())
            share[i, k] = nin / P.n
            if nin >= MIN_IN and P.n - nin >= MIN_IN:
                valid[i, k] = True
                Dm[i, k] = YU[m].mean(axis=0, dtype=np.float64) - (tot - YU[m].sum(axis=0, dtype=np.float64)) / (P.n - nin)
                Fd[i, k] = P.det[m].mean(axis=0) - P.det[~m].mean(axis=0)
    return valid, share, Dm, Fd


def signature(Dm, valid, disc, k, exclude_pos):
    d = [i for i in disc if valid[i, k]]
    diffs = Dm[d, k, :].astype(np.float64)
    sd = diffs.std(axis=0, ddof=1)
    t = np.where(sd > 0, diffs.mean(axis=0) / np.where(sd > 0, sd, 1.0) * np.sqrt(len(d)), 0.0)
    if exclude_pos is not None and exclude_pos >= 0:
        t[exclude_pos] = -np.inf
    return np.argsort(-t)[:SIG_N]


# %% the splits
def make_split(n_patients, s, strata=None):
    """Seeded split into discovery and validation halves; with strata, each group is split evenly (odd one to validation)."""
    rng = np.random.default_rng(1000 + s)
    if strata is None:
        perm = rng.permutation(n_patients)
        return np.sort(perm[:n_patients // 2]), np.sort(perm[n_patients // 2:])
    disc, val = [], []
    for g in np.unique(strata):
        idx = rng.permutation(np.where(strata == g)[0])
        h = len(idx) // 2
        disc += list(idx[:h])
        val += list(idx[h:])
    return np.sort(disc), np.sort(val)


def run_splits(patients, clusters, valid, share, Dm, Fd, expo, uni_pos, n_splits, strata=None):
    Np, G = len(patients), len(expo)
    cache = {}
    rows = []
    for s in range(n_splits):
        disc, val = make_split(Np, s, strata)
        for g in range(G):
            best, best_score = None, -np.inf
            for k in range(len(clusters)):
                d = [i for i in disc if valid[i, k]]
                if len(d) >= MIN_PATIENTS:
                    v = np.nanmean(Fd[d, k, g])
                    if v > best_score:
                        best, best_score = k, v
            rec = {"split": s, "gene": expo[g], "cluster": np.nan, "disc_diff": np.nan, "val_diff": np.nan,
                   "val_share_ok": np.nan, "val_auc": np.nan, "val_auc_share": np.nan}
            if best is not None:
                rec["cluster"], rec["disc_diff"] = clusters[best], best_score
                vv = [i for i in val if valid[i, best]]
                if len(vv) >= MIN_PATIENTS:
                    rec["val_diff"] = float(np.nanmean(Fd[vv, best, g]))
                rec["val_share_ok"] = float(np.mean(share[val, best] >= SHARE_MIN))
                top = signature(Dm, valid, disc, best, uni_pos[g])
                key = (s, best, tuple(top))
                if key not in cache:
                    cache[key] = {}
                    for i in val:
                        P = patients[i]
                        sc = ((P.Y[:, top] - P.mu[top]) * P.inv_sd[top]).mean(axis=1)
                        cache[key][i] = sc - P.D3 @ (P.pinv @ sc)
                aucs = []
                for i in val:
                    e = patients[i].det[:, g]
                    if e.sum() >= MIN_IN and (~e).sum() >= MIN_IN:
                        aucs.append(m3.auc(cache[key][i], e))
                if len(aucs) >= MIN_PATIENTS:
                    rec["val_auc"], rec["val_auc_share"] = float(np.mean(aucs)), float(np.mean(np.array(aucs) > 0.5))
            rows.append(rec)
    return pd.DataFrame(rows)


def evaluate(res, n_splits):
    out = []
    for s in range(n_splits):
        d = res[res["split"] == s].reset_index(drop=True)
        g0, comp = d.iloc[0], d.iloc[1:]
        n_comp = len(comp)
        c0 = bool(pd.notna(g0["val_diff"]) and g0["val_diff"] > 0 and pd.notna(g0["val_auc"]) and g0["val_auc"] > 0.5
                  and g0["val_auc_share"] >= C1_SHARE)
        c1 = bool(g0["val_share_ok"] >= C1_SHARE) if pd.notna(g0["val_share_ok"]) else False
        beats_diff = int((g0["val_diff"] > comp["val_diff"]).sum()) if pd.notna(g0["val_diff"]) else 0
        c2 = bool(beats_diff / n_comp >= RANK_SHARE)
        beats_auc = int((g0["val_auc"] > comp["val_auc"]).sum()) if pd.notna(g0["val_auc"]) else 0
        c3 = bool(pd.notna(g0["val_auc"]) and g0["val_auc"] > 0.5 and beats_auc / n_comp >= RANK_SHARE
                  and g0["val_auc_share"] >= C1_SHARE)
        out.append({"split": s, "igf2bp2_cluster": g0["cluster"], "val_share_ok": g0["val_share_ok"],
                    "val_diff": g0["val_diff"], "comp_beaten_diff": beats_diff, "val_auc": g0["val_auc"],
                    "comp_beaten_auc": beats_auc, "n_comparison": n_comp, "C0": c0, "C1": c1, "C2": c2, "C3": c3})
    return pd.DataFrame(out)


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--population", choices=["untreated", "all_eligible"], default="untreated")
    ap.add_argument("--variant", choices=["main", "S1"], default="main")
    ap.add_argument("--resolution", default="r0.3")
    ap.add_argument("--n_splits", type=int, default=N_SPLITS)
    args = ap.parse_args()
    out = Path(args.outdir)
    tag = f"{args.population}_{args.variant}_{args.resolution}" + ("" if args.n_splits == N_SPLITS else f"_n{args.n_splits}")

    a, genes, X, counts, obs, counts_t = m11.load_cache(args.cache)
    total_counts = obs["total_counts"].to_numpy(dtype=float)
    expo = [g.strip() for g in (out / f"exposure_genes_{args.population}.txt").read_text().splitlines() if g.strip()]
    gl = list(genes)
    expo_idx = np.array([gl.index(g) for g in expo])
    tab = pd.read_csv(out / f"clusters_{args.population}_{args.variant}.tsv", sep="\t", index_col=0)
    if args.resolution not in tab.columns:
        print(f"{args.resolution} is not in the cluster table (variant not formed); nothing to test")
        return
    labels = tab[args.resolution].astype(str)
    keep = obs.index.isin(labels.index)
    pids = m11.population_pids(obs, counts_t, args.population)
    patients, universe = prepare(obs, X, counts, counts_t, total_counts, genes, pids, keep, labels, expo_idx)
    treat = obs.assign(pid=obs["pid"].astype(str)).groupby("pid")["treatment_status"].first()
    strata = np.array([0 if treat[P.pid] == "Untreated" else 1 for P in patients]) if args.population == "all_eligible" else None
    uni_map = {int(g): i for i, g in enumerate(universe)}
    uni_pos = [uni_map.get(int(g), -1) for g in expo_idx]
    clusters = sorted(set(np.concatenate([P.cl for P in patients])))
    print(f"{tag}: {len(patients)} patients, {len(clusters)} clusters, {len(expo)} exposure genes, "
          f"{sum(P.n for P in patients)} nuclei")

    valid, share, Dm, Fd = summarize(patients, clusters, len(expo))
    res = run_splits(patients, clusters, valid, share, Dm, Fd, expo, uni_pos, args.n_splits, strata)
    res.round(4).to_csv(out / f"state_genes_{tag}.tsv", sep="\t", index=False)
    ev = evaluate(res, args.n_splits)
    ev.round(4).to_csv(out / f"state_splits_{tag}.tsv", sep="\t", index=False)

    # composition of every cluster (descriptive)
    k_igf = ev["igf2bp2_cluster"].mode().iloc[0] if ev["igf2bp2_cluster"].notna().any() else None
    comp_rows = []
    for k, c in enumerate(clusters):
        n_c = sum(int((P.cl == c).sum()) for P in patients)
        sh = np.array([(P.cl == c).mean() for P in patients])
        per_patient_n = np.array([(P.cl == c).sum() for P in patients])
        comp_rows.append({"cluster": c, "nuclei": n_c, "patients_with_5pct": int((sh >= SHARE_MIN).sum()),
                          "max_patient_share_of_cluster": float(per_patient_n.max() / max(n_c, 1)),
                          "times_igf2bp2_cluster": int((ev["igf2bp2_cluster"] == c).sum())})
    pd.DataFrame(comp_rows).round(4).to_csv(out / f"state_clusters_{tag}.tsv", sep="\t", index=False)

    summary = {"population": args.population, "variant": args.variant, "resolution": args.resolution,
               "patients": len(patients), "splits": args.n_splits, "C0_rate": ev["C0"].mean(), "C1_rate": ev["C1"].mean(),
               "C2_rate": ev["C2"].mean(), "C3_rate": ev["C3"].mean()}
    for c in ("C0", "C1", "C2", "C3"):
        summary[c + "_splits_passing"] = int(ev[c].sum())
        summary[c + "_holds"] = bool(summary[c + "_rate"] >= SPLIT_SHARE)
    summary["state_supported"] = bool(summary["C1_holds"] and summary["C2_holds"] and summary["C3_holds"])
    if summary["state_supported"] and summary["C0_holds"]:
        summary["outcome"] = "A: the group repeats and IGF2BP2 is special (strong support)"
    elif summary["state_supported"]:
        summary["outcome"] = "inconsistent (supported but C0 fails): inspect the per-split table"
    elif summary["C0_holds"] and summary["C1_holds"]:
        summary["outcome"] = "B: the group repeats but IGF2BP2 is not special (IGF2BP2 is one marker among many)"
    else:
        summary["outcome"] = "C: the group does not repeat (no support in this dataset)"
    summary["modal_igf2bp2_cluster"] = k_igf
    summary["mean_val_auc_igf2bp2"] = float(ev["val_auc"].mean())
    summary["mean_comparison_beaten_auc"] = float(ev["comp_beaten_auc"].mean())

    if summary["state_supported"]:
        ks = clusters.index(k_igf)
        use = [P for i, P in enumerate(patients) if valid[i, ks]]
        for P in use:
            P.exposure["state"] = (P.cl == k_igf)
        s4, per4 = m8.run_e2(use, "state", np.random.default_rng(0), len(universe))
        per4.round(4).to_csv(out / f"state_c4_patients_{tag}.tsv", sep="\t", index=False)
        summary.update({"C4_patients": s4["n_patients"], "C4_mean_effect": s4["mean_effect"], "C4_ci_low": s4["ci_low"],
                        "C4_ci_high": s4["ci_high"], "C4_share_positive": s4["share_positive"],
                        "C4_median_p": s4["median_p_vs_random"], "C4_holds": s4["meets_criterion"]})
    else:
        summary["C4_holds"] = "not evaluated (state not supported)"
    pd.Series(summary).to_csv(out / f"state_summary_{tag}.tsv", sep="\t", header=False)

    pd.set_option("display.width", 200)
    print("\n== clusters ==")
    print(pd.DataFrame(comp_rows).round(3).to_string(index=False))
    print("\n== per split (IGF2BP2) ==")
    print(ev.round(3).to_string(index=False))
    print("\n== summary ==")
    for k, v in summary.items():
        print(f"  {k}: {round(v, 4) if isinstance(v, float) else v}")


if __name__ == "__main__":
    main()