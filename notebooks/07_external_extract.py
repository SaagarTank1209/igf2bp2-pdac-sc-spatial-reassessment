# %% [markdown]
# # 07 External data (GSE202051): extract cancer nuclei and check eligibility
# Reads the dataset authors' single-nucleus file, checks that the three cancer-cell label columns agree and that
# the counts layer holds integers, lists patients with their treatment status, and applies the eligibility rule
# locked in docs/input_manifest.md. It writes the cancer nuclei to a cache file for the later analysis.
# Only the NUMBER of nuclei with IGF2BP2 counts is read here, never any association.
#
# Run: python notebooks/07_external_extract.py --h5ad FILE --outdir DIR --cache_dir DIR

# %% imports
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

# %% locked settings (see the external check plan in docs/input_manifest.md)
TARGET = "IGF2BP2"
LABEL_COLUMN, LABEL_VALUE = "Level 1 Annotation", "Epithelial (malignant)"
OTHER_LABELS = [("new_celltypes", "Epithelial-Malignant"), ("Level 2 Annotation", "Malignant")]
MIN_DETECTED = 30       # cancer nuclei with IGF2BP2 counts, per patient
MIN_UNDETECTED = 30     # cancer nuclei without IGF2BP2 counts, per patient
MIN_UNTREATED_ELIGIBLE = 8
KEEP_OBS = ["pid", "sampleid", "treatment_status", "new_treatment", "response", "total_counts", "n_genes",
            "Level 1 Annotation", "Level 3 Annotation", "new_celltypes"]


def to_dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)


def eligibility(obs_mal, detected):
    """One row per patient: cancer nuclei, detected and undetected counts, treatment, eligibility."""
    d = obs_mal[["pid", "sampleid", "treatment_status"]].copy()
    d["detected"] = detected
    rows = []
    for pid, g in d.groupby("pid", observed=True):
        n, k = len(g), int(g["detected"].sum())
        rows.append({"pid": pid, "n_samples": g["sampleid"].nunique(),
                     "treatment_status": "/".join(sorted(map(str, g["treatment_status"].unique()))),
                     "n_cancer_nuclei": n, "n_detected": k, "n_undetected": n - k,
                     "eligible": bool(k >= MIN_DETECTED and (n - k) >= MIN_UNDETECTED)})
    return pd.DataFrame(rows)


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cache_dir", required=True)
    args = ap.parse_args()
    out, cache = Path(args.outdir), Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(args.h5ad)
    print("loaded", a.n_obs, "nuclei x", a.n_vars, "genes")
    obs = a.obs
    mal = (obs[LABEL_COLUMN] == LABEL_VALUE).to_numpy()
    print(f"cancer nuclei by {LABEL_COLUMN} == {LABEL_VALUE}: {int(mal.sum())}")
    for col, val in OTHER_LABELS:
        other = (obs[col] == val).to_numpy()
        print(f"  {col} == {val}: {int(other.sum())} nuclei, agrees with the primary label for "
              f"{int((other == mal).sum())} of {len(mal)} nuclei")

    cnt = a.layers["counts"]
    rows = np.random.default_rng(0).choice(a.n_obs, 2000, replace=False)
    vals = to_dense(cnt[rows]).ravel()
    vals = vals[vals != 0]
    print(f"counts layer: max {vals.max():.4g}, integer-valued among nonzeros of 2000 random nuclei: "
          f"{bool(np.allclose(vals, np.round(vals)))}")

    j = list(a.var_names).index(TARGET)
    detected = to_dense(cnt[:, j]).ravel() > 0
    elig = eligibility(obs.loc[mal], detected[mal])
    elig.to_csv(out / "eligibility.tsv", sep="\t", index=False)
    multi = int((elig["treatment_status"].str.contains("/")).sum())
    print(f"\npatients with at least one cancer nucleus: {len(elig)}; patients with more than one treatment status: {multi}")
    pd.set_option("display.width", 200)
    print(elig.groupby("treatment_status").agg(patients=("pid", "size"), eligible=("eligible", "sum"),
                                              nuclei=("n_cancer_nuclei", "sum")).to_string())
    untreated = elig[elig["treatment_status"] == "Untreated"]
    n_el = int(untreated["eligible"].sum())
    print(f"\nuntreated patients: {len(untreated)}; eligible untreated patients: {n_el} (threshold {MIN_UNTREATED_ELIGIBLE})")
    print("external check feasible:", n_el >= MIN_UNTREATED_ELIGIBLE)
    print("\nper patient (untreated):")
    print(untreated.sort_values("n_cancer_nuclei", ascending=False).to_string(index=False))

    try:
        a.raw = None
    except Exception:
        pass
    sub = a[mal].copy()
    sub.obs = sub.obs[[c for c in KEEP_OBS if c in sub.obs.columns]].copy()
    for key in list(sub.obsm.keys()):
        del sub.obsm[key]
    sub.write_h5ad(cache / "GSE202051_cancer_nuclei.h5ad", compression="gzip")
    print("\nwrote", cache / "GSE202051_cancer_nuclei.h5ad", "with", sub.n_obs, "nuclei")


if __name__ == "__main__":
    main()