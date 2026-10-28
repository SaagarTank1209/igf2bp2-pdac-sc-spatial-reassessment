# %% [markdown]
# # 14 Post-hoc sensitivity of the spatial result: without array A3, and dropping one tumor at a time
# Reads array_tests.tsv of step 04 and repeats the tumor level aggregation of step 04 (mean of the array effects per
# tumor, 95 percent t interval over tumors, median p against random gene sets over arrays). Post hoc: the locked
# analysis uses all nine arrays. Writes leave_one_out.tsv next to the input.
#
# Run: python notebooks/14_leave_one_out.py --results docs/results/04_spatial

# %% imports
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist

OUTCOMES = ["C1_stress", "H1_lipid", "H2_macrophage", "S1_lam", "S2_exhaustion"]


def summarize(df):
    per = df.groupby("tumor")["beta"].mean()
    n = len(per)
    half = t_dist.ppf(0.975, n - 1) * per.std(ddof=1) / np.sqrt(n)
    p = df["p_vs_random"].dropna()
    return {"tumors": n, "tumors_positive": int((per > 0).sum()), "mean": per.mean(),
            "ci_low": per.mean() - half, "ci_high": per.mean() + half,
            "median_p_vs_random": p.median() if len(p) else np.nan}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    args = ap.parse_args()
    folder = Path(args.results)
    a = pd.read_csv(folder / "array_tests.tsv", sep="\t")
    rows = []
    for o in OUTCOMES:
        sub = a[a["outcome"] == o]
        if sub.empty:
            continue
        analyses = [("all arrays (locked)", sub), ("without A3", sub[sub["array"] != "A3"])]
        analyses += [(f"without tumor {t}", sub[sub["tumor"] != t]) for t in sorted(sub["tumor"].unique())]
        for name, df in analyses:
            rows.append({"outcome": o, "analysis": name, **summarize(df)})
    out = pd.DataFrame(rows)
    out.round(4).to_csv(folder / "leave_one_out.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(out.round(3).to_string(index=False))
    for o in ("H1_lipid", "H2_macrophage"):
        d = out[out["outcome"] == o]
        a3 = d[d["analysis"] == "without A3"].iloc[0]
        one = d[d["analysis"].str.startswith("without tumor")]
        print(f"\n{o}: without A3 mean {a3['mean']:.3f} (interval {a3['ci_low']:.3f} to {a3['ci_high']:.3f}); "
              f"dropping one tumor: means {one['mean'].min():.3f} to {one['mean'].max():.3f}; "
              f"intervals that exclude zero: {int(((one['ci_low'] > 0) | (one['ci_high'] < 0)).sum())} of {len(one)}")
    print("\nwrote", folder / "leave_one_out.tsv")


if __name__ == "__main__":
    main()
