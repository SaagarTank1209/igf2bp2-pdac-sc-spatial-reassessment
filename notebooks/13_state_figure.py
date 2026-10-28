# %% [markdown]
# # 13 Figure 7: the state extension (GSE202051 cancer nuclei)
# Reads the per-split tables of step 12 and draws, for the primary run, how many of the 100 comparison genes IGF2BP2
# beat in each of the 20 splits (against the 95 required), and, for every run found, how many splits pass C0 to C3
# (against the 14 required).
#
# Run: python notebooks/13_state_figure.py --results DIR --outdir DIR     (DIR: docs/results/11_state and docs/figures)

# %% imports
import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out_name", default="fig7_state_extension.png")
    ap.add_argument("--title", default="State extension on GSE202051 cancer nuclei: the group repeats, but IGF2BP2 is not special")
    ap.add_argument("--left_title", default="Primary run: IGF2BP2 never beats 95 of the 100 matched genes")
    ap.add_argument("--open_label", default="open: split where IGF2BP2 picked a different cluster")
    ap.add_argument("--patients", type=int, default=18)
    ap.add_argument("--noun", default="untreated patients")
    ap.add_argument("--pop_label", default="untreated")
    a = ap.parse_args()
    RUNS = [(f"{a.pop_label}, res 0.3 (primary)", "untreated_main_r0.3"), (f"{a.pop_label}, res 0.2", "untreated_main_r0.2"),
            (f"{a.pop_label}, res 0.5", "untreated_main_r0.5"), ("26 patients, res 0.3", "all_eligible_main_r0.3")]
    res, out = Path(a.results), Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    prim = pd.read_csv(res / "state_splits_untreated_main_r0.3.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    x = prim["split"].to_numpy()
    main_cl = prim["igf2bp2_cluster"].mode().iloc[0]
    is_main = (prim["igf2bp2_cluster"] == main_cl).to_numpy()
    ax.axhline(95, color="k", ls="--", lw=1)
    ax.text(x.max() + 0.3, 95, "required: 95", va="bottom", ha="right", fontsize=8)
    ax.scatter(x[is_main] - 0.12, prim.loc[is_main, "comp_beaten_diff"], color="tab:blue", label="enrichment (C2)", zorder=3)
    ax.scatter(x[is_main] + 0.12, prim.loc[is_main, "comp_beaten_auc"], color="tab:orange", label="signature AUC (C3)", zorder=3)
    ax.scatter(x[~is_main] - 0.12, prim.loc[~is_main, "comp_beaten_diff"], facecolors="none", edgecolors="tab:blue", zorder=3)
    ax.scatter(x[~is_main] + 0.12, prim.loc[~is_main, "comp_beaten_auc"], facecolors="none", edgecolors="tab:orange",
               label=a.open_label, zorder=3)
    ax.set_ylim(0, 100)
    ax.set_xticks(x)
    ax.set_xlabel(f"split (random halves of the {a.patients} {a.noun})")
    ax.set_ylabel("comparison genes beaten by IGF2BP2 (of 100)")
    ax.set_title(a.left_title)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=8)

    ax = axes[1]
    runs = []
    for label, tag in RUNS:
        f = res / f"state_splits_{tag}.tsv"
        if f.exists():
            d = pd.read_csv(f, sep="\t")
            runs.append((label, [int(d[c].sum()) for c in ("C0", "C1", "C2", "C3")], len(d)))
    w = 0.8 / len(runs)
    for i, (label, counts, n) in enumerate(runs):
        xs = np.arange(4) + (i - (len(runs) - 1) / 2) * w
        ax.bar(xs, counts, w, label=label)
        for xv, cv in zip(xs, counts):
            ax.text(xv, cv + 0.2, str(cv), ha="center", va="bottom", fontsize=7)
    ax.axhline(14, color="k", ls="--", lw=1)
    ax.text(3.45, 14.3, "needed: 14 of 20", ha="right", fontsize=8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["C0\ngroup repeats", "C1\nspread over\npatients", "C2\nenrichment\nspecific", "C3\nsignature\nspecific"])
    ax.set_ylim(0, 21)
    ax.set_ylabel("splits passing (of 20)")
    ax.set_title("Splits passing each criterion")
    ax.legend(fontsize=7, loc="upper right")
    fig.suptitle(a.title, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out / a.out_name, dpi=150, bbox_inches="tight")
    print("wrote", out / a.out_name, "| runs found:", [r[0] for r in runs])


if __name__ == "__main__":
    main()
