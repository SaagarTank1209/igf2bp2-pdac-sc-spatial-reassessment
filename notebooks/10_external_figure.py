# %% [markdown]
# # 10 Figure 6: the external check in one picture
# Reads the tables of steps 08 and 09 (untreated patients, primary exposure) and draws four panels:
# IGF2BP2 against detection-matched random genes, per-patient AUCs, the robustness variants, and the lipid effect.
#
# Run: python notebooks/10_external_figure.py --results DIR --outdir DIR
# DIR holds the folders 08_external and 09_external (docs/results in the repository, or the results folder on scratch).

# %% imports
import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read(res, folder, name):
    return pd.read_csv(Path(res) / folder / name, sep="\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    var = read(a.results, "09_external", "e1_variants.tsv")
    rnd = read(a.results, "09_external", "e1_matched_random_exposures.tsv")
    e1 = read(a.results, "08_external", "e1_patients.tsv")
    e2 = read(a.results, "08_external", "e2_patients.tsv")
    e1 = e1[(e1["population"] == "untreated") & (e1["exposure"] == "primary")].sort_values("auc")
    e2 = e2[(e2["population"] == "untreated") & (e2["exposure"] == "primary")].sort_values("beta")
    v0 = var.iloc[0]
    igf = float(v0["mean_auc"])
    n_at_or_above = int((rnd["mean_auc"] >= igf).sum())

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    ax.hist(rnd["mean_auc"], bins=16, color="lightgrey", edgecolor="grey")
    ax.axvline(igf, color="tab:red", lw=2, label=f"IGF2BP2 ({igf:.3f})")
    ax.axvline(float(v0["null_mean"]), color="k", ls="--", lw=1, label=f"shuffled-label null ({float(v0['null_mean']):.3f})")
    ax.set_xlabel("mean held-out AUC of the program")
    ax.set_ylabel(f"random genes (n = {len(rnd)})")
    ax.set_title(f"IGF2BP2 versus detection-matched random genes\n{n_at_or_above} of {len(rnd)} reach or exceed IGF2BP2")
    ax.legend()

    ax = axes[0, 1]
    x = np.arange(len(e1))
    ax.scatter(x, e1["auc"], color="tab:blue", zorder=3)
    ax.axhline(0.5, color="k", lw=0.8)
    ax.axhline(igf, color="tab:red", ls="--", lw=1)
    ax.set_xticks([])
    ax.set_xlabel("untreated patients, sorted")
    ax.set_ylabel("held-out AUC of the IGF2BP2 program")
    ax.set_title(f"E1 per patient: {int((e1['auc'] > 0.5).sum())} of {len(e1)} above 0.5 (dashed: mean)")

    ax = axes[1, 0]
    short = {"V0": "V0 baseline", "V1": "V1 + detected genes", "V2": "V2 + lncRNA control", "V3": "V3 both"}
    names = [short.get(str(v)[:2], str(v)[:30]) for v in var["variant"]]
    y = np.arange(len(var))[::-1]
    ax.axvspan(rnd["mean_auc"].mean() - rnd["mean_auc"].std(), rnd["mean_auc"].mean() + rnd["mean_auc"].std(),
               color="lightgrey", alpha=0.6, label="random genes, mean ± SD")
    ax.errorbar(var["mean_auc"], y, xerr=[var["mean_auc"] - var["ci_low"], var["ci_high"] - var["mean_auc"]],
                fmt="o", color="tab:blue", capsize=4, label="IGF2BP2, 95% interval over patients")
    ax.axvline(0.5, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("mean held-out AUC")
    ax.set_title("Robustness variants (post-hoc): all overlap the random genes")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=8)

    ax = axes[1, 1]
    x = np.arange(len(e2))
    ax.scatter(x, e2["beta"], color="tab:green", zorder=3)
    m = float(e2["beta"].mean())
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(m, color="tab:green", ls="--", lw=1)
    ax.set_xticks([])
    ax.set_xlabel("untreated patients, sorted")
    ax.set_ylabel("lipid score effect (SD) of IGF2BP2 detection")
    ax.set_title(f"E2 lipid score: {int((e2['beta'] > 0).sum())} of {len(e2)} positive, mean {m:.3f} (dashed)")

    fig.suptitle("External check (GSE202051, untreated patients): an association that is not specific to IGF2BP2, and no lipid link",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "fig6_external_check.png", dpi=150, bbox_inches="tight")
    print("wrote", out / "fig6_external_check.png")


if __name__ == "__main__":
    main()