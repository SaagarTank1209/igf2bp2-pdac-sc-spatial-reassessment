# %% [markdown]
# # 06 Figures
# Builds the report figures from the outputs of steps 01 to 05. Nothing is recomputed here except
# class-level summaries of tables that earlier steps wrote.
#
# Run: python notebooks/06_figures.py --results DIR --labels FILE --outdir DIR
# DIR holds the folders 02_annotation, 03_state, 04_spatial and 05_posthoc.

# %% imports
import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt

TUMORS = ["A", "B", "D", "E", "G"]


def read(res, folder, name, **kw):
    return pd.read_csv(Path(res) / folder / name, sep="\t", **kw)


# %% figure 1: single-cell composition and IGF2BP2 detection
def fig_cell_types(res, labels, out):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    groups = ["candidate cancer", "other epithelial", "all other cells"]
    det = {}
    for ax, tumor in zip(axes[:2], ("A", "B")):
        cl = read(res, "02_annotation", f"clusters_{tumor}.tsv", index_col=0)["leiden"].astype(int)
        cls = cl.map(labels[tumor])
        counts = cls.value_counts().sort_values()
        ax.barh(counts.index, counts.values, color="tab:blue")
        for y, v in enumerate(counts.values):
            ax.text(v, y, f" {v}", va="center", fontsize=8)
        ax.set_title(f"tumor {tumor}: {len(cl)} cells (provisional labels)")
        ax.set_xlabel("cells")
        frac = read(res, "02_annotation", f"target_gene_fraction_{tumor}.tsv", index_col=0).loc["IGF2BP2"]
        frac.index = frac.index.astype(int)
        sizes = cl.value_counts()
        def grp(c):
            c = labels[tumor].get(c)
            return groups[0] if c == "malignant" else groups[1] if c in ("epithelial", "tuft") else groups[2]
        tot, num = dict.fromkeys(groups, 0.0), dict.fromkeys(groups, 0.0)
        for k, n in sizes.items():
            tot[grp(k)] += n
            num[grp(k)] += n * frac[k]
        det[tumor] = [num[g] / tot[g] if tot[g] else np.nan for g in groups]
    x = np.arange(len(groups))
    axes[2].bar(x - 0.2, det["A"], 0.4, label="tumor A", color="tab:blue")
    axes[2].bar(x + 0.2, det["B"], 0.4, label="tumor B", color="tab:orange")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(groups, rotation=15)
    axes[2].set_ylabel("fraction of cells with IGF2BP2 counts")
    axes[2].set_title("IGF2BP2 detection by cell group")
    axes[2].legend()
    fig.tight_layout()
    fig.savefig(Path(out) / "fig1_single_cell.png", dpi=150)
    plt.close(fig)


# %% figure 2: the reproducibility gate
def fig_gate(res, out):
    st = read(res, "03_state", "state_stats.tsv")
    rp = read(res, "03_state", "reproducibility.tsv").set_index("statistic")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1.1, 1]})
    both = st["top_A"] & st["top_B"]
    ax1.scatter(st["t_A"], st["t_B"], s=4, c="lightgrey")
    ax1.scatter(st.loc[st["top_A"] & ~both, "t_A"], st.loc[st["top_A"] & ~both, "t_B"], s=10, c="tab:blue", label="top 100 in A")
    ax1.scatter(st.loc[st["top_B"] & ~both, "t_A"], st.loc[st["top_B"] & ~both, "t_B"], s=10, c="tab:orange", label="top 100 in B")
    ax1.scatter(st.loc[both, "t_A"], st.loc[both, "t_B"], s=16, c="tab:red", label="top in both")
    ax1.axhline(0, c="k", lw=0.5)
    ax1.axvline(0, c="k", lw=0.5)
    ax1.set_xlabel("t statistic, IGF2BP2 detected vs not, tumor A")
    ax1.set_ylabel("t statistic, tumor B")
    ax1.set_title(f"{len(st)} genes; Spearman rho = {rp.loc['spearman_t', 'observed']:.2f}")
    ax1.legend()
    order = ["top_overlap", "auc_A_to_B", "auc_B_to_A", "spearman_t"]
    z = [(rp.loc[s, "observed"] - rp.loc[s, "null_mean"]) / rp.loc[s, "null_sd"] for s in order]
    gate_passed = all(z[order.index(s)] > 1.645 for s in ("top_overlap", "auc_A_to_B", "auc_B_to_A"))
    ax2.barh(order, z, color="tab:grey")
    ax2.axvline(1.645, ls="--", c="k", lw=0.8)
    for i, s in enumerate(order):
        ax2.text(max(z[i], 0) + 0.05, i, f"p = {rp.loc[s, 'p_upper']:.3f}", va="center", fontsize=9)
    ax2.set_xlabel("observed minus null mean, in null SDs\n(dashed line: one-sided 5 percent)")
    ax2.set_title("Reproducibility gate: " + ("PASSED" if gate_passed else "NOT passed")
                  + "\nneeds overlap and both AUCs beyond the dashed line")
    ax2.set_xlim(min(-1, min(z) - 0.5), max(3, max(z) + 1.5))
    fig.tight_layout()
    fig.savefig(Path(out) / "fig2_reproducibility_gate.png", dpi=150)
    plt.close(fig)


# %% figure 3: spot maps of two arrays chosen by a fixed rule
def fig_spatial_maps(res, out):
    sp = read(res, "04_spatial", "spot_scores.tsv")
    sp["IGF2BP2_detected"] = sp["IGF2BP2_detected"].astype(bool)
    n_exposed = sp.groupby("array")["IGF2BP2_detected"].sum()
    chosen = [n_exposed.idxmax(), n_exposed.idxmin()]   # most and fewest IGF2BP2-detected spots
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for row, arr in zip(axes, chosen):
        d = sp[sp["array"] == arr]
        for ax, col, title in zip(row[:2], ("epithelial", "lipid"), ("epithelial score", "lipid score")):
            s = ax.scatter(d["x"], d["y"], c=d[col], s=22, cmap="viridis")
            fig.colorbar(s, ax=ax, fraction=0.046)
            ax.set_title(f"array {arr}: {title}")
        ax = row[2]
        ax.scatter(d["x"], d["y"], s=22, c="lightgrey")
        e = d[d["IGF2BP2_detected"]]
        ax.scatter(e["x"], e["y"], s=30, c="tab:red")
        ax.set_title(f"array {arr}: IGF2BP2 detected in {len(e)} of {len(d)} spots")
        for a in row:
            a.set_aspect("equal")
            a.invert_yaxis()
            a.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            a.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.suptitle(f"Arrays shown: {chosen[0]} (most detected spots) and {chosen[1]} (fewest), chosen by rule, not by appearance",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(Path(out) / "fig3_spatial_maps.png", dpi=150)
    plt.close(fig)


# %% figure 4: tumor-level effects, locked and bias-corrected
def forest(ax, rows, title):
    for i, (name, per_tumor, mean, lo, hi) in enumerate(rows):
        ax.scatter(per_tumor, np.full(len(per_tumor), i), s=26, c="tab:blue", alpha=0.7, zorder=2)
        ax.errorbar([mean], [i], xerr=[[mean - lo], [hi - mean]], fmt="D", c="k", capsize=4, zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    ax.axvline(0, c="k", lw=0.5)
    ax.set_title(title)


def fig_effects(res, out):
    at = read(res, "04_spatial", "array_tests.tsv")
    ts = read(res, "04_spatial", "tumor_summary.tsv").set_index("outcome")
    bs = read(res, "05_posthoc", "random_set_bias_arrays.tsv")
    bsum = read(res, "05_posthoc", "random_set_bias_summary.tsv").set_index("model")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    rows = [(o, at[at["outcome"] == o].groupby("tumor")["beta"].mean().to_numpy(),
             ts.loc[o, "mean_beta"], ts.loc[o, "ci_low"], ts.loc[o, "ci_high"])
            for o in ["C1_stress", "H1_lipid", "H2_macrophage", "S1_lam", "S2_exhaustion"] if o in ts.index]
    forest(ax1, rows, "Locked analysis: adjusted effect per tumor (dots)\nmean and 95% interval over 5 tumors (diamond, bars)")
    ax1.set_xlabel("effect in SD of the outcome per unit of exposure")
    rows = [(m, bs[bs["model"] == m].groupby("tumor")["z_vs_random"].mean().to_numpy(),
             bsum.loc[m, "z_mean"], bsum.loc[m, "z_ci_low"], bsum.loc[m, "z_ci_high"])
            for m in ["C1_stress", "H1_lipid", "H2_macrophage"] if m in bsum.index]
    forest(ax2, rows, "Post-hoc: bias-corrected effect per tumor\n(real t versus matched random gene sets)")
    ax2.set_xlabel("z of real t relative to random-set t")
    fig.tight_layout()
    fig.savefig(Path(out) / "fig4_effects.png", dpi=150)
    plt.close(fig)


# %% figure 5: post-hoc diagnostics
def fig_diagnostics(res, out):
    sp = read(res, "05_posthoc", "spike_in_summary.tsv")
    sh = read(res, "05_posthoc", "split_half_summary.tsv")
    bs = read(res, "05_posthoc", "random_set_bias_arrays.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    ax.fill_between(sp["delta"], sp["ci_low"], sp["ci_high"], color="tab:blue", alpha=0.2)
    ax.plot(sp["delta"], sp["mean_beta"], c="tab:blue")
    met = sp["meets_criterion"].astype(bool)
    ax.scatter(sp.loc[met, "delta"], sp.loc[met, "mean_beta"], c="tab:green", s=40, zorder=3, label="criterion met")
    ax.scatter(sp.loc[~met, "delta"], sp.loc[~met, "mean_beta"], facecolors="none", edgecolors="tab:red", s=40, zorder=3, label="not met")
    ax.axhline(0, c="k", lw=0.5)
    ax.set_xlabel("shift added to the lipid score of IGF2BP2-detected spots (SD)")
    ax.set_ylabel("tumor-level effect (SD), 95% interval")
    ax.set_title("Spike-in: when would the locked criterion be met?")
    ax.legend()
    ax = axes[1]
    xs = np.arange(len(sh))
    ax.errorbar(xs, sh["mean"], yerr=[sh["mean"] - sh["ci_low"], sh["ci_high"] - sh["mean"]], fmt="o", c="k", capsize=4)
    ax.set_xticks(xs)
    ax.set_xticklabels(sh["signature"])
    ax.axhline(0, c="k", lw=0.5)
    ax.set_ylabel("effect of one random half on the other (SD)")
    ax.set_title("Split-half plumbing control (post-hoc)")
    ax = axes[2]
    piv = bs.pivot_table(index="model", columns="array", values="rand_mean_t")
    piv = piv[[c for c in ["A1", "A2", "A3", "B1", "B2", "B3", "D", "E", "G"] if c in piv.columns]]
    im = ax.imshow(piv.to_numpy(), cmap="coolwarm", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{piv.iloc[i, j]:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title("Mean t of matched random gene sets (0 = no bias)")
    fig.tight_layout()
    fig.savefig(Path(out) / "fig5_posthoc_diagnostics.png", dpi=150)
    plt.close(fig)


NOTES = """Figures written by notebooks/06_figures.py
fig1_single_cell.png: cluster sizes by provisional class for tumors A and B, and IGF2BP2 detection by cell group. Labels are provisional and have no copy-number support.
fig2_reproducibility_gate.png: t statistics of IGF2BP2-detected versus undetected candidate cancer cells in tumor A against tumor B, and the gate statistics against the shuffled-label null.
fig3_spatial_maps.png: spot maps of the arrays with the most and the fewest IGF2BP2-detected spots (chosen by rule, not appearance).
fig4_effects.png: tumor-level adjusted effects from the locked analysis, and the bias-corrected values (post-hoc).
fig5_posthoc_diagnostics.png: spike-in detectability, split-half plumbing control, and random-set bias per array. All post-hoc.
"""


# %% main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    with open(a.labels) as fh:
        labels = yaml.safe_load(fh)
    fig_cell_types(a.results, labels, out)
    fig_gate(a.results, out)
    fig_spatial_maps(a.results, out)
    fig_effects(a.results, out)
    fig_diagnostics(a.results, out)
    (out / "figure_notes.txt").write_text(NOTES)
    print("figures written to", out, "->", ", ".join(sorted(p.name for p in out.glob("*.png"))))


if __name__ == "__main__":
    main()