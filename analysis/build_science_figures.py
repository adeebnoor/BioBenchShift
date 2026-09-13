#!/usr/bin/env python3
"""Build the four main Science-style figures from frozen source-data rows.

The script intentionally reads figures/source_data_main.csv instead of embedding
headline values in plotting code. Generated files are presentation artifacts;
the frozen result files named in source_data_main.csv remain the scientific
sources of record.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "figures" / "source_data_main.csv"
OUT = ROOT / "figures" / "generated"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SRC)


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig1():
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))
    ax = axes[0]
    d = df[(df.figure == 1) & (df.analysis == "degree_only_auc")]
    order = ["DTI", "HuRI_PPI", "Compound_Disease", "Disease_Gene"]
    x = np.arange(len(order)); w = 0.36
    conv = [float(d[(d.series == "Conventional") & (d.x == k)].value.iloc[0]) for k in order]
    neut = [float(d[(d.series == "Neutralized") & (d.x == k)].value.iloc[0]) for k in order]
    ax.bar(x-w/2, conv, w, label="Conventional")
    ax.bar(x+w/2, neut, w, label="Structure-neutralized")
    ax.axhline(0.5, ls="--", lw=1)
    ax.set_ylim(0.45, 1.02); ax.set_ylabel("Degree-only AUROC")
    ax.set_xticks(x, ["DTI", "HuRI PPI", "Compound–disease", "Disease–gene"], rotation=15, ha="right")
    ax.legend(frameon=False); ax.set_title("a  Structural observability can look predictive", loc="left", fontweight="bold")

    ax = axes[1]
    d = df[(df.figure == 1) & (df.analysis == "learned_dti_auc")]
    models = ["SVD", "NeuralMF", "LightGCN"]; x = np.arange(len(models))
    conv = [float(d[(d.series == "Conventional") & (d.x == k)].value.iloc[0]) for k in models]
    neut = [float(d[(d.series == "Neutralized") & (d.x == k)].value.iloc[0]) for k in models]
    ax.bar(x-w/2, conv, w, label="Conventional")
    ax.bar(x+w/2, neut, w, label="Structure-neutralized")
    ax.set_ylim(0.80, 1.01); ax.set_ylabel("DTI AUROC")
    ax.set_xticks(x, models); ax.set_title("b  Learned models depend on it unequally", loc="left", fontweight="bold")
    for s in ["top", "right"]:
        for a in axes: a.spines[s].set_visible(False)
    fig.suptitle("Figure 1 | A benchmark can reward structural observability", y=1.03, fontweight="bold")
    fig.tight_layout(); save(fig, "Figure1_structural_observability")


def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    for ax, analysis, title in [
        (axes[0], "dti_turnover", "a  Drug–target hypotheses"),
        (axes[1], "ppi_turnover", "b  Protein-interaction hypotheses")]:
        d = df[(df.figure == 2) & (df.analysis == analysis)]
        ks = sorted(d.x.astype(int).unique())
        for series, lab in [("Cross_selected", "Different selected models"),
                            ("Within_NeuralMF_LOO", "Within NeuralMF"),
                            ("Within_SVD_LOO", "Within SVD")]:
            dd = d[d.series == series].sort_values("x")
            ax.plot(dd.x.astype(int), dd.value, marker="o", label=lab)
        ax.set_xscale("log"); ax.set_ylim(0, 1.04); ax.set_ylabel("Hypothesis turnover (HT@k)"); ax.set_xlabel("Top-k")
        ax.set_title(title, loc="left", fontweight="bold"); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Figure 2 | Benchmark choice changes the model and redirects hypothesis identity", y=1.03, fontweight="bold")
    fig.tight_layout(); save(fig, "Figure2_hypothesis_identity")


def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    ax = axes[0]
    d = df[(df.figure == 3) & (df.analysis == "biogrid_future_hits")]
    for series in ["SVD", "NeuralMF", "LightGCN"]:
        dd = d[d.series == series].sort_values("x")
        ax.plot(dd.x.astype(int), dd.value, marker="o", label=series)
    ax.set_xscale("log"); ax.set_xlabel("Top-k historical candidates"); ax.set_ylabel("Later-supported BioGRID relations")
    ax.set_title("a  Same 70,041,100-pair historical universe", loc="left", fontweight="bold"); ax.legend(frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax = axes[1]
    d = df[(df.figure == 3) & (df.analysis == "biogrid_recall_delta")].sort_values("x")
    x = np.arange(len(d)); y = d.value.to_numpy(float); lo = d.lower.to_numpy(float); hi = d.upper.to_numpy(float)
    ax.errorbar(x, y, yerr=[y-lo, hi-y], fmt="o", capsize=3)
    ax.axhline(0, ls="--", lw=1); ax.set_xticks(x, [f"{int(v):,}" for v in d.x]); ax.set_xlabel("Top-k")
    ax.set_ylabel("Recall difference: SVD − NeuralMF")
    ax.set_title("b  Paired bootstrap uncertainty", loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.suptitle("Figure 3 | Historical benchmark choice changes later evidence recovery", y=1.03, fontweight="bold")
    fig.tight_layout(); save(fig, "Figure3_biogrid_future_evidence")


def fig4():
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    ax = axes[0]
    d = df[(df.figure == 4) & (df.analysis == "chembl_primary_hits")]
    ks = [100, 500, 1000]; x = np.arange(len(ks)); w = 0.25
    for i, series in enumerate(["SVD", "NeuralMF", "LightGCN"]):
        vals = [float(d[(d.series == series) & (d.x.astype(int) == k)].value.iloc[0]) for k in ks]
        ax.bar(x + (i-1)*w, vals, w, label=series)
    ax.set_xticks(x, [str(k) for k in ks]); ax.set_xlabel("Top-k"); ax.set_ylabel("ChEMBL-supported relations")
    ax.set_title("a  Independent evidence at the experimental frontier", loc="left", fontweight="bold"); ax.legend(frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax = axes[1]
    d = df[(df.figure == 4) & (df.analysis == "chembl_broad_hits")]
    for series in ["SVD", "NeuralMF"]:
        dd = d[d.series == series].sort_values("x")
        ax.plot(dd.x.astype(int), dd.value, marker="o", label=series)
    ax.set_xscale("log"); ax.set_xlabel("Top-k"); ax.set_ylabel("ChEMBL-supported relations")
    ax.set_title("b  Broad-cutoff crossover is retained", loc="left", fontweight="bold"); ax.legend(frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.suptitle("Figure 4 | Independent evidence differs at the drug–target experimental frontier", y=1.03, fontweight="bold")
    fig.tight_layout(); save(fig, "Figure4_chembl_frontier")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print(f"Wrote figures to {OUT}")
