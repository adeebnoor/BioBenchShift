#!/usr/bin/env python3
"""Complete the Science A1/A2 replay while retaining execution drift.

The frozen GraphBAN metrics are used only as a reproducibility comparator.
If the new execution differs, it is reported as a distinct replay rather than
tuned, pooled, or attached to the frozen AUROCs. Queue turnover is always
paired with AUROC recomputed from the SAME newly fitted model.
"""
from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

import run_science_queue_resolution_a1 as core

PROTOCOL = "BioBenchShift-Science-A1-A2-secondary-20260915-v2-replay"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--drug-map", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--drug-features", required=True)
    ap.add_argument("--protein-features", required=True)
    ap.add_argument("--expected-panel", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seeds", default="0,1,2,3,4")
    a = ap.parse_args()

    E = core.gb.parse_edges(a.input)
    dd = pd.read_csv(a.drug_map, dtype=str).fillna("")
    gd = pd.read_csv(a.gene_map, dtype=str).fillna("")
    dd = dd[dd.smiles.astype(bool)].drop_duplicates("drug")
    gd = gd[gd.sequence.astype(bool)].drop_duplicates("gene")
    dm, gm = set(dd.drug), set(gd.gene)
    mapped = [e for e in E if e[0] in dm and e[1] in gm]
    L = sorted({x for x, _ in mapped})
    R = sorted({y for _, y in mapped})
    baseXd = np.load(a.drug_features).astype("float32")
    baseXp = np.load(a.protein_features).astype("float32")
    if baseXd.shape[0] != len(L) or baseXp.shape[0] != len(R):
        raise RuntimeError("feature/order mismatch")

    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]
    expected = pd.read_csv(a.expected_panel)
    torch.set_num_threads(4)
    device = torch.device("cpu")

    prefix = Path(a.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    eval_rows, pair_rows, split_rows, top_rows = [], [], [], []
    candidate_hash = None

    for seed in seeds:
        print(f"A1 replay seed {seed}...", flush=True)
        cand_ids, scores, erows, smeta = core.fit_seed(E, mapped, L, R, baseXd, baseXp, seed, device)
        eval_rows.extend(erows)
        split_rows.append(smeta)
        if candidate_hash is None:
            candidate_hash = smeta["candidate_universe_hash"]
        elif candidate_hash != smeta["candidate_universe_hash"]:
            raise RuntimeError("candidate universe differs across split seeds")

        # Audit top-1000 identities.
        for model in core.MODELS:
            ids1000 = core.top_ids(scores[model], cand_ids, 1000)
            pos = {int(cid): i + 1 for i, cid in enumerate(ids1000.tolist())}
            score_at = {int(cid): float(sc) for cid, sc in zip(cand_ids.tolist(), scores[model].tolist())}
            for cid in ids1000:
                cid = int(cid)
                top_rows.append({
                    "seed": seed,
                    "model": model,
                    "candidate_id": cid,
                    "left_index": cid // len(R),
                    "right_index": cid % len(R),
                    "rank": pos[cid],
                    "score": score_at[cid],
                })

        emap = {(r["model"], r["regime"]): r for r in erows}
        for i, ma in enumerate(core.MODELS):
            for mb in core.MODELS[i + 1:]:
                hts = {k: core.ht(scores[ma], scores[mb], cand_ids, k) for k in core.CUTOFFS}
                for regime in ("random", "matched"):
                    aa = float(emap[(ma, regime)]["auroc"])
                    ab = float(emap[(mb, regime)]["auroc"])
                    pair_rows.append({
                        "seed": seed,
                        "regime": regime,
                        "model_a": ma,
                        "model_b": mb,
                        "auroc_a": aa,
                        "auroc_b": ab,
                        "delta_auroc_signed_a_minus_b": aa - ab,
                        "delta_auroc_abs": abs(aa - ab),
                        "ht_100": hts[100],
                        "ht_500": hts[500],
                        "ht_1000": hts[1000],
                    })
        del scores

    eval_df = pd.DataFrame(eval_rows)
    pair_df = pd.DataFrame(pair_rows)
    split_df = pd.DataFrame(split_rows)
    top_df = pd.DataFrame(top_rows)

    # Reproducibility comparison only; never used to substitute frozen values.
    fidelity = []
    for _, r in expected.iterrows():
        seed = int(r.seed)
        for model in core.MODELS:
            for regime in ("random", "matched"):
                got = float(eval_df.query("seed == @seed and model == @model and regime == @regime").iloc[0].auroc)
                exp = float(r[f"{model}_auc_{regime}"])
                fidelity.append({
                    "seed": seed,
                    "model": model,
                    "regime": regime,
                    "frozen_auroc": exp,
                    "replay_auroc": got,
                    "replay_minus_frozen": got - exp,
                    "abs_diff": abs(got - exp),
                })
    fidelity_df = pd.DataFrame(fidelity)

    # Baseline models should replay closely; GraphBAN drift is retained explicitly.
    by_model = fidelity_df.groupby("model").abs_diff.agg(["mean", "max"]).reset_index()
    replay_status = "exact_within_5e-7" if fidelity_df.abs_diff.max() <= 5e-7 else "distinct_replay_execution"

    eval_df.to_csv(str(prefix) + "_evaluation_metrics.csv", index=False)
    pair_df.to_csv(str(prefix) + "_pairwise_resolution.csv", index=False)
    split_df.to_csv(str(prefix) + "_split_provenance.csv", index=False)
    fidelity_df.to_csv(str(prefix) + "_frozen_vs_replay.csv", index=False)
    top_df.to_csv(str(prefix) + "_top1000_identities.csv.gz", index=False, compression="gzip")

    nom = pair_df[(pair_df.model_a == "neuralmf") & (pair_df.model_b == "graphban")].copy()
    if len(nom) != 10:
        raise RuntimeError(f"nominated pair rows {len(nom)} != 10")

    per_seed = []
    for seed in seeds:
        q = nom[(nom.seed == seed) & (nom.regime == "random")].iloc[0]
        m = nom[(nom.seed == seed) & (nom.regime == "matched")].iloc[0]
        per_seed.append({
            "seed": seed,
            "ht100": float(q.ht_100),
            "ht500": float(q.ht_500),
            "ht1000": float(q.ht_1000),
            "abs_delta_auroc_random": float(q.delta_auroc_abs),
            "abs_delta_auroc_matched": float(m.delta_auroc_abs),
        })

    summary = {
        "protocol": PROTOCOL,
        "analysis_status": "completed_secondary_replay",
        "frozen_comparison_status": replay_status,
        "candidate_universe_size": int(split_df.n_candidates.iloc[0]),
        "candidate_universe_hash": candidate_hash,
        "nominated_pair": ["graphban", "neuralmf"],
        "nominated_pair_was_chosen_after_existing_auroc_table": True,
        "graphban_neuralmf": {
            "per_seed": per_seed,
            "mean_ht100": float(nom[nom.regime == "random"].ht_100.mean()),
            "mean_ht500": float(nom[nom.regime == "random"].ht_500.mean()),
            "mean_ht1000": float(nom[nom.regime == "random"].ht_1000.mean()),
            "random_mean_abs_delta_auroc": float(nom[nom.regime == "random"].delta_auroc_abs.mean()),
            "matched_mean_abs_delta_auroc": float(nom[nom.regime == "matched"].delta_auroc_abs.mean()),
            "random_delta_range": [float(nom[nom.regime == "random"].delta_auroc_abs.min()), float(nom[nom.regime == "random"].delta_auroc_abs.max())],
            "matched_delta_range": [float(nom[nom.regime == "matched"].delta_auroc_abs.min()), float(nom[nom.regime == "matched"].delta_auroc_abs.max())],
        },
        "fidelity_by_model": by_model.to_dict("records"),
        "interpretation_boundary": (
            "Queue turnover and AUROC gaps are paired within this replay execution. "
            "The frozen GraphBAN AUROCs are not substituted into the replay. A small AUROC gap is descriptive and not an equivalence test."
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "torch": torch.__version__,
            "device": str(device),
            "torch_threads": torch.get_num_threads(),
        },
        "input_hashes": {
            "targetdecagon": core.sha256_file(a.input),
            "drug_map": core.sha256_file(a.drug_map),
            "gene_map": core.sha256_file(a.gene_map),
            "drug_features": core.sha256_file(a.drug_features),
            "protein_features": core.sha256_file(a.protein_features),
            "frozen_panel": core.sha256_file(a.expected_panel),
        },
    }
    Path(str(prefix) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # A2 descriptive figure. Same queue HT is plotted under both evaluation regimes;
    # each x-coordinate is the AUROC gap recomputed for that regime.
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8.2, 5.6))
    ax = fig.add_subplot(111)
    for _, r in pair_df.iterrows():
        marker = "o" if r.regime == "random" else "x"
        ax.scatter(r.delta_auroc_abs, r.ht_100, marker=marker, alpha=0.72)
    ax.set_xlabel("Absolute paired AUROC difference")
    ax.set_ylabel("Hypothesis turnover at top 100")
    ax.set_title("Same-fit model-score separation and queue turnover")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(str(prefix) + "_A2_resolution.png", dpi=220)
    plt.close(fig)

    lines = [
        "# Science A1/A2 same-fit queue-resolution replay",
        "",
        f"Protocol: `{PROTOCOL}`",
        f"Candidate universe: **{summary['candidate_universe_size']:,}** mapped unknown drug-target pairs.",
        f"Frozen-vs-replay status: **{replay_status}**. Frozen metrics are retained separately and are not substituted into this analysis.",
        "",
        "## GraphBAN-NeuralMF nominated comparison",
        "",
        "| Seed | HT@100 | HT@500 | HT@1000 | |ΔAUROC| random | |ΔAUROC| matched |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for r in per_seed:
        lines.append(
            f"| {r['seed']} | {r['ht100']:.3f} | {r['ht500']:.3f} | {r['ht1000']:.3f} | {r['abs_delta_auroc_random']:.6f} | {r['abs_delta_auroc_matched']:.6f} |"
        )
    lines += [
        "",
        f"Mean HT@100 = **{summary['graphban_neuralmf']['mean_ht100']:.3f}**; mean HT@500 = **{summary['graphban_neuralmf']['mean_ht500']:.3f}**; mean HT@1000 = **{summary['graphban_neuralmf']['mean_ht1000']:.3f}**.",
        f"Mean absolute AUROC gap = **{summary['graphban_neuralmf']['random_mean_abs_delta_auroc']:.6f}** under random controls and **{summary['graphban_neuralmf']['matched_mean_abs_delta_auroc']:.6f}** under matched controls.",
        "",
        "All six model pairs are retained in the pairwise CSV. These are descriptive same-fit results. They do not establish statistical equivalence or a universal benchmark-resolution threshold.",
    ]
    Path(str(prefix) + ".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
