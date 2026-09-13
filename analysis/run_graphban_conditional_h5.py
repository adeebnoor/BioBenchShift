#!/usr/bin/env python3
"""Conditional GraphBAN-style full-graph H5 analysis.

Protocol: GRAPHBAN_CONDITIONAL_H5_CHEMBL_AMENDMENT_20260913.md
The script is written before the full-graph GraphBAN ranking outcome is inspected.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler

import run_chembl_dti_external as h6
from run_graphban_targetdecagon_clean import (
    make_graph,
    model_scores,
    parse_edges,
    train_model,
)

KS = (100, 500, 1000)
GRAPHBAN_SEEDS = (0, 1, 2)
COMPARATORS = ("NeuralMF", "SVD", "LightGCN")


def top_ids(score: np.ndarray, mask: np.ndarray, k: int) -> np.ndarray:
    ids = np.flatnonzero(mask.ravel())
    kk = min(int(k), len(ids))
    if kk <= 0:
        return np.asarray([], dtype=np.int64)
    vals = score.ravel()[ids]
    take = np.argpartition(vals, -kk)[-kk:]
    sel = ids[take]
    return sel[np.argsort(score.ravel()[sel])[::-1]]


def overlap_stats(a: np.ndarray, b: np.ndarray, k: int) -> dict:
    A, B = set(map(int, a)), set(map(int, b))
    inter = len(A & B)
    denom = min(int(k), len(A), len(B))
    union = len(A | B)
    return {
        "shared": inter,
        "turnover": float(1 - inter / denom) if denom else float("nan"),
        "jaccard": float(inter / union) if union else float("nan"),
    }


def fit_graphban_full(mapped, drugs, genes, Xd0, Xp0, seed, device):
    # Complete mapped graph is the downstream ranking training graph; this is
    # not a held-out performance experiment.
    Xd = StandardScaler().fit_transform(Xd0).astype("float32")
    Xp = StandardScaler().fit_transform(Xp0).astype("float32")
    data, li, ri = make_graph(drugs, genes, mapped, Xd, Xp, device)
    model = train_model(
        data,
        li,
        ri,
        mapped,
        set(mapped),
        drugs,
        genes,
        int(seed),
        device,
        epochs=20,
    )
    mask = np.ones((len(drugs), len(genes)), dtype=bool)
    for a, b in mapped:
        mask[li[a], ri[b]] = False
    flat = np.full(mask.size, -np.inf, dtype=np.float32)
    ids = np.flatnonzero(mask.ravel())
    nr = len(genes)
    candidates = [(drugs[int(x // nr)], genes[int(x % nr)]) for x in ids]
    flat[ids] = model_scores(model, data, candidates, li, ri, device).astype("float32")
    return flat.reshape(mask.shape), mask


def instability_rows(seed_scores, mask):
    rows = []
    seed_tops = {s: {k: top_ids(sc, mask, k) for k in KS} for s, sc in seed_scores.items()}
    for a, b in combinations(sorted(seed_scores), 2):
        for k in KS:
            st = overlap_stats(seed_tops[a][k], seed_tops[b][k], k)
            rows.append({"kind": "seed_pair", "a": str(a), "b": str(b), "k": k, **st})

    seeds = sorted(seed_scores)
    loo_scores = {}
    for left_out in seeds:
        keep = [seed_scores[s] for s in seeds if s != left_out]
        loo_scores[left_out] = np.mean(np.stack(keep), axis=0)
    loo_tops = {s: {k: top_ids(sc, mask, k) for k in KS} for s, sc in loo_scores.items()}
    for a, b in combinations(seeds, 2):
        for k in KS:
            st = overlap_stats(loo_tops[a][k], loo_tops[b][k], k)
            rows.append({"kind": "loo_pair", "a": f"drop_{a}", "b": f"drop_{b}", "k": k, **st})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--drug-map", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--drug-features", required=True)
    ap.add_argument("--protein-features", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--comparator-inits", type=int, default=5)
    args = ap.parse_args()

    edges = parse_edges(args.input)
    dd = pd.read_csv(args.drug_map, dtype=str).fillna("")
    gd = pd.read_csv(args.gene_map, dtype=str).fillna("")
    dd = dd[dd.smiles.astype(bool)].drop_duplicates("drug")
    gd = gd[gd.sequence.astype(bool)].drop_duplicates("gene")
    dm = set(dd.drug.astype(str))
    gm = set(gd.gene.astype(str))
    mapped = sorted(e for e in edges if e[0] in dm and e[1] in gm)
    drugs = sorted({a for a, _ in mapped})
    genes = sorted({b for _, b in mapped})

    # Frozen mapping preflight values.
    if len(drugs) != 284 or len(genes) != 3607 or len(mapped) != 18631:
        raise RuntimeError(
            f"GraphBAN mapping fingerprint drift: drugs={len(drugs)}, genes={len(genes)}, edges={len(mapped)}"
        )

    Xd = np.load(args.drug_features).astype("float32")
    Xp = np.load(args.protein_features).astype("float32")
    if Xd.shape[0] != len(drugs) or Xp.shape[0] != len(genes):
        raise RuntimeError(
            f"feature row mismatch: drug={Xd.shape[0]}/{len(drugs)} protein={Xp.shape[0]}/{len(genes)}"
        )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device", device)
    seed_scores = {}
    mask = None
    for seed in GRAPHBAN_SEEDS:
        print("full-graph GraphBAN seed", seed, flush=True)
        sc, m = fit_graphban_full(mapped, drugs, genes, Xd, Xp, seed, device)
        if mask is None:
            mask = m
        elif not np.array_equal(mask, m):
            raise RuntimeError("candidate mask changed across GraphBAN seeds")
        seed_scores[seed] = sc
    assert mask is not None
    ensemble = np.mean(np.stack([seed_scores[s] for s in GRAPHBAN_SEEDS]), axis=0).astype("float32")

    # Comparator full-graph ensembles use the existing frozen H6 definitions,
    # then are restricted to the identical GraphBAN-mappable candidate universe.
    all_drugs, all_genes, comparator_full = h6.fit_ensembles(edges, args.comparator_inits)
    adi = {x: i for i, x in enumerate(all_drugs)}
    agi = {x: i for i, x in enumerate(all_genes)}
    drows = [adi[x] for x in drugs]
    grows = [agi[x] for x in genes]
    comparator = {name: sc[np.ix_(drows, grows)] for name, sc in comparator_full.items()}

    comparison_rows = []
    gb_top = {k: top_ids(ensemble, mask, k) for k in KS}
    for name in COMPARATORS:
        sc = comparator[name]
        rho = float(spearmanr(ensemble[mask], sc[mask]).statistic)
        for k in KS:
            ct = top_ids(sc, mask, k)
            st = overlap_stats(gb_top[k], ct, k)
            comparison_rows.append({
                "model_a": "GraphBAN-style",
                "model_b": name,
                "k": k,
                "candidate_count": int(mask.sum()),
                "spearman_all_candidates": rho,
                **st,
            })
    comp_df = pd.DataFrame(comparison_rows)
    inst_df = instability_rows(seed_scores, mask)

    # Transparent top-1000 identities for all four models.
    nr = len(genes)
    top_rows = []
    score_sets = {"GraphBAN-style": ensemble, **comparator}
    for name, sc in score_sets.items():
        ids = top_ids(sc, mask, 1000)
        for rank, x in enumerate(ids, 1):
            i, j = int(x // nr), int(x % nr)
            top_rows.append({
                "model": name,
                "rank": rank,
                "drug": drugs[i],
                "gene": genes[j],
                "score": float(sc[i, j]),
            })
    top_df = pd.DataFrame(top_rows)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    comp_df.to_csv(str(prefix) + "_comparisons.csv", index=False)
    inst_df.to_csv(str(prefix) + "_instability.csv", index=False)
    top_df.to_csv(str(prefix) + "_top1000_pairs.csv", index=False)
    np.savez_compressed(
        str(prefix) + "_scores.npz",
        drugs=np.asarray(drugs, dtype="U32"),
        genes=np.asarray(genes, dtype="U32"),
        mask=mask,
        graphban_ensemble=ensemble,
        graphban_seed0=seed_scores[0],
        graphban_seed1=seed_scores[1],
        graphban_seed2=seed_scores[2],
        neuralmf=comparator["NeuralMF"].astype("float32"),
        svd=comparator["SVD"].astype("float32"),
        lightgcn=comparator["LightGCN"].astype("float32"),
    )

    inst_summary = {}
    for kind, q in inst_df.groupby("kind"):
        inst_summary[kind] = {}
        for k, z in q.groupby("k"):
            inst_summary[kind][str(int(k))] = {
                "turnover_mean": float(z.turnover.mean()),
                "turnover_max": float(z.turnover.max()),
                "jaccard_mean": float(z.jaccard.mean()),
            }
    summary = {
        "protocol": "GRAPHBAN_CONDITIONAL_H5_CHEMBL_AMENDMENT_20260913.md",
        "baseline_commit": "1298e4c83ff4fa488a784d52219d20e706118a9e",
        "graphban_seeds": list(GRAPHBAN_SEEDS),
        "comparator_initializations": args.comparator_inits,
        "mapped_drugs": len(drugs),
        "mapped_genes": len(genes),
        "mapped_positive_edges": len(mapped),
        "common_unknown_candidate_count": int(mask.sum()),
        "within_graphban_instability": inst_summary,
        "interpretation_boundary": "Queue turnover measures scientific-decision identity, not biological correctness. GraphBAN-style is not a verbatim reproduction of full published GraphBAN.",
    }
    Path(str(prefix) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# Conditional GraphBAN-style H5 — common-universe hypothesis identity",
        "",
        f"Common unknown candidate universe: **{int(mask.sum()):,}** pairs across {len(drugs)} drugs and {len(genes)} sequence-mapped genes.",
        "",
        "| Comparator | Spearman (all candidates) | HT@100 | HT@500 | HT@1000 |",
        "|---|---:|---:|---:|---:|",
    ]
    for name in COMPARATORS:
        z = comp_df[comp_df.model_b == name].set_index("k")
        lines.append(
            f"| {name} | {z.spearman_all_candidates.iloc[0]:.3f} | "
            f"{z.loc[100,'turnover']:.3f} | {z.loc[500,'turnover']:.3f} | {z.loc[1000,'turnover']:.3f} |"
        )
    lines += ["", "Within-GraphBAN instability controls:"]
    for kind in ("seed_pair", "loo_pair"):
        q = inst_df[inst_df.kind == kind]
        for k in KS:
            z = q[q.k == k]
            lines.append(f"- {kind} HT@{k}: mean {z.turnover.mean():.3f}, max {z.turnover.max():.3f}.")
    lines += [
        "",
        "All comparator queues use the identical GraphBAN-mappable candidate universe. No external evidence source enters ranking. Cross-model turnover is not evidence that one list is biologically correct.",
    ]
    md = "\n".join(lines) + "\n"
    Path(str(prefix) + ".md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
