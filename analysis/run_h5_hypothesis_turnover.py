#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import run_biosnap_dti_gate2_models as dti
import run_hetionet_gate2_models as het
import run_lightgcn_gate2 as lgcn


def split_edges(edges, seed):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(edges))
    nt = max(1, int(round(0.2 * len(edges))))
    test = [edges[i] for i in idx[:nt]]
    train = [edges[i] for i in idx[nt:]]
    return train, test


def mask_candidates(edges, train, left, right):
    li = {x: i for i, x in enumerate(left)}
    ri = {x: i for i, x in enumerate(right)}
    dl = {x: 0 for x in left}
    dr = {x: 0 for x in right}
    for a, b in train:
        dl[a] += 1
        dr[b] += 1
    seen_l = np.asarray([dl[x] > 0 for x in left], dtype=bool)
    seen_r = np.asarray([dr[x] > 0 for x in right], dtype=bool)
    mask = seen_l[:, None] & seen_r[None, :]
    for a, b in edges:
        mask[li[a], ri[b]] = False
    return mask


def mf_matrix(model):
    li, ri, U, V, bu, bv = model
    return U @ V.T + bu[:, None] + bv[None, :]


def svd_matrix(model):
    li, ri, U, V = model
    return U @ V.T


def lgcn_matrix(model, left, right):
    li, ri, z = model
    L = np.vstack([z[li[x]] for x in left])
    R = np.vstack([z[ri[x]] for x in right])
    return L @ R.T


def top_flat(score, mask, k):
    flat_score = score.ravel()
    flat_mask = mask.ravel()
    ids = np.flatnonzero(flat_mask)
    if len(ids) < k:
        k = len(ids)
    vals = flat_score[ids]
    if k == 0:
        return np.array([], dtype=np.int64)
    take = np.argpartition(vals, -k)[-k:]
    selected = ids[take]
    order = np.argsort(flat_score[selected])[::-1]
    return selected[order]


def rank_stats(score_a, score_b, mask, ks=(100, 500, 1000)):
    x = score_a[mask]
    y = score_b[mask]
    rho = float(spearmanr(x, y).statistic)
    out = {"candidate_count": int(mask.sum()), "spearman_all_candidates": rho}
    for k in ks:
        A = top_flat(score_a, mask, k)
        B = top_flat(score_b, mask, k)
        sa, sb = set(A.tolist()), set(B.tolist())
        inter = len(sa & sb)
        kk = min(k, len(sa), len(sb))
        union = len(sa | sb)
        out[f"k{k}_overlap_count"] = inter
        out[f"k{k}_turnover"] = float(1 - inter / kk) if kk else float("nan")
        out[f"k{k}_jaccard"] = float(inter / union) if union else float("nan")
    return out


def bootstrap_seed_ci(df, col, B=10000, seed=20260913):
    vals = np.asarray(df[col], dtype=float)
    rng = np.random.default_rng(seed)
    means = np.empty(B, dtype=float)
    for b in range(B):
        means[b] = rng.choice(vals, size=len(vals), replace=True).mean()
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def run_dti(path, seeds):
    edges = dti.parse(path)
    left = sorted({a for a, _ in edges})
    right = sorted({b for _, b in edges})
    blocked = set(edges)
    rows = []
    for seed in seeds:
        train, _ = split_edges(edges, seed)
        mask = mask_candidates(edges, train, left, right)
        # Gate-2 mean winner under conventional evaluation: NeuralMF.
        # Gate-2 mean winner after degree matching: SVD.
        mf = dti.fit_logistic_mf(train, left, right, blocked, np.random.default_rng(seed + 1000))
        svd = dti.fit_svd(train, left, right, 32, seed)
        s_random_winner = mf_matrix(mf)
        s_neutral_winner = svd_matrix(svd)
        r = rank_stats(s_random_winner, s_neutral_winner, mask)
        r.update({"family": "DTI", "seed": seed, "random_winner": "NeuralMF", "neutral_winner": "SVD"})
        rows.append(r)
    return rows


def run_dag(path, seeds):
    edges = het.parse_edges(path, "DaG", "Disease::", "Gene::")
    left = sorted({a for a, _ in edges})
    right = sorted({b for _, b in edges})
    blocked = set(edges)
    rows = []
    for seed in seeds:
        train, _ = split_edges(edges, seed)
        mask = mask_candidates(edges, train, left, right)
        # Gate-2 mean winner under conventional evaluation: NeuralMF.
        # Gate-2 mean winner after degree matching: LightGCN.
        mf = het.fit_logistic_mf(train, left, right, blocked, np.random.default_rng(seed + 10000))
        graph = lgcn.fit_lightgcn(train, left, right, blocked, seed)
        s_random_winner = mf_matrix(mf)
        s_neutral_winner = lgcn_matrix(graph, left, right)
        r = rank_stats(s_random_winner, s_neutral_winner, mask)
        r.update({"family": "Disease-Gene", "seed": seed, "random_winner": "NeuralMF", "neutral_winner": "LightGCN"})
        rows.append(r)
    return rows


def summarize(df):
    summary = {}
    for family, q in df.groupby("family"):
        z = {
            "replicates": int(len(q)),
            "random_winner": q.random_winner.iloc[0],
            "neutral_winner": q.neutral_winner.iloc[0],
            "candidate_count_mean": float(q.candidate_count.mean()),
            "spearman_mean": float(q.spearman_all_candidates.mean()),
            "spearman_sd": float(q.spearman_all_candidates.std(ddof=1)),
        }
        for k in (100, 500, 1000):
            col = f"k{k}_turnover"
            lo, hi = bootstrap_seed_ci(q, col)
            z[f"ht{k}_mean"] = float(q[col].mean())
            z[f"ht{k}_sd"] = float(q[col].std(ddof=1))
            z[f"ht{k}_bootstrap95"] = [lo, hi]
            z[f"jaccard{k}_mean"] = float(q[f"k{k}_jaccard"].mean())
        summary[family] = z
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dti", required=True)
    ap.add_argument("--hetionet", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--replicates", type=int, default=5)
    args = ap.parse_args()

    seeds = list(range(args.replicates))
    rows = run_dti(args.dti, seeds) + run_dag(args.hetionet, seeds)
    df = pd.DataFrame(rows)
    p = Path(args.out_prefix)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(p) + "_replicates.csv", index=False)
    s = summarize(df)
    Path(str(p) + "_summary.json").write_text(json.dumps(s, indent=2) + "\n")

    lines = [
        "# H5 pilot — benchmark-induced hypothesis-selection instability",
        "",
        "The candidate universe is identical within each seed. Known positives and endpoints absent from training are excluded before either winner ranks candidates.",
        "",
        "| Family | Conventional winner | Neutralized winner | Candidate pairs | Spearman | HT@100 | HT@500 | HT@1000 |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for fam in ["DTI", "Disease-Gene"]:
        z = s[fam]
        lines.append(
            f"| {fam} | {z['random_winner']} | {z['neutral_winner']} | {z['candidate_count_mean']:.0f} | "
            f"{z['spearman_mean']:.3f} | {z['ht100_mean']:.3f} | {z['ht500_mean']:.3f} | {z['ht1000_mean']:.3f} |"
        )
    lines += [
        "",
        "HT@k = 1 - overlap/k. Values are means across frozen split seeds. The JSON records percentile bootstrap 95% intervals over split-level estimates.",
        "",
        "Interpretation boundary: turnover is scientific-decision instability, not evidence that either candidate list is biologically correct. Temporal/independent validation is required to assign external validity.",
    ]
    md = "\n".join(lines) + "\n"
    Path(str(p) + ".md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
