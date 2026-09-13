#!/usr/bin/env python3
"""Expansion D: prespecified matching-method robustness on BioSNAP DTI.

Rules are frozen in MATCHING_ROBUSTNESS_AMENDMENT_20260913.md. Model fits are
shared across matching methods within each seed; only the evaluation controls
change. The original log2-bin rule remains the primary method.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.metrics import roc_auc_score

import run_biosnap_dti_gate2_models as dti
import run_lightgcn_gate2 as lgcn

METHODS = ("log2_bins", "degree_deciles", "nearest_caliper")


def auc(pos, neg, scorer):
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    s = np.r_[scorer(pos), scorer(neg)]
    return float(roc_auc_score(y, s))


def all_nonedges(L, R, blocked):
    return [(a, b) for a in L for b in R if (a, b) not in blocked]


def pool_match(rng, positives, candidates, signature):
    pools = defaultdict(list)
    for e in candidates:
        pools[signature(e)].append(e)
    for v in pools.values():
        rng.shuffle(v)
    used = Counter()
    mp, mn = [], []
    for i in rng.permutation(len(positives)):
        p = positives[int(i)]
        sig = signature(p)
        j = used[sig]
        pool = pools.get(sig, [])
        if j < len(pool):
            mp.append(p)
            mn.append(pool[j])
            used[sig] += 1
    return mp, mn


def quantile_mapper(values, q=10):
    z = np.asarray([math.log1p(float(x)) for x in values], float)
    cuts = np.unique(np.quantile(z, np.linspace(0, 1, q + 1)))
    inner = cuts[1:-1]
    def f(x):
        return int(np.searchsorted(inner, math.log1p(float(x)), side="right"))
    return f, [float(x) for x in cuts]


def nearest_match(positives, candidates, dl, dr):
    left_vals = np.asarray([math.log1p(float(v)) for v in dl.values()] or [0.0])
    right_vals = np.asarray([math.log1p(float(v)) for v in dr.values()] or [0.0])
    mu = np.array([left_vals.mean(), right_vals.mean()], float)
    sd = np.array([left_vals.std(ddof=0), right_vals.std(ddof=0)], float)
    sd[sd == 0] = 1.0

    def feat(e):
        return (np.array([math.log1p(dl.get(e[0], 0)), math.log1p(dr.get(e[1], 0))], float) - mu) / sd

    X = np.vstack([feat(e) for e in candidates])
    tree = cKDTree(X)
    used = set()
    mp, mn = [], []
    radius = math.sqrt(2.0) * 0.25 + 1e-12
    for p in positives:
        fp = feat(p)
        ids = tree.query_ball_point(fp, r=radius)
        eligible = []
        for j in ids:
            if j in used:
                continue
            delta = np.abs(X[j] - fp)
            if np.all(delta <= 0.25 + 1e-12):
                eligible.append((float(np.dot(delta, delta)), j))
        if not eligible:
            continue
        eligible.sort(key=lambda z: (z[0], z[1]))
        j = eligible[0][1]
        used.add(j)
        mp.append(p)
        mn.append(candidates[j])
    return mp, mn


def one_seed(E, seed):
    rng = np.random.default_rng(seed)
    ix = rng.permutation(len(E))
    nt = round(0.2 * len(E))
    test = [E[i] for i in ix[:nt]]
    train = [E[i] for i in ix[nt:]]
    blocked = set(E)
    L = sorted({a for a, _ in E})
    R = sorted({b for _, b in E})
    dl = Counter(a for a, _ in train)
    dr = Counter(b for _, b in train)
    seen = [e for e in test if dl[e[0]] > 0 and dr[e[1]] > 0]
    random_neg = dti.sample_random(rng, len(seen), L, R, blocked)

    svd = dti.fit_svd(train, L, R, 32, seed)
    mf = dti.fit_logistic_mf(train, L, R, blocked, np.random.default_rng(seed + 1000))
    gcn = lgcn.fit_lightgcn(train, L, R, blocked, seed)
    scorers = {
        "SVD": lambda z: dti.score_svd(z, svd),
        "NeuralMF": lambda z: dti.score_mf(z, mf),
        "LightGCN": lambda z: lgcn.score_edges(z, gcn),
    }
    conv = {m: auc(seen, random_neg, fn) for m, fn in scorers.items()}

    candidates = all_nonedges(L, R, blocked)

    sig_log2 = lambda e: (dti.dbin(dl.get(e[0], 0)), dti.dbin(dr.get(e[1], 0)))
    mp1, mn1 = pool_match(np.random.default_rng(seed + 20000), seen, candidates, sig_log2)

    ql, left_cuts = quantile_mapper([dl.get(x, 0) for x in L], 10)
    qr, right_cuts = quantile_mapper([dr.get(x, 0) for x in R], 10)
    sig_q = lambda e: (ql(dl.get(e[0], 0)), qr(dr.get(e[1], 0)))
    mp2, mn2 = pool_match(np.random.default_rng(seed + 30000), seen, candidates, sig_q)

    mp3, mn3 = nearest_match(seen, candidates, dl, dr)
    pairs = {
        "log2_bins": (mp1, mn1),
        "degree_deciles": (mp2, mn2),
        "nearest_caliper": (mp3, mn3),
    }
    rows = []
    for method, (mp, mn) in pairs.items():
        coverage = len(mp) / len(seen) if seen else 0.0
        for model, fn in scorers.items():
            rows.append({
                "seed": seed,
                "method": method,
                "model": model,
                "n_seen": len(seen),
                "n_matched": len(mp),
                "matching_coverage": coverage,
                "auc_conventional": conv[model],
                "auc_neutralized": auc(mp, mn, fn) if mp else float("nan"),
            })
    meta = {"seed": seed, "left_decile_cuts_log1p": left_cuts, "right_decile_cuts_log1p": right_cuts}
    return rows, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seeds", default="0,1,2,3,4,5,6,7,8,9")
    a = ap.parse_args()
    E = dti.parse(a.input)
    allrows, metas = [], []
    for seed in [int(x) for x in a.seeds.split(",") if x.strip()]:
        rows, meta = one_seed(E, seed)
        allrows.extend(rows)
        metas.append(meta)
    df = pd.DataFrame(allrows)
    p = Path(a.out_prefix)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(p) + "_replicates.csv", index=False)

    winner_rows = []
    for seed in sorted(df.seed.unique()):
        dseed = df[df.seed == seed]
        convtab = dseed.drop_duplicates("model").set_index("model").auc_conventional
        convwinner = str(convtab.idxmax())
        for method in METHODS:
            z = dseed[dseed.method == method].set_index("model")
            winner_rows.append({
                "seed": int(seed),
                "method": method,
                "conventional_winner": convwinner,
                "neutralized_winner": str(z.auc_neutralized.idxmax()),
                "winner_changed": convwinner != str(z.auc_neutralized.idxmax()),
                "matching_coverage": float(z.matching_coverage.iloc[0]),
            })
    winners = pd.DataFrame(winner_rows)
    winners.to_csv(str(p) + "_winners.csv", index=False)

    summary = {
        "protocol": "MATCHING_ROBUSTNESS_AMENDMENT_20260913.md",
        "baseline_commit": "1298e4c83ff4fa488a784d52219d20e706118a9e",
        "seeds": sorted(int(x) for x in df.seed.unique()),
        "methods": {},
    }
    for method in METHODS:
        z = df[df.method == method]
        w = winners[winners.method == method]
        summary["methods"][method] = {
            "matching_coverage_mean": float(w.matching_coverage.mean()),
            "winner_change_fraction": float(w.winner_changed.mean()),
            "neutralized_winner_counts": {str(k): int(v) for k, v in w.neutralized_winner.value_counts().items()},
            "mean_auc": {m: {
                "conventional": float(z[z.model == m].auc_conventional.mean()),
                "neutralized": float(z[z.model == m].auc_neutralized.mean()),
            } for m in ["SVD", "NeuralMF", "LightGCN"]},
        }
    Path(str(p) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    Path(str(p) + "_bin_metadata.json").write_text(json.dumps(metas, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
