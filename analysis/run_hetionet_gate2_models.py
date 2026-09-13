#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import roc_auc_score


def dbin(d: int) -> int:
    return 0 if d <= 0 else int(math.floor(math.log2(d))) + 1


def parse_edges(path: str, metaedge: str, left_prefix: str, right_prefix: str):
    d = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    d.columns = [c.lower() for c in d.columns]
    q = d[d.metaedge.eq(metaedge)]
    edges = sorted({
        (str(a).strip(), str(b).strip())
        for a, b in zip(q.source, q.target)
        if str(a).startswith(left_prefix) and str(b).startswith(right_prefix)
    })
    if len(edges) < 100:
        raise ValueError(f"{metaedge}: only {len(edges)} edges parsed")
    return edges


def sample_random_nonedges(rng, n, left, right, blocked):
    out = set()
    limit = max(200000, n * 300)
    k = 0
    while len(out) < n and k < limit:
        k += 1
        e = (left[int(rng.integers(len(left)))], right[int(rng.integers(len(right)))])
        if e in blocked or e in out:
            continue
        out.add(e)
    if len(out) < n:
        raise RuntimeError(f"sampled only {len(out)}/{n} random nonedges")
    return list(out)


def targeted_match(rng, positives, left, right, blocked, dl, dr):
    lb, rb = defaultdict(list), defaultdict(list)
    for x in left:
        lb[dbin(dl.get(x, 0))].append(x)
    for x in right:
        rb[dbin(dr.get(x, 0))].append(x)
    used, mp, mn = set(), [], []
    for i in rng.permutation(len(positives)):
        p = positives[int(i)]
        A = lb.get(dbin(dl.get(p[0], 0)), [])
        B = rb.get(dbin(dr.get(p[1], 0)), [])
        found = None
        for _ in range(1200):
            if not A or not B:
                break
            e = (A[int(rng.integers(len(A)))], B[int(rng.integers(len(B)))])
            if e not in blocked and e not in used:
                found = e
                break
        if found is not None:
            used.add(found)
            mp.append(p)
            mn.append(found)
    return mp, mn


def sigmoid(x):
    return 1 / (1 + np.exp(-np.clip(x, -30, 30)))


def fit_svd(train, left, right, k=32, seed=0):
    li = {x: i for i, x in enumerate(left)}
    ri = {x: i for i, x in enumerate(right)}
    rr = [li[a] for a, _ in train]
    cc = [ri[b] for _, b in train]
    X = sparse.csr_matrix((np.ones(len(rr)), (rr, cc)), shape=(len(left), len(right)))
    ncomp = max(1, min(k, min(X.shape) - 1))
    svd = TruncatedSVD(n_components=ncomp, random_state=seed)
    U = svd.fit_transform(X)
    V = svd.components_.T
    return li, ri, U, V


def score_svd(edges, model):
    li, ri, U, V = model
    return np.array([np.dot(U[li[a]], V[ri[b]]) for a, b in edges])


def fit_logistic_mf(train, left, right, blocked, rng, k=24, epochs=30, lr=0.04, reg=2e-4):
    li = {x: i for i, x in enumerate(left)}
    ri = {x: i for i, x in enumerate(right)}
    U = rng.normal(0, 0.08, (len(left), k))
    V = rng.normal(0, 0.08, (len(right), k))
    bu = np.zeros(len(left))
    bv = np.zeros(len(right))
    pos = np.array([(li[a], ri[b]) for a, b in train], dtype=int)
    batch = min(1024, max(64, len(pos)))
    for _ in range(epochs):
        rng.shuffle(pos)
        for s in range(0, len(pos), batch):
            pp = pos[s:s + batch]
            n = len(pp)
            neg = []
            while len(neg) < n:
                a = left[int(rng.integers(len(left)))]
                b = right[int(rng.integers(len(right)))]
                if (a, b) not in blocked:
                    neg.append((li[a], ri[b]))
            nn = np.asarray(neg, dtype=int)
            ii = np.r_[pp[:, 0], nn[:, 0]]
            jj = np.r_[pp[:, 1], nn[:, 1]]
            y = np.r_[np.ones(n), np.zeros(n)]
            pred = sigmoid(np.sum(U[ii] * V[jj], axis=1) + bu[ii] + bv[jj])
            g = pred - y
            gu = g[:, None] * V[jj] + reg * U[ii]
            gv = g[:, None] * U[ii] + reg * V[jj]
            np.add.at(U, ii, -lr * gu)
            np.add.at(V, jj, -lr * gv)
            np.add.at(bu, ii, -lr * g)
            np.add.at(bv, jj, -lr * g)
    return li, ri, U, V, bu, bv


def score_mf(edges, model):
    li, ri, U, V, bu, bv = model
    return np.array([np.dot(U[li[a]], V[ri[b]]) + bu[li[a]] + bv[ri[b]] for a, b in edges])


def eval_auc(pos, neg, scorer):
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    s = np.r_[scorer(pos), scorer(neg)]
    return float(roc_auc_score(y, s))


def one(edges, seed):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(edges))
    nt = max(1, int(round(0.2 * len(edges))))
    test = [edges[i] for i in idx[:nt]]
    train = [edges[i] for i in idx[nt:]]
    blocked = set(edges)
    left = sorted({a for a, _ in edges})
    right = sorted({b for _, b in edges})
    dl = Counter(a for a, _ in train)
    dr = Counter(b for _, b in train)
    seen = [e for e in test if dl.get(e[0], 0) > 0 and dr.get(e[1], 0) > 0]
    rnd = sample_random_nonedges(rng, len(seen), left, right, blocked)
    mp, mn = targeted_match(rng, seen, left, right, blocked, dl, dr)

    svd = fit_svd(train, left, right, seed=seed)
    mf = fit_logistic_mf(train, left, right, blocked, np.random.default_rng(seed + 10000))

    sr = eval_auc(seen, rnd, lambda z: score_svd(z, svd))
    sm = eval_auc(mp, mn, lambda z: score_svd(z, svd))
    nr = eval_auc(seen, rnd, lambda z: score_mf(z, mf))
    nm = eval_auc(mp, mn, lambda z: score_mf(z, mf))
    return {
        "seed": seed,
        "n_test_seen": len(seen),
        "seen_fraction": len(seen) / len(test),
        "n_matched": len(mp),
        "match_fraction": len(mp) / len(seen),
        "svd_auc_random": sr,
        "svd_auc_matched": sm,
        "svd_drop": sr - sm,
        "neuralmf_auc_random": nr,
        "neuralmf_auc_matched": nm,
        "neuralmf_drop": nr - nm,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metaedge", required=True)
    ap.add_argument("--left-prefix", required=True)
    ap.add_argument("--right-prefix", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--replicates", type=int, default=5)
    a = ap.parse_args()

    edges = parse_edges(a.input, a.metaedge, a.left_prefix, a.right_prefix)
    out = pd.DataFrame([one(edges, s) for s in range(a.replicates)])
    p = Path(a.out_prefix)
    p.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(str(p) + "_replicates.csv", index=False)

    summary = {
        "benchmark": a.label,
        "metaedge": a.metaedge,
        "edge_count": len(edges),
        "replicates": a.replicates,
        "match_fraction_mean": float(out.match_fraction.mean()),
        "seen_fraction_mean": float(out.seen_fraction.mean()),
    }
    for m in ["svd", "neuralmf"]:
        for x in ["auc_random", "auc_matched", "drop"]:
            summary[f"{m}_{x}_mean"] = float(out[f"{m}_{x}"].mean())
            summary[f"{m}_{x}_sd"] = float(out[f"{m}_{x}"].std(ddof=1))
    Path(str(p) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    md = f"""# Gate 2 pilot — {a.label}\n\n| Model | Random-negative AUC | Degree-matched AUC | Mean drop |\n|---|---:|---:|---:|\n| Truncated-SVD latent factors | {summary['svd_auc_random_mean']:.3f} ± {summary['svd_auc_random_sd']:.3f} | {summary['svd_auc_matched_mean']:.3f} ± {summary['svd_auc_matched_sd']:.3f} | {summary['svd_drop_mean']:.3f} |\n| Logistic neural matrix factorization | {summary['neuralmf_auc_random_mean']:.3f} ± {summary['neuralmf_auc_random_sd']:.3f} | {summary['neuralmf_auc_matched_mean']:.3f} ± {summary['neuralmf_auc_matched_sd']:.3f} | {summary['neuralmf_drop_mean']:.3f} |\n\nMean matching coverage: **{summary['match_fraction_mean']:.3f}**; mean seen-endpoint test fraction: **{summary['seen_fraction_mean']:.3f}**. These are controlled structural latent models, not the final SOTA panel.\n"""
    Path(str(p) + ".md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
