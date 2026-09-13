#!/usr/bin/env python3
"""Generic structural-shortcut audit for pairwise biomedical benchmarks.

This is a model-free diagnostic. It asks whether endpoint popularity/degree alone
can discriminate positive edges from a supplied negative/counter-evidence set or
from randomly sampled unlabelled pairs, and how that discrimination changes when
endpoint degree is balanced.

It does not establish biological validity, clinical safety, or causal mechanism.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def norm_pair(a: str, b: str, directed: bool) -> tuple[str, str]:
    a, b = str(a).strip(), str(b).strip()
    if directed:
        return a, b
    return tuple(sorted((a, b)))


def auc_mann_whitney(pos_scores, neg_scores) -> float:
    p = np.asarray(pos_scores, dtype=float)
    n = np.asarray(neg_scores, dtype=float)
    if len(p) == 0 or len(n) == 0:
        return float("nan")
    ranks = pd.Series(np.concatenate([p, n])).rank(method="average").to_numpy()
    u = ranks[: len(p)].sum() - len(p) * (len(p) + 1) / 2
    return float(u / (len(p) * len(n)))


def average_precision(pos_scores, neg_scores) -> float:
    scores = np.concatenate([np.asarray(pos_scores, float), np.asarray(neg_scores, float)])
    labels = np.concatenate([np.ones(len(pos_scores), dtype=int), np.zeros(len(neg_scores), dtype=int)])
    if labels.sum() == 0:
        return float("nan")
    # Deterministic tie handling at score thresholds.
    order = np.argsort(-scores, kind="mergesort")
    scores, labels = scores[order], labels[order]
    tp = fp = 0
    ap_num = 0.0
    i = 0
    total_pos = int(labels.sum())
    while i < len(scores):
        j = i
        group_tp = group_fp = 0
        while j < len(scores) and scores[j] == scores[i]:
            if labels[j]:
                group_tp += 1
            else:
                group_fp += 1
            j += 1
        tp += group_tp
        fp += group_fp
        if group_tp:
            ap_num += group_tp * (tp / (tp + fp))
        i = j
    return float(ap_num / total_pos)


def degree_bin(d: int) -> int:
    # Prespecified log2-style bins: 0, 1, 2-3, 4-7, 8-15, ...
    if d <= 0:
        return 0
    return int(math.floor(math.log2(d))) + 1


def load_pairs(path: str, src_col: str, dst_col: str, directed: bool) -> list[tuple[str, str]]:
    d = pd.read_csv(path)
    missing = [c for c in (src_col, dst_col) if c not in d.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}; available={list(d.columns)}")
    out = set()
    for a, b in zip(d[src_col], d[dst_col]):
        if pd.isna(a) or pd.isna(b):
            continue
        p = norm_pair(a, b, directed)
        if p[0] == p[1]:
            continue
        out.add(p)
    return sorted(out)


def build_degrees(positives, bipartite: bool):
    if bipartite:
        left = Counter(a for a, _ in positives)
        right = Counter(b for _, b in positives)
        return left, right
    deg = Counter()
    for a, b in positives:
        deg[a] += 1
        deg[b] += 1
    return deg, deg


def score_pair(p, left_deg, right_deg, method: str) -> float:
    da = float(left_deg.get(p[0], 0))
    db = float(right_deg.get(p[1], 0))
    if method == "log_product":
        return float(np.log1p(da) * np.log1p(db))
    if method == "product":
        return da * db
    if method == "sum":
        return da + db
    if method == "min":
        return min(da, db)
    raise ValueError(method)


def signature(p, left_deg, right_deg, bipartite: bool):
    a = degree_bin(int(left_deg.get(p[0], 0)))
    b = degree_bin(int(right_deg.get(p[1], 0)))
    return (a, b) if bipartite else tuple(sorted((a, b)))


def sample_random_unlabelled(
    rng,
    n: int,
    left_vocab: list[str],
    right_vocab: list[str],
    blocked: set[tuple[str, str]],
    directed: bool,
    bipartite: bool,
) -> list[tuple[str, str]]:
    out = set()
    max_attempts = max(10000, n * 200)
    attempts = 0
    while len(out) < n and attempts < max_attempts:
        attempts += 1
        a = left_vocab[int(rng.integers(len(left_vocab)))]
        b = right_vocab[int(rng.integers(len(right_vocab)))]
        if not bipartite and a == b:
            continue
        p = norm_pair(a, b, directed or bipartite)
        if p in blocked or p in out:
            continue
        out.add(p)
    if len(out) < n:
        raise RuntimeError(
            f"Could sample only {len(out)} unlabelled pairs of requested {n}; "
            "graph may be too dense or vocab too small."
        )
    return sorted(out)


def degree_match(rng, positives, negatives, left_deg, right_deg, bipartite: bool):
    pools = defaultdict(list)
    for q in negatives:
        pools[signature(q, left_deg, right_deg, bipartite)].append(q)
    for pool in pools.values():
        rng.shuffle(pool)

    used = Counter()
    mp, mn = [], []
    for idx in rng.permutation(len(positives)):
        p = positives[int(idx)]
        s = signature(p, left_deg, right_deg, bipartite)
        pool = pools.get(s, [])
        u = used[s]
        if u >= len(pool):
            continue
        mp.append(p)
        mn.append(pool[u])
        used[s] += 1
    return mp, mn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positives", required=True)
    ap.add_argument("--positive-src", default="src")
    ap.add_argument("--positive-dst", default="dst")
    ap.add_argument("--negatives")
    ap.add_argument("--negative-src", default="src")
    ap.add_argument("--negative-dst", default="dst")
    ap.add_argument("--bipartite", action="store_true")
    ap.add_argument("--directed", action="store_true")
    ap.add_argument("--replicates", type=int, default=20)
    ap.add_argument("--random-pool-multiplier", type=int, default=5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    if args.bipartite and args.directed:
        raise ValueError("Use --bipartite for role-specific endpoints; do not also set --directed.")

    positives = load_pairs(
        args.positives, args.positive_src, args.positive_dst, directed=(args.directed or args.bipartite)
    )
    posset = set(positives)
    if not positives:
        raise ValueError("No positive pairs loaded.")

    left_deg, right_deg = build_degrees(positives, args.bipartite)
    if args.bipartite:
        left_vocab = sorted(left_deg)
        right_vocab = sorted(right_deg)
    else:
        vocab = sorted(set(left_deg) | set(right_deg))
        left_vocab = right_vocab = vocab

    curated = None
    overlap_excluded = 0
    if args.negatives:
        raw_neg = load_pairs(
            args.negatives,
            args.negative_src,
            args.negative_dst,
            directed=(args.directed or args.bipartite),
        )
        overlap_excluded = len(set(raw_neg) & posset)
        curated = sorted(set(raw_neg) - posset)
        if not curated:
            raise ValueError("No negative pairs remain after positive-overlap exclusion.")

    methods = ["log_product", "product", "sum", "min"]
    rows = []
    for seed in range(args.replicates):
        rng = np.random.default_rng(seed)
        n_eval = min(len(positives), len(curated)) if curated is not None else len(positives)
        if n_eval < 2:
            raise ValueError("Need at least two evaluable positives/negatives.")

        p_eval = [positives[i] for i in rng.choice(len(positives), size=n_eval, replace=False)]
        blocked = set(posset)
        if curated is not None:
            blocked |= set(curated)

        random_neg = sample_random_unlabelled(
            rng, n_eval, left_vocab, right_vocab, blocked, directed=args.directed, bipartite=args.bipartite
        )

        if curated is not None:
            match_source = curated
        else:
            pool_n = min(
                max(n_eval * args.random_pool_multiplier, n_eval),
                max(n_eval, len(left_vocab) * len(right_vocab) - len(posset) - 1),
            )
            match_source = sample_random_unlabelled(
                rng, pool_n, left_vocab, right_vocab, blocked, directed=args.directed, bipartite=args.bipartite
            )

        mp, mn = degree_match(rng, p_eval, match_source, left_deg, right_deg, args.bipartite)

        for method in methods:
            p_scores = [score_pair(x, left_deg, right_deg, method) for x in p_eval]
            r_scores = [score_pair(x, left_deg, right_deg, method) for x in random_neg]
            rec = {
                "seed": seed,
                "score": method,
                "n_positive": len(p_eval),
                "n_random": len(random_neg),
                "n_degree_matched": len(mp),
                "matched_fraction": len(mp) / len(p_eval),
                "auc_random": auc_mann_whitney(p_scores, r_scores),
                "ap_random": average_precision(p_scores, r_scores),
            }
            if curated is not None:
                c_scores = [score_pair(x, left_deg, right_deg, method) for x in curated]
                rec["auc_curated"] = auc_mann_whitney(p_scores, c_scores)
                rec["ap_curated"] = average_precision(p_scores, c_scores)
            if mp:
                mp_scores = [score_pair(x, left_deg, right_deg, method) for x in mp]
                mn_scores = [score_pair(x, left_deg, right_deg, method) for x in mn]
                rec["auc_degree_matched"] = auc_mann_whitney(mp_scores, mn_scores)
                rec["ap_degree_matched"] = average_precision(mp_scores, mn_scores)
            rows.append(rec)

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    summary = {
        "positive_pairs": len(positives),
        "positive_left_nodes": len(left_vocab),
        "positive_right_nodes": len(right_vocab),
        "bipartite": bool(args.bipartite),
        "directed": bool(args.directed),
        "curated_negative_pairs": None if curated is None else len(curated),
        "positive_negative_overlap_excluded": overlap_excluded,
        "replicates": int(args.replicates),
        "degree_bin_rule": "0; then floor(log2(degree))+1",
        "results": {},
        "note": (
            "Model-free structural diagnostic only. High performance indicates that endpoint "
            "popularity/degree contains label information under the evaluated benchmark construction; "
            "it does not establish or refute biological validity."
        ),
    }
    for method, g in out.groupby("score"):
        stats = {}
        for col in [
            "auc_random", "ap_random", "auc_curated", "ap_curated",
            "auc_degree_matched", "ap_degree_matched", "matched_fraction",
        ]:
            if col in g and g[col].notna().any():
                stats[col + "_mean"] = float(g[col].mean())
                stats[col + "_sd"] = float(g[col].std(ddof=1)) if len(g) > 1 else 0.0
        summary["results"][method] = stats

    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
