#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def dbin(d: int) -> int:
    return 0 if d <= 0 else int(math.floor(math.log2(d))) + 1


def load_edges(args):
    if args.format == "biosnap":
        d = pd.read_csv(args.input, compression="infer", comment="#", header=None,
                        names=["left", "right"], sep=",", dtype=str).dropna()
        edges = sorted({(str(a).strip(), str(b).strip()) for a, b in zip(d.left, d.right)})
    else:
        d = pd.read_csv(args.input, sep="\t", compression="infer", dtype=str)
        d.columns = [c.lower() for c in d.columns]
        q = d[d.metaedge.eq(args.metaedge)]
        edges = sorted({
            (str(a).strip(), str(b).strip())
            for a, b in zip(q.source, q.target)
            if str(a).startswith(args.left_prefix) and str(b).startswith(args.right_prefix)
        })
    if len(edges) < 100:
        raise ValueError(f"Only {len(edges)} edges parsed")
    return edges


def random_nonedges(rng, n, left, right, blocked):
    out = set()
    limit = max(200000, n * 300)
    k = 0
    while len(out) < n and k < limit:
        k += 1
        e = (left[int(rng.integers(len(left)))], right[int(rng.integers(len(right)))])
        if e not in blocked and e not in out:
            out.add(e)
    if len(out) < n:
        raise RuntimeError(f"Only {len(out)}/{n} random nonedges")
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


def fit_lightgcn(train, left, right, blocked, seed, dim=32, layers=2, epochs=35, lr=0.03, reg=1e-5):
    import torch

    torch.manual_seed(seed)
    rng = np.random.default_rng(seed + 12345)
    li = {x: i for i, x in enumerate(left)}
    ri = {x: i + len(left) for i, x in enumerate(right)}
    n_nodes = len(left) + len(right)

    # Symmetric normalized adjacency built from TRAINING EDGES ONLY.
    src, dst = [], []
    deg = np.zeros(n_nodes, dtype=np.float32)
    for a, b in train:
        i, j = li[a], ri[b]
        src += [i, j]
        dst += [j, i]
        deg[i] += 1
        deg[j] += 1
    vals = []
    for i, j in zip(src, dst):
        vals.append(1.0 / math.sqrt(max(deg[i], 1.0) * max(deg[j], 1.0)))
    idx = torch.tensor([src, dst], dtype=torch.long)
    val = torch.tensor(vals, dtype=torch.float32)
    A = torch.sparse_coo_tensor(idx, val, (n_nodes, n_nodes)).coalesce()

    emb0 = torch.nn.Embedding(n_nodes, dim)
    torch.nn.init.normal_(emb0.weight, std=0.08)
    opt = torch.optim.Adam(emb0.parameters(), lr=lr, weight_decay=reg)
    pos_i = torch.tensor([li[a] for a, _ in train], dtype=torch.long)
    pos_j = torch.tensor([ri[b] for _, b in train], dtype=torch.long)

    def propagated():
        z0 = emb0.weight
        zs = [z0]
        z = z0
        for _ in range(layers):
            z = torch.sparse.mm(A, z)
            zs.append(z)
        return torch.stack(zs, dim=0).mean(dim=0)

    for _ in range(epochs):
        neg_a, neg_b = [], []
        while len(neg_a) < len(train):
            a = left[int(rng.integers(len(left)))]
            b = right[int(rng.integers(len(right)))]
            if (a, b) not in blocked:
                neg_a.append(li[a])
                neg_b.append(ri[b])
        ni = torch.tensor(neg_a, dtype=torch.long)
        nj = torch.tensor(neg_b, dtype=torch.long)
        z = propagated()
        ps = (z[pos_i] * z[pos_j]).sum(dim=1)
        ns = (z[ni] * z[nj]).sum(dim=1)
        loss = torch.nn.functional.softplus(-ps).mean() + torch.nn.functional.softplus(ns).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()

    with torch.no_grad():
        z = propagated().cpu().numpy()
    return li, ri, z


def score_edges(edges, model):
    li, ri, z = model
    return np.asarray([float(np.dot(z[li[a]], z[ri[b]])) for a, b in edges])


def auc(pos, neg, model):
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    s = np.r_[score_edges(pos, model), score_edges(neg, model)]
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
    rnd = random_nonedges(rng, len(seen), left, right, blocked)
    mp, mn = targeted_match(rng, seen, left, right, blocked, dl, dr)
    model = fit_lightgcn(train, left, right, blocked, seed)
    ar = auc(seen, rnd, model)
    am = auc(mp, mn, model)
    return {
        "seed": seed,
        "n_seen": len(seen),
        "seen_fraction": len(seen) / len(test),
        "n_matched": len(mp),
        "match_fraction": len(mp) / len(seen),
        "lightgcn_auc_random": ar,
        "lightgcn_auc_matched": am,
        "lightgcn_drop": ar - am,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--format", choices=["biosnap", "hetionet"], required=True)
    ap.add_argument("--metaedge")
    ap.add_argument("--left-prefix", default="")
    ap.add_argument("--right-prefix", default="")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--replicates", type=int, default=3)
    args = ap.parse_args()

    edges = load_edges(args)
    out = pd.DataFrame([one(edges, s) for s in range(args.replicates)])
    p = Path(args.out_prefix)
    p.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(str(p) + "_replicates.csv", index=False)
    summary = {
        "benchmark": args.label,
        "edge_count": len(edges),
        "replicates": args.replicates,
        "lightgcn_auc_random_mean": float(out.lightgcn_auc_random.mean()),
        "lightgcn_auc_random_sd": float(out.lightgcn_auc_random.std(ddof=1)),
        "lightgcn_auc_matched_mean": float(out.lightgcn_auc_matched.mean()),
        "lightgcn_auc_matched_sd": float(out.lightgcn_auc_matched.std(ddof=1)),
        "lightgcn_drop_mean": float(out.lightgcn_drop.mean()),
        "lightgcn_drop_sd": float(out.lightgcn_drop.std(ddof=1)),
        "match_fraction_mean": float(out.match_fraction.mean()),
        "seen_fraction_mean": float(out.seen_fraction.mean()),
        "model": "2-layer LightGCN-style message passing, embeddings trained with random unlabelled negatives",
    }
    Path(str(p) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    md = f"""# Gate 2 — LightGCN on {args.label}\n\n- Random-negative AUC: **{summary['lightgcn_auc_random_mean']:.3f} ± {summary['lightgcn_auc_random_sd']:.3f}**\n- Degree-matched AUC: **{summary['lightgcn_auc_matched_mean']:.3f} ± {summary['lightgcn_auc_matched_sd']:.3f}**\n- Mean AUC drop: **{summary['lightgcn_drop_mean']:.3f}**\n- Matching coverage: **{summary['match_fraction_mean']:.3f}**\n- Seen-endpoint test fraction: **{summary['seen_fraction_mean']:.3f}**\n\nAdjacency and degree are constructed from training edges only; held-out positive edges are excluded from message passing.\n"""
    Path(str(p) + ".md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
