#!/usr/bin/env python3
from __future__ import annotations

import argparse, json, math
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def dbin(d: int) -> int:
    return 0 if d <= 0 else int(math.floor(math.log2(d))) + 1


def degree_score(df, ddrug, dprot):
    return np.asarray([
        math.log1p(ddrug.get(s, 0)) + math.log1p(dprot.get(p, 0))
        for s, p in zip(df.SMILES.astype(str), df.Protein.astype(str))
    ], dtype=float)


def matched_indices(df, ddrug, dprot, rng):
    pos = defaultdict(list); neg = defaultdict(list)
    for i, r in df.reset_index(drop=True).iterrows():
        key = (dbin(ddrug.get(str(r.SMILES), 0)), dbin(dprot.get(str(r.Protein), 0)))
        (pos if int(r.Y) == 1 else neg)[key].append(i)
    out = []
    npos = 0
    for key in sorted(set(pos) & set(neg)):
        n = min(len(pos[key]), len(neg[key]))
        if n <= 0:
            continue
        pp = rng.choice(pos[key], size=n, replace=False)
        nn = rng.choice(neg[key], size=n, replace=False)
        out.extend(pp.tolist()); out.extend(nn.tolist()); npos += n
    return np.asarray(sorted(out), dtype=int), npos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train', required=True)
    ap.add_argument('--test', required=True)
    ap.add_argument('--pred', required=True)
    ap.add_argument('--seed-label', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--match-reps', type=int, default=100)
    a = ap.parse_args()

    tr = pd.read_csv(a.train)
    te = pd.read_csv(a.test).reset_index(drop=True)
    pr = pd.read_csv(a.pred).reset_index(drop=True)
    req = {'SMILES','Protein','Y'}
    if not req.issubset(tr.columns) or not req.issubset(te.columns):
        raise ValueError('GraphBAN CSVs must contain SMILES, Protein, Y')
    if len(pr) != len(te) or not {'pred','target'}.issubset(pr.columns):
        raise ValueError('Prediction CSV is not aligned with test rows')
    if not np.array_equal(pr.target.astype(int).to_numpy(), te.Y.astype(int).to_numpy()):
        raise ValueError('Prediction targets do not match test Y row order')

    trp = tr[tr.Y.astype(int) == 1]
    ddrug = Counter(trp.SMILES.astype(str))
    dprot = Counter(trp.Protein.astype(str))
    y = te.Y.astype(int).to_numpy()
    score = pr.pred.astype(float).to_numpy()
    dscore = degree_score(te, ddrug, dprot)

    conv_auc = float(roc_auc_score(y, score))
    conv_deg_auc = float(roc_auc_score(y, dscore))
    npos_total = int((y == 1).sum())

    rows = []
    for r in range(a.match_reps):
        idx, npos = matched_indices(te, ddrug, dprot, np.random.default_rng(910000 + r))
        if npos == 0:
            raise RuntimeError('No degree-bin matched test examples')
        yy = y[idx]
        rows.append({
            'replicate': r,
            'n_pairs_per_class': npos,
            'positive_coverage': npos / npos_total,
            'graphban_auc_matched': float(roc_auc_score(yy, score[idx])),
            'degree_only_auc_matched': float(roc_auc_score(yy, dscore[idx])),
        })
    R = pd.DataFrame(rows)
    p = Path(a.out_prefix); p.parent.mkdir(parents=True, exist_ok=True)
    R.to_csv(str(p) + '_matching_replicates.csv', index=False)
    s = {
        'seed': str(a.seed_label),
        'n_train_rows': len(tr),
        'n_train_positives': len(trp),
        'n_test_rows': len(te),
        'n_test_positives': npos_total,
        'graphban_auc_conventional': conv_auc,
        'degree_only_auc_conventional': conv_deg_auc,
        'match_replicates': a.match_reps,
        'matching_positive_coverage_mean': float(R.positive_coverage.mean()),
        'graphban_auc_matched_mean': float(R.graphban_auc_matched.mean()),
        'graphban_auc_matched_sd': float(R.graphban_auc_matched.std(ddof=1)),
        'graphban_auc_matched_min': float(R.graphban_auc_matched.min()),
        'graphban_auc_matched_max': float(R.graphban_auc_matched.max()),
        'graphban_auc_drop_mean': float(conv_auc - R.graphban_auc_matched.mean()),
        'degree_only_auc_matched_mean': float(R.degree_only_auc_matched.mean()),
        'degree_only_auc_drop_mean': float(conv_deg_auc - R.degree_only_auc_matched.mean()),
    }
    Path(str(p) + '_summary.json').write_text(json.dumps(s, indent=2) + '\n')
    md = f"""# GraphBAN BioSNAP structural audit — seed {a.seed_label}\n\n- Conventional GraphBAN AUROC: **{conv_auc:.3f}**\n- Degree-matched GraphBAN AUROC: **{s['graphban_auc_matched_mean']:.3f} ± {s['graphban_auc_matched_sd']:.3f}**\n- Mean GraphBAN AUROC change: **{s['graphban_auc_drop_mean']:.3f}**\n- Matching positive coverage: **{s['matching_positive_coverage_mean']:.3f}**\n- Degree-only AUROC: **{conv_deg_auc:.3f} conventional → {s['degree_only_auc_matched_mean']:.3f} matched**\n\nThe checkpoint and official BioSNAP split come from the pinned upstream GraphBAN repository. Matching uses positive training degrees only and 100 deterministic within-stratum replicates.\n"""
    Path(str(p) + '.md').write_text(md)
    print(md)


if __name__ == '__main__':
    main()
