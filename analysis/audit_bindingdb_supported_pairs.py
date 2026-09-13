#!/usr/bin/env python3
"""Transparent pair-level audit for the preregistered BindingDB expansion.

This script adds no new selection rule or inferential endpoint. It enumerates
ALL <=1000 nM BindingDB-supported BioSNAP candidate pairs and records their
ranks under every frozen ranking model so biological examples cannot be chosen
post hoc by perceived interest.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

import run_biosnap_dti_gate2_models as dti
import run_bindingdb_external as bdb


def rank_supported(ranked, support_flat):
    need = set(int(x) for x in support_flat)
    out = {}
    for i, x in enumerate(ranked, 1):
        xx = int(x)
        if xx in need:
            out[xx] = i
            if len(out) == len(need):
                break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--biosnap', required=True)
    ap.add_argument('--bindingdb', required=True)
    ap.add_argument('--gene-map', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--inits', type=int, default=5)
    a = ap.parse_args()

    edges = dti.parse(a.biosnap)
    support, strict, _ = bdb.load_bindingdb(Path(a.bindingdb), edges, Path(a.gene_map))
    L, R, li, ri, mask, rankings = bdb.full_rankings(edges, a.inits)
    sf = bdb.flat_support(support, li, ri, len(R))
    inv = {m: rank_supported(r, sf) for m, r in rankings.items()}

    rows = []
    for pair, rec in sorted(support.items()):
        if pair[0] not in li or pair[1] not in ri:
            continue
        flat = li[pair[0]] * len(R) + ri[pair[1]]
        row = {
            'drug': pair[0],
            'gene': pair[1],
            **rec,
            'strict_le_100nM': pair in strict,
            'published_2018_or_later': pd.to_datetime(rec['publication_date'], errors='coerce') >= pd.Timestamp('2018-01-01'),
        }
        for m in ['SVD', 'NeuralMF', 'LightGCN']:
            rank = inv.get(m, {}).get(flat)
            row[f'{m}_rank'] = rank
            row[f'{m}_top1000'] = bool(rank is not None and rank <= 1000)
            row[f'{m}_top5000'] = bool(rank is not None and rank <= 5000)
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(['SVD_rank', 'NeuralMF_rank'], na_position='last')
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(out.to_string(index=False))


if __name__ == '__main__':
    main()
