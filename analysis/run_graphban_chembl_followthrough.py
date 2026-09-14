#!/usr/bin/env python3
"""Frozen fair-universe GraphBAN -> ChEMBL follow-through.

Protocol: protocols/GRAPHBAN_CHEMBL_FOLLOWTHROUGH_20260914.md
No ChEMBL information enters model fitting or model selection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler

import run_chembl_dti_external as ext
import run_graphban_targetdecagon_clean as gb


def graphban_ensemble(mapped, L, R, Xd, Xp, ninit, device):
    blocked = set(mapped)
    data, li, ri = gb.make_graph(L, R, mapped, Xd, Xp, device)
    candidates = [(a, b) for a in L for b in R if (a, b) not in blocked]
    flat = np.asarray([li[a] * len(R) + ri[b] for a, b in candidates], dtype=np.int64)
    score_runs = []
    for i in range(ninit):
        seed = 91000 + i
        model = gb.train_model(data, li, ri, mapped, blocked, L, R, seed, device, 20)
        s = gb.model_scores(model, data, candidates, li, ri, device, batch=8192).astype('float32')
        score_runs.append(s)
        print(f'GraphBAN full-graph fit {i+1}/{ninit}: seed={seed}, candidates={len(candidates)}')
    S = np.vstack(score_runs)
    mean = S.mean(axis=0)
    ranked = flat[np.argsort(mean)[::-1]]
    loo = {}
    for k in (100, 500, 1000):
        top = set(ranked[:k].tolist())
        vals = []
        for i in range(ninit):
            m = np.delete(S, i, axis=0).mean(axis=0)
            r = flat[np.argsort(m)[::-1]]
            vals.append(1.0 - len(top & set(r[:k].tolist())) / k)
        loo[str(k)] = {
            'mean': float(np.mean(vals)),
            'max': float(np.max(vals)),
            'values': [float(x) for x in vals],
        }
    return ranked, loo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--biosnap', required=True)
    ap.add_argument('--drug-map', required=True)
    ap.add_argument('--gene-map', required=True)
    ap.add_argument('--drug-features', required=True)
    ap.add_argument('--protein-features', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--inits', type=int, default=5)
    args = ap.parse_args()

    status = ext.get_json('https://www.ebi.ac.uk/chembl/api/data/status.json')
    version = status.get('chembl_db_version') or status.get('version')
    if version != 'ChEMBL_37':
        raise RuntimeError(f'Frozen protocol requires ChEMBL_37, got {version!r}')

    edges = gb.parse_edges(args.biosnap)
    dd = pd.read_csv(args.drug_map, dtype=str).fillna('')
    gd = pd.read_csv(args.gene_map, dtype=str).fillna('')
    dd = dd[dd.smiles.astype(bool)].drop_duplicates('drug')
    gd = gd[gd.sequence.astype(bool)].drop_duplicates('gene')
    dm = set(dd.drug)
    gm = set(gd.gene)
    mapped = sorted([e for e in edges if e[0] in dm and e[1] in gm])
    coverage = len(mapped) / len(edges)
    if coverage < 0.90:
        raise RuntimeError(f'GraphBAN mapping gate unexpectedly failed: {coverage:.4f}')

    L = sorted({a for a, _ in mapped})
    R = sorted({b for _, b in mapped})
    li = {x: i for i, x in enumerate(L)}
    ri = {x: i for i, x in enumerate(R)}
    blocked = set(mapped)

    Xd0 = np.load(args.drug_features).astype('float32')
    Xp0 = np.load(args.protein_features).astype('float32')
    if Xd0.shape[0] != len(L) or Xp0.shape[0] != len(R):
        raise RuntimeError(f'feature rows mismatch: drugs {Xd0.shape[0]}/{len(L)}, proteins {Xp0.shape[0]}/{len(R)}')
    Xd = StandardScaler().fit_transform(Xd0).astype('float32')
    Xp = StandardScaler().fit_transform(Xp0).astype('float32')

    # Baseline ensembles use the already frozen five-fit external-validation definitions,
    # now on exactly the GraphBAN-mapped relation universe.
    L0, R0, baseline_scores = ext.fit_ensembles(mapped, args.inits)
    if L0 != L or R0 != R:
        raise RuntimeError('baseline and GraphBAN canonical endpoint orders differ')

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    graphban_ranked, graphban_loo = graphban_ensemble(mapped, L, R, Xd, Xp, args.inits, device)

    # Re-open the independently frozen evidence source without changing its rule.
    chembl_drug_map = ext.map_pubchem_to_chembl(L)
    _, gene_targets = ext.target_maps(R)
    support6, support7, support_exact, evidence_meta, raw = ext.evidence_sets(edges, chembl_drug_map, gene_targets, 6.0)

    mask = np.ones((len(L), len(R)), dtype=bool)
    for a, b in blocked:
        mask[li[a], ri[b]] = False
    eligible_n = int(mask.sum())
    sf6 = ext.flat_support(support6, li, ri, len(R))
    sf7 = ext.flat_support(support7, li, ri, len(R))
    sfexact = ext.flat_support(support_exact, li, ri, len(R))

    rankings = {}
    for name, sc in baseline_scores.items():
        rankings[name] = ext.ranked_indices(sc, mask)
    rankings['GraphBAN'] = graphban_ranked

    primary_rows = []
    strict_rows = []
    exact_rows = []
    for name, ranked in rankings.items():
        primary_rows += ext.metrics_for_rank(name, ranked, sf6, eligible_n)
        for r in ext.metrics_for_rank(name, ranked, sf7, eligible_n):
            r['sensitivity'] = 'pchembl_ge_7'
            strict_rows.append(r)
        for r in ext.metrics_for_rank(name, ranked, sfexact, eligible_n):
            r['sensitivity'] = 'exact_relation_pchembl_ge_6'
            exact_rows.append(r)

    p = Path(args.out_prefix)
    p.parent.mkdir(parents=True, exist_ok=True)
    primary = pd.DataFrame(primary_rows)
    sensitivity = pd.DataFrame(strict_rows + exact_rows)
    primary.to_csv(str(p) + '_primary_table.csv', index=False)
    sensitivity.to_csv(str(p) + '_sensitivity_table.csv', index=False)
    raw.to_csv(str(p) + '_eligible_activity_records.csv', index=False)

    graphban_primary = primary[primary.model == 'GraphBAN'].copy()
    graphban_strict = sensitivity[(sensitivity.model == 'GraphBAN') & (sensitivity.sensitivity == 'pchembl_ge_7')].copy()
    summary = {
        'benchmark': 'BioSNAP TargetDecagon DTI - fair GraphBAN-mapped universe',
        'chembl_version': version,
        'edge_count_full': len(edges),
        'edge_count_mapped': len(mapped),
        'mapping_fraction': coverage,
        'candidate_universe_size': eligible_n,
        'ensemble_fits_per_model': args.inits,
        'selected_model_conventional': 'GraphBAN',
        'selected_model_neutralized': 'GraphBAN',
        'benchmark_induced_model_switch': False,
        'graphban_training_seeds': [91000 + i for i in range(args.inits)],
        'graphban_loo_turnover': graphban_loo,
        'evidence_meta': evidence_meta,
        'supported_pairs_primary_in_ranked_universe': len(sf6),
        'supported_pairs_strict_in_ranked_universe': len(sf7),
        'graphban_primary': graphban_primary[['k','hits','precision_lower_bound','recall_supported','enrichment_uniform']].to_dict('records'),
        'graphban_strict': graphban_strict[['k','hits','precision_lower_bound','recall_supported','enrichment_uniform']].to_dict('records'),
    }
    Path(str(p) + '_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    piv = primary.pivot(index='model', columns='k', values='hits')
    lines = [
        '# Fair-universe GraphBAN -> ChEMBL 37 follow-through',
        '',
        f'Mapped TargetDecagon coverage: **{coverage:.1%}**. Candidate universe: **{eligible_n:,}** unknown mapped drug-target pairs.',
        f'ChEMBL primary supported candidates in this universe: **{len(sf6)}**; strict pChEMBL >= 7: **{len(sf7)}**.',
        '',
        'GraphBAN is the mean-AUROC winner under both conventional and structure-neutralized evaluation in the frozen fair four-model panel. Therefore the contemporary DTI panel has **no benchmark-induced model-identity switch**; the same GraphBAN ranking is the relevant external-evidence ranking for both selection regimes.',
        '',
        '## Primary ChEMBL support (pChEMBL >= 6)',
        '',
        '| Model | @100 | @500 | @1,000 | @5,000 | @10,000 | @50,000 |',
        '|---|---:|---:|---:|---:|---:|---:|',
    ]
    for model in ['SVD','NeuralMF','LightGCN','GraphBAN']:
        row = piv.loc[model]
        vals = [int(row.get(k, 0)) for k in (100,500,1000,5000,10000,50000)]
        lines.append('| ' + model + ' | ' + ' | '.join(map(str, vals)) + ' |')
    lines += [
        '',
        '## GraphBAN ranking stability',
        '',
        f"Leave-one-fit-out turnover relative to the five-fit ensemble: HT@100 mean **{graphban_loo['100']['mean']:.3f}**, HT@500 **{graphban_loo['500']['mean']:.3f}**, HT@1000 **{graphban_loo['1000']['mean']:.3f}**.",
        '',
        '**Interpretation boundary.** The full contemporary DTI panel does not support a benchmark-induced winner reversal: GraphBAN remains selected under both evaluation regimes despite a material performance drop after structural neutralization. ChEMBL therefore characterizes the independently supported content of the common GraphBAN discovery frontier rather than a benchmark-induced DTI switch. The PPI experiment remains the clean full chain from benchmark rule to winner reversal, hypothesis turnover, and later-evidence consequence.',
    ]
    Path(str(p) + '.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
