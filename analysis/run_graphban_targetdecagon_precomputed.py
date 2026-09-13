#!/usr/bin/env python3
"""Run the frozen leakage-free GraphBAN S6 from precomputed frozen features.

Only feature materialization is parallelized. Split construction, scaling,
training, negative sampling, matching, metrics and model architecture are
imported unchanged from run_graphban_targetdecagon_clean.py.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from run_graphban_targetdecagon_clean import parse_edges, one_seed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--drug-map", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--drug-features", required=True)
    ap.add_argument("--protein-features", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seeds", default="0,1,2")
    args = ap.parse_args()

    edges = parse_edges(args.input)
    dd = pd.read_csv(args.drug_map, dtype=str).fillna("")
    gd = pd.read_csv(args.gene_map, dtype=str).fillna("")
    dd = dd[dd.smiles.astype(bool)].drop_duplicates("drug")
    gd = gd[gd.sequence.astype(bool)].drop_duplicates("gene")
    dm = dict(zip(dd.drug, dd.smiles))
    gm = dict(zip(gd.gene, gd.sequence))
    mapped = [e for e in edges if e[0] in dm and e[1] in gm]
    coverage = len(mapped) / len(edges)
    if coverage < 0.90:
        raise RuntimeError(f"predeclared GraphBAN mapping gate failed: {coverage:.4f}")

    drugs = sorted({a for a, _ in mapped})
    genes = sorted({b for _, b in mapped})
    drug_features = np.load(args.drug_features).astype("float32")
    protein_features = np.load(args.protein_features).astype("float32")
    if drug_features.shape[0] != len(drugs):
        raise RuntimeError(
            f"drug feature rows {drug_features.shape[0]} != canonical drugs {len(drugs)}"
        )
    if protein_features.shape[0] != len(genes):
        raise RuntimeError(
            f"protein feature rows {protein_features.shape[0]} != canonical genes {len(genes)}"
        )

    di = {x: i for i, x in enumerate(drugs)}
    gi = {x: i for i, x in enumerate(genes)}
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device", device)

    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    rows = [
        one_seed(
            mapped,
            edges,
            drugs,
            genes,
            drug_features,
            protein_features,
            di,
            gi,
            seed,
            device,
        )
        for seed in seeds
    ]
    out = pd.DataFrame(rows)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(str(prefix) + "_replicates.csv", index=False)

    summary = {
        "benchmark": "BioSNAP TargetDecagon DTI",
        "architecture": "leakage-free GraphBAN-style transductive GraphSAGE + ChemBERTa + ESM-1b",
        "feature_execution": "deterministic precomputed parallel features; scientific protocol unchanged",
        "mapping_edge_fraction": coverage,
        "seeds": [int(x) for x in out.seed],
        "auc_random_mean": float(out.auc_random.mean()),
        "auc_random_sd": float(out.auc_random.std(ddof=1)),
        "auc_matched_mean": float(out.auc_matched.mean()),
        "auc_matched_sd": float(out.auc_matched.std(ddof=1)),
        "auc_drop_mean": float(out.auc_drop.mean()),
        "auprc_random_mean": float(out.auprc_random.mean()),
        "auprc_matched_mean": float(out.auprc_matched.mean()),
        "auprc_drop_mean": float(out.auprc_drop.mean()),
        "match_fraction_mean": float(out.match_fraction.mean()),
        "heldout_edges_in_message_passing": False,
        "epochs": 20,
        "hidden_dim": 256,
    }
    Path(str(prefix) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    md = f"""# Leakage-free contemporary GraphBAN challenge — TargetDecagon DTI

Mapping edge coverage: **{coverage:.1%}**. Held-out positives in message passing: **no**.
Feature execution: deterministic parallel precomputation using the frozen ChemBERTa and ESM-1b definitions; downstream S6 scientific protocol unchanged.

| Metric | Conventional | Structure-neutralized | Difference |
|---|---:|---:|---:|
| AUROC | {summary['auc_random_mean']:.3f} ± {summary['auc_random_sd']:.3f} | {summary['auc_matched_mean']:.3f} ± {summary['auc_matched_sd']:.3f} | {summary['auc_drop_mean']:+.3f} conventional-minus-neutralized |
| AUPRC | {summary['auprc_random_mean']:.3f} | {summary['auprc_matched_mean']:.3f} | {summary['auprc_drop_mean']:+.3f} |

Mean matching coverage: **{summary['match_fraction_mean']:.3f}** across seeds {summary['seeds']}.

This is the primary S6 contemporary-model challenge defined before outcome inspection. The public GraphBAN transductive test-graph construction is not used here because held-out positive edges must not enter message passing. By the frozen S6 rule, either a sensitive or comparatively robust GraphBAN result closes the contemporary-model gate if the execution is valid and the result is incorporated without cherry-picking.
"""
    Path(str(prefix) + ".md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
