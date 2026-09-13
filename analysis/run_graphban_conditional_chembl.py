#!/usr/bin/env python3
"""Conditional four-model ChEMBL 37 validation on a common GraphBAN-mappable universe.

Protocol: GRAPHBAN_CONDITIONAL_H5_CHEMBL_AMENDMENT_20260913.md
This script hard-fails on drift from the frozen H6 ChEMBL mapping/evidence fingerprint.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_chembl_dti_external as h6

EXPECTED = {
    "mapped_drugs": 264,
    "total_drugs": 284,
    "mapped_genes": 1897,
    "total_genes": 3648,
    "mapped_candidate_universe": 486560,
    "n_mapped_molecules": 264,
    "n_mapped_targets_raw": 1876,
    "n_eligible_single_protein_human_targets": 1790,
    "n_activities_pchembl_ge_6_on_mapped_targets": 9509,
    "n_good_binding_conf9_assays": 4636,
    "n_supported_candidate_pairs_pchembl_ge_6": 66,
    "n_supported_candidate_pairs_pchembl_ge_7": 26,
    "n_supported_candidate_pairs_exact_relation_ge_6": 66,
}
MODELS = ("GraphBAN-style", "NeuralMF", "SVD", "LightGCN")
SCORE_KEYS = {
    "GraphBAN-style": "graphban_ensemble",
    "NeuralMF": "neuralmf",
    "SVD": "svd",
    "LightGCN": "lightgcn",
}


def fingerprint(edges, drug_map, gene_targets, emeta):
    drugs = sorted({a for a, _ in edges})
    genes = sorted({b for _, b in edges})
    li = {x: i for i, x in enumerate(drugs)}
    ri = {x: i for i, x in enumerate(genes)}
    mask = np.ones((len(drugs), len(genes)), dtype=bool)
    for a, b in edges:
        mask[li[a], ri[b]] = False
    mapped_drugs = {d for d in drugs if drug_map.get(d)}
    mapped_genes = {g for g in genes if gene_targets.get(g)}
    mm = np.zeros_like(mask)
    for d in mapped_drugs:
        mm[li[d], :] = True
    col = np.zeros(len(genes), dtype=bool)
    col[[ri[g] for g in mapped_genes]] = True
    mm &= col[None, :]
    mm &= mask
    return {
        "mapped_drugs": len(mapped_drugs),
        "total_drugs": len(drugs),
        "mapped_genes": len(mapped_genes),
        "total_genes": len(genes),
        "mapped_candidate_universe": int(mm.sum()),
        **emeta,
    }


def metrics_table(scores, mask, support_flat):
    rows = []
    for name in MODELS:
        ranked = h6.ranked_indices(scores[name], mask)
        rows.extend(h6.metrics_for_rank(name, ranked, support_flat, int(mask.sum())))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--biosnap", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    edges = h6.dti.parse(args.biosnap)
    all_drugs = sorted({a for a, _ in edges})
    all_genes = sorted({b for _, b in edges})

    status = h6.get_json("https://www.ebi.ac.uk/chembl/api/data/status.json")
    version = status.get("chembl_db_version") or status.get("version")
    if version != "ChEMBL_37":
        raise RuntimeError(f"frozen protocol requires ChEMBL_37, got {version!r}")

    drug_map = h6.map_pubchem_to_chembl(all_drugs)
    gene_acc, gene_targets = h6.target_maps(all_genes)
    S6, S7, Sexact, emeta, raw = h6.evidence_sets(edges, drug_map, gene_targets, 6.0)
    observed = fingerprint(edges, drug_map, gene_targets, emeta)
    drift = {k: {"expected": v, "observed": observed.get(k)} for k, v in EXPECTED.items() if observed.get(k) != v}

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(prefix) + "_fingerprint.json").write_text(
        json.dumps({"chembl_version": version, "expected": EXPECTED, "observed": observed, "drift": drift}, indent=2) + "\n"
    )
    if drift:
        raise RuntimeError(f"frozen ChEMBL fingerprint drift; conditional analysis stopped: {drift}")

    z = np.load(args.scores, allow_pickle=False)
    drugs = [str(x) for x in z["drugs"]]
    genes = [str(x) for x in z["genes"]]
    mask = z["mask"].astype(bool)
    if len(drugs) != 284 or len(genes) != 3607 or int(mask.size - mask.sum()) != 18631:
        raise RuntimeError("GraphBAN common-universe fingerprint drift")
    scores = {name: z[SCORE_KEYS[name]].astype(float) for name in MODELS}
    for name, sc in scores.items():
        if sc.shape != mask.shape:
            raise RuntimeError(f"score shape mismatch for {name}: {sc.shape} vs {mask.shape}")

    li = {x: i for i, x in enumerate(drugs)}
    ri = {x: i for i, x in enumerate(genes)}
    sf6 = h6.flat_support(S6, li, ri, len(genes))
    sf7 = h6.flat_support(S7, li, ri, len(genes))
    sfexact = h6.flat_support(Sexact, li, ri, len(genes))

    # Common GraphBAN-mappable universe for all four models.
    primary = metrics_table(scores, mask, sf6)
    primary.to_csv(str(prefix) + "_common_universe_primary.csv", index=False)

    # Common GraphBAN + ChEMBL-mapped universe for all four models.
    mapped_drugs = {d for d in drugs if drug_map.get(d)}
    mapped_genes = {g for g in genes if gene_targets.get(g)}
    mapped_mask = np.zeros_like(mask)
    for d in mapped_drugs:
        mapped_mask[li[d], :] = True
    col = np.zeros(len(genes), dtype=bool)
    col[[ri[g] for g in mapped_genes]] = True
    mapped_mask &= col[None, :]
    mapped_mask &= mask
    mapped = metrics_table(scores, mapped_mask, sf6)
    mapped.to_csv(str(prefix) + "_common_mapped_primary.csv", index=False)

    sens_rows = []
    for label, support in (("pchembl_ge_7", sf7), ("exact_relation_pchembl_ge_6", sfexact)):
        for name in MODELS:
            ranked = h6.ranked_indices(scores[name], mask)
            for row in h6.metrics_for_rank(name, ranked, support, int(mask.sum())):
                row["sensitivity"] = label
                sens_rows.append(row)
    pd.DataFrame(sens_rows).to_csv(str(prefix) + "_common_universe_sensitivity.csv", index=False)
    raw.to_csv(str(prefix) + "_eligible_activity_records.csv", index=False)

    summary = {
        "protocol": "GRAPHBAN_CONDITIONAL_H5_CHEMBL_AMENDMENT_20260913.md",
        "baseline_commit": "1298e4c83ff4fa488a784d52219d20e706118a9e",
        "chembl_version": version,
        "frozen_fingerprint_match": True,
        "common_graphban_drugs": len(drugs),
        "common_graphban_genes": len(genes),
        "common_unknown_candidate_universe": int(mask.sum()),
        "common_chembl_mapped_drugs": len(mapped_drugs),
        "common_chembl_mapped_genes": len(mapped_genes),
        "common_chembl_mapped_candidate_universe": int(mapped_mask.sum()),
        "supported_pairs_pchembl_ge_6_total_frozen": len(S6),
        "supported_pairs_pchembl_ge_6_in_common_universe": len(sf6),
        "supported_pairs_pchembl_ge_7_in_common_universe": len(sf7),
        "supported_pairs_exact_relation_ge_6_in_common_universe": len(sfexact),
        "interpretation_boundary": "Fixed-K support yield is compared only on identical GraphBAN-mappable candidate universes. Unlabelled ChEMBL pairs are not biological negatives.",
    }
    Path(str(prefix) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# Conditional GraphBAN-style external validation — ChEMBL 37",
        "",
        "Frozen ChEMBL fingerprint: **MATCH**.",
        f"Common GraphBAN-mappable unknown universe: **{int(mask.sum()):,}** pairs; supported pChEMBL>=6 pairs represented: **{len(sf6)}**.",
        "",
        "| Model | K | ChEMBL-supported hits | Precision lower bound | Recall of common supported pairs | Enrichment vs uniform |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, r in primary.iterrows():
        lines.append(
            f"| {r.model} | {int(r.k):,} | {int(r.hits)} | {r.precision_lower_bound:.5f} | {r.recall_supported:.5f} | {r.enrichment_uniform:.2f}× |"
        )
    lines += [
        "",
        "The historical frozen H6 table is not overwritten. This conditional table re-ranks all four models on the identical GraphBAN-mappable universe. Unlabelled ChEMBL pairs remain unobserved, not negatives.",
    ]
    md = "\n".join(lines) + "\n"
    Path(str(prefix) + ".md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
