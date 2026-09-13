#!/usr/bin/env python3
"""Expansion B: functional consequences of benchmark-selected DTI queues.

Protocol: BIOLOGICAL_CONSEQUENCE_AMENDMENT_20260913.md
No pathway term is manually selected or suppressed after outcome inspection.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.spatial.distance import jensenshannon

import run_biosnap_dti_gate2_models as dti
import run_bindingdb_external as bdb

KS = (100, 500, 1000)
MODELS = ("SVD", "NeuralMF")
GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


def target_sets(edges, inits=5):
    L, R, li, ri, mask, rankings = bdb.full_rankings(edges, inits)
    out = {}
    for model in MODELS:
        ranked = rankings[model]
        for k in KS:
            genes = []
            for x in ranked[:k]:
                _, g = bdb.pair_from_flat(int(x), L, R)
                genes.append(str(g))
            out[(model, k)] = sorted(set(genes), key=lambda x: (len(x), x))
    return R, out


def gost(query, background, tries=6):
    payload = {
        "organism": "hsapiens",
        "query": [str(x) for x in query],
        "sources": ["REAC", "GO:BP"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": [str(x) for x in background],
        "numeric_ns": "ENTREZGENE_ACC",
        "ordered": False,
        "no_evidences": True,
        "output": "json",
    }
    err = None
    for i in range(tries):
        try:
            r = requests.post(GPROFILER, json=payload, timeout=120,
                              headers={"User-Agent": "BioBenchShift-Science-expansion/1.0"})
            r.raise_for_status()
            j = r.json()
            if "result" not in j:
                raise RuntimeError(f"Unexpected g:Profiler response keys: {list(j)}")
            return j
        except Exception as e:
            err = e
            time.sleep(min(2 ** i, 20))
    raise RuntimeError(f"g:Profiler request failed: {err}")


def result_frame(j, model, k):
    rows = []
    for r in j.get("result", []):
        if r.get("source") not in {"REAC", "GO:BP"}:
            continue
        p = float(r.get("p_value", 1.0))
        significant = bool(r.get("significant", p < 0.05)) and p < 0.05
        if not significant:
            continue
        rows.append({
            "model": model,
            "k": k,
            "source": r.get("source"),
            "term_id": r.get("native"),
            "term_name": r.get("name"),
            "adjusted_p": p,
            "query_size": r.get("query_size"),
            "term_size": r.get("term_size"),
            "intersection_size": r.get("intersection_size"),
            "effective_domain_size": r.get("effective_domain_size"),
        })
    return pd.DataFrame(rows)


def js_from_terms(a: pd.DataFrame, b: pd.DataFrame):
    pa = {str(r.term_id): -math.log10(max(float(r.adjusted_p), 1e-300)) for _, r in a.iterrows()}
    pb = {str(r.term_id): -math.log10(max(float(r.adjusted_p), 1e-300)) for _, r in b.iterrows()}
    terms = sorted(set(pa) | set(pb))
    if not terms:
        return None
    x = np.asarray([pa.get(t, 0.0) for t in terms], float)
    y = np.asarray([pb.get(t, 0.0) for t in terms], float)
    if x.sum() == 0 and y.sum() == 0:
        return None
    if x.sum() == 0:
        x = np.ones_like(x) / len(x)
    else:
        x /= x.sum()
    if y.sum() == 0:
        y = np.ones_like(y) / len(y)
    else:
        y /= y.sum()
    return float(jensenshannon(x, y, base=2.0) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--biosnap", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--inits", type=int, default=5)
    a = ap.parse_args()

    edges = dti.parse(a.biosnap)
    background, sets = target_sets(edges, a.inits)
    p = Path(a.out_prefix)
    p.parent.mkdir(parents=True, exist_ok=True)

    all_enrich = []
    api_meta = {}
    for model in MODELS:
        for k in KS:
            genes = sets[(model, k)]
            j = gost(genes, background)
            api_meta[f"{model}_{k}"] = j.get("meta", {})
            df = result_frame(j, model, k)
            if not df.empty:
                all_enrich.append(df)
    enrich = pd.concat(all_enrich, ignore_index=True) if all_enrich else pd.DataFrame(
        columns=["model", "k", "source", "term_id", "term_name", "adjusted_p", "query_size", "term_size", "intersection_size", "effective_domain_size"]
    )
    enrich.to_csv(str(p) + "_all_significant_terms.csv", index=False)

    target_rows = []
    for model in MODELS:
        for k in KS:
            for g in sets[(model, k)]:
                target_rows.append({"model": model, "k": k, "gene": g})
    pd.DataFrame(target_rows).to_csv(str(p) + "_target_sets.csv", index=False)

    summary_rows = []
    for k in KS:
        ga, gb = set(sets[("SVD", k)]), set(sets[("NeuralMF", k)])
        ta = enrich[(enrich.model == "SVD") & (enrich.k == k)]
        tb = enrich[(enrich.model == "NeuralMF") & (enrich.k == k)]
        A, B = set(ta.term_id.astype(str)), set(tb.term_id.astype(str))
        min_n = min(len(A), len(B))
        summary_rows.append({
            "k": k,
            "svd_unique_targets": len(ga),
            "neuralmf_unique_targets": len(gb),
            "target_jaccard": len(ga & gb) / len(ga | gb) if ga | gb else float("nan"),
            "svd_significant_terms": len(A),
            "neuralmf_significant_terms": len(B),
            "shared_significant_terms": len(A & B),
            "term_overlap_coefficient": len(A & B) / min_n if min_n else (1.0 if not A and not B else 0.0),
            "svd_unique_significant_terms": len(A - B),
            "neuralmf_unique_significant_terms": len(B - A),
            "js_divergence_enrichment_profile": js_from_terms(ta, tb),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(str(p) + "_summary.csv", index=False)

    meta = {
        "protocol": "BIOLOGICAL_CONSEQUENCE_AMENDMENT_20260913.md",
        "baseline_commit": "1298e4c83ff4fa488a784d52219d20e706118a9e",
        "models": list(MODELS),
        "k": list(KS),
        "background_targets": len(background),
        "gprofiler_endpoint": GPROFILER,
        "gprofiler_request": {
            "organism": "hsapiens",
            "sources": ["REAC", "GO:BP"],
            "significance_threshold_method": "fdr",
            "user_threshold": 0.05,
            "domain_scope": "custom",
            "numeric_ns": "ENTREZGENE_ACC",
            "ordered": False,
        },
        "api_metadata": api_meta,
        "interpretation_boundary": "Enrichment differences indicate different biological programs emphasized by the queues; they do not validate individual predicted pairs or establish causal correctness.",
    }
    Path(str(p) + "_metadata.json").write_text(json.dumps(meta, indent=2, default=str) + "\n")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
