#!/usr/bin/env python3
"""Science-expansion external DTI validation using BindingDB curated articles.

Protocol: protocols/SCIENCE_EXPANSION_PROTOCOL_20260913.md
Baseline scientific state: commit 1298e4c83ff4fa488a784d52219d20e706118a9e

The primary input must be the BindingDB 2026-09 file containing only records
curated from articles by BindingDB, not the BindingDB subset imported from
ChEMBL. The script never treats absent BindingDB evidence as a negative label.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import run_biosnap_dti_gate2_models as dti
import run_chembl_dti_external as ext

KS = (100, 500, 1000, 5000, 10000, 50000)
NULL_KS = (100, 500, 1000)
N_NULL = 20000
LATER_DATE = pd.Timestamp("2018-01-01")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cid_norm(x) -> str | None:
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return None
    if s.upper().startswith("CID"):
        s = s[3:]
    m = re.search(r"\d+", s)
    if not m:
        return None
    try:
        return str(int(m.group(0)))
    except Exception:
        return None


def parse_multi_cids(x) -> set[str]:
    s = str(x)
    out = set()
    for token in re.findall(r"\d+", s):
        try:
            out.add(str(int(token)))
        except Exception:
            pass
    return out


def find_col(columns, exact=None, regex=None, required=True):
    cols = list(columns)
    if exact:
        for c in cols:
            if c == exact:
                return c
    if regex:
        rx = re.compile(regex, re.I)
        for c in cols:
            if rx.search(c):
                return c
    if required:
        raise KeyError(f"Required BindingDB column not found: exact={exact!r} regex={regex!r}")
    return None


def parse_measurement(cell):
    """Return (relation, value_nM) or None; conservative relation parsing."""
    if cell is None:
        return None
    s = str(cell).strip().replace(",", "")
    if not s or s.lower() in {"nan", "na", "n/a", "none", "null"}:
        return None
    m = re.search(r"(<=|>=|<|>|=|≤|≥|~)?\s*([0-9]*\.?[0-9]+(?:[eE][+-]?\d+)?)", s)
    if not m:
        return None
    rel = m.group(1) or "="
    try:
        val = float(m.group(2))
    except Exception:
        return None
    if not np.isfinite(val) or val <= 0:
        return None
    return rel, val


def qualifies(rel, val, threshold):
    # Values reported as lower bounds (>, >=) do not establish affinity <= threshold.
    if rel in {">", ">=", "≥"}:
        return False
    # Approximate values are retained as quantitative evidence but separately visible.
    return val <= threshold


def load_bindingdb(path: Path, edges, gene_map_path: Path):
    header = pd.read_csv(path, sep="\t", nrows=0, low_memory=False).columns
    c_cid = find_col(header, exact="PubChem CID")
    c_org = find_col(header, exact="Target Source Organism According to Curator or DataSource")
    c_chains = find_col(header, regex=r"^Number of Protein Chains in Target")
    c_sp = find_col(header, regex=r"^UniProt \(SwissProt\) Primary ID of Target Chain 1")
    c_lig = find_col(header, exact="BindingDB Ligand Name", required=False)
    c_tgt = find_col(header, exact="Target Name", required=False)
    c_pub = find_col(header, exact="Date of publication")
    c_pmid = find_col(header, exact="PMID", required=False)
    c_doi = find_col(header, exact="Article DOI", required=False)
    meas_cols = [c for c in ["Ki (nM)", "Kd (nM)", "IC50 (nM)", "EC50 (nM)"] if c in header]
    if not meas_cols:
        raise RuntimeError("No expected quantitative affinity/activity columns found in BindingDB TSV")

    usecols = [c_cid, c_org, c_chains, c_sp, c_pub] + meas_cols
    for c in [c_lig, c_tgt, c_pmid, c_doi]:
        if c:
            usecols.append(c)
    df = pd.read_csv(path, sep="\t", usecols=list(dict.fromkeys(usecols)), dtype=str,
                     low_memory=False, on_bad_lines="skip").fillna("")

    drugs = sorted({a for a, _ in edges})
    cid_to_drug = {}
    for d in drugs:
        c = cid_norm(d)
        if c:
            cid_to_drug[c] = d

    gm = pd.read_csv(gene_map_path, dtype=str).fillna("")
    if not {"gene", "accession"}.issubset(gm.columns):
        raise RuntimeError("gene-map must contain gene and accession columns")
    acc_to_genes = defaultdict(set)
    for _, r in gm.iterrows():
        g = str(r["gene"]).strip()
        a = str(r["accession"]).strip()
        if g and a:
            acc_to_genes[a].add(g)

    blocked = set(edges)
    best1000 = {}
    best100 = {}
    mapped_rows = 0
    eligible_rows = 0

    for _, r in df.iterrows():
        org = str(r[c_org]).strip().lower()
        if "homo sapiens" not in org and org not in {"human", "9606"}:
            continue
        try:
            chains = int(float(str(r[c_chains]).strip()))
        except Exception:
            continue
        if chains != 1:
            continue
        acc = str(r[c_sp]).strip().split()[0] if str(r[c_sp]).strip() else ""
        genes = acc_to_genes.get(acc, set())
        if not genes:
            continue
        drugs_here = {cid_to_drug[c] for c in parse_multi_cids(r[c_cid]) if c in cid_to_drug}
        if not drugs_here:
            continue
        mapped_rows += 1

        vals = []
        for c in meas_cols:
            z = parse_measurement(r[c])
            if z:
                vals.append((c.split()[0], z[0], z[1]))
        if not vals:
            continue
        eligible_rows += 1
        vals.sort(key=lambda x: x[2])
        best = vals[0]
        pubdate = pd.to_datetime(str(r[c_pub]).strip(), errors="coerce")
        if pd.isna(pubdate):
            continue

        record = {
            "drug_name": str(r[c_lig]).strip() if c_lig else "",
            "target_name": str(r[c_tgt]).strip() if c_tgt else "",
            "uniprot": acc,
            "measurement_type": best[0],
            "relation": best[1],
            "affinity_nM": float(best[2]),
            "publication_date": pubdate.strftime("%Y-%m-%d"),
            "pmid": str(r[c_pmid]).strip() if c_pmid else "",
            "doi": str(r[c_doi]).strip() if c_doi else "",
        }
        for d in drugs_here:
            for g in genes:
                pair = (d, g)
                if pair in blocked:
                    continue
                if qualifies(best[1], best[2], 1000.0):
                    old = best1000.get(pair)
                    if old is None or record["affinity_nM"] < old["affinity_nM"]:
                        best1000[pair] = dict(record)
                if qualifies(best[1], best[2], 100.0):
                    old = best100.get(pair)
                    if old is None or record["affinity_nM"] < old["affinity_nM"]:
                        best100[pair] = dict(record)

    meta = {
        "bindingdb_rows": int(len(df)),
        "mapped_human_single_protein_rows": int(mapped_rows),
        "mapped_rows_with_quantitative_measurement": int(eligible_rows),
        "supported_candidate_pairs_le_1000nM": int(len(best1000)),
        "supported_candidate_pairs_le_100nM": int(len(best100)),
    }
    return best1000, best100, meta


def full_rankings(edges, inits=5):
    L, R, scores = ext.fit_ensembles(edges, ninit=inits)
    li = {x: i for i, x in enumerate(L)}
    ri = {x: i for i, x in enumerate(R)}
    mask = np.ones((len(L), len(R)), dtype=bool)
    for a, b in edges:
        mask[li[a], ri[b]] = False
    rankings = {}
    for name, sc in scores.items():
        flat = sc.ravel()
        idx = np.flatnonzero(mask.ravel())
        rankings[name] = idx[np.argsort(flat[idx])[::-1]]
    return L, R, li, ri, mask, rankings


def flat_support(support, li, ri, nr):
    return {li[a] * nr + ri[b] for a, b in support if a in li and b in ri}


def pair_from_flat(x, L, R):
    nr = len(R)
    return L[int(x) // nr], R[int(x) % nr]


def metrics_table(rankings, support_flat, eligible_n):
    rows = []
    total = len(support_flat)
    base = total / eligible_n if eligible_n else float("nan")
    for name, ranked in rankings.items():
        for k in KS:
            kk = min(k, len(ranked))
            hits = sum(int(int(x) in support_flat) for x in ranked[:kk])
            rows.append({
                "model": name,
                "k": kk,
                "hits": hits,
                "precision_lower_bound": hits / kk if kk else float("nan"),
                "recall_supported": hits / total if total else float("nan"),
                "enrichment_uniform": (hits / kk) / base if kk and base > 0 else float("nan"),
            })
    return pd.DataFrame(rows)


def degree_null(edges, L, R, mask, rankings, support_flat):
    dl = Counter(a for a, _ in edges)
    dr = Counter(b for _, b in edges)
    nr = len(R)
    pools_n = Counter()
    pools_m = Counter()
    for i, a in enumerate(L):
        da = dti.dbin(dl.get(a, 0))
        js = np.flatnonzero(mask[i])
        for j in js:
            sig = (da, dti.dbin(dr.get(R[int(j)], 0)))
            flat = i * nr + int(j)
            pools_n[sig] += 1
            if flat in support_flat:
                pools_m[sig] += 1

    rng = np.random.default_rng(20260913)
    rows = []
    for model in [m for m in ["SVD", "NeuralMF"] if m in rankings]:
        ranked = rankings[model]
        for k in NULL_KS:
            selected = ranked[:k]
            sig_n = Counter()
            obs = 0
            for flat in selected:
                a, b = pair_from_flat(int(flat), L, R)
                sig_n[(dti.dbin(dl.get(a, 0)), dti.dbin(dr.get(b, 0)))] += 1
                obs += int(int(flat) in support_flat)
            draws = np.zeros(N_NULL, dtype=int)
            for sig, nsel in sig_n.items():
                N = pools_n[sig]
                M = pools_m[sig]
                if N <= 0 or nsel <= 0:
                    continue
                nsel = min(nsel, N)
                draws += rng.hypergeometric(M, N - M, nsel, size=N_NULL)
            mu = float(draws.mean())
            sd = float(draws.std(ddof=1))
            dev = abs(obs - mu)
            p2 = float((1 + np.sum(np.abs(draws - mu) >= dev)) / (N_NULL + 1))
            rows.append({
                "model": model,
                "k": k,
                "observed_supported_hits": int(obs),
                "degree_matched_null_mean": mu,
                "degree_matched_null_sd": sd,
                "z_vs_degree_matched_null": (obs - mu) / sd if sd > 0 else float("nan"),
                "empirical_two_sided_p": p2,
                "null_draws": N_NULL,
            })
    return pd.DataFrame(rows)


def case_studies(rankings, support, L, R):
    if "SVD" not in rankings or "NeuralMF" not in rankings:
        return pd.DataFrame()
    sv = [int(x) for x in rankings["SVD"][:1000]]
    nm = [int(x) for x in rankings["NeuralMF"][:1000]]
    svset, nmset = set(sv), set(nm)
    rank_sv = {x: i + 1 for i, x in enumerate(sv)}
    rank_nm = {x: i + 1 for i, x in enumerate(nm)}
    rows = []
    for direction, own, other, ownrank, othername in [
        ("SVD_only_top1000", sv, nmset, rank_sv, "NeuralMF"),
        ("NeuralMF_only_top1000", nm, svset, rank_nm, "SVD"),
    ]:
        cand = []
        for flat in own:
            if flat in other:
                continue
            pair = pair_from_flat(flat, L, R)
            rec = support.get(pair)
            if not rec:
                continue
            dt = pd.to_datetime(rec["publication_date"], errors="coerce")
            if pd.isna(dt) or dt < LATER_DATE:
                continue
            item = {
                "direction": direction,
                "drug": pair[0],
                "gene": pair[1],
                "supporting_model_rank": int(ownrank[flat]),
                "competing_model": othername,
                "competing_model_top1000": False,
                **rec,
            }
            cand.append(item)
        cand.sort(key=lambda x: (x["supporting_model_rank"], x["affinity_nM"], x["drug"], x["gene"]))
        rows.extend(cand[:10])
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--biosnap", required=True)
    ap.add_argument("--bindingdb", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--inits", type=int, default=5)
    args = ap.parse_args()

    biosnap = Path(args.biosnap)
    bdb = Path(args.bindingdb)
    gmap = Path(args.gene_map)
    edges = dti.parse(str(biosnap))

    support1000, support100, bmeta = load_bindingdb(bdb, edges, gmap)
    L, R, li, ri, mask, rankings = full_rankings(edges, args.inits)
    sf1000 = flat_support(support1000, li, ri, len(R))
    sf100 = flat_support(support100, li, ri, len(R))

    primary = metrics_table(rankings, sf1000, int(mask.sum()))
    strict = metrics_table(rankings, sf100, int(mask.sum()))
    null = degree_null(edges, L, R, mask, rankings, sf1000)
    cases = case_studies(rankings, support1000, L, R)

    p = Path(args.out_prefix)
    p.parent.mkdir(parents=True, exist_ok=True)
    primary.to_csv(str(p) + "_primary_table.csv", index=False)
    strict.to_csv(str(p) + "_strict100nM_table.csv", index=False)
    null.to_csv(str(p) + "_degree_matched_null.csv", index=False)
    cases.to_csv(str(p) + "_later_case_studies.csv", index=False)

    comp = {}
    for k in KS:
        z = primary[primary.k == min(k, int(mask.sum()))].set_index("model")
        if "SVD" in z.index and "NeuralMF" in z.index:
            comp[str(k)] = {
                "svd_hits": int(z.loc["SVD", "hits"]),
                "neuralmf_hits": int(z.loc["NeuralMF", "hits"]),
                "difference": int(z.loc["SVD", "hits"] - z.loc["NeuralMF", "hits"]),
            }

    summary = {
        "protocol": "SCIENCE_EXPANSION_PROTOCOL_20260913.md",
        "baseline_commit": "1298e4c83ff4fa488a784d52219d20e706118a9e",
        "bindingdb_file": bdb.name,
        "bindingdb_sha256": sha256(bdb),
        "biosnap_sha256": sha256(biosnap),
        "ensemble_initializations": args.inits,
        "candidate_universe": int(mask.sum()),
        "bindingdb": bmeta,
        "confirmatory_SVD_minus_NeuralMF_hits": comp,
        "named_later_case_studies": int(len(cases)),
        "interpretation_boundary": "BindingDB absence is not a negative label; curated biochemical evidence may share literature provenance with other databases despite excluding the BindingDB ChEMBL-import subset.",
    }
    Path(str(p) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# BindingDB curated-article external validation",
        "",
        f"BindingDB input SHA-256: `{summary['bindingdb_sha256']}`",
        f"BioSNAP input SHA-256: `{summary['biosnap_sha256']}`",
        f"Primary externally supported candidate pairs (<=1,000 nM): **{len(sf1000):,}**.",
        f"Strict <=100 nM supported candidate pairs: **{len(sf100):,}**.",
        "",
        "## Primary ranking endpoint",
        "",
        "| Model | K | BindingDB-supported hits | Precision lower bound | Recall supported | Enrichment vs uniform |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, x in primary.iterrows():
        lines.append(f"| {x.model} | {int(x.k):,} | {int(x.hits)} | {x.precision_lower_bound:.5f} | {x.recall_supported:.5f} | {x.enrichment_uniform:.2f}x |")
    lines += ["", "## Degree/popularity-matched null", "", "| Model | K | Observed | Matched-null mean | z | empirical two-sided p |", "|---|---:|---:|---:|---:|---:|"]
    for _, x in null.iterrows():
        lines.append(f"| {x.model} | {int(x.k):,} | {int(x.observed_supported_hits)} | {x.degree_matched_null_mean:.2f} | {x.z_vs_degree_matched_null:.2f} | {x.empirical_two_sided_p:.5f} |")
    lines += ["", f"Mechanically selected later-published discordant case-study rows: **{len(cases)}**.", "", "No absent BindingDB pair is interpreted as a negative interaction."]
    Path(str(p) + ".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
