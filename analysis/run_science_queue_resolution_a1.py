#!/usr/bin/env python3
"""Secondary A1/A2 queue-resolution replay for the Science upgrade.

This analysis replays the frozen five TargetDecagon split seeds and compares
candidate queues from the SAME fitted models whose held-out AUROCs are
recomputed on the original conventional and degree-matched controls.

It is explicitly secondary: the GraphBAN--NeuralMF pair was nominated after
seeing the existing AUROC table. All six model pairs are therefore reported.
No equivalence claim is made from a small AUROC gap.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler

import run_graphban_targetdecagon_clean as gb
import run_dti_fair_model_selection_panel as panel

PROTOCOL = "BioBenchShift-Science-A1-A2-secondary-20260915-v1"
MODELS = ("svd", "neuralmf", "lightgcn", "graphban")
CUTOFFS = (100, 500, 1000)


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_json(obj) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(raw).hexdigest()


def top_ids(scores: np.ndarray, ids: np.ndarray, k: int) -> np.ndarray:
    scores = np.asarray(scores)
    ids = np.asarray(ids)
    cut = np.partition(scores, len(scores) - k)[len(scores) - k]
    above = np.flatnonzero(scores > cut)
    equal = np.flatnonzero(scores == cut)
    if len(above) < k:
        equal = equal[np.argsort(ids[equal], kind="stable")][: k - len(above)]
    chosen = np.r_[above, equal]
    chosen = chosen[np.lexsort((ids[chosen], -scores[chosen].astype(float)))]
    return ids[chosen]


def ht(a: np.ndarray, b: np.ndarray, ids: np.ndarray, k: int) -> float:
    A = top_ids(a, ids, k)
    B = top_ids(b, ids, k)
    return float(1.0 - len(np.intersect1d(A, B)) / k)


def encode_edges(edges, li, ri, nR):
    return np.asarray([li[a] * nR + ri[b] for a, b in edges], dtype=np.int64)


def candidate_universe(L, R, blocked):
    nL, nR = len(L), len(R)
    li = {x: i for i, x in enumerate(L)}
    ri = {x: i for i, x in enumerate(R)}
    blocked_flat = np.fromiter((li[a] * nR + ri[b] for a, b in blocked), dtype=np.int64)
    all_ids = np.arange(nL * nR, dtype=np.int64)
    mask = np.ones(nL * nR, dtype=bool)
    mask[blocked_flat] = False
    ids = all_ids[mask]
    return ids, ids // nR, ids % nR, li, ri


def chunk_dot(A, B, left_idx, right_idx, extra=None, batch=100000):
    out = np.empty(len(left_idx), dtype=np.float64)
    for s in range(0, len(out), batch):
        t = min(len(out), s + batch)
        val = np.einsum("ij,ij->i", A[left_idx[s:t]], B[right_idx[s:t]])
        if extra is not None:
            val = val + extra[0][left_idx[s:t]] + extra[1][right_idx[s:t]]
        out[s:t] = val
    return out


def graphban_flat_scores(model, data, left_idx, right_idx, device, batch=32768):
    model.eval()
    out = np.empty(len(left_idx), dtype=np.float64)
    with torch.no_grad():
        for s in range(0, len(out), batch):
            t = min(len(out), s + batch)
            ix = torch.tensor(
                np.vstack([left_idx[s:t], right_idx[s:t]]),
                dtype=torch.long,
                device=device,
            )
            out[s:t] = torch.sigmoid(model(data.x_dict, data.edge_index_dict, ix)).cpu().numpy()
    return out


def binary_metrics(pos_scores, neg_scores):
    y = np.r_[np.ones(len(pos_scores)), np.zeros(len(neg_scores))]
    s = np.r_[pos_scores, neg_scores]
    return float(roc_auc_score(y, s)), float(average_precision_score(y, s))


def fit_seed(E_all, E_mapped, L, R, baseXd, baseXp, seed, device):
    rng = np.random.default_rng(seed)
    ix = rng.permutation(len(E_all))
    nt = round(0.2 * len(E_all))
    test_all = {E_all[i] for i in ix[:nt]}
    train_all = {E_all[i] for i in ix[nt:]}
    mapped = set(E_mapped)
    train = sorted(train_all & mapped)
    test = sorted(test_all & mapped)
    blocked = set(E_mapped)

    dl = Counter(a for a, b in train)
    dr = Counter(b for a, b in train)
    seen = [e for e in test if dl[e[0]] > 0 and dr[e[1]] > 0]
    erng = np.random.default_rng(seed + 40000)
    rn = gb.sample_random(erng, len(seen), L, R, blocked)
    mp, mn = gb.degree_matched(erng, seen, L, R, blocked, dl, dr)

    # Baselines: unchanged frozen panel fits.
    svd = panel.fit_svd(train, L, R, 32, seed)
    nmf = panel.fit_logistic_mf(train, L, R, blocked, np.random.default_rng(seed + 1000))
    lg = panel.fit_lightgcn(train, L, R, blocked, seed)

    # GraphBAN: unchanged frozen S6 feature scaling and fit.
    di = {x: i for i, x in enumerate(L)}
    gi = {x: i for i, x in enumerate(R)}
    td = {a for a, b in train}
    tg = {b for a, b in train}
    dsc = StandardScaler().fit(baseXd[[di[x] for x in sorted(td)]])
    psc = StandardScaler().fit(baseXp[[gi[x] for x in sorted(tg)]])
    Xd = dsc.transform(baseXd).astype("float32")
    Xp = psc.transform(baseXp).astype("float32")
    data, gli, gri = gb.make_graph(L, R, train, Xd, Xp, device)
    gmodel = gb.train_model(data, gli, gri, train, blocked, L, R, seed, device, 20)

    cand_ids, cL, cR, li, ri = candidate_universe(L, R, blocked)
    if len(cand_ids) != 1005757:
        raise RuntimeError(f"candidate universe drift: {len(cand_ids)} != 1005757")

    # Full candidate scores, one fitted model per split.
    _, _, Us, Vs = svd
    svd_scores = chunk_dot(Us, Vs, cL, cR)
    _, _, Un, Vn, bun, bvn = nmf
    nmf_scores = chunk_dot(Un, Vn, cL, cR, (bun, bvn))
    _, _, z = lg
    lg_scores = chunk_dot(z, z, cL, cR + len(L))
    gb_scores = graphban_flat_scores(gmodel, data, cL, cR, device)
    scores = {
        "svd": svd_scores,
        "neuralmf": nmf_scores,
        "lightgcn": lg_scores,
        "graphban": gb_scores,
    }

    # Recompute held-out evaluation metrics from these exact same fits.
    eval_rows = []
    for regime, pos, neg in (("random", seen, rn), ("matched", mp, mn)):
        for name in MODELS:
            if name == "svd":
                ps, ns = panel.score_svd(pos, svd), panel.score_svd(neg, svd)
            elif name == "neuralmf":
                ps, ns = panel.score_mf(pos, nmf), panel.score_mf(neg, nmf)
            elif name == "lightgcn":
                ps, ns = panel.score_lg(pos, lg), panel.score_lg(neg, lg)
            else:
                ps = gb.model_scores(gmodel, data, pos, gli, gri, device)
                ns = gb.model_scores(gmodel, data, neg, gli, gri, device)
            auc, auprc = binary_metrics(ps, ns)
            eval_rows.append(
                {
                    "seed": seed,
                    "regime": regime,
                    "model": name,
                    "auroc": auc,
                    "auprc": auprc,
                    "n_pos": len(pos),
                    "n_neg": len(neg),
                }
            )

    training_hash = sha256_json(train)
    universe_hash = hashlib.sha256(cand_ids.tobytes()).hexdigest()
    split_meta = {
        "seed": seed,
        "n_train": len(train),
        "n_test_mapped": len(test),
        "n_seen": len(seen),
        "n_matched": len(mp),
        "match_fraction": len(mp) / len(seen),
        "training_hash": training_hash,
        "candidate_universe_hash": universe_hash,
        "n_candidates": len(cand_ids),
    }
    return cand_ids, scores, eval_rows, split_meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--drug-map", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--drug-features", required=True)
    ap.add_argument("--protein-features", required=True)
    ap.add_argument("--expected-panel", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seeds", default="0,1,2,3,4")
    args = ap.parse_args()

    E = gb.parse_edges(args.input)
    dd = pd.read_csv(args.drug_map, dtype=str).fillna("")
    gd = pd.read_csv(args.gene_map, dtype=str).fillna("")
    dd = dd[dd.smiles.astype(bool)].drop_duplicates("drug")
    gd = gd[gd.sequence.astype(bool)].drop_duplicates("gene")
    dm = set(dd.drug)
    gm = set(gd.gene)
    mapped = [e for e in E if e[0] in dm and e[1] in gm]
    L = sorted({a for a, b in mapped})
    R = sorted({b for a, b in mapped})
    dd = dd.set_index("drug").loc[L].reset_index()
    gd = gd.set_index("gene").loc[R].reset_index()
    baseXd = np.load(args.drug_features).astype("float32")
    baseXp = np.load(args.protein_features).astype("float32")
    if baseXd.shape[0] != len(L) or baseXp.shape[0] != len(R):
        raise RuntimeError(
            f"feature/order mismatch: drug {baseXd.shape[0]}/{len(L)}, protein {baseXp.shape[0]}/{len(R)}"
        )

    seeds = [int(x) for x in args.seeds.split(",") if x.strip()]
    expected = pd.read_csv(args.expected_panel)
    device = torch.device("cpu")
    torch.set_num_threads(4)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    eval_rows, pair_rows, split_rows, top_rows = [], [], [], []
    candidate_hash = None

    for seed in seeds:
        print(f"A1 seed {seed}...", flush=True)
        cand_ids, scores, erows, smeta = fit_seed(E, mapped, L, R, baseXd, baseXp, seed, device)
        eval_rows.extend(erows)
        split_rows.append(smeta)
        uh = smeta["candidate_universe_hash"]
        if candidate_hash is None:
            candidate_hash = uh
        elif candidate_hash != uh:
            raise RuntimeError("candidate universe differs across split seeds")

        # Save top-1000 identities for audit without committing million-score arrays.
        for model in MODELS:
            ids1000 = top_ids(scores[model], cand_ids, 1000)
            rank = {int(x): i + 1 for i, x in enumerate(ids1000.tolist())}
            s_lookup = dict(zip(cand_ids.tolist(), scores[model].tolist()))
            for cid in ids1000:
                top_rows.append(
                    {
                        "seed": seed,
                        "model": model,
                        "candidate_id": int(cid),
                        "left_index": int(cid // len(R)),
                        "right_index": int(cid % len(R)),
                        "rank": rank[int(cid)],
                        "score": float(s_lookup[int(cid)]),
                    }
                )

        emap = {(r["model"], r["regime"]): r for r in erows}
        for i, a in enumerate(MODELS):
            for b in MODELS[i + 1 :]:
                rho = float(spearmanr(scores[a], scores[b]).statistic)
                for regime in ("random", "matched"):
                    da = emap[(a, regime)]["auroc"]
                    db = emap[(b, regime)]["auroc"]
                    row = {
                        "seed": seed,
                        "regime": regime,
                        "model_a": a,
                        "model_b": b,
                        "delta_auroc_signed_a_minus_b": float(da - db),
                        "delta_auroc_abs": float(abs(da - db)),
                        "candidate_score_spearman": rho,
                    }
                    for k in CUTOFFS:
                        row[f"ht_{k}"] = ht(scores[a], scores[b], cand_ids, k)
                    pair_rows.append(row)

        del scores

    eval_df = pd.DataFrame(eval_rows)
    pair_df = pd.DataFrame(pair_rows)
    split_df = pd.DataFrame(split_rows)
    top_df = pd.DataFrame(top_rows)

    # Hard fidelity gate against frozen five-seed table S17b.
    mismatches = []
    for _, r in expected.iterrows():
        seed = int(r.seed)
        for model in MODELS:
            for regime in ("random", "matched"):
                col = f"{model}_auc_{regime}"
                got = float(eval_df.query("seed == @seed and model == @model and regime == @regime").iloc[0].auroc)
                exp = float(r[col])
                diff = abs(got - exp)
                if diff > 5e-7:
                    mismatches.append({"seed": seed, "model": model, "regime": regime, "expected": exp, "replay": got, "abs_diff": diff})
    if mismatches:
        Path(str(prefix) + "_FIDELITY_FAILURE.json").write_text(json.dumps(mismatches, indent=2) + "\n")
        raise RuntimeError(f"replay fidelity gate failed: {mismatches[:3]}")

    eval_df.to_csv(str(prefix) + "_evaluation_metrics.csv", index=False)
    pair_df.to_csv(str(prefix) + "_pairwise_resolution.csv", index=False)
    split_df.to_csv(str(prefix) + "_split_provenance.csv", index=False)
    top_df.to_csv(str(prefix) + "_top1000_identities.csv.gz", index=False, compression="gzip")

    nom = pair_df[(pair_df.model_a == "neuralmf") & (pair_df.model_b == "graphban")].copy()
    if nom.empty:
        nom = pair_df[(pair_df.model_a == "graphban") & (pair_df.model_b == "neuralmf")].copy()
    summary = {
        "protocol": PROTOCOL,
        "secondary_analysis": True,
        "nominated_pair": ["graphban", "neuralmf"],
        "n_seeds": len(seeds),
        "candidate_universe_size": int(split_df.n_candidates.iloc[0]),
        "candidate_universe_hash": candidate_hash,
        "mapping_fraction": len(mapped) / len(E),
        "graphban_neuralmf": {
            regime: {
                "mean_abs_delta_auroc": float(nom[nom.regime == regime].delta_auroc_abs.mean()),
                "median_abs_delta_auroc": float(nom[nom.regime == regime].delta_auroc_abs.median()),
                "mean_ht100": float(nom[nom.regime == regime].ht_100.mean()),
                "range_ht100": [float(nom[nom.regime == regime].ht_100.min()), float(nom[nom.regime == regime].ht_100.max())],
                "mean_ht500": float(nom[nom.regime == regime].ht_500.mean()),
                "mean_ht1000": float(nom[nom.regime == regime].ht_1000.mean()),
            }
            for regime in ("random", "matched")
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "torch": torch.__version__,
            "device": str(device),
            "torch_threads": torch.get_num_threads(),
        },
        "input_hashes": {
            "targetdecagon": sha256_file(args.input),
            "drug_map": sha256_file(args.drug_map),
            "gene_map": sha256_file(args.gene_map),
            "drug_features": sha256_file(args.drug_features),
            "protein_features": sha256_file(args.protein_features),
            "expected_panel": sha256_file(args.expected_panel),
        },
        "interpretation": "Descriptive same-fit queue comparison. A small AUROC gap is not an equivalence test; HT quantifies queue identity difference only.",
    }
    Path(str(prefix) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # A2 descriptive plot: absolute paired AUROC difference versus HT@100.
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8.2, 5.6))
    ax = fig.add_subplot(111)
    for _, r in pair_df.iterrows():
        marker = "o" if r.regime == "random" else "x"
        ax.scatter(r.delta_auroc_abs, r.ht_100, marker=marker, alpha=0.75)
    ax.set_xlabel("Absolute paired AUROC difference")
    ax.set_ylabel("Hypothesis turnover at top 100")
    ax.set_title("Same-fit model-score separation and queue turnover")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(str(prefix) + "_A2_resolution.png", dpi=220)
    plt.close(fig)

    lines = [
        "# Science A1/A2 queue-resolution replay",
        "",
        f"Protocol: `{PROTOCOL}`",
        f"Five split seeds; common candidate universe: **{summary['candidate_universe_size']:,}** mapped unknown drug-target pairs.",
        "Replay AUROCs passed the frozen five-seed fidelity gate for all four models under both evaluation regimes.",
        "",
        "## Nominated GraphBAN-NeuralMF comparison",
        "",
        "| Regime | Mean |ΔAUROC| | Mean HT@100 | Mean HT@500 | Mean HT@1000 |",
        "|---|---:|---:|---:|---:|",
    ]
    for regime in ("random", "matched"):
        x = summary["graphban_neuralmf"][regime]
        lines.append(
            f"| {regime} | {x['mean_abs_delta_auroc']:.6f} | {x['mean_ht100']:.3f} | {x['mean_ht500']:.3f} | {x['mean_ht1000']:.3f} |"
        )
    lines += [
        "",
        "All six model pairs and all five seeds are retained in the CSV. These are descriptive results, not a statistical-equivalence test.",
    ]
    Path(str(prefix) + ".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
