#!/usr/bin/env python3
"""Deterministic feature extraction for the frozen leakage-free GraphBAN S6.

This script changes execution only. Feature definitions match
run_graphban_targetdecagon_clean.py. Protein chunks are indexed against the
canonical sorted mapped-protein endpoint list and can be computed independently.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from run_graphban_targetdecagon_clean import parse_edges


def canonical_endpoints(input_path: str, drug_map: str, gene_map: str):
    edges = parse_edges(input_path)
    dd = pd.read_csv(drug_map, dtype=str).fillna("")
    gd = pd.read_csv(gene_map, dtype=str).fillna("")
    dd = dd[dd.smiles.astype(bool)].drop_duplicates("drug")
    gd = gd[gd.sequence.astype(bool)].drop_duplicates("gene")
    dm = dict(zip(dd.drug, dd.smiles))
    gm = dict(zip(gd.gene, gd.sequence))
    mapped = [e for e in edges if e[0] in dm and e[1] in gm]
    if len(mapped) / len(edges) < 0.90:
        raise RuntimeError("Frozen GraphBAN mapping gate failed during feature extraction")
    drugs = sorted({a for a, _ in mapped})
    genes = sorted({b for _, b in mapped})
    dd = dd.set_index("drug").loc[drugs].reset_index()
    gd = gd.set_index("gene").loc[genes].reset_index()
    return drugs, genes, dd, gd


def extract_drugs(drug_df: pd.DataFrame, device: torch.device):
    from transformers import AutoTokenizer, RobertaModel

    tok = AutoTokenizer.from_pretrained("DeepChem/ChemBERTa-77M-MTR")
    model = RobertaModel.from_pretrained(
        "DeepChem/ChemBERTa-77M-MTR", add_pooling_layer=True
    ).eval().to(device)
    smiles = drug_df.smiles.astype(str).tolist()
    chunks = []
    for start in range(0, len(smiles), 16):
        batch = smiles[start : start + 16]
        enc = tok(
            batch,
            return_tensors="pt",
            padding="max_length",
            max_length=290,
            truncation=True,
        ).to(device)
        with torch.no_grad():
            feat = model(**enc).last_hidden_state[:, 0, :]
        chunks.append(feat.cpu().numpy().astype("float32"))
    return np.vstack(chunks)


def extract_protein_chunk(
    gene_df: pd.DataFrame,
    indices: np.ndarray,
    device: torch.device,
):
    import esm

    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    model = model.eval().to(device)
    converter = alphabet.get_batch_converter()
    seqs = gene_df.sequence.fillna("").astype(str).tolist()
    features = []

    # Same representation rule as the frozen serial implementation.
    # Work is sorted by sequence length only within this chunk to reduce padding;
    # canonical indices are retained and restored during combination.
    work = [(int(i), seqs[int(i)][:1022]) for i in indices]
    work.sort(key=lambda x: len(x[1]))
    ordered_indices = []
    for start in range(0, len(work), 4):
        batch_work = work[start : start + 4]
        data = [(f"p{i}", seq) for i, seq in batch_work]
        _, _, tokens = converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            rep = model(tokens, repr_layers=[33], return_contacts=False)[
                "representations"
            ][33]
        for j, (idx, seq) in enumerate(batch_work):
            features.append(
                rep[j, 1 : len(seq) + 1].mean(0).cpu().numpy().astype("float32")
            )
            ordered_indices.append(idx)
    return np.asarray(ordered_indices, dtype=np.int64), np.vstack(features)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--drug-map", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--kind", choices=["drug", "protein"], required=True)
    ap.add_argument("--chunk-index", type=int, default=0)
    ap.add_argument("--num-chunks", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    torch.set_num_threads(max(1, min(4, torch.get_num_threads())))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    drugs, genes, dd, gd = canonical_endpoints(
        args.input, args.drug_map, args.gene_map
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.kind == "drug":
        feat = extract_drugs(dd, device)
        if feat.shape[0] != len(drugs):
            raise RuntimeError("Drug feature row count does not match canonical endpoint order")
        np.save(out, feat)
        Path(str(out) + ".order.txt").write_text("\n".join(drugs) + "\n")
        print(f"drug features: {feat.shape}")
        return

    if args.num_chunks < 1 or not (0 <= args.chunk_index < args.num_chunks):
        raise ValueError("Invalid chunk specification")
    idx = np.arange(len(genes), dtype=np.int64)
    idx = idx[idx % args.num_chunks == args.chunk_index]
    ordered_idx, feat = extract_protein_chunk(gd, idx, device)
    if len(ordered_idx) != len(idx):
        raise RuntimeError("Protein chunk did not return all assigned endpoints")
    np.savez_compressed(out, indices=ordered_idx, features=feat)
    print(
        f"protein chunk {args.chunk_index}/{args.num_chunks}: "
        f"{len(ordered_idx)} endpoints, features {feat.shape}"
    )


if __name__ == "__main__":
    main()
