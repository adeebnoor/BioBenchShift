#!/usr/bin/env python3
"""Combine deterministic GraphBAN ESM-1b feature chunks into canonical order."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunks", nargs="+", required=True)
    ap.add_argument("--n-rows", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = None
    seen = np.zeros(args.n_rows, dtype=bool)
    for name in sorted(args.chunks):
        z = np.load(name)
        idx = z["indices"].astype(np.int64)
        feat = z["features"].astype("float32")
        if len(idx) != feat.shape[0]:
            raise RuntimeError(f"Index/feature mismatch in {name}")
        if rows is None:
            rows = np.empty((args.n_rows, feat.shape[1]), dtype=np.float32)
        elif rows.shape[1] != feat.shape[1]:
            raise RuntimeError("Feature dimension differs across chunks")
        if np.any(idx < 0) or np.any(idx >= args.n_rows):
            raise RuntimeError(f"Out-of-range canonical index in {name}")
        if np.any(seen[idx]):
            raise RuntimeError(f"Duplicate canonical index in {name}")
        rows[idx] = feat
        seen[idx] = True

    if rows is None or not seen.all():
        missing = np.where(~seen)[0]
        raise RuntimeError(f"Missing canonical protein features: {missing[:20].tolist()}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.save(out, rows)
    print(f"combined features: {rows.shape}; every canonical row present exactly once")


if __name__ == "__main__":
    main()
