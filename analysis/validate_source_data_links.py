#!/usr/bin/env python3
"""Validate that every manuscript source-data row points to an existing frozen result.

This is an integrity check only; it does not recompute scientific results.
"""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "figures" / "source_data_main.csv"

with SOURCE.open(newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

if not rows:
    raise SystemExit("figures/source_data_main.csv contains no data rows")

missing = []
empty = []
for i, row in enumerate(rows, start=2):
    ref = (row.get("source_result") or "").strip()
    if not ref:
        empty.append(i)
        continue
    p = ROOT / ref
    if not p.is_file():
        missing.append((i, ref))

if empty or missing:
    if empty:
        print("Rows with empty source_result:", empty)
    for line, ref in missing:
        print(f"Missing source_result at CSV line {line}: {ref}")
    raise SystemExit(1)

print(f"PASS: {len(rows)} source-data rows reference existing frozen result files.")
