# GraphBAN architecture-fidelity note — 13 September 2026

**Frozen before the S6 outcome was available.**

This note prevents overstatement of the contemporary-model challenge. The official GraphBAN paper and repository describe a broader framework with graph knowledge distillation, teacher/student components, bilinear attention (BAN), and domain-adaptation machinery, with different components used across transductive and inductive modes (Hadipour et al., Nature Communications 16, 2541, 2025; DOI: 10.1038/s41467-025-57536-9; official repository: `HamidHadipour/GraphBAN`).

The primary leakage-free S6 implementation in this project is therefore labeled **GraphBAN-style leakage-free adaptation**, not a verbatim or exact reproduction of the published full GraphBAN framework.

## What S6 preserves

S6 preserves the aspects frozen in `CONTEMPORARY_MODEL_CHALLENGE.md`, `GRAPHBAN_PROTOCOL_AMENDMENT_20260913.md`, and `GRAPHBAN_EXECUTION_PARALLELIZATION_20260913.md`:

- the same frozen TargetDecagon relation universe used by the project;
- mapped PubChem/SMILES and UniProt protein-sequence attributes;
- feature-rich ChemBERTa and ESM-1b node representations;
- a heterogeneous graph neural encoder over the training bipartite graph;
- a learned nonlinear edge decoder;
- no held-out positive edge in message passing;
- identical trained scores evaluated under conventional and prespecified degree-matched controls;
- frozen seeds, optimization, matching rule, and reporting metrics.

## What S6 does not claim

S6 is **not** described as:

- a faithful reproduction of every teacher/student, BAN, knowledge-distillation, or domain-adaptation component in the published GraphBAN framework;
- a reproduction of the published GraphBAN BioSNAP performance number;
- evidence about the original authors' implementation quality;
- a substitute for the published model's inductive/cross-domain experiments.

## Why the adaptation remains informative

The purpose of S6 is a prespecified stress test of whether the benchmark-framing effect survives a substantially richer, contemporary feature-aware graph model on the project's identical frozen relation universe. The exact architecture-fidelity limitation is retained in all reporting.

The S6 result must be reported regardless of direction. A favorable result cannot upgrade the label from “GraphBAN-style adaptation” to “GraphBAN reproduction,” and an unfavorable result cannot trigger replacement or relabeling.
