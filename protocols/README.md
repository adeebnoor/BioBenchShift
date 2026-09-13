# Protocol and provenance record

This directory preserves analysis decisions that were frozen before the corresponding outcomes were inspected, plus timestamped amendments made for technical validity.

## Core records

- `SCIENCE_PROJECT.md` — **historical original prespecified analysis plan**. The filename is retained verbatim for provenance; its journal-target wording reflects the project state at the time it was frozen and should not be read as repository branding or a publication claim.
- `RESULTS_CHECKPOINT_20260913.md` — outcome checkpoint kept separate from the original plan to avoid retrospective rewriting.
- `BENCHMARK_MANIFEST.md` — historical dataset-acquisition and admission rules.
- `CHECKSUMS.sha256` — frozen integrity records.
- `CONTEMPORARY_MODEL_CHALLENGE.md` — prespecified stronger-model robustness gate.
- `GRAPHBAN_AUDIT_PROTOCOL.md` and `GRAPHBAN_PROTOCOL_AMENDMENT_20260913.md` — audit and leakage-free amendment for the contemporary challenge.
- `EXTERNAL_DTI_CHEMBL_PROTOCOL.md` — independent ChEMBL evidence rule frozen before model-specific recovery was inspected.
- `EVIDENCE_STATE_OVERLAP_PROTOCOL.md` — supportive evidence-state sensitivity analysis and balance rules.

## Interpretation rule

The original plan is not rewritten to match later outcomes. Where evidence changed the interpretation—most notably the disease–gene stability boundary—the result is recorded in checkpoint/result files rather than by weakening the original criterion.

Editorial strategy, reviewer suggestions, cover letters, submission-system metadata, and literature-positioning notes are intentionally excluded from the public reproducibility repository.