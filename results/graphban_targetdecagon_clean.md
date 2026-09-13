# Leakage-free contemporary GraphBAN challenge — TargetDecagon DTI

Mapping edge coverage: **99.7%**. Held-out positives in message passing: **no**.
Feature execution: deterministic parallel precomputation using the frozen ChemBERTa and ESM-1b definitions; downstream S6 scientific protocol unchanged.

| Metric | Conventional | Structure-neutralized | Difference |
|---|---:|---:|---:|
| AUROC | 0.997 ± 0.000 | 0.928 ± 0.014 | +0.069 conventional-minus-neutralized |
| AUPRC | 0.998 | 0.943 | +0.055 |

Mean matching coverage: **0.588** across seeds [0, 1, 2].

This is the primary S6 contemporary-model challenge defined before outcome inspection. The public GraphBAN transductive test-graph construction is not used here because held-out positive edges must not enter message passing. By the frozen S6 rule, either a sensitive or comparatively robust GraphBAN result closes the contemporary-model gate if the execution is valid and the result is incorporated without cherry-picking.
