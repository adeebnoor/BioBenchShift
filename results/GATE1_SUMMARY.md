# Gate 1 — cross-domain structural-inflation summary

## Decision

**PASS.** The prespecified project gate required reproducible structural inflation in at least 3 of 4 independent external relation families. All 4 tested families passed.

| Relation family | Benchmark | Edges | Degree-only AUC: random negatives | Degree-only AUC: degree-matched | Mean AUC inflation | Matching coverage | Gate |
|---|---|---:|---:|---:|---:|---:|---|
| Drug–target | BioSNAP TargetDecagon | 18,690 | 0.983 ± 0.001 | 0.620 ± 0.012 | 0.364 | 0.599 | PASS |
| Protein–protein | HuRI | 52,068 unique non-self edges | 0.928 ± 0.002 | 0.513 ± 0.001 | 0.415 | 1.000 | PASS |
| Drug–disease treatment | Hetionet CtD | 755 | 0.871 ± 0.019 | 0.519 ± 0.010 | 0.352 | 1.000 | PASS |
| Disease–gene association | Hetionet DaG | 12,623 | 0.876 ± 0.006 | 0.513 ± 0.002 | 0.363 | 1.000 | PASS |

The Anti-DDI seed analysis independently motivated the project: a popularity-only score at ATC5 class level achieved approximately 0.90 AUC under conventional comparison and approximately 0.50 after degree matching. It is retained as the evidence-aware DDI case, not counted among the four external Gate-1 families above.

## What Gate 1 establishes

Across four distinct biomedical relation types, a score using only endpoint observation frequency / degree can achieve apparently strong discrimination under conventional random-negative evaluation. Matching test negatives to the endpoint-degree structure removes most of that advantage. This establishes a reproducible **benchmark structural-inflation phenomenon** across domains.

## What Gate 1 does NOT establish

Gate 1 does not show that modern biomedical AI models are useless, that all reported performance is caused by degree, or that degree matching is a complete evaluation solution. It is a model-free diagnostic. The project proceeds to learned-model sensitivity, model-ranking stability, evidence-aware negative states, and independent / temporal validation.

## Project implication

The Science-level question is therefore no longer whether degree bias exists. Prior literature already establishes versions of that problem. The remaining high-value question is whether conventional biomedical AI evaluation systematically confounds **structural observability with biological generalization**, and whether a combined structure-neutralized + evidence-aware evaluation better predicts performance on genuinely new biology.
