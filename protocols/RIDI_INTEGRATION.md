# RIDI-inspired layer for the Science structural-shortcuts project

## Purpose

RIDI is used here as a **measurement grammar for scientific consequence**, not as a second RIDI/Nature paper.

The structural-shortcuts project asks whether benchmark construction rewards structural opportunity rather than pair-specific biology. The RIDI-inspired extension asks a downstream question:

> **When benchmark construction changes which model appears best, does it also change which biological hypotheses are prioritized for discovery?**

This layer turns model-rank reversal into a decision-identity result on a fixed biological candidate universe.

## Hard non-overlap boundary with the RIDI/Nature manuscript

The Science project must **not** reuse as primary evidence:

- the RIDI RAG experiment,
- EPSS or COMPAS experiments,
- RTX-KG2 representation/canonicalization contrasts from the RIDI manuscript,
- ReVerb45K open-domain replication,
- the RIDI allocation-identity theorem as the principal novelty,
- the same identification → consequence → deployment → remedy narrative.

Those remain part of the separate RIDI/Nature program.

This paper instead introduces a biomedical-benchmark-specific application: **evaluation design can alter model selection and thereby alter the identity of prioritized biological hypotheses.**

## New object: hypothesis-selection identity

For each biomedical relation family, freeze a candidate universe U of biologically admissible, currently unconfirmed pairs.

Let:

- B_rand = conventional/random-negative benchmark,
- B_neutral = structure-neutralized benchmark,
- M_rand = model selected by B_rand,
- M_neutral = model selected by B_neutral,
- TopK(M, U) = the k highest-scoring candidate relations in the same frozen universe U.

Define downstream hypothesis turnover:

`HT@k = 1 - |TopK(M_rand,U) ∩ TopK(M_neutral,U)| / k`

and overlap-adjusted rank disagreement using Jaccard, Spearman/Kendall on the shared candidate universe, and rank-biased overlap where appropriate.

The key causal discipline is that U is identical for both selected models. The only upstream change is the evaluation rule that determines which model wins.

## H5 — benchmark-induced scientific decision instability

**Pre-specified hypothesis:** In at least two biomedical relation families where conventional and structure-neutralized evaluation select different winning models, the resulting winners will prioritize materially different top-k biological hypotheses on the same frozen candidate universe.

Primary descriptive thresholds for the pilot:

- report HT@100, HT@500, HT@1000,
- report bootstrap confidence intervals over candidate resampling,
- report exact overlap counts and rank correlation,
- do not declare a universal cutoff for "material" instability until the temporal validation is examined.

## Strongest external-validity test

Temporal validation converts H5 from descriptive instability into scientific consequence.

For a historical cutoff t0:

1. train all models using only information available at or before t0;
2. construct B_rand(t0) and B_neutral(t0);
3. select M_rand and M_neutral without seeing post-t0 data;
4. rank the same frozen candidate universe U(t0);
5. identify relations first supported after t0 in an independent later snapshot;
6. compare prospective discovery yield of the two selected models.

Primary quantities:

- Precision@k for post-t0 discoveries,
- Recall@k where the future-positive universe is defined,
- enrichment over random eligible candidates,
- time-to-discovery distribution if dates are sufficiently resolved,
- paired bootstrap or permutation tests on the same future-positive set.

The decisive Science-level result would be:

> conventional evaluation selects a model with superior headline benchmark performance, but structure-neutralized evaluation selects a different model whose prioritized hypotheses are validated more often by later independent evidence.

## RIDI-style audit sufficiency test

A secondary experiment will search for model pairs/seeds with nearly identical conventional AUROC/AUPRC but substantially different:

- structure-neutralized performance,
- hypothesis-selection identity,
- temporal discovery yield.

This is not presented as the RIDI theorem. It is a biomedical benchmark sufficiency test: **headline predictive metrics may be insufficient to identify the scientific decisions induced by a model.**

## Figure integration

### Figure 4 — model selection changes
Show random-vs-neutralized model rank reversals across relation families.

### Figure 5 — hypothesis-selection identity
For each relation family with a rank reversal, show the top-k overlap/turnover between the model selected conventionally and the model selected after neutralization on the same frozen candidate universe.

### Figure 6 — future discovery
Historical freeze → model selection under each benchmark → fixed candidate ranking → later independent discoveries. Compare prospective discovery yield.

### Figure 7 — proposed reporting standard
A biomedical relation-prediction claim should report:

1. structural null,
2. structure-neutralized evaluation,
3. model-rank stability,
4. hypothesis-selection stability,
5. temporal/independent validation.

## Interpretation boundary

Hypothesis turnover is not automatically error. Different models may encode genuinely different biological hypotheses. The scientific claim requires the external-validity layer to determine which evaluation regime better selects models that generalize to later evidence.

The paper therefore separates:

- **structural inflation** — a benchmark property,
- **model sensitivity** — a model-by-benchmark property,
- **model-rank reversal** — a model-selection property,
- **hypothesis-selection turnover** — a scientific-decision property,
- **future discovery yield** — an external-validity property.

This decomposition is the intended RIDI contribution to the Science project.