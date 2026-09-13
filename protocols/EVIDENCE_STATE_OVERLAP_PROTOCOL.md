# Frozen evidence-state overlap protocol — Anti-DDI sensitivity analysis

**Frozen:** 2026-09-13, before inspecting overlap-weighted learned-model outcomes.  
**Role in Science project:** extended/sensitivity evidence separating negative evidence state from structural exposure.  
**Claim boundary:** this analysis does not establish clinical non-interaction and does not treat persistent database absence as a true negative.

## Scientific question

The earlier Anti-DDI pilot showed that learned-model discrimination differed when randomly unobserved DDI pairs were replaced by curated T1/T2 counter-evidence pairs, but exact endpoint-degree matching covered only a minority of held-out positives and forced nearest matching left residual structural imbalance.

The confirmatory sensitivity question is therefore:

> When curated counter-evidence pairs and randomly unobserved pairs are standardized to their common structural-exposure distribution, do learned-model AUCs against the same frozen positive test set remain different?

The aim is to reduce structural-distribution confounding without forcing one-to-one matches outside common support.

## Frozen inputs

Positive reference:
- frozen `goldd2_atc5_positive_reference.csv`;
- 8,094 unique unordered ATC5 positive class pairs over 437 nodes;
- required SHA-256: `1f96a2c3197db78ea7fe31444efec2e4c4cfffb19a5f1fed7bd2d0ca56166018`.

Curated evidence-state reference:
- `antiddi_v3_benchmark.csv` from the validated Anti-DDI v3.0.1 Supplementary Code and Data archive;
- only `T1_wellpowered` and `T2_moderate` records;
- project each record to all valid unordered ATC level-5 class pairs from `atc_a × atc_b`;
- remove self-pairs and pairs overlapping the frozen positive reference;
- for this analysis, retain only curated pairs whose two ATC5 nodes occur in the frozen positive vocabulary.

The frozen positive CSV itself is not required to be committed to the public Science branch; byte-level provenance is enforced by checksum.

## Frozen split and training rules

- Replicates: **10** deterministic seeds, `0..9`.
- For each seed, split the 8,094 positive pairs into 80% train / 20% test using a seeded permutation.
- Degrees and every structural covariate are computed from **training positives only**.
- The primary positive test set contains held-out pairs whose endpoints were both seen in training.
- All known frozen positives and all curated Anti-DDI ATC5 pairs are excluded from random-unlabelled training/evaluation controls.
- Training unknowns are always random unlabelled nonedges; curated counter-evidence does not enter model fitting.

## Learned-model panel

Use the same three structural relation-model families already employed in the project:

1. truncated SVD latent factors, rank 32;
2. logistic NeuralMF, latent dimension 32, using the already frozen Gate-2 optimization convention;
3. two-layer LightGCN-style message passing, using the already frozen Gate-2 convention.

No model-specific tuning is permitted after evidence-state outcomes are observed.

## Random evaluation pool

For each seed, generate a random-unlabelled evaluation pool from the same 437-node positive vocabulary, excluding:
- every frozen positive pair;
- every curated T1/T2 pair;
- duplicates and self-pairs.

Pool size is fixed at **20 × the number of eligible curated pairs**, or the complete eligible nonedge universe if smaller.

## Structural covariates

For each unordered negative candidate `(u,v)`, compute training-positive endpoint degrees and define:

- `lo = min(log1p(deg(u)), log1p(deg(v)))`;
- `hi = max(log1p(deg(u)), log1p(deg(v)))`;
- `sum = lo + hi`;
- `gap = hi - lo`;
- `product = lo × hi`.

These five covariates are frozen before outcome inspection.

## Propensity standardization

Fit a logistic-regression classifier to the combined curated and random negative pools, where:
- outcome 1 = curated counter-evidence pair;
- outcome 0 = random-unlabelled pair;
- predictors = standardized five structural covariates above;
- regularization = scikit-learn default L2 logistic regression with `max_iter=2000` and deterministic random state equal to the replicate seed.

Define estimated propensity `p = P(curated | structure)`.

### Common support

Primary overlap analysis retains negative candidates with **0.05 ≤ p ≤ 0.95**. This threshold is frozen before outcome inspection and will not be tightened to improve results.

### Overlap weights

Within common support:
- curated pair weight = `1 - p`;
- random pair weight = `p`.

Weights are normalized within each negative group only for reporting; weighted AUC is invariant to a common scale within a group.

This targets the structural-overlap population rather than forcing poor one-to-one matching.

## Balance and interpretability gates

For each structural covariate report standardized mean differences (SMDs) between curated and random negative pools:
- before weighting;
- after support restriction and overlap weighting.

Also report effective sample size (ESS): `(Σw)^2 / Σw²` for each negative group.

The evidence-state-standardized contrast is considered **structurally interpretable** for a seed only if:
- maximum absolute weighted SMD across the five covariates is **< 0.10**;
- weighted ESS is **≥100** in both curated and random groups.

If these criteria fail systematically, the analysis remains a boundary result and no isolated evidence-state claim is promoted.

## Weighted AUC definition

The held-out positive test set is identical for the curated and random comparisons within a seed and receives unit weight.

For a model score `s`, weighted AUC against a weighted negative set is the weighted Mann–Whitney probability:

`P[s(positive) > s(negative)] + 0.5 P[tie]`,

where negative observations contribute their overlap weights and positives contribute equal unit weight.

Report per model and seed:
- ordinary AUC vs random-unlabelled pool;
- ordinary AUC vs curated counter-evidence;
- overlap-standardized AUC vs random-unlabelled pool;
- overlap-standardized AUC vs curated counter-evidence;
- standardized difference `AUC_random_overlap - AUC_curated_overlap`;
- maximum weighted SMD and both ESS values.

## Uncertainty

Primary uncertainty for each model is paired across the 10 frozen seeds:
- mean standardized AUC difference;
- sample SD;
- percentile 95% bootstrap CI over seeds using **20,000** bootstrap resamples with fixed bootstrap RNG seed 20260913.

The analysis is supportive only if balance/ESS gates are met. The CI is descriptive sensitivity inference over the frozen split seeds, not a claim of a biological population-sampling mechanism.

## Interpretation rules

### If balanced overlap-standardized curated pairs remain harder
Allowed:
> Evidence state contributes additional evaluation difficulty beyond the measured endpoint-degree exposure captured by the prespecified structural covariates.

Do not claim all structural confounding is removed.

### If the contrast disappears after balancing
Allowed:
> The earlier random-versus-curated difference was largely explained by measured structural exposure, reinforcing structure as the primary benchmark dimension.

### If overlap/balance is inadequate
Allowed:
> Curated counter-evidence occupies a structurally distinct support region, so evidence state cannot be isolated cleanly with the available reference; evidence state and structural exposure remain separate reporting dimensions.

All three outcomes are publishable boundary information and must be retained.