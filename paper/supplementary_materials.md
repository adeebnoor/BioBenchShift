# Supplementary Materials for “Benchmark design redirects biomedical discovery”

**Status:** PRE-SUBMISSION working draft.  
The final version will be frozen only after the leakage-free contemporary-model challenge completes and all file hashes/source-data links are updated.

## Supplementary overview

The main paper reports a four-step scientific-decision chain:

1. conventional sampled-unknown evaluation can reward structural observability;
2. learned model families depend on that signal unequally enough to change model selection;
3. benchmark-induced winner changes can redirect the identity of top biological hypotheses beyond ordinary training instability;
4. selected models can differ in later or independent evidence recovery on the same frozen candidate universe.

This Supplementary Materials file contains the detailed protocols, negative controls, stability analyses, provenance and boundaries deliberately removed from the four-figure main-text arc.

---

## Supplementary Methods

### S1. General design and outcome-free protocol freezing

For each analysis, the split rule, candidate universe, structural control, model-selection criterion and primary cutoffs were fixed before the corresponding held-out or future-evidence outcome was opened. `SCIENCE_PROJECT.md` preserves the original gates; outcome-dependent interpretation is appended in `RESULTS_CHECKPOINT_20260913.md` rather than rewriting the original hypotheses.

All structural quantities used for evaluation controls were derived from training relations only. Held-out positive edges were not allowed to determine endpoint degree in matched evaluation and were excluded from message-passing adjacency in the leakage-free contemporary graph-model challenge.

### S2. Structural-only audit across relation families

The model-free structural diagnostic assigns a score using only endpoint degree/popularity from the training graph. It intentionally contains no pair-specific molecular, sequence, disease or pharmacological feature.

The audit was applied to four external relation families:

- BioSNAP TargetDecagon drug–target interaction (DTI);
- HuRI protein–protein interaction (PPI);
- Hetionet compound–treats–disease (CtD);
- Hetionet disease–associates–gene (DaG).

Conventional evaluation contrasted held-out positives with randomly sampled unobserved pairs. Structure-neutralized evaluation paired positives with unobserved pairs matched on endpoint-degree bins under the prespecified matching procedure. Matching coverage was reported rather than silently discarding unmatched positives.

Headline degree-only AUROCs changed from 0.983→0.620 (DTI), 0.928→0.513 (HuRI PPI), 0.871→0.519 (CtD), and 0.876→0.513 (DaG). The interpretation is benchmark vulnerability, not a claim that learned models use degree alone.

### S3. Learned model families

Controlled learned-model analyses used truncated SVD latent factors, NeuralMF and LightGCN-style message passing under the same frozen relation-family splits. Model hyperparameters and random seeds are contained in the executable analysis files and result manifests.

The purpose was not to establish a new state of the art. It was to test whether plausible model families depend unequally on benchmark structure enough to change the model selected for downstream discovery.

In DTI, conventional evaluation ranked NeuralMF first whereas structure-neutralized evaluation ranked SVD first. Historical PPI produced the same NeuralMF→SVD reversal. CtD served as a non-reversal control because LightGCN remained first under both evaluation regimes.

### S4. Hypothesis-turnover definition

For a fixed biological candidate universe U and top-k sets selected by models M1 and M2, hypothesis turnover is:

**HT@k = 1 − |TopK(M1,U) ∩ TopK(M2,U)| / k.**

The candidate universe is identical within each comparison. Only the benchmark rule used to select the winning model changes.

To distinguish benchmark-induced model identity from training randomness, independently initialized model fits were averaged into score ensembles. Leave-one-initialization-out (LOO) ensembles quantified within-family ranking instability.

#### DTI

Cross-selected five-fit ensembles had HT@100 = 1.000 and HT@500 = 1.000. Within-family LOO HT@100 was 0.085 for NeuralMF and 0.042 for SVD; at top 500 it was 0.005 and 0.033.

#### PPI

On the complete 70,041,100-pair historical non-edge universe, NeuralMF-versus-SVD ensemble turnover was 0.990, 0.986 and 0.976 at top 100, 500 and 1,000. Corresponding LOO turnover was 0.137/0.139/0.130 for NeuralMF and 0.050/0.035/0.038 for SVD.

#### Disease–gene boundary

Single NeuralMF fits in the disease–gene analysis were intrinsically unstable; therefore the initial cross-model turnover could not be cleanly attributed to benchmark-induced model selection. Five-fit ensembling reduced but did not erase this limitation. Disease–gene is retained as a boundary/sensitivity analysis rather than a headline replication.

### S5. Historical BioGRID later-evidence experiment

BioGRID MV-Physical release 5.0.250 was frozen as the historical state and release 5.0.261 as the later evidence state. Download URLs and SHA-256 checksums are recorded in the repository.

The historical human network contained 93,146 unique interactions among 11,844 proteins. Model selection was performed using only the historical snapshot. The conventionally evaluated winner was NeuralMF; the degree-matched winner was SVD.

For the non-circular consequence analysis, each model was trained on the complete historical graph and ranked the identical 70,041,100 unordered protein pairs absent from the historical snapshot. No degree matching, sampled future negatives or future-control construction was used in this endpoint.

The later release added 5,635 relations for which both endpoints were already represented historically. SVD recovered 7/18/26/87/159/522 later-supported relations at top 100/500/1,000/5,000/10,000/50,000. NeuralMF recovered 0/1/5/22/42/181; LightGCN recovered 0/4/5/30/66/241.

Later additions are treated as evidence accumulation, not an unbiased biological truth set.

### S6. BioGRID uncertainty

The 5,635 later-added relations were resampled with replacement in 20,000 paired bootstrap replicates. At each replicate and cutoff, recovery was recomputed for the same selected-model rankings. SVD-minus-NeuralMF recall differences remained positive at every prespecified reported cutoff.

Reported differences (95% bootstrap interval):

- top 100: 0.0012 (0.0004, 0.0023);
- top 1,000: 0.0037 (0.0020, 0.0057);
- top 10,000: 0.0208 (0.0163, 0.0252);
- top 50,000: 0.0605 (0.0531, 0.0680).

### S7. Independent ChEMBL 37 DTI evidence

The ChEMBL evidence protocol was frozen before outcome inspection. BioSNAP alone determined model fitting and model selection. ChEMBL was used only as a separately curated evidence source.

Primary evidence required:

- Homo sapiens;
- direct SINGLE PROTEIN target;
- binding assay;
- target confidence score 9;
- pChEMBL ≥ 6.

pChEMBL ≥ 7 was frozen as a stricter sensitivity threshold.

Mapped coverage included 264/284 BioSNAP drugs (93.0%) and 1,897/3,648 genes (52.0%). Sixty-six mapped unknown relations met the primary evidence definition and 26 met the stricter threshold.

Primary evidence recovery at top 100/500/1,000 was 3/5/6 for SVD and 0/0/2 for NeuralMF. Under pChEMBL ≥ 7 the corresponding counts were 2/3/3 and 0/0/2. At broad cutoffs the ranking crossed: at top 10,000 NeuralMF recovered 19 primary-evidence relations versus 15 for SVD, and at top 50,000 the counts were 35 versus 23.

The broad-cutoff crossover is a required boundary, not an inconvenient result to omit.

### S8. Anti-DDI evidence-state sensitivity

Anti-DDI is not a principal evidence source for the paper. It is used to test the additional proposition that “random unknown” and “curated counter-evidence” are not interchangeable negative states.

A frozen 8,094-pair ATC5 positive reference was combined with the T1/T2 Anti-DDI class-pair projection. Exact matching had limited coverage and nearest-degree matching left material residual imbalance, so the final sensitivity analysis used overlap weighting on prespecified structural covariates with explicit balance and effective-sample-size diagnostics.

After overlap standardization, mean random-minus-curated AUROC differences were approximately +0.048 for NeuralMF, +0.031 for LightGCN and +0.008 for SVD. However, the strict prespecified balance gate passed in only 5/10 seeds. Therefore the paper states only that a residual model-dependent evidence-state signal is compatible with the data within common structural support; it does not claim causal isolation of evidence state.

### S9. Contemporary GraphBAN challenge

GraphBAN was selected as a contemporary DTI architecture because it is a recent feature-rich graph method evaluated on BioSNAP. A reproducibility audit identified that the upstream transductive evaluation code constructs validation/test graph adjacency using positive edges from those held-out splits. That upstream evaluation is therefore not used as clean primary evidence in this paper.

A protocol amendment was frozen before the clean result was opened:

- frozen TargetDecagon relation set;
- ChemBERTa drug features and ESM-1b protein features following the public GraphBAN feature logic;
- train-positive graph only for message passing;
- held-out positives excluded from adjacency;
- same conventional versus training-degree-matched evaluation logic as the controlled DTI analyses;
- prespecified mapping gate ≥90% positive-edge feature coverage.

The mapping gate passed at 18,631/18,690 positive edges (99.7%). The final model outcome will be inserted regardless of direction. Structural sensitivity will be reported as architectural extension; comparative robustness will be reported as a model-specific boundary.

### S10. Software, data provenance and reproducibility

The branch contains:

- analysis scripts;
- workflow definitions;
- release and database identifiers;
- checksum/provenance files;
- seed-level replicate outputs;
- summary outputs;
- main-figure source data;
- source-to-result manifest files.

Before submission, the exact manuscript commit will be tagged and archived under a persistent DOI. External datasets remain under their original terms and are referenced rather than redistributed where redistribution is inappropriate.

---

## Extended Data plan

### Extended Data Fig. 1 — Full structural-only replicates

Seed-level conventional and degree-matched AUROC distributions for DTI, HuRI PPI, CtD and DaG, with matching coverage.

### Extended Data Fig. 2 — Full learned-model sensitivity

SVD, NeuralMF and LightGCN conventional versus neutralized evaluation across DTI, CtD and DaG; include CtD non-reversal control explicitly.

### Extended Data Fig. 3 — Disease–gene instability boundary

Single-fit and five-fit ensemble turnover compared with within-family instability.

### Extended Data Fig. 4 — DTI turnover across split seeds

HT@100/500/1,000 across original frozen split seeds plus ensemble-control comparison.

### Extended Data Fig. 5 — PPI complete-universe turnover

Cross-selected versus LOO turnover at top 100/500/1,000.

### Extended Data Fig. 6 — BioGRID full cutoff and enrichment curves

Later-supported hits, precision lower bounds, recall and enrichment relative to uniform selection.

### Extended Data Fig. 7 — BioGRID sensitivity

Degree-stratified later-evidence recovery and any additional frozen sensitivity analyses.

### Extended Data Fig. 8 — ChEMBL mapping and broad-cutoff boundary

Mapping flow, primary pChEMBL ≥6 results, pChEMBL ≥7 sensitivity, and full cutoff curves showing crossover.

### Extended Data Fig. 9 — Anti-DDI overlap weighting

Raw and weighted covariate balance, ESS, per-seed random-versus-curated performance and 5/10 strict balance-gate result.

### Extended Data Fig. 10 — Contemporary architecture audit

Upstream GraphBAN protocol audit, mapping gate and final leakage-free challenge result after freezing.

---

## Supplementary tables

- **Table S1:** datasets, relation families, versions/releases, entity counts and provenance.
- **Table S2:** all model hyperparameters and seed rules.
- **Table S3:** degree-only replicate results and matching coverage.
- **Table S4:** learned-model conventional/neutralized results.
- **Table S5:** DTI and PPI hypothesis-turnover results.
- **Table S6:** BioGRID later-evidence counts, recall, precision lower bounds and enrichment.
- **Table S7:** BioGRID paired bootstrap estimates.
- **Table S8:** ChEMBL evidence definitions, mappings and cutoff results.
- **Table S9:** Anti-DDI overlap-weighting balance and ESS diagnostics.
- **Table S10:** contemporary GraphBAN challenge configuration and result.
- **Table S11:** file hashes, database release identifiers and software environment.

---

## Claims and interpretation guardrails

The Supplementary Materials preserve the same claim limits as the main paper. The analyses do not establish that all biomedical AI models are structurally biased, that degree matching is uniquely correct, that SVD is biologically superior in general, that unobserved pairs are negatives, that later BioGRID additions are unbiased truth, or that curated Anti-DDI candidates are clinically safe. The central inference is narrower: **benchmark design can change model selection and thereby change the scientific hypotheses and evidence concentration produced by a biomedical AI discovery workflow.**
