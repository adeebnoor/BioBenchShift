# GraphBAN -> ChEMBL follow-through

**Frozen:** 14 September 2026, before inspecting any GraphBAN ChEMBL recovery outcome.

## Trigger

The prespecified contemporary-model protocol required an external-evidence follow-through if GraphBAN became the selected DTI model under either frozen evaluation regime. The fair five-seed four-model panel subsequently selected GraphBAN by mean AUROC under both conventional and structure-neutralized evaluation. This follow-through therefore completes a predeclared branch of the analysis; it is not an outcome-driven addition.

## Scientific question

Does the contemporary GraphBAN winner concentrate independently supported ChEMBL 37 drug-target relations near the top of the same frozen DTI candidate universe, and how does that external-evidence profile compare descriptively with SVD, NeuralMF and LightGCN when all four models are trained and ranked on the identical GraphBAN-mapped TargetDecagon universe?

Because GraphBAN is selected under both evaluation regimes, there is **no benchmark-induced DTI model-identity switch in the full contemporary panel**. Consequently this analysis is an external-evidence characterization and robustness boundary, not a new winner-turnover test.

## Frozen universe and models

- TargetDecagon input and GraphBAN feature mappings are unchanged.
- Candidate universe: every mapped drug x mapped target pair not present among known mapped TargetDecagon positives.
- All four models rank exactly this same candidate universe.
- Baseline families: SVD, NeuralMF, LightGCN, using the existing five-fit ensemble definitions from the frozen ChEMBL analysis.
- Contemporary family: leakage-free GraphBAN-style architecture from S6, using the already frozen ChemBERTa and ESM-1b node features.
- GraphBAN is trained on the complete mapped known-positive graph because the external-evidence analysis ranks unknown pairs after model selection, matching the full-graph training convention used by the existing ChEMBL ensemble analysis.
- GraphBAN ensemble: five independently initialized fits with fixed training seeds 91000-91004; architecture, optimizer, learning rate, negative-sampling rule, hidden dimension and 20 epochs are unchanged from S6.
- Node feature scaling is fit once on the complete mapped node set before the five GraphBAN fits; no ChEMBL information enters feature scaling, training or model selection.

## Frozen ChEMBL evidence rule

No evidence definition is changed from the preregistered ChEMBL 37 analysis:

- ChEMBL 37 only;
- Homo sapiens;
- direct SINGLE PROTEIN targets;
- binding assays;
- confidence score 9;
- primary support: pChEMBL >= 6;
- strict sensitivity: pChEMBL >= 7;
- known TargetDecagon positives excluded before external-evidence counting.

The prespecified cutoffs remain K = 100, 500, 1,000, 5,000, 10,000 and 50,000.

## Reporting rule

For every model and K, report supported-pair count, precision lower bound, supported-pair recall and enrichment over uniform ranking on the identical candidate universe. Report GraphBAN regardless of whether its external-evidence profile is favorable or unfavorable. Do not use ChEMBL results to alter model selection, thresholds, architecture, seeds, or the manuscript's stated boundary that GraphBAN remains the DTI winner under both benchmark regimes.
