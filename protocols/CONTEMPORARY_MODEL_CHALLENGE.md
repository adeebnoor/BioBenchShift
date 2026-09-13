# Frozen contemporary-model challenge — GraphBAN

**Protocol freeze:** 2026-09-13  
**Branch:** `main`  
**Primary architecture:** GraphBAN (Nature Communications, 2025)  
**Purpose:** test whether the benchmark-structure findings survive a contemporary feature-rich compound–protein interaction model rather than only latent-factor / LightGCN-style baselines.

## Why GraphBAN

GraphBAN is a contemporary CPI model with molecular and protein features plus bipartite graph information, published in Nature Communications in 2025. Its paper explicitly benchmarks BioSNAP, BindingDB and KIBA and evaluates transductive and inductive settings. It is therefore a more relevant challenge to the DTI conclusion than adding an arbitrary generic GNN.

Important dataset note: the GraphBAN paper's processed BioSNAP benchmark is not assumed to be identical to the frozen Stanford TargetDecagon positive-only network used in this project. The paper reports 13,741 BioSNAP drug–target pairs, whereas our frozen TargetDecagon source contains 18,690 edges. We will therefore **adapt the published GraphBAN architecture to our already frozen TargetDecagon relation universe** rather than substitute GraphBAN's processed benchmark and call it the same experiment.

## Non-negotiable anti-cherry-picking rule

This protocol is frozen before inspecting GraphBAN performance under the project's random versus structure-neutralized evaluation regimes.

We will not:

- change the primary DTI split because GraphBAN performs poorly;
- tune separate hyperparameters for conventional and neutralized evaluation;
- change the degree-matching rule after seeing GraphBAN results;
- replace GraphBAN with another contemporary model solely because its result is more favorable;
- drop GraphBAN from reporting if it is structurally robust and weakens the universality claim.

## Frozen data universe

Primary data are the same frozen BioSNAP TargetDecagon positive relations used in H1/H2/H5.

- Drugs: frozen PubChem-CID nodes.
- Targets: frozen GeneID nodes.
- Positives: frozen TargetDecagon edge set.
- Positive train/test split: use the same project split convention and seeds already used for Gate 2; no ChEMBL information enters this stage.
- Held-out positives are excluded from graph construction/message passing.

## Feature mapping

GraphBAN requires compound and protein attributes.

### Compound features

- PubChem CID -> canonical/isomeric SMILES through the official PubChem PUG REST service.
- RDKit processing follows the GraphBAN implementation as closely as possible.

### Protein features

- GeneID -> human UniProtKB accession through the official UniProt ID-mapping service.
- UniProtKB protein sequence is retrieved from the official UniProt service.
- If multiple eligible human proteins map to a GeneID, use the reviewed canonical mapping when unambiguous; mapping multiplicity is recorded.

### Mapping gate

The primary GraphBAN challenge is declared interpretable only if both:

- >=90% of positive edges have valid mapped compound and target features; and
- >=90% of held-out seen-endpoint positive edges used in evaluation have valid mapped features.

If this gate fails, GraphBAN is labeled **mapping-limited/inconclusive**, not selectively replaced.

## Architecture and hyperparameters

- Start from the authors' public GraphBAN code associated with the 2025 Nature Communications paper.
- Primary mode: transductive GraphBAN, because the main structural-bias comparison concerns seen endpoints and already excludes cold-start positives from the primary matched test.
- Use the authors' published/default BioSNAP transductive hyperparameters where compatible with our mapped TargetDecagon feature universe.
- Any engineering change needed solely to accept our identifiers/data format is documented in a patch manifest.
- No hyperparameter may be optimized separately for random-negative versus degree-matched evaluation.
- Training randomness: at least 3 frozen seeds; 5 if computationally feasible.

## Training negatives

To preserve comparability with the controlled model ladder, GraphBAN training uses random unlabelled nonedges drawn from the training node universe while excluding all known TargetDecagon positives.

The same frozen training-negative protocol is used regardless of how the model will later be evaluated.

Training unknowns are not called biological negatives in the manuscript.

## Primary evaluation

For each frozen test seed, GraphBAN is evaluated on:

1. **conventional random-unlabelled controls**, using the same project sampling rule;
2. **structure-neutralized degree-matched controls**, using the same prespecified log2 endpoint-degree-bin rule and training-edge degrees only.

Report for each regime:

- AUROC;
- AUPRC;
- number of evaluable held-out positives;
- matching coverage;
- seed-wise values and mean ± SD / confidence interval.

## Primary contemporary-model question

> Does a contemporary feature-rich GraphBAN model retain a materially different apparent performance under conventional versus structure-neutralized evaluation on the frozen TargetDecagon task?

A large drop supports generalization of benchmark sensitivity beyond the simple model ladder.

A small or absent drop is equally important: it bounds the phenomenon and argues that richer molecular/protein features can reduce susceptibility.

## Model-selection consequence

After GraphBAN results are frozen, add it to the same model-selection panel with SVD, NeuralMF and LightGCN.

We then report, without changing the evaluation rule:

- which model is selected under conventional evaluation;
- which model is selected under structure-neutralized evaluation;
- whether the existing NeuralMF -> SVD winner reversal persists, disappears, or becomes a different reversal.

This outcome is reported regardless of direction.

## Hypothesis-identity consequence

Only if GraphBAN becomes the winner under at least one frozen evaluation regime do we promote it into the H5 scientific-decision comparison.

In that case:

- all compared winners rank the identical frozen unknown DTI candidate universe;
- use ensemble scores across frozen initializations;
- report HT@100/500/1000 and score-rank correlation;
- compare cross-selected-model turnover with within-GraphBAN initialization/ensemble instability.

If GraphBAN wins neither regime, its role remains a robustness challenge rather than an artificially inserted Top-K comparison.

## ChEMBL external evidence

ChEMBL 37 remains completely external to training and model selection.

If GraphBAN becomes a selected winner under a frozen evaluation regime, its full-universe ranking is added to the already frozen ChEMBL external-validation protocol **without altering ChEMBL evidence definitions, thresholds, mappings or K values**.

If it is not a selected winner, GraphBAN may be shown descriptively but does not replace the prespecified SVD-versus-NeuralMF confirmatory contrast.

## Secondary stress tests

Only after the primary GraphBAN challenge is frozen may we consider:

- GraphBAN inductive/cold-start evaluation;
- BioPathNet / NBFNet-style relational path reasoning on a compatible heterogeneous graph;
- TAPB/DrugBAN-style prior-debiasing comparison.

These are secondary analyses and cannot replace an unfavorable primary GraphBAN result.

## Interpretation rules

### If GraphBAN is strongly benchmark-sensitive

Allowed conclusion:

> The evaluation effect is not confined to simple latent or message-passing baselines and persists in a contemporary feature-rich CPI architecture.

### If GraphBAN is comparatively robust

Allowed conclusion:

> Structural benchmark sensitivity is model-dependent; a contemporary feature-rich architecture is more robust, while benchmark design still changes the ranking and scientific decisions among susceptible model families.

This weakens universality but strengthens mechanistic precision.

### If mapping/training is inconclusive

Report the failure transparently and do not count S6 as passed.

## Science send criterion S6

S6 is passed when a contemporary architecture is successfully evaluated under the frozen dual-regime protocol and the result — positive or negative — is fully incorporated into the claim boundary. Passing S6 does **not** require GraphBAN to support our preferred direction; it requires that the central paper survive the strongest fair contemporary challenge.
