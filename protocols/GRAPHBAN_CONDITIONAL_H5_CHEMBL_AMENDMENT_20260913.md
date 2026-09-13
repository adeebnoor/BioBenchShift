# Conditional GraphBAN-style H5 + ChEMBL amendment — 13 September 2026

**Frozen before generating GraphBAN-style full-graph rankings or inspecting their H5/ChEMBL outcomes.**

Baseline submission commit: `1298e4c83ff4fa488a784d52219d20e706118a9e`.

This amendment implements the conditional branch already required by `CONTEMPORARY_MODEL_CHALLENGE.md`. The gate was activated mechanically because the frozen GraphBAN-style S6 mean AUROC is the numerical maximum in at least one frozen evaluation regime when compared on seeds 0–2; see `results/science_expansion/GRAPHBAN_MODEL_SELECTION_GATE_20260913.md`.

## 1. Architecture and feature lock

Use the exact frozen GraphBAN-style implementation and feature definitions already used by S6:

- mapped TargetDecagon drug/gene universe defined by the committed mapping files;
- `DeepChem/ChemBERTa-77M-MTR` drug features and ESM-1b layer-33 protein features;
- three-layer heterogeneous GraphSAGE encoder, hidden dimension 256;
- four-layer nonlinear edge MLP decoder with dropout 0.2;
- Adam learning rate 0.001, 20 epochs;
- frozen GraphBAN seeds `0,1,2`;
- random training negatives sampled only from pairs absent from the frozen positive graph.

No architecture, feature model, epoch count, hidden dimension, learning rate, or seed may be changed after outcome inspection.

## 2. Full-graph ranking fit

The S6 held-out protocol remains the only performance-evaluation experiment. For downstream scientific-decision ranking only, refit the same frozen architecture separately at seeds 0, 1 and 2 on **all mapped positive TargetDecagon edges**. Because there is no held-out performance estimate in this stage, all mapped positive edges may be used for message passing and training. Fit feature standardization on the endpoints present in the complete mapped training graph and apply it to all mapped nodes.

For each seed, score every mapped drug–gene pair not present in the frozen positive graph. The GraphBAN-style ensemble score is the arithmetic mean of the three seed-level scores. No ChEMBL, BindingDB, pathway, BioGRID-future, or Open Targets information enters training or ranking.

## 3. Identical candidate universe for H5

H5 comparisons are restricted to the exact GraphBAN-mappable unknown universe: mapped drugs × mapped genes minus all frozen known positive TargetDecagon edges. Comparator models (NeuralMF, SVD and LightGCN) are fitted/ranked using their already-frozen full-graph ensemble definitions, but their score matrices are restricted to this identical candidate universe before Top-K extraction.

Primary K values are unchanged: `100, 500, 1000`.

For GraphBAN-style versus each comparator, report without selection:

1. Top-K hypothesis turnover `1 - |A∩B|/K`;
2. Jaccard overlap;
3. Spearman correlation of full candidate-universe scores (equivalently rank correlation because Spearman is computed over scores);
4. the number of shared Top-K pairs.

The comparison to the previously selected NeuralMF and SVD queues is primary for continuity; LightGCN is retained as a prespecified secondary comparator. No comparator may be dropped after results are known.

## 4. Within-GraphBAN instability control

Using the three frozen full-graph GraphBAN fits, compute pairwise seed-to-seed Top-K turnover for all three seed pairs at K=100,500,1000. Also construct three leave-one-seed-out two-model averages and compute pairwise turnover among those LOO ensembles. Report mean and maximum instability.

A cross-model queue difference is interpreted as model/benchmark consequential only when it is materially larger than the corresponding within-GraphBAN instability. No numerical threshold is invented after seeing results; all raw values are reported.

## 5. ChEMBL 37 conditional external validation

Use the unchanged protocol and evidence definitions from `EXTERNAL_DTI_CHEMBL_PROTOCOL.md` / `analysis/run_chembl_dti_external.py`:

- ChEMBL database version must equal `ChEMBL_37`;
- human `SINGLE PROTEIN` targets only;
- binding assays (`assay_type=B`) with confidence score 9;
- valid data-validity comment;
- primary support threshold pChEMBL >= 6;
- sensitivity support pChEMBL >= 7 and exact-relation pChEMBL >= 6;
- K = `100, 500, 1000, 5000, 10000, 50000`;
- known BioSNAP TargetDecagon positives excluded;
- unlabelled ChEMBL pairs are not negatives.

Before any GraphBAN hit count is interpreted, the live ChEMBL mapping/evidence reconstruction must reproduce the frozen baseline fingerprint exactly:

- mapped drugs = `264` of `284`;
- mapped genes = `1897` of `3648`;
- mapped candidate universe = `486560`;
- mapped molecules = `264`;
- raw mapped targets = `1876`;
- eligible human single-protein targets = `1790`;
- eligible pChEMBL>=6 activities on mapped targets = `9509`;
- confidence-9 binding assays = `4636`;
- supported unknown candidate pairs pChEMBL>=6 = `66`;
- supported pairs pChEMBL>=7 = `26`;
- exact-relation pChEMBL>=6 supported pairs = `66`.

If any fingerprint count differs, the conditional ChEMBL run must stop and report snapshot drift. It must not substitute the new counts or continue to a favorable result.

For a matching fingerprint, add `GraphBAN-style` to the existing full-universe primary table and mapped-universe table using exactly the same evidence sets and K values. The existing NeuralMF/SVD/LightGCN rows remain unchanged and are not recomputed for narrative convenience.

## 6. Interpretation lock

This conditional analysis is not designed to restore the original NeuralMF→SVD winner reversal. It asks whether a contemporary feature-rich selected model produces a materially different scientific queue and where independently supported ChEMBL pairs occur in that queue.

Permitted conclusions include:

- GraphBAN-style is sensitive to benchmark construction even if its model identity remains unchanged;
- GraphBAN-style changes or preserves the identity of early hypotheses relative to prior selected models;
- ChEMBL support appears earlier, later, similarly, or not at all in the GraphBAN-style ranking.

Prohibited conclusions include:

- calling ChEMBL-unlabelled pairs negatives;
- claiming GraphBAN-style is the full published GraphBAN reproduction;
- claiming universal winner reversal across model classes or matching definitions;
- altering the model, candidate universe, K values, evidence thresholds, seeds or feature representations after seeing outcomes.

All outcomes, including null or contradictory outcomes, must be retained in the expansion evidence matrix and final manuscript decision.
