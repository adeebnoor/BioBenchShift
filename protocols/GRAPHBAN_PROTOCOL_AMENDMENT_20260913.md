# GraphBAN protocol amendment — leakage-risk safeguard

**Frozen:** 2026-09-13, before any GraphBAN benchmark-sensitivity outcome was inspected.  
**Applies to:** `CONTEMPORARY_MODEL_CHALLENGE.md` / Science send gate S6.

## Why this amendment was necessary

During implementation of the predeclared GraphBAN challenge, inspection of the authors' public transductive code showed that the graph used to score a validation or test split is constructed from the positive edges of that same split. In the public `load_edge_csv` routine, rows with `Y == 1` become `edge_index`; `run_model` then separately constructs `data_val.edge_index` from validation positives and `data_test.edge_index` from test positives before evaluating edge labels.

Because message passing can therefore receive information about held-out positive edges, the upstream transductive prediction path is not acceptable as the primary leakage-free S6 test for this project. This is a protocol-integrity decision, not a response to an unfavorable GraphBAN result: it was made while the first upstream checkpoint inference was still running and before its random-versus-matched AUC was available.

## Status of the upstream-checkpoint run

The authors' pinned BioSNAP seed-12 checkpoint run may be retained as an **upstream-protocol diagnostic** only. Its output, if completed, cannot by itself pass S6 and cannot be used to claim leakage-free generalization of the structural-benchmark effect.

## Primary S6 test remains the predeclared clean adaptation

The primary contemporary-model test remains the experiment already specified in `CONTEMPORARY_MODEL_CHALLENGE.md`:

- frozen Stanford TargetDecagon DTI positives used elsewhere in this project;
- GraphBAN transductive architecture adapted to the frozen node universe;
- molecular features from ChemBERTa and protein features from ESM-1b, following the public GraphBAN implementation;
- **message-passing graph contains training positives only**;
- held-out positives never enter graph construction, feature fitting, negative construction, or model selection;
- one trained model per frozen split is evaluated against both conventional random-unlabelled controls and the existing log2 degree-matched controls;
- degrees are computed from training positives only;
- the same trained scores are used for both evaluation regimes;
- no ChEMBL information enters training or evaluation;
- mapping coverage and all exclusions are reported.

## Frozen implementation choices

Before outcome inspection, the clean adaptation fixes the following choices:

- positive split: the same deterministic 80/20 seed convention used by the existing TargetDecagon Gate-2 code;
- primary seeds: 0, 1, 2;
- seen-endpoint held-out positives are the primary evaluable positives;
- GraphBAN-style encoder: three heterogeneous GraphSAGE layers, hidden dimension 256;
- decoder: the public GraphBAN transductive four-layer edge MLP over concatenated drug/protein node embeddings;
- initial drug feature: `DeepChem/ChemBERTa-77M-MTR` CLS embedding;
- initial target feature: ESM-1b (`esm1b_t33_650M_UR50S`) mean residue embedding, truncated at the model's 1022-residue input limit as in the public code;
- training controls: random unlabelled nonedges excluding all known TargetDecagon positives;
- fixed optimization: Adam, learning rate 0.001, balanced positive/control BCE objective, 20 epochs, no regime-specific tuning;
- feature embeddings are computed once per mapped node and reused across the three seeds;
- primary evaluation reports AUROC and AUPRC under conventional and degree-matched controls, matching coverage, and seed-wise values.

The 20-epoch choice is a fixed engineering compromise for the reproducible challenge and is not tuned after observing outcomes.

## Mapping rule

- PubChem CID → canonical SMILES through PubChem PUG REST.
- GeneID → Homo sapiens UniProtKB through the official UniProt ID-mapping service.
- Prefer reviewed human UniProt entries when available; otherwise use a unique human mapping, with multiplicity recorded.
- The original mapping gate remains: at least 90% of all positive edges and at least 90% of evaluable seen-endpoint held-out positives must have both feature mappings.

If the mapping gate fails, S6 is mapping-limited/inconclusive.

## Interpretation

A large random-to-matched drop supports benchmark sensitivity beyond the simpler model ladder. A small or absent drop is equally reportable and bounds the claim. The S6 criterion is successful evaluation of a contemporary architecture under the leakage-free frozen dual-regime protocol, not a favorable direction.

## Reproducibility boundary

The public GraphBAN implementation and pinned upstream commit remain cited as architectural provenance. This project's clean adaptation is not represented as a verbatim reproduction of the published transductive benchmark because the held-out-edge construction is intentionally corrected to satisfy the leakage-free design required here.