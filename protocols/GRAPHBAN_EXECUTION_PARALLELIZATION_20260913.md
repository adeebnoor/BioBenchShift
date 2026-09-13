# GraphBAN S6 execution-only parallelization amendment

**Frozen:** 13 September 2026, before observing any leakage-free GraphBAN S6 performance outcome.

## Reason

The original leakage-free S6 run uses the prespecified GraphBAN-style feature definition but computes ESM-1b representations for 3,607 mapped proteins serially on a CPU GitHub Actions runner. Feature extraction, rather than model training or evaluation, is the wall-clock bottleneck.

This amendment changes **execution only**. It does not alter the scientific experiment.

## Unchanged scientific specification

The following remain exactly as frozen in `CONTEMPORARY_MODEL_CHALLENGE.md` and `GRAPHBAN_PROTOCOL_AMENDMENT_20260913.md`:

- frozen TargetDecagon positive relation universe;
- mapping gate and mapped endpoint set;
- PubChem/SMILES and human UniProt sequence mappings;
- ChemBERTa checkpoint: `DeepChem/ChemBERTa-77M-MTR`;
- ChemBERTa representation: first-token hidden representation under the existing tokenizer/max-length rule;
- ESM checkpoint: `esm1b_t33_650M_UR50S`;
- ESM representation: layer 33, mean over residue representations, sequences truncated at 1,022 residues;
- canonical sorted drug and protein endpoint ordering;
- positive train/test split convention;
- seeds 0, 1 and 2;
- exclusion of held-out positive edges from message passing;
- training-negative sampling rule;
- GraphSAGE-style heterogeneous encoder and decoder;
- hidden dimension, optimizer, learning rate, epochs and batch rules;
- conventional random-unlabelled evaluation;
- structure-neutralized log2 degree-bin matching using training-edge degrees only;
- AUROC/AUPRC metrics and matching-coverage reporting;
- interpretation rule: report the result regardless of direction.

## Parallel execution

Protein endpoints are partitioned deterministically by their index in the canonical sorted mapped-protein list. Independent workers compute the same frozen ESM-1b representation on disjoint index subsets. Each output stores both canonical indices and feature vectors. A deterministic combine step requires every canonical index exactly once and reconstructs the original endpoint order before any scaling, graph construction, training or evaluation.

Drug features are computed once using the same ChemBERTa definition and canonical sorted drug order.

Only after feature arrays are reassembled does the original S6 training/evaluation logic run. No biological relation, endpoint, score, split, seed, threshold, matching rule or model-selection criterion is inspected or changed during feature parallelization.

## Reproducibility check

The final S6 summary must record that precomputed features were generated under this amendment. If both the original serial run and this parallelized run complete, they are treated as replicate implementations of the same frozen experiment. Material disagreement triggers investigation; results are not selected based on which direction is more favorable.