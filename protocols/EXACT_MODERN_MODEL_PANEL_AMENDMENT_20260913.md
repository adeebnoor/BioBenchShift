# Exact Modern-Model Panel Amendment — Science Expansion

**Frozen:** 13 September 2026  
**Prerequisite:** `PROSPECTIVE_DTI_VALIDATION_PROTOCOL_20260913.md`  
**Immutable baseline:** `1298e4c83ff4fa488a784d52219d20e706118a9e`

## Rationale

The previously frozen GraphBAN-style analysis is a legitimate leakage-free architecture stress test, but it is not a verbatim implementation of the published GraphBAN model. For the Science expansion, this amendment opens an **exact published-code model panel** so that benchmark susceptibility can be tested across contemporary DTI architectures rather than inferred from one adapted model.

No prospective wet-lab outcome is available at the time of this amendment. The panel is frozen before any new exact-model benchmark result is used to change manuscript claims.

## Primary contemporary architectures

The primary 2025 panel consists of public author implementations pinned to exact repository commits observed on 13 September 2026:

1. **GraphBAN** — `HamidHadipour/GraphBAN` @ `2dded31a46edfb18d685fb1bb361a95e99eaae49`.
2. **GS-DTI** — `purvavideha/GSDTI` @ `868efff7d7b678f3773a59acfa1fda11f1c294f6`.
3. **SP-DTI** — `Steven51516/SP-DTI` @ `3f487a26eeb117711d9e9111ee7c44828fa194f4`.
4. **DTIAM** — `CSUBioGroup/DTIAM` @ `61d68b8148a86b2aea657a3e9d707a491de70ad8`.

These models were selected because they represent materially different contemporary DTI strategies: inductive graph/attribute learning, graph/sequence large-model features, structure/subpocket-aware modeling, and self-supervised drug/protein representation learning.

## Secondary published reference architectures

Two earlier but widely used published implementations are included as secondary reference points:

5. **DrugBAN** — `peizhenbai/DrugBAN` @ `9923f8c99959e00263103ff9ac61ba0eaccc8e02`.
6. **MolTrans** — `kexinhuang12345/MolTrans` @ `47ac16b8c158b080ba6cdaec74cd7aa9c1332b73`.

The original SVD, NeuralMF, and LightGCN ladder remains the frozen baseline and is not retroactively redefined.

## Fidelity rule

For each external architecture:

- preserve the model-defining representation and interaction modules;
- use the authors' public implementation rather than rewriting the architecture from memory;
- record the exact repository commit, environment, pretrained asset identifiers, and any compatibility patch;
- any compatibility patch must be mechanical (path, dependency, input-adapter, device, serialization) unless explicitly disclosed as a scientific modification;
- if a scientific modification is required to place the model on the frozen BioSNAP universe, label the run `adapted` and do not present it as an exact reproduction.

## Common benchmark contract

All models are evaluated against the same frozen BioSNAP DTI positive universe and the same evaluation intervention.

For each split seed:

1. use the same train/test positive partition;
2. exclude held-out positives from training and all graph/message-passing structures;
3. construct conventional random-unlabelled controls from the allowed candidate universe;
4. construct the original joint endpoint-degree-matched neutralized controls without changing the matching definition after model results are seen;
5. record matching coverage;
6. report AUROC and AUPRC as benchmark metrics;
7. preserve raw prediction probabilities/scores for all evaluated pairs;
8. where technically feasible, generate full-candidate scores for top-queue analysis.

## Replication

Target 10 seeds for each model. If the published implementation has deterministic frozen feature extraction plus a stochastic downstream learner, repeat the stochastic learner 10 times. If an architecture is prohibitively expensive, a minimum of 5 completed seeds is allowed only with the computational limitation disclosed before winner comparison.

## Model-selection rule

Within the **primary contemporary panel**, models are ranked separately under conventional and neutralized evaluation by mean AUROC, with mean AUPRC as the frozen tie-breaker.

The secondary DrugBAN/MolTrans models are reported but do not change the primary 2025 winner unless a primary model cannot be executed faithfully. Any replacement must be documented before opening candidate-queue turnover.

## Decision-relevant robustness outputs

For every primary model report:

- conventional AUROC/AUPRC;
- neutralized AUROC/AUPRC;
- absolute and relative performance drop;
- rank under each benchmark;
- matching coverage;
- seed-level uncertainty;
- model-class descriptors needed to interpret susceptibility.

The main scientific question is not whether every modern model reverses. It is whether **benchmark susceptibility is moderated by model class** and whether benchmark choice changes the model a scientist would rationally select.

## Modern winner-change gate

Only if conventional and neutralized evaluation select different primary contemporary winners is a full modern hypothesis-queue turnover analysis opened.

If the winner changes:

- train/freeze a 5-member ensemble for each selected architecture on one fixed split;
- score the identical assay-eligible candidate universe;
- report HT@100, HT@500, HT@1,000, rank-biased overlap, and within-model ensemble instability;
- open the secondary prospective wet-lab panel defined in the prospective protocol (25 candidates per selected policy).

If the winner does not change:

- retain the null winner result;
- report model-specific benchmark sensitivity;
- do not manufacture a post hoc model contrast;
- the primary SVD-vs-NeuralMF prospective experiment remains valid because it tests the already frozen causal decision chain rather than universal modern-model reversal.

## Interpretation boundary

A robust feature-rich model is evidence against a universal reversal claim, not a failure of the project. The intended Science-level inference is conditional: benchmark construction can alter scientific prioritization, while the magnitude of that effect depends on model class, neutralization rule, and decision scale.

## Reproducibility outputs

Create one machine-readable manifest with:

- model name;
- publication DOI;
- repository URL;
- pinned commit SHA;
- environment lock/checksum;
- pretrained asset identifiers/checksums;
- exact command;
- split seed;
- benchmark condition;
- metrics;
- output prediction checksum;
- fidelity status (`exact`, `mechanical_patch`, `adapted`, `failed`).

All failures remain part of the record.
