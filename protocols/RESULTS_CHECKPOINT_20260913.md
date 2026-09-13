# Science structural-shortcuts project — results checkpoint

**Freeze date:** 2026-09-13

This file records outcomes separately from `SCIENCE_PROJECT.md`. The original hypotheses and go/no-go gates are intentionally left unchanged after outcome inspection.

## H1 / Gate 1 — cross-domain structural inflation

Status: **PASS (4/4 external relation families)**.

Degree/popularity-only structural null under conventional versus degree-matched evaluation:

| Relation family | Conventional AUC | Degree-matched AUC | Difference | Matching coverage |
|---|---:|---:|---:|---:|
| BioSNAP DTI | 0.983 | 0.620 | 0.364 | 0.599 |
| HuRI PPI | 0.928 | 0.513 | 0.415 | 1.000 |
| Hetionet Compound–treats–Disease | 0.871 | 0.519 | 0.352 | 1.000 |
| Hetionet Disease–associates–Gene | 0.876 | 0.513 | 0.363 | 1.000 |

Anti-DDI remains the motivating/seed observation and is not counted among these four external families.

## H2 / Gate 2 — learned-model sensitivity

Status: **SUPPORTED, BUT ORIGINAL STRICT GATE NOT YET DECLARED PASS**.

Across the current controlled model ladder, NeuralMF is strongly sensitive in DTI, CtD and DaG. LightGCN shows larger sensitivity in DTI and DaG but only a small change in CtD; SVD is comparatively stable. This model-by-task heterogeneity is retained rather than collapsed into a universal claim.

Selected mean AUC changes:

- DTI: NeuralMF 0.997 → 0.908; LightGCN 0.989 → 0.883; SVD 0.950 → 0.914.
- CtD: NeuralMF 0.879 → 0.560; LightGCN 0.902 → 0.882; SVD 0.717 → 0.705.
- DaG: NeuralMF 0.863 → 0.640; LightGCN 0.794 → 0.721; SVD 0.620 → 0.604.

The prespecified Gate-2 wording required at least two learned model families with meaningful sensitivity in at least three relation families. We do **not** retroactively weaken that wording.

## H3 / Gate 3a — model-rank instability

Status: **PASS as pilot**.

- DTI conventional winner: NeuralMF; degree-matched winner: SVD.
- Disease–gene conventional winner: NeuralMF; degree-matched winner: LightGCN.
- Compound–disease: LightGCN remains winner under both evaluations and serves as a useful non-reversal control.

## H5 — RIDI-inspired hypothesis-selection identity

Status: **DTI ROBUST; DISEASE–GENE BOUNDARY / INCONCLUSIVE FOR CAUSAL ATTRIBUTION**.

The same frozen candidate universe is ranked by the model selected under conventional evaluation and the model selected under structure-neutralized evaluation.

Initial turnover estimates:

| Family | Candidate pairs | Score Spearman | HT@100 | HT@500 | HT@1000 |
|---|---:|---:|---:|---:|---:|
| DTI | ~817,292 | 0.091 | 1.000 | 0.9996 | 0.9406 |
| Disease–gene | ~614,361 | 0.258 | 0.942 | 0.8988 | 0.8760 |

A fixed-split, fixed-candidate-universe stability control then varied only model initialization/training randomness:

| Family | Comparison | HT@100 | HT@500 | HT@1000 |
|---|---|---:|---:|---:|
| DTI | within NeuralMF | 0.307 | 0.036 | 0.212 |
| DTI | within SVD | 0.138 | 0.136 | 0.155 |
| DTI | conventional winner vs neutralized winner | 1.000 | 0.999 | 0.945 |
| Disease–gene | within NeuralMF | 0.934 | 0.877 | 0.849 |
| Disease–gene | within LightGCN | 0.684 | 0.658 | 0.642 |
| Disease–gene | conventional winner vs neutralized winner | 0.944 | 0.898 | 0.872 |

Interpretation: **DTI provides a clean benchmark-induced scientific-decision result** because cross-selected-model turnover greatly exceeds ordinary within-model initialization turnover. In disease–gene prediction, NeuralMF itself is highly unstable at top-k, so the initial cross-model turnover cannot be attributed mainly to benchmark-induced model selection. This is retained as a boundary result rather than counted as a second robust H5 replication. A confirmatory ensemble-ranking analysis and PPI hypothesis-identity analysis are warranted.

## H4 — temporal external-validity pilot

Status: **PASS as pilot**.

Frozen BioGRID comparison:

- historical snapshot: 5.0.250,
- later snapshot: 5.0.261,
- historical human multi-validated physical network: 93,146 unique edges / 11,844 nodes,
- later-added edges with both endpoints already present historically: 5,635.

Internal historical evaluation:

| Model | Random AUC | Degree-matched AUC | Future AUC against degree-matched persistent-unobserved controls |
|---|---:|---:|---:|
| NeuralMF | 0.931 | 0.552 | 0.549 |
| LightGCN | 0.908 | 0.595 | 0.589 |
| SVD | 0.877 | 0.766 | 0.703 |

Conventional evaluation selected NeuralMF. Structure-neutralized evaluation selected SVD. On later-added edges versus degree-matched persistent-unobserved controls, the neutralized-selected model exceeded the conventional-selected model by **0.154 AUC**.

### H4b — non-circular prospective discovery yield

Status: **PASS and stronger than the matched-control H4 pilot**.

To remove the concern that the future test reused the same degree-balancing principle as model selection, all three models were retrained on the complete historical network and used to rank the **same complete historical non-edge universe**. No degree matching, negative sampling, or future-control construction was used. The frozen universe contained **70,041,100** unknown human protein pairs; **5,635** were subsequently present as BioGRID MV-Physical relations in release 5.0.261.

| Model | Future hits @100 | @1,000 | @10,000 | @50,000 |
|---|---:|---:|---:|---:|
| SVD — neutralized-selected | 7 | 26 | 159 | 522 |
| NeuralMF — conventional-selected | 0 | 5 | 42 | 181 |
| LightGCN | 0 | 5 | 66 | 241 |

For SVD, the corresponding enrichment above a uniform draw from the full candidate universe was **870.1× at 100**, **323.2× at 1,000**, **197.6× at 10,000**, and **129.8× at 50,000**. At 50,000 hypotheses, SVD recovered 9.26% of all later-added closed-world edges versus 3.21% for NeuralMF and 4.28% for LightGCN.

This is the strongest current external-validity result because the future-yield endpoint is not created by the same degree-matching intervention used to select the model. It remains a single PPI temporal domain, and BioGRID additions reflect research and curation processes rather than an unbiased sample of biological truth.

## Evidence-aware Anti-DDI axis

Status: **NEXT HARDENING STEP**.

The frozen supplementary archive contains:

- 8,094 GoldD2-derived ATC5 positive class pairs,
- 538 higher-support Anti-DDI candidate records,
- explicit evidence tiers and contradiction safeguards.

Independent re-execution supports a performance shift when random missing pairs are replaced by curated counter-evidence pairs, but the evidence-state sets differ materially in structural exposure and exact matching coverage is low. Therefore the final analysis will not present forced one-to-one matching as if it removed confounding. Evidence state and structural exposure will be reported as separate axes with explicit balance diagnostics.

## Prior-art boundary

The project does not claim discovery of degree/rich-node/prior bias or invention of temporal biomedical hypothesis benchmarking. Closely related prior work already includes:

- PNAS 2025 bias-aware PPI evaluation showing rich-node bias and temporal BioGRID effects;
- Nature Communications 2025 target-prior bias and causal debiasing in DTI;
- Bioinformatics 2026 biomedical-KG leakage/evaluation work comparing random/cold-start tests with independent evidence;
- Dyport 2024 dynamic/temporal biomedical hypothesis-generation benchmarking.

The intended new contribution is the joint chain:

**cross-domain structural inflation → model-rank reversal → stable hypothesis-selection change where identifiable → non-circular later-evidence consequence**, with evidence-state-aware negative controls where available.

## Current claim ceiling

The current evidence supports a strong claim that benchmark construction can redirect model selection, with robust downstream hypothesis-priority change demonstrated in DTI and one non-circular temporal PPI experiment showing substantially higher later-evidence yield for the structure-neutralized-selected model. A Science-level general claim still requires at least one independent temporal/external replication outside PPI, stronger contemporary models, formal uncertainty for prospective yield, and a second stable hypothesis-identity replication. It does **not** justify a claim that structure-neutralized evaluation universally improves prospective discovery.