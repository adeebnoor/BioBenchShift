# Reproducibility guide

This document maps each manuscript-level claim to executable code and frozen outputs. It is intentionally separate from journal correspondence and editorial materials.

## Environment

Core analyses use Python 3.11.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The contemporary GraphBAN robustness challenge requires the additional packages in `requirements-graphban.txt` and is computationally heavier because it materializes ChemBERTa and ESM-1b features.

## Claim-to-code map

| Claim | Code | Frozen evidence |
|---|---|---|
| Structural observability inflates conventional evaluation across four biomedical relation families | `analysis/run_biosnap_dti_gate1.py`; `analysis/run_huri_ppi_gate1.py`; `analysis/run_hetionet_ctd_gate1.py`; `analysis/run_hetionet_bipartite_gate1.py` | `results/biosnap_dti_gate1*`; `results/huri_ppi_gate1*`; `results/hetionet_*_gate1*`; `results/GATE1_SUMMARY.md` |
| Learned models respond unequally to structural neutralization and can reverse rank | `analysis/run_biosnap_dti_gate2_models.py`; `analysis/run_hetionet_gate2_models.py`; `analysis/run_lightgcn_gate2.py` | `results/*gate2*`; `results/*lightgcn*`; `results/GATE2_RANK_REVERSAL_PILOT.md` |
| Benchmark-selected models produce different prioritized hypotheses beyond initialization noise | `analysis/run_h5_hypothesis_turnover.py`; `analysis/run_h5_stability_control.py`; `analysis/run_h5_ensemble_confirmatory.py`; `analysis/run_biogrid_h5_identity.py` | `results/h5_*`; `results/h5b_*`; `results/h5c_*` |
| Neutralized-selected historical PPI model recovers more later-added BioGRID relations | `analysis/run_biogrid_temporal_h4.py`; `analysis/run_biogrid_future_yield.py` | `results/h4_biogrid_temporal*`; `results/h4b_biogrid_future_yield*` |
| BioGRID future-yield difference is robust to paired uncertainty | `analysis/run_biogrid_future_uncertainty.py` | `results/h4c_biogrid_future_uncertainty*` |
| Independent ChEMBL evidence is more concentrated at the earliest DTI priorities after structural neutralization, with a broad-cutoff boundary | `analysis/run_chembl_dti_external.py` | `results/h6_chembl_dti_external*` |
| Evidence-state sensitivity is model-dependent and incompletely balanced | evidence-state overlap analysis archived in results/protocols | `results/antiddi_evidence_overlap*`; `protocols/EVIDENCE_STATE_OVERLAP_PROTOCOL.md` |
| Contemporary feature-rich architecture remains benchmark-sensitive under a leakage-free protocol | `analysis/run_graphban_targetdecagon_clean.py`; `analysis/run_graphban_targetdecagon_precomputed.py` | `results/graphban_targetdecagon_mapping*`; `results/graphban_targetdecagon_clean*`; `protocols/GRAPHBAN_PROTOCOL_AMENDMENT_20260913.md`; `protocols/GRAPHBAN_EXECUTION_PARALLELIZATION_20260913.md` |

## Frozen headline values

### Structural nulls

- BioSNAP DTI degree-only AUROC: random **0.983 ± 0.001**, structure-matched **0.620 ± 0.012**.
- HuRI PPI: **0.928 ± 0.002 → 0.513 ± 0.001**.
- Hetionet compound–disease: **0.871 ± 0.019 → 0.519 ± 0.010**.
- Hetionet disease–gene: **0.876 ± 0.006 → 0.513 ± 0.002**.

### Hypothesis identity

- DTI selected-model ensemble turnover at top 100: **1.000**.
- Historical PPI selected-model ensemble turnover at top 100: **0.990**.
- These cross-model differences materially exceed within-model initialization instability; see frozen H5b/H5c outputs.

### Later BioGRID evidence

Historical candidate universe: **70,041,100** unknown protein pairs. Later-added relations with both endpoints historically observed: **5,635**.

At top 50,000 predictions:

- structure-neutralized-selected SVD: **522** later-added relations;
- conventional-selected NeuralMF: **181** later-added relations;
- recall difference: **0.0605**, paired-bootstrap 95% CI **[0.0531, 0.0680]**.

### Independent ChEMBL 37 evidence

Under the frozen primary evidence rule (human, direct single-protein binding, confidence score 9, pChEMBL ≥6), **66** supported BioSNAP candidate relations were identified. At top 100 / 500 / 1,000 priorities, the structure-neutralized-selected SVD recovered **3 / 5 / 6** versus **0 / 0 / 2** for NeuralMF. The advantage does not persist at all broad cutoffs and that crossover is part of the reported boundary.

### Leakage-free contemporary GraphBAN challenge

The prespecified contemporary-model completion used a GraphBAN-style heterogeneous GraphSAGE encoder with ChemBERTa drug features and ESM-1b protein features. Mapping retained **18,631 / 18,690 (99.7%)** TargetDecagon positive edges. Held-out positive edges never entered message passing. Across seeds 0, 1 and 2:

- conventional random-unlabelled AUROC: **0.9972 ± 0.0002**;
- structure-neutralized AUROC: **0.9278 ± 0.0141**;
- mean AUROC difference: **0.0694**;
- conventional AUPRC: **0.9976**;
- structure-neutralized AUPRC: **0.9425**;
- mean AUPRC difference: **0.0551**;
- mean matching coverage: **0.5880**.

The model therefore retained high absolute performance while showing a material evaluation-regime shift. This supports architectural generality of benchmark sensitivity without supporting a claim that the contemporary model contains no biological signal.

## Figure rebuild

All four main figures are generated from one public source-data table:

```bash
python analysis/build_science_figures.py
```

Expected outputs are written to `figures/generated/`.

## Integrity rules

- Do not retune thresholds using BioGRID later-release or ChEMBL evidence.
- Do not replace negative/boundary results after observing them.
- Preserve the original prespecified gate document and timestamped amendments.
- Report the exact repository commit and source-data release identifiers.
- Do not interpret unobserved relation pairs as experimentally verified negatives.