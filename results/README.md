# Frozen results index

This directory stores immutable or versioned outputs produced by the study analyses. Manuscript claims should point to these files rather than to ad hoc notebook state.

## Structural audits

- `GATE1_SUMMARY.md` — cross-domain structural-null summary.
- `biosnap_dti_gate1*` — drug–target structural audit.
- `huri_ppi_gate1*` — HuRI protein–protein audit.
- `hetionet_ctd_gate1*` — compound–disease audit.
- `hetionet_dag_gate1*` — disease–gene audit.

## Learned models and rank changes

- `GATE2_RANK_REVERSAL_PILOT.md`
- `biosnap_dti_gate2*`
- `biosnap_dti_lightgcn*`
- `hetionet_*_gate2*`
- `hetionet_*_lightgcn*`

## Hypothesis identity

- `h5_hypothesis_turnover*` — original turnover analysis.
- `h5_stability_control*` — within-model training-noise control.
- `h5b_ensemble_confirmatory*` — ensemble confirmatory DTI/disease–gene analysis.
- `h5c_biogrid_ppi_identity*` — ensemble PPI hypothesis-identity analysis.

## Temporal and independent evidence

- `h4_biogrid_temporal*` — historical model selection and initial temporal evaluation.
- `h4b_biogrid_future_yield*` — non-circular complete-universe later-evidence analysis.
- `h4c_biogrid_future_uncertainty*` — paired uncertainty for later-evidence differences.
- `h6_chembl_dti_external*` — independent ChEMBL 37 DTI evidence.

## Evidence-state sensitivity

- `antiddi_evidence_overlap*` — supportive overlap-standardized sensitivity. Because the prespecified balance criterion passed in only 5/10 seeds, this result is a boundary/sensitivity analysis rather than a principal causal claim.

## Contemporary-model challenge

- `graphban_targetdecagon_mapping*` — frozen mapping gate and feature coverage.
- Final S6 outputs are added only after the leakage-free prespecified run completes. No interim performance outcome is inferred from runtime status.

## Input integrity

Files ending in `.sha256` record frozen inputs or acquisition boundaries used by the corresponding analyses. See `../DATA_PROVENANCE.md` and `../protocols/CHECKSUMS.sha256`.