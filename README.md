# BioBenchShift

## Benchmark design can reshape biomedical hypothesis prioritization

[![Scientific Integrity CI](https://github.com/adeebnoor/BioBenchShift/actions/workflows/ci.yml/badge.svg)](https://github.com/adeebnoor/BioBenchShift/actions/workflows/ci.yml)

BioBenchShift is the companion reproducibility repository for a study of how biomedical benchmark controls can affect model selection and computational hypothesis prioritization. The contribution is conditional: changing the control population can change the selected model and its ranked hypotheses, but this does not occur under every matching rule or model menu.

> benchmark control definition -> selected model -> ranked hypotheses -> later recorded evidence

This repository does not establish that these rankings were used to choose actual laboratory experiments, that an algorithm is generally superior, or that later database additions measure unbiased discovery yield.

## Evidence and boundaries

- Endpoint degree alone discriminated observed from randomly sampled unknown pairs across four relation families. Degree matching reduced that signal. In DTI, AUROC remained 0.620 rather than reaching chance; the other three reported matched AUROCs were 0.513, 0.519 and 0.513.
- The original log2 degree-bin control changed the selected model within the SVD/NeuralMF/LightGCN comparisons in DTI and historical PPI. In the DTI matching-rule sensitivity, log2 matching reversed the winner in 7/10 split seeds, whereas degree-decile and caliper matching did not reverse it in those ten seeds.
- For the original SVD-versus-NeuralMF contrasts, independently initialized ensembles differed on 100% of DTI top-100 and 99% of PPI top-100 hypotheses. These values must not be attributed to a GraphBAN-versus-NeuralMF comparison or a different candidate universe.
- On the same 70,041,100-pair historical BioGRID candidate universe, SVD recovered 522 of 5,635 later-added interactions within its top 50,000, compared with 181 for NeuralMF. Later curation, historical degree and assay visibility constrain the biological interpretation; a popularity-matched temporal null is a separate analysis, not already established by these hit counts.
- The completed fair four-model comparison across five split seeds retained GraphBAN as the mean-AUROC winner under both control definitions: 0.997122 conventionally and 0.927247 after log2 matching. GraphBAN won under both rules in four seeds; the remaining seed selected NeuralMF conventionally and SVD after matching. There was no aggregate winner reversal in this expanded panel.
- On the common 1,005,757-pair mapped DTI universe, the sequential and parallel GraphBAN ChEMBL follow-through executions differed numerically. Both outputs remain visible. Neither provides a benchmark-induced model switch because the aggregate selector chose GraphBAN under both rules.

A small difference in average AUROC is not, by itself, evidence of statistical equivalence. Neither the original queue-turnover values nor the aggregate score tables establish disjoint queues between GraphBAN and NeuralMF. That proposed secondary comparison requires identity-aligned candidate predictions and evaluation predictions from the corresponding fitted models. Full-history ranking artifacts must not be relabelled as the five held-out split fits.

## Immutable evidence states

| State | Commit |
|---|---|
| Baseline analyses | `1298e4c83ff4fa488a784d52219d20e706118a9e` |
| Additional matching, case and biological-program analyses | `d58c816ccbe30226247087009ff0d22ef44c608a` |
| Completed modern-model panel and both ChEMBL executions | `66b5c050739bf912f531f0de354925824122eedc` |

The expansion and main branches are separate evidence lineages. A branch name is not an immutable archive. Documentation updates do not alter the frozen results above.

## Repository layout

- `analysis/`: executable analysis scripts.
- `results/`: frozen replicate-level and summary outputs.
- `figures/`: figure files and source data; use the source and evidence lineage appropriate to the manuscript version.
- `protocols/`: analysis specifications, amendments, checksums and provenance.
- `REPRODUCIBILITY.md` and `DATA_PROVENANCE.md`: commands, inputs, resource releases and interpretation boundaries.

## Reproduction entry points

| Analysis | Script or output |
|---|---|
| Structural-only audits | `analysis/run_biosnap_dti_gate1.py`, `analysis/run_huri_ppi_gate1.py`, `analysis/run_hetionet_ctd_gate1.py`, `analysis/run_hetionet_bipartite_gate1.py` |
| Baseline models | `analysis/run_biosnap_dti_gate2_models.py`, `analysis/run_hetionet_gate2_models.py`, `analysis/run_lightgcn_gate2.py` |
| Original hypothesis-identity controls | `analysis/run_h5_ensemble_confirmatory.py`, `analysis/run_biogrid_h5_identity.py` |
| BioGRID later-evidence results | `analysis/run_biogrid_future_yield.py`, `analysis/run_biogrid_future_uncertainty.py` |
| Feature-rich model challenge | `analysis/run_graphban_targetdecagon_clean.py`, `analysis/run_graphban_targetdecagon_precomputed.py` |
| Completed fair comparison | `results/dti_fair_model_selection_panel_5seed_replicates.csv`, `results/dti_fair_model_selection_panel_5seed_summary.json` |
| Both ChEMBL executions | `results/graphban_chembl_followthrough_primary_table.csv`, `results/graphban_chembl_followthrough_parallel_primary_table.csv` |

The fair-panel script imports saved GraphBAN metrics; it does not train GraphBAN or preserve its per-pair predictions. A complete score-export replay therefore needs both model-training entry points, not a four-line addition to the fair-panel summary writer alone.

## Reproducibility policy

1. Retain frozen protocols and distinguish later, secondary analysis plans from original prospective specifications.
2. Keep negative, incomplete-balance and execution-sensitive results visible.
3. Keep external evidence separate from the original model-selection rule.
4. Record candidate identity, mapping, training split, model-fit identity, software environment and prediction-file hashes when comparing queues.
5. Do not infer statistical equivalence from a nonsignificant difference or from overlapping split variation.
6. Do not claim submission readiness, independent validation or prospective discovery yield from a passing software check alone.

The independent [ANTI-DDI](https://github.com/adeebnoor/ANTI-DDI) repository remains the resource of record for Anti-DDI. It is used here only in a secondary evidence-state sensitivity analysis. Submission correspondence and reviewer suggestions are not part of this public repository.

## Citation and license

Until an archival release and article DOI are available, cite the exact repository commit used. `CITATION.cff` and `.zenodo.json` provide machine-readable metadata; consult the immutable evidence state appropriate to the analysis.

Original code is MIT-licensed. Original documentation, figure source data and frozen outputs are CC BY 4.0 unless otherwise stated. Third-party resources retain their original terms; see `DATA_PROVENANCE.md`.
