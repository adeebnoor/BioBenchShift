# BioBenchShift

## Benchmark design can reshape biomedical hypothesis prioritization

[![Scientific Integrity CI](https://github.com/adeebnoor/BioBenchShift/actions/workflows/ci.yml/badge.svg)](https://github.com/adeebnoor/BioBenchShift/actions/workflows/ci.yml)

BioBenchShift is the companion reproducibility repository for a study of how biomedical benchmark controls and score resolution can affect model selection and computational hypothesis prioritization. The contribution is conditional: changing the control population can change the selected model and its ranked hypotheses, but this does not occur under every matching rule or model menu.

> benchmark control definition -> model scores / selection -> ranked hypotheses -> later recorded evidence

This repository does not establish that these rankings were used to choose actual laboratory experiments, that an algorithm is generally superior, or that later database additions measure unbiased discovery yield.

## Evidence and boundaries

- Endpoint degree alone discriminated observed from randomly sampled unknown pairs across four relation families. Degree matching reduced that signal. In DTI, AUROC remained 0.620 rather than reaching chance; the other three reported matched AUROCs were 0.513, 0.519 and 0.513.
- The original log2 degree-bin control changed the selected model within the SVD/NeuralMF/LightGCN comparisons in DTI and historical PPI. In the DTI matching-rule sensitivity, log2 matching reversed the winner in 7/10 split seeds, whereas degree-decile and caliper matching did not reverse it in those ten seeds.
- For the original SVD-versus-NeuralMF contrasts, independently initialized ensembles differed on 100% of DTI top-100 and 99% of PPI top-100 hypotheses.
- On the same 70,041,100-pair historical BioGRID candidate universe, SVD recovered 522 of 5,635 later-added interactions within its top 50,000, compared with 181 for NeuralMF. Later curation, historical degree and assay visibility constrain the biological interpretation; a popularity-matched temporal null is a separate analysis, not already established by these hit counts.
- The completed fair four-model comparison across five split seeds retained GraphBAN as the mean-AUROC winner under both control definitions: 0.997122 conventionally and 0.927247 after log2 matching. GraphBAN won under both rules in four seeds; the remaining seed selected NeuralMF conventionally and SVD after matching. There was no aggregate winner reversal in this expanded panel.
- A later, explicitly secondary same-fit replay asked whether the close conventional scores of GraphBAN and NeuralMF resolve the same candidate queue. On the common 1,005,757-pair mapped DTI universe, their top-100 queues had zero overlap in all five split seeds (HT@100 = 1.000); mean HT@500 was 0.9632 and mean HT@1,000 was 0.8638. The mean absolute AUROC gap in those same replayed fits was 0.000939 under conventional controls and 0.032213 under log2-matched controls. This is descriptive score closeness, not a statistical-equivalence test.
- The secondary replay reproduced the frozen split-level AUROCs for SVD, NeuralMF and LightGCN to numerical precision but not for GraphBAN. The new GraphBAN execution is therefore kept separate from the frozen five-seed model-selection panel; frozen GraphBAN AUROCs are not substituted into the secondary queue analysis.
- On the common 1,005,757-pair mapped DTI universe, the sequential and parallel GraphBAN ChEMBL follow-through executions also differed numerically. Both outputs remain visible. Neither provides a benchmark-induced model switch because the aggregate frozen selector chose GraphBAN under both rules.

A small AUROC difference is not, by itself, evidence of statistical equivalence. The modern GraphBAN-NeuralMF queue result is a secondary analysis: the pair was nominated after the existing performance table was known, all six model pairs were retained in the analysis output, and the observed DTI relationship is not presented as a universal benchmark-resolution threshold.

## Immutable evidence states

| State | Commit |
|---|---|
| Baseline analyses | `1298e4c83ff4fa488a784d52219d20e706118a9e` |
| Additional matching, case and biological-program analyses | `d58c816ccbe30226247087009ff0d22ef44c608a` |
| Completed modern-model panel and both ChEMBL executions | `66b5c050739bf912f531f0de354925824122eedc` |
| Secondary A1/A2 same-fit replay | `71bfb714de8cea73fd84bc763bcb4875a9ae7c43` on `science-a1-20260915` |

The secondary replay does not replace the prespecified evidence states. Branch names are not immutable archives; the commit hashes above are the reproducibility anchors.

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
| Secondary same-fit queue resolution | commit `71bfb714de8cea73fd84bc763bcb4875a9ae7c43`: `analysis/run_science_queue_resolution_a1_replay.py`, `results/science_queue_resolution_a1_*` |

The original fair-panel script imports saved GraphBAN metrics; it does not preserve per-pair GraphBAN predictions. The secondary replay therefore trains and evaluates each fitted model and ranks the common candidate universe within the same execution, while auditing the replay against the frozen performance table.

## Reproducibility policy

1. Retain frozen protocols and distinguish later, secondary analyses from original prospective specifications.
2. Keep negative, incomplete-balance and execution-sensitive results visible.
3. Keep external evidence separate from the original model-selection rule.
4. Record candidate identity, mapping, training split, model-fit identity, software environment and prediction-file hashes when comparing queues.
5. Do not infer statistical equivalence from a small or nonsignificant performance difference.
6. Do not combine metrics from one model execution with rankings from a different execution.
7. Do not claim submission readiness, independent validation or prospective discovery yield from a passing software check alone.

The independent [ANTI-DDI](https://github.com/adeebnoor/ANTI-DDI) repository remains the resource of record for Anti-DDI. It is used here only in a secondary evidence-state sensitivity analysis. Submission correspondence and reviewer suggestions are not part of this public repository.

## Citation and license

Until an archival release and article DOI are available, cite the exact repository commit used. `CITATION.cff` and `.zenodo.json` provide machine-readable metadata; consult the evidence state appropriate to the analysis.

Original code is MIT-licensed. Original documentation, figure source data and frozen outputs are CC BY 4.0 unless otherwise stated. Third-party resources retain their original terms; see `DATA_PROVENANCE.md`.
