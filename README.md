# BioBenchShift

## Benchmark design redirects biomedical discovery

[![Scientific Integrity CI](https://github.com/adeebnoor/BioBenchShift/actions/workflows/ci.yml/badge.svg)](https://github.com/adeebnoor/BioBenchShift/actions/workflows/ci.yml)

**BioBenchShift is the companion reproducibility repository for the manuscript _Benchmark design redirects biomedical discovery_.**

Biomedical AI benchmarks do more than measure models: they can change which model wins, which biological hypotheses are prioritized, and which hypotheses reach the front of the experimental queue.

> **benchmark construction → model selection → hypothesis identity → later / independent evidence**

### Main findings

- Structural observability alone strongly separates observed relations from randomly sampled unknown relations across drug–target, protein–protein, compound–disease and disease–gene benchmarks, and largely collapses after structural neutralization.
- Benchmark construction changes model ranking in drug–target and historical protein-interaction analyses.
- Stability-controlled selected-model ensembles disagree on **100% of DTI top-100** and **99% of PPI top-100** hypotheses, far beyond within-model initialization variability.
- On the same frozen **70,041,100-pair** historical BioGRID candidate universe, the structure-neutralized-selected model recovers **522** of **5,635** later-added interactions in its top 50,000 predictions, compared with **181** for the conventional winner.
- In independent ChEMBL 37 evidence, support is more concentrated among the earliest drug–target priorities selected after structural neutralization; the broad-cutoff crossover is retained as an explicit boundary rather than hidden.

## Repository layout

```text
BioBenchShift/
├── analysis/       # executable analyses used in the study
├── figures/        # figure source data and reproducibly generated figures
├── results/        # frozen seed-level and summary outputs
├── protocols/      # prespecified analyses, amendments, checksums and provenance
├── CITATION.cff
├── DATA_PROVENANCE.md
├── REPRODUCIBILITY.md
├── requirements.txt
└── LICENSE
```

Submission correspondence, reviewer suggestions and journal-system metadata are intentionally **not** part of the public reproducibility repository.

## Quick start

```bash
git clone https://github.com/adeebnoor/BioBenchShift.git
cd BioBenchShift
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python analysis/build_science_figures.py
```

The command above rebuilds the four manuscript figures from `figures/source_data_main.csv`. The full analysis ladder uses external public datasets and, for the contemporary GraphBAN robustness analysis, an additional environment documented in `requirements-graphban.txt`.

## Reproduce the headline evidence

| Scientific question | Primary script | Frozen output |
|---|---|---|
| Does structural observability inflate conventional evaluation? | `analysis/run_biosnap_dti_gate1.py`, `analysis/run_huri_ppi_gate1.py`, `analysis/run_hetionet_ctd_gate1.py`, `analysis/run_hetionet_bipartite_gate1.py` | `results/*gate1*` |
| Do learned-model rankings change? | `analysis/run_biosnap_dti_gate2_models.py`, `analysis/run_hetionet_gate2_models.py`, `analysis/run_lightgcn_gate2.py` | `results/*gate2*`, `results/*lightgcn*` |
| Does benchmark-induced model selection change hypothesis identity? | `analysis/run_h5_ensemble_confirmatory.py`, `analysis/run_biogrid_h5_identity.py` | `results/h5b_*`, `results/h5c_*` |
| Do selected models differ on later evidence? | `analysis/run_biogrid_future_yield.py`, `analysis/run_biogrid_future_uncertainty.py` | `results/h4b_*`, `results/h4c_*` |
| Does the pattern transfer to an independent DTI evidence source? | `analysis/run_chembl_dti_external.py` | `results/h6_chembl_*` |
| Does the conclusion survive a contemporary architecture challenge? | `analysis/run_graphban_targetdecagon_clean.py` | frozen S6 output when the prespecified run completes |

For exact release identifiers, hashes, expected outputs, and commands, see [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) and [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md).

## Reproducibility policy

1. **Protocols are frozen before outcome inspection.** Amendments are timestamped and retained.
2. **Negative and boundary results remain visible.** Disease–gene instability, the ChEMBL broad-cutoff crossover, and incomplete Anti-DDI overlap balance are not removed from the record.
3. **External evidence is not used for model selection.** BioGRID later-release evidence and ChEMBL 37 evidence are opened only after the relevant model-selection rule is fixed.
4. **Third-party datasets are not silently repackaged.** Release identifiers and checksums are provided where redistribution is restricted or unnecessary.
5. **The same candidate universe is used when testing hypothesis identity.** Only the benchmark rule used to select the model changes.

## Scope and boundaries

BioBenchShift does **not** claim that degree bias is newly discovered, that biomedical AI generally learns no biology, that structural neutralization is universally optimal, that one model family is universally superior, or that later database additions are unbiased biological truth. The tested contribution is narrower and falsifiable: **benchmark design can operate upstream of biomedical discovery by changing model selection and the identity of the hypotheses prioritized for follow-up.**

The independent [`ANTI-DDI`](https://github.com/adeebnoor/ANTI-DDI) repository remains the resource of record for Anti-DDI v3.0.1. BioBenchShift uses that resource only for a supportive evidence-state sensitivity analysis; the principal evidence in this repository is the cross-domain benchmark, hypothesis-identity, temporal BioGRID, and independent ChEMBL analyses.

## Status

The core cross-domain, hypothesis-identity, BioGRID temporal, uncertainty, and ChEMBL analyses are frozen. One prespecified leakage-free contemporary-model robustness run is retained as an explicit completion gate; its outcome will be reported regardless of direction.

## Citation

Until an archival release and article DOI are available, cite the exact repository commit used. `CITATION.cff` and `.zenodo.json` provide machine-readable metadata.

## License

Original BioBenchShift code is MIT-licensed. Original documentation, figure source data, and frozen outputs are CC BY 4.0 unless a file states otherwise. Third-party datasets and identifiers remain subject to their original terms; see `DATA_PROVENANCE.md`.