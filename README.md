# BioBenchShift

## Benchmark design redirects biomedical discovery

**Pre-submission reproducibility repository for a planned _Science_ Research Article.**

BioBenchShift tests a simple but consequential question: **can the benchmark used to select a biomedical AI model change which model wins, which biological hypotheses are prioritized, and which later or independent evidence is concentrated among those priorities?**

> **benchmark construction → model identity → hypothesis identity → later / independent evidence**

### Headline findings

- Structural observability alone strongly discriminates observed from randomly sampled unknown relations across DTI, PPI, compound–disease and disease–gene benchmarks, and largely collapses after structural neutralization.
- Benchmark choice reverses the selected model in drug–target and historical protein-interaction analyses.
- The selected-model ensembles disagree on **100% of DTI top-100** and **99% of PPI top-100** hypotheses, far beyond within-model instability.
- On the same frozen **70,041,100-pair** historical BioGRID universe, the structure-neutralized-selected model recovers **522** later-added interactions in its top 50,000 predictions versus **181** for the conventional winner.
- In independent ChEMBL 37 evidence, the structure-neutralized-selected model concentrates more supported DTI relations at the earliest experimental cutoffs, while the broad-cutoff crossover is retained as an explicit boundary.

### Repository map

- `paper/` — manuscript, cover letter, supplementary materials and submission metadata.
- `analysis/` — executable analyses for structural audits, learned models, hypothesis turnover, BioGRID temporal evidence, ChEMBL validation and contemporary-model robustness.
- `results/` — frozen seed-level and summary outputs.
- `figures/` — figure source data and generated main figures.
- `protocols/` — prespecified gates, amendments, provenance and novelty boundaries.
- `.github/workflows/` — reproducible execution and submission-package builds.

### Scientific boundary

This work does **not** claim that degree bias is new, that biomedical AI generally learns no biology, that structural neutralization is universally optimal, that SVD is universally superior, or that later database additions are unbiased biological truth. The tested contribution is narrower: **benchmark design can act upstream of scientific discovery by changing model selection and the identity of the hypotheses that enter the experimental queue.**

### Related work

The independent [`ANTI-DDI`](https://github.com/adeebnoor/ANTI-DDI) repository remains the resource of record for the Anti-DDI v3.0.1 evidence-state dataset. BioBenchShift uses Anti-DDI only as a supportive evidence-state boundary; the Science manuscript's principal evidence is cross-domain DTI/PPI/CtD/DaG benchmarking, fixed-universe hypothesis turnover, temporal BioGRID evidence and independent ChEMBL evidence.

### Status

Core evidence is frozen. The prespecified leakage-free contemporary GraphBAN challenge is being completed and will be reported regardless of direction before submission. This repository is a pre-submission research package and is not a statement of peer review, acceptance or publication.

### Citation

Until the final archival release is minted, cite the exact BioBenchShift commit used. A `CITATION.cff` and Zenodo metadata file are included in this repository.