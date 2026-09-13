# Manuscript companion

## Benchmark design redirects biomedical discovery

**Author:** Adeeb Noor  
**Affiliation:** Department of Information Technology, Faculty of Computing and Information Technology, King Abdulaziz University, Jeddah, Saudi Arabia  
**Status:** pre-submission research manuscript; repository evidence is frozen except for the explicitly prespecified contemporary-model completion gate.

## Abstract

Biomedical AI increasingly determines which molecular relations are prioritized for experiments, yet models are often selected on benchmarks that contrast known relations with sampled unknown pairs. Across four relation families, structural observability alone produced strong discrimination and affected learned models unequally, reversing the selected model in drug–target and protein-interaction analyses. These reversals changed 100% and 99% of top-100 hypotheses, respectively, far beyond within-model instability. In a frozen historical BioGRID network, the model selected after structural neutralization recovered 522 of 5,635 later-added interactions among its top 50,000 predictions versus 181 for the conventional winner. Independent ChEMBL evidence was also more concentrated at the earliest drug–target cutoffs, although the advantage did not persist at broad cutoffs. Thus benchmark design can redirect biomedical discovery by changing which model wins and which biology is tested first.

## Public repository / submission-package boundary

This repository contains the scientific reproducibility record: code, protocols, frozen results, source data, generated figures, provenance, checksums and citation metadata.

The following journal-submission materials are deliberately maintained outside the public repository:

- formatted manuscript submission file;
- supplementary-materials submission file;
- cover letter;
- reviewer suggestions and conflict checks;
- submission-system metadata;
- editorial-positioning notes.

This separation follows the role of a publication companion repository: enable independent inspection and reproduction of the scientific evidence without publishing private editorial correspondence.

## Evidence map

See:

- `REPRODUCIBILITY.md` for claim → code → frozen-output mapping;
- `DATA_PROVENANCE.md` for data releases, hashes and external-evidence boundaries;
- `results/README.md` for the result index;
- `protocols/` for prespecified rules and timestamped amendments;
- `figures/source_data_main.csv` for main-figure source data.