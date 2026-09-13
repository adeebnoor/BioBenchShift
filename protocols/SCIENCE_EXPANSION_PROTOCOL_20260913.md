# Science Expansion Protocol — 13 September 2026

## Purpose

This branch extends, but does not overwrite, the frozen submission state at commit `1298e4c83ff4fa488a784d52219d20e706118a9e`. The objective is to test whether the benchmark-to-model-to-hypothesis-to-evidence chain survives stronger external, biological, and contemporary-model challenges suitable for a Science-only submission strategy.

No result from the analyses below has been inspected before this protocol is committed. Any null, mixed, or adverse result will be retained and reported.

## Baseline that remains immutable

The following frozen findings are treated as prior results, not re-estimated to optimize the expansion:

- structural observability audit across DTI, PPI, compound–disease, and disease–gene relations;
- controlled model ladder using SVD, NeuralMF, and LightGCN;
- stability-controlled DTI and PPI hypothesis turnover;
- BioGRID 5.0.250 -> 5.0.261 later-evidence recovery;
- ChEMBL 37 independent DTI evidence;
- Anti-DDI boundary analysis;
- leakage-free GraphBAN-style benchmark-sensitivity challenge.

## Expansion A — independent direct DTI evidence from BindingDB

### Source freeze

Use the dated **BindingDB 2026-09 release** and, for the primary analysis, only the file explicitly described by BindingDB as **data curated from articles by BindingDB** (`BindingDB_BindingDB_Articles_202609_tsv.zip`). Do not use the BindingDB subset drawn from ChEMBL in the primary analysis. Record the exact file checksum.

### Mapping

- BioSNAP drug identifiers: STITCH/PubChem CID -> BindingDB `PubChem CID`.
- BioSNAP gene identifiers: NCBI Entrez Gene -> reviewed human UniProt accession -> BindingDB SwissProt target accession.
- Restrict the primary evidence set to Homo sapiens, single-protein targets when these fields are available.
- Exclude all BioSNAP-positive drug–target pairs from the external candidate evidence set.

### Primary evidence definition

A candidate pair is externally supported when at least one directly curated BindingDB article record satisfies all of the following:

1. candidate pair was not positive in the frozen BioSNAP benchmark;
2. direct target mapping succeeds for both endpoints;
3. at least one quantitative binding measurement is reported as Ki, Kd, IC50, or EC50;
4. the strongest qualifying measurement is <= 1,000 nM after relation-aware parsing;
5. the publication date is available.

A stricter sensitivity threshold of <= 100 nM is frozen before outcome inspection.

### Temporal case-study gate

For the named case-study table, require publication date >= 2018-01-01 so that the supporting article post-dates the benchmark era used in the frozen study. Case studies are selected mechanically rather than manually:

- enumerate all externally supported pairs among the top 1,000 candidates of each benchmark-selected model;
- identify pairs present in one selected model's top 1,000 but absent from the competing selected model's top 1,000;
- order by the supporting model rank, then by stronger affinity, then lexical pair identifier;
- report up to the first 10 pairs per direction, with drug name, target name, rank, affinity, publication year, PMID/DOI when available.

No pair may be substituted because it appears more biologically interesting.

### Primary endpoints

For SVD, NeuralMF, and LightGCN, report externally supported hits, precision lower bound, and recall of the externally supported candidate set at K = 100, 500, 1,000, 5,000, 10,000, and 50,000.

The confirmatory contrast remains SVD minus NeuralMF because these are the structure-neutralized-selected and conventional-selected models in the frozen DTI ladder.

### Degree/popularity control

For K = 100, 500, and 1,000, compare the observed BindingDB-supported hit count against 20,000 random candidate sets matched to the selected list on the joint log2 endpoint-degree bins computed from the frozen BioSNAP training graph. Report the empirical two-sided tail probability and the standardized enrichment over the matched null. This tests whether external support concentration is explainable by endpoint popularity alone.

## Expansion B — biological consequence of different hypothesis queues

The goal is not to assert that pathway enrichment validates a drug–target pair. It tests whether the benchmark-selected queues direct attention toward different biological programs.

### Target-set definitions

For each selected model, derive unique target-gene sets from Top-100, Top-500, and Top-1,000 DTI candidate pairs. Use the full mapped BioSNAP target set as the background universe.

### Pathway analysis

Run Reactome and Gene Ontology Biological Process enrichment with Benjamini–Hochberg FDR correction using a reproducible API or package version. Report all terms passing FDR < 0.05 and do not curate terms post hoc.

Compare the two selected-model enrichment profiles using:

- number of significant terms unique to each queue;
- overlap coefficient of significant terms;
- Jensen–Shannon divergence of normalized enrichment-score profiles over the union of tested terms.

### Disease relevance analysis

Use one frozen Open Targets Platform release (prefer 26.06 if available) to annotate target–disease evidence. Report the distribution of high-level disease areas represented by targets in each queue and quantify divergence between selected-model queues. Open Targets is interpretive context, not pair-level validation.

## Expansion C — contemporary model ladder

A single GraphBAN model cannot by itself produce a benchmark-induced hypothesis-turnover claim because turnover in this study is downstream of **model selection**. A modern turnover test therefore requires a ladder of multiple modern architectures evaluated under the same frozen benchmark rules.

### Architectures

Primary modern ladder:

1. existing leakage-free GraphBAN-style model with ChemBERTa + ESM-1b features;
2. GS-DTI (2025), using the public implementation and its graph/ESM-based representations;
3. SP-DTI (2025), using the public BioSNAP-compatible implementation and frozen processed structural inputs when available.

DTIAM (2025) is a prespecified replacement only if one of the primary external implementations cannot be executed reproducibly on the frozen BioSNAP universe without altering its defining architecture.

### Common evaluation contract

- same frozen BioSNAP positive edge universe;
- same train/test split seeds for all architectures;
- held-out positives excluded from message passing or training graph construction;
- same conventional random-unlabelled controls and same structure-neutralized endpoint-degree matching rule;
- 10 independent seeds where computationally feasible; if an external architecture is deterministically pretrained and only the downstream head is stochastic, repeat the stochastic portion 10 times;
- report AUROC, AUPRC, precision@K, recall@K, and top-tail calibration summary for K = 100, 500, 1,000;
- report matching coverage for every structure-neutralized evaluation.

### Modern winner-reversal gate

Rank architectures by the prespecified primary metric under conventional evaluation and under structure-neutralized evaluation. Only if the winning architecture changes is a modern hypothesis-turnover analysis opened.

If the winner changes:

- fit 5-member score ensembles for the two selected architectures on a fixed split;
- rank the identical candidate universe;
- report HT@100, HT@500, HT@1,000 and rank-biased overlap;
- compare cross-selected turnover with within-architecture leave-one-initialization-out turnover.

If the winner does not change, report the null result and do not manufacture a turnover contrast.

## Expansion D — matching-method robustness

Repeat the principal DTI structural audit and model-selection comparison under three prespecified neutralization variants:

1. original joint log2 degree-bin matching;
2. finer joint degree quantile matching with fixed quantile edges learned from training endpoints;
3. nearest-neighbor matching on standardized log1p endpoint degrees with a fixed caliper.

For each method report coverage, model ranking, and whether the conventional-to-neutralized winner identity changes. No method is designated primary after results are observed; the original frozen method remains the original primary analysis.

## Expansion E — explicit limitations and decision rule

The expanded manuscript must include a main-text `Limitations` subsection stating at minimum:

- computational/retrospective evidence is not a substitute for prospective wet-lab validation;
- database additions reflect research and curation processes;
- external biochemical databases can share literature provenance despite source-level exclusions;
- pathway and disease-area analyses demonstrate different biological prioritization, not causal correctness of individual pairs;
- the modern architecture ladder can return a null winner-reversal result;
- structural neutralization is a diagnostic intervention, not a universally optimal benchmark.

## Wet-lab validation boundary

No experimental result will be claimed without a real experimental collaborator, assay protocol, raw measurements, and contributor approval. If a wet-lab collaboration becomes available, the assay panel will be selected **before assay outcomes are known** from mechanically defined discordant candidates generated by Expansion A/C. A proposed default panel is 20 candidates: 10 high-ranked by the structure-neutralized-selected model and 10 high-ranked by the conventional-selected model, matched as closely as possible on drug and target popularity and tested blind to model identity.

## Authorship boundary

No collaborator will be added solely to strengthen perceived credibility. Authorship requires a genuine qualifying contribution and explicit approval of the submitted manuscript.

## Stop / go rule for Science-only strategy

The expanded package will be considered materially stronger for Science only if at least two of the following three evidence classes survive without claim inflation:

1. independent direct BindingDB evidence supports materially different yield or named later-supported candidates after degree/popularity control;
2. biological interpretation shows reproducible, materially different pathway/disease prioritization between benchmark-selected queues;
3. the contemporary model ladder shows a benchmark-induced winner change with stability-controlled hypothesis turnover, or an equally consequential benchmark sensitivity that changes the decision-relevant top queue under a prespecified alternative decision rule.

Failure of these gates is not hidden; it becomes evidence that the strongest defensible destination is a top specialist or broad multidisciplinary journal rather than Science.
