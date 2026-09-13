# Biological consequence amendment — 13 September 2026

This amendment refines Expansion B before any pathway- or disease-prioritization outcome has been generated or inspected.

## Ranked queue definitions

Use the same five-initialization SVD and NeuralMF full-candidate ensembles used in the frozen DTI external-validation pipeline. Exclude all known BioSNAP positive edges. For each model, define Top-100, Top-500 and Top-1,000 candidate-pair queues and take the unique NCBI Entrez target-gene IDs represented in each queue. The background universe is the full set of BioSNAP target genes, not the human genome.

## Functional enrichment

Use the g:Profiler public API for organism `hsapiens`, numeric namespace `ENTREZGENE_ACC`, sources `REAC` and `GO:BP`, custom background equal to the full BioSNAP target-gene universe, ordered-query mode disabled, and false-discovery-rate correction. Retain every returned term with adjusted p < 0.05. Do not manually delete, merge, or prioritize terms after seeing results.

For each K separately, compare the SVD and NeuralMF target queues by:

1. number of unique target genes and target-gene Jaccard overlap;
2. number of significant Reactome/GO:BP terms per model;
3. significant-term overlap coefficient `|A∩B| / min(|A|,|B|)` when both sets are nonempty;
4. number of significant terms unique to each queue;
5. Jensen–Shannon divergence between normalized `-log10(adjusted p)` profiles over the union of significant terms, with absent terms assigned zero mass. If both profiles have zero total mass, divergence is undefined and reported as null.

All returned enrichment rows and all summary values are retained.

## Disease-area context

Disease-area analysis, if executed, will use the Open Targets Platform 26.06 release. It is interpretive context and cannot validate an individual drug–target pair. Because target-to-disease aggregation choices materially affect this endpoint, no disease-area result will be generated until a separate explicit aggregation amendment is committed.

## Interpretation boundary

Functional enrichment is evidence that the benchmark-selected queues emphasize different biological programs only if the prespecified profile-divergence metrics differ. It is not evidence that one queue is biologically correct, that an enriched pathway is causally involved, or that any individual predicted drug–target pair is valid.
