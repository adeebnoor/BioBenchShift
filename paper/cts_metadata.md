# Science CTS metadata — pre-submission sheet

**Journal:** Science  
**Article type:** Research Article  
**Manuscript title:** Benchmark design redirects biomedical discovery  
**Corresponding author:** Adeeb Noor  
**Affiliation:** Department of Information Technology, Faculty of Computing and Information Technology, King Abdulaziz University, Jeddah, Saudi Arabia  
**Email:** arnoor@kau.edu.sa  
**ORCID:** 0000-0002-8251-1853

## Abstract

Biomedical AI increasingly determines which molecular relations are prioritized for experiments, yet models are often selected on benchmarks that contrast known relations with sampled unknown pairs. Across four relation families, structural observability alone produced strong discrimination and affected learned models unequally, reversing the selected model in drug–target and protein-interaction analyses. These reversals changed 100% and 99% of top-100 hypotheses, respectively, far beyond within-model instability. In a frozen historical BioGRID network, the model selected after structural neutralization recovered 522 of 5,635 later-added interactions among its top 50,000 predictions versus 181 for the conventional winner. Independent ChEMBL evidence was also more concentrated at the earliest drug–target cutoffs, although the advantage did not persist at broad cutoffs. Thus benchmark design can redirect biomedical discovery by changing which model wins and which biology is tested first.

## Suggested subject-area positioning

Primary positioning should emphasize **AI as an upstream scientific-decision system**, not a specialist DTI benchmark.

Recommended descriptors to use where the CTS vocabulary permits:

- artificial intelligence / machine learning;
- computational biology / bioinformatics;
- systems biology / network biology;
- drug discovery / pharmacology;
- scientific methodology / reproducibility.

Avoid making “drug–drug interaction” the primary subject because Anti-DDI is only a supportive boundary analysis in this paper.

## Keywords

- biomedical artificial intelligence
- benchmark design
- model selection
- hypothesis prioritization
- link prediction
- network biology
- drug–target interaction
- protein–protein interaction
- structural bias
- reproducibility

## One-sentence editor pitch

A benchmark can change the biomedical AI model that wins, the hypotheses that reach the front of the experimental queue, and the later or independent evidence concentrated among those hypotheses.

## Data and code availability — submission wording

All analysis code, frozen protocols, seed-level outputs, figure source data, release identifiers and checksum/provenance files supporting the manuscript are maintained at the public GitHub repository `adeebnoor/BioBenchShift`, branch `main`. Before submission, the exact commit used for the manuscript will be frozen as a release and archived under a persistent DOI. Third-party datasets are referenced by release/version and remain subject to their original access and licensing terms. BioGRID historical and later-release identifiers, ChEMBL 37 evidence criteria and frozen input hashes are recorded in the repository provenance files.

## Related-manuscript disclosure

This manuscript is scientifically distinct from two related lines of work:

1. **Anti-DDI v3.0.1** is an evidence-state resource for drug non-interaction research. It motivates the distinction between unknown and curated counter-evidence but is not the scope of this manuscript. The Science manuscript adds cross-domain relation families, learned-model comparisons, benchmark-induced winner reversals, fixed-universe hypothesis-identity analyses, historical BioGRID later-evidence testing and independent ChEMBL drug–target evidence.
2. A separate **allocation-identity / RIDI** manuscript studies identity changes under finite-capacity allocation in non-biomedical primary examples. Its EPSS/COMPAS evidence, theorem and primary argument are not reused as evidence here. The present manuscript independently defines hypothesis turnover to measure biomedical scientific-decision identity after benchmark-induced model selection.

No result should be represented as novel here if it is already a principal result in either related manuscript.

## Funding

**Action before submission:** enter only verified grant/award identifiers. Do not infer or invent funding. If this work had no dedicated external award, use the institutionally correct no-specific-funding wording permitted by CTS.

## Competing interests

Current working statement: **The author declares no competing interests.**  
**Action before submission:** verify against Science/AAAS definitions and any institutional/IP disclosures connected to the work.

## Human / animal subjects

No new human-participant or animal experiments are reported in the present computational study. External biomedical databases are analyzed under their published data-access terms. **Action before submission:** confirm whether any underlying controlled-access data or institutional review statement is required; current headline analyses use public benchmark/database resources.

## Reviewer strategy

Target independent reviewers who collectively cover:

1. graph/link-prediction benchmark validity and degree bias;
2. network biology and temporal PPI evaluation;
3. drug–target interaction modeling and biomedical ML bias;
4. biomedical knowledge-graph/generalization evaluation;
5. AI-for-scientific-discovery methodology.

Prefer reviewers able to evaluate the **decision-consequence chain**, not only one model family. Avoid close collaborators, recent coauthors, institutional conflicts, personal conflicts, and anyone involved in the related Anti-DDI/RIDI manuscripts.

Potential literature-derived candidate pool to vet for conflicts before entering CTS:

- Sadamori Kojaku / Yong-Yeol Ahn — generic link-prediction degree bias;
- Mehmet Koyutürk — network biology bias-aware evaluation;
- Anaïs Baudot or Galadriel Brière — biomedical KG leakage/generalization;
- Kasper Lage or Fatma-Elzahraa Eid — systematic auditing of biological ML;
- Ilya Safro — temporal biomedical hypothesis benchmarking.

These are **candidate domains/names, not a final reviewer list**; conflict, availability and current contact details must be checked immediately before submission.

## Possible excluded-reviewer logic

Use exclusions only for concrete conflict or documented competitive/confidentiality reasons. Do not exclude authors merely because they published close prior art. Record the reason privately before entering any exclusion.

## Files to upload in CTS

Required/primary package:

1. `Science_Research_Article_PRE_SUBMISSION.docx` — replace PRE-SUBMISSION label at final freeze.
2. Main figures integrated in manuscript or uploaded as requested by the live CTS instructions.
3. `Science_Supplementary_Materials_PRE_SUBMISSION.docx`.
4. `Science_Cover_Letter.docx`.
5. Any source-data or supplementary files requested by the current submission interface.

Before final submit, verify live CTS instructions rather than relying on this working sheet for technical limits.
