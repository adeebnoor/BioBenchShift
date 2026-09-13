# Science submission strategy — 13 September 2026

**Target:** *Science* (AAAS flagship)  
**Article type:** Research Article  
**Project branch:** `main`  
**Manuscript positioning:** AI for scientific discovery / biomedical computation / network biology — **not** a DTI benchmark paper and **not** an Anti-DDI paper.

## 1. Editorial objective

The first objective is to survive *Science*'s editorial/Board screening, not merely to satisfy a specialist reviewer. The paper must be legible as a broad scientific result from the title, abstract, first figure and cover letter:

> **The benchmark used to select a biomedical AI model can redirect what biology is tested first.**

The paper should therefore sell a scientific consequence, not a metric correction. Degree/rich-node/prior bias is prior art. The new contribution is the measured chain:

**benchmark construction → model identity → hypothesis identity → later/independent evidence at the discovery frontier.**

## 2. Current submission mechanism

The AAAS Science Content Tracking System (CTS) is the active manuscript-management portal: `https://cts.sciencemag.org/`. It supports new submissions, continuing an incomplete submission and tracking a submitted manuscript.

The current CTS submission tutorial exposes the practical sequence:

1. **New Submission** and select the correct journal. The tutorial warns that the journal cannot be changed later within that submission.
2. Read the submission requirements and accept the AAAS terms and conditions.
3. Enter the author grid and identify the first author.
4. Enter manuscript metadata.
5. Enter funding information (including an explicit no-funding/not-listed path where applicable).
6. Select subject area(s).
7. Add suggested reviewers and excluded reviewers if desired; the tutorial marks this step as optional.
8. Upload documents. The manuscript itself is a required upload and the CTS tutorial accepts it as **.docx or .pdf**.
9. Review all required fields and submit.

We should therefore prepare the portal metadata, subject-area positioning, reviewer/exclusion list and disclosure language **before** opening the final submission rather than making these decisions ad hoc in CTS.

Science consolidated original research into the **Research Article** category in 2023, discontinuing the shorter Reports category. The editorial description at that change indicated a typical Research Article of roughly five printed pages, commonly about **2,000–3,000 words, 3–5 figures and ~50 references**. Treat these as design targets rather than immutable technical limits until the live author-instructions page is checked at final upload.

Current author-instructions page: `https://www.science.org/content/page/instructions-preparing-initial-manuscript`  
Current editorial-policies page: `https://www.science.org/content/page/science-journals-editorial-policies`

The initial manuscript should be built as one coherent main file with figures/tables integrated for editorial reading, with supplementary material separate. Data/code/software provenance must be explicit and persistent. For this computational paper, the public repository should be frozen to a release/commit and preferably archived to Zenodo before submission.

## 3. What the editorial screen must understand in <2 minutes

### Title
Do not lead with “degree matching,” “link prediction,” “negative sampling,” “RIDI,” or a model name. If “benchmark” remains in the title, it must immediately resolve to the discovery consequence.

Preferred title:

> **Benchmark design redirects biomedical discovery**

Strong alternatives:
- **Benchmark choice redirects biomedical discovery**
- **Benchmark structure redirects biomedical AI discovery**

### Abstract
The abstract should contain only five moves:
1. biomedical AI increasingly determines what is tested;
2. sampled-unknown benchmarks can encode structural opportunity;
3. this changes which model wins;
4. winner changes redirect top hypotheses beyond training noise;
5. temporal BioGRID and independent ChEMBL evidence show a consequence at the discovery frontier.

Do not enumerate every AUC or every domain in the abstract. Use only the decisive numbers.

### Figure 1
Figure 1 must communicate the complete problem visually without specialist background:

**same biology + different benchmark rule → different apparent winner → different hypotheses**

Then use the cross-domain degree-only result as evidence that the vulnerability is not a one-dataset curiosity.

## 4. Lessons from recent Science papers

### Huang et al., 2026 — “Autonomous biomedical research with an artificial intelligence agent”
*Science* 393, eadz4351 (20 Aug 2026), DOI `10.1126/science.adz4351`.

Editorial pattern:
- title states the scientific capability, not the internal product name;
- abstract starts from a broad biomedical research bottleneck;
- introduces the system only after the problem;
- benchmarks across heterogeneous tasks;
- then moves to real-world cases;
- closes on scientific discovery, not benchmark performance.

**Implication for us:** the paper should open with the fact that benchmark-selected models decide which hypotheses receive experiments. “Structural bias” is the mechanism, not the headline.

### King et al., 2026 — “Generative design of bacteriophages with genome language models”
*Science* 393, eaec2657 (6 Aug 2026), DOI `10.1126/science.aec2657`.

Editorial pattern:
- broad biological problem;
- one crisp computational advance;
- direct experimentally grounded consequence;
- closing statement generalizes the scientific capability.

**Implication for us:** BioGRID temporal evidence and ChEMBL independent evidence must function as the consequence layer. The manuscript cannot end at “AUC changed.”

### Beaglehole et al., 2026 — “Toward universal steering and monitoring of AI models”
*Science* 391, 787–792 (19 Feb 2026), DOI `10.1126/science.aea6792`.

Editorial pattern:
- a general conceptual problem about model representations;
- a method capable of testing the concept;
- demonstrations across many concepts/settings;
- practical consequence for monitoring/capability.

**Implication for us:** methodological AI work can fit flagship *Science* when it exposes a general scientific problem and demonstrates consequences beyond a single benchmark.

### Current issue check
The 10 Sep 2026 issue (Vol. 393, Issue 6816) continues the same editorial style: short, causal/consequence-led titles such as “Growth rate overrides the benefit of extended growing season from boreal to semiarid conifers” and mechanism-driven biomedical titles. The title makes the result visible before methods are read.

## 5. Main-paper architecture locked for v0.6

### Figure 1 — A benchmark can reward structural opportunity
Combine:
- conceptual schematic;
- four relation-family degree-only conventional vs neutralized AUCs;
- learned-model sensitivity as a compact panel.

**Editorial message:** the benchmark is an active measurement choice, not a passive scoreboard.

### Figure 2 — Benchmark choice changes model and hypothesis identity
Main evidence:
- DTI winner reversal;
- PPI winner reversal;
- DTI/PPI ensemble HT against within-family LOO noise.

Move to Extended Data:
- CtD non-reversal control;
- disease–gene single-fit instability and ensemble boundary;
- secondary cutoffs.

**Editorial message:** this is where a metric artifact becomes a scientific-decision artifact.

### Figure 3 — Historical benchmark choice changes later evidence recovery
BioGRID historical freeze → fixed complete 70,041,100-pair candidate universe → later release.

Show at a glance:
- model selected by conventional vs neutralized benchmark;
- later-supported hits across prespecified K;
- paired bootstrap uncertainty;
- explicit statement that later database additions are evidence, not biological truth.

**Editorial message:** the benchmark choice made before future evidence existed changes what future evidence is recovered.

### Figure 4 — Independent evidence changes at the experimental frontier
ChEMBL 37 is never used in model selection/fitting/tuning.

Show:
- independent mapping/evidence definition;
- top-100/500/1000 result prominently;
- broad-cutoff crossover as a boundary, not hidden.

**Editorial message:** the effect replicates outside the original data lineage at the small candidate sets that experiments can realistically test first.

### Extended Data / Supplementary
- GraphBAN leakage-free contemporary-model challenge (once frozen result is complete);
- original GraphBAN evaluation-leakage audit and protocol amendment;
- CtD non-reversal;
- disease–gene instability/ensemble analysis;
- Anti-DDI evidence-state overlap-weighted sensitivity and failed full balance gate;
- matching diagnostics, alternative structural definitions, full seed tables;
- BioGRID sensitivity analyses and all K;
- ChEMBL pChEMBL≥7 and broad-cutoff results;
- full reproducibility/provenance manifests.

## 6. Claims allowed in the main paper

Allowed:
- sampled-unknown benchmark construction can encode structural observability;
- learned models respond differently enough to change model selection;
- in DTI and PPI, selection reversal changes top hypotheses far beyond within-family initialization noise;
- in the historical BioGRID experiment, the neutralized-selected model recovered more later-supported interactions across prespecified cutoffs;
- in independent ChEMBL evidence, the neutralized-selected DTI model concentrated more support at the earliest discovery cutoffs;
- benchmark design can therefore become part of the discovery process.

Not allowed:
- degree bias is newly discovered;
- structural neutralization is universally superior;
- SVD is the biologically best model;
- persistent unknowns are negatives;
- later BioGRID additions are unbiased biological truth;
- ChEMBL validates superiority at all cutoffs;
- disease–gene turnover is cleanly attributable to the benchmark without reporting its intrinsic instability;
- Anti-DDI evidence state has been causally isolated;
- official GraphBAN transductive results are leakage-free.

## 7. Cover-letter logic

Keep the cover letter to three short substantive paragraphs.

**Paragraph 1 — broad problem**  
Biomedical AI increasingly determines which molecular hypotheses are tested, but model selection is usually treated as if the benchmark were a passive measuring instrument.

**Paragraph 2 — decisive evidence**  
Across multiple biomedical relation families, structural opportunity produces apparent predictive signal; learned models depend on it unequally; changing the benchmark reverses selected models; the reversals change 99–100% of top hypotheses in stability-controlled DTI/PPI analyses; and temporal BioGRID plus independent ChEMBL evidence show different evidence concentration among those priorities.

**Paragraph 3 — why Science**  
The result is not another negative-sampling correction. It is a general scientific-decision consequence: evaluation design can determine what biology enters the experimental queue. Emphasize the cross-domain result, prospective/independent evidence and transparent negative controls.

Do not discuss journal prestige, career importance, or claim the manuscript is “perfect for Science.”

## 8. Submission package

Before CTS submission, require:
- final Research Article manuscript in the current Science template;
- 4 main figures optimized for single-page editorial reading;
- Supplementary Materials with methods, extended analyses and Extended Data figures/tables;
- frozen GitHub commit/release and preferably Zenodo DOI;
- source data for every plotted value;
- model/data provenance table including BioSNAP, HuRI, Hetionet, BioGRID releases and ChEMBL 37;
- explicit code/data availability statement;
- contribution statement using CRediT terminology where appropriate;
- competing-interests statement;
- related-manuscript disclosure/firewall for the separate Nature/RIDI and Anti-DDI work;
- suggested reviewers selected for biomedical network science + AI evaluation, avoiding close collaborators and direct conflicts;
- concise cover letter;
- prewritten CTS metadata: title, author order, funding, subject areas and reviewer/exclusion information.

Preprints are compatible with AAAS policy when placed in recognized preprint repositories such as bioRxiv/arXiv/ChemRxiv/medRxiv; do not distribute other manuscript versions broadly without checking the live policy.

## 9. Science-specific red-team screen before Send

The manuscript does **not** go to CTS unless all answers are YES:

1. Can a broad scientist explain the result after title + abstract + Fig. 1?
2. Does the central claim remain important if the reader already knows degree/rich-node bias?
3. Is there a consequence beyond benchmark scores? **Yes: hypothesis identity + later/independent evidence.**
4. Is the consequence demonstrated in more than one biological setting? **Yes: DTI + PPI; BioGRID temporal + ChEMBL external.**
5. Are negative controls/boundaries visible rather than buried? **CtD, disease–gene stability boundary, broad-cutoff ChEMBL crossover, Anti-DDI balance limitation.**
6. Does the paper avoid claiming that neutralization is universally optimal?
7. Is the contemporary-model challenge complete and interpreted without leakage?
8. Can every main number be regenerated from a frozen public artifact?
9. Are the Nature/RIDI and Anti-DDI manuscripts clearly non-overlapping in question, evidence and contribution?
10. Are main text and figures compressed to a general-science narrative rather than a benchmarking catalog?

## 10. Current decision

**HOLD — one substantive pre-submission gate remains: the frozen leakage-free contemporary-model challenge.**

If the contemporary model is structurally sensitive, it broadens the architectural evidence. If it is robust, it becomes an important boundary showing that susceptibility is model-specific rather than universal. Either outcome is publishable evidence if frozen and reported without post-hoc criterion changes.

Once S6 is frozen, update the manuscript once, run the final red-team screen, freeze the repository, prepare the four-figure package and submit through CTS.