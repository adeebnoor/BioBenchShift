# Alignment audit against recent Science papers — 13 September 2026

## Scope

This audit asks a narrow editorial question: **does the current manuscript read like a flagship Science Research Article, rather than a specialist biomedical-AI benchmark paper?**

The comparison prioritizes recent *Science* original research at the AI/biology interface and the latest issue available before 13 September 2026.

## Most recent directly relevant anchor found

### Townley et al. — De novo design of RNA pseudoknots with deep learning

*Science* 393, 931–937 (27 Aug 2026). DOI: `10.1126/science.aeg6829`.

Why this is the strongest current editorial comparator:
- AI is central but the headline is a **biological design result**, not an algorithm name;
- the abstract opens with the biological bottleneck;
- the claim is tested in **blind challenges**;
- the computational result is followed by orthogonal experimental evidence: chemical mapping, compensatory mutagenesis and cryo-EM;
- the closing sentence generalizes what biological design may now be possible.

Editorial lesson for our paper:

> Do not sell “better evaluation.” Sell a discovery-control result that is opened against evidence the selection process did not see.

Our nearest analogue to a blind challenge is not another random split. It is the frozen decision chain:
- select under historical/independent-free information;
- freeze the candidate universe;
- then open BioGRID future additions or ChEMBL external evidence.

This is why BioGRID and ChEMBL belong in the main four-figure arc.

## Other recent anchors

### Huang et al. — Autonomous biomedical research with an artificial intelligence agent

*Science* 393, eadz4351 (20 Aug 2026). DOI: `10.1126/science.adz4351`.

Pattern:
1. broad biomedical bottleneck;
2. general AI capability;
3. systematic heterogeneous benchmarking;
4. real-world scientific cases;
5. discovery-level implication.

Lesson: heterogeneous benchmarking earns space only because it supports a broader scientific capability. For us, the cross-domain benchmark experiments establish the mechanism; the discovery consequence is the headline.

### King et al. — Generative design of bacteriophages with genome language models

*Science* 393, eaec2657 (6 Aug 2026). DOI: `10.1126/science.aec2657`.

Pattern:
1. broad biological problem;
2. one direct advance;
3. quantitative experimental outcomes;
4. orthogonal structural/functional confirmation;
5. general biological implication.

Lesson: avoid ending on model metrics. A Science reader should remember what changed biologically/scientifically.

### Beaglehole et al. — Toward universal steering and monitoring of AI models

*Science* 391, 787–792 (19 Feb 2026). DOI: `10.1126/science.aea6792`.

Pattern:
1. general AI problem;
2. controlled representation method;
3. broad demonstrations;
4. practical safety/capability consequence.

Lesson: a methodological AI paper can fit *Science* without wet-lab experiments when the problem is general and the consequence is broad, controlled and empirically deep.

### Yang et al. — Structural ontogeny of protein-protein interactions

*Science* 391, eadx6931 (12 Feb 2026). DOI: `10.1126/science.adx6931`.

Pattern:
- starts from a broad protein-interaction question;
- uses computational/machine-learning analysis as one part of a mechanistic biological argument;
- closes with a biological interpretation, not a method claim.

Lesson: because our paper contains PPI analysis, we must avoid making PPI/network prediction itself the endpoint. PPI is one replication of a broader discovery-selection phenomenon.

## Latest-issue style check

The latest issue checked before this audit is *Science* vol. 393, issue 6816 (10 Sep 2026). Its research titles remain strongly result-forward, for example:
- “Growth rate overrides the benefit of extended growing season from boreal to semiarid conifers”;
- “Lean adipocyte oxylipin signaling restrains breast cancer through ferroptosis”;
- “Mechanism of membrane perforation in rotavirus cell entry.”

These titles state the result/mechanism rather than the assay, benchmark or software.

**Implication:** `Benchmark design redirects biomedical discovery` is substantially better aligned than a technical title built around link prediction, degree matching or negative sampling.

## Abstract-density audit

Approximate word counts of recent *Science* abstracts:
- Huang/BIOMNI: ~123 words;
- King/phage design: ~120 words;
- Beaglehole/AI monitoring: ~128 words;
- Townley/RNA pseudoknots: similarly compact and result-led.

The first v0.6 abstract was ~172 words and therefore too dense relative to these comparators. The preferred replacement is ~134 words:

> Biomedical AI increasingly determines which molecular relations are prioritized for experiments, yet models are often selected on benchmarks that contrast known relations with sampled unknown pairs. Across four relation families, we found that structural observability alone produced strong discrimination and affected learned models unequally, reversing the selected model in drug–target and protein-interaction analyses. These reversals changed 100% and 99% of top-100 hypotheses, respectively, far beyond within-model instability. In a frozen historical BioGRID network, the structure-neutralized-selected model recovered 522 of 5,635 later-added interactions among its top 50,000 predictions versus 181 for the conventional winner. Independent ChEMBL evidence was also more concentrated at the earliest drug–target cutoffs, although the advantage did not persist at broad cutoffs. Thus benchmark design can redirect biomedical discovery by changing which model wins and which biology is tested first.

This version should replace the longer v0.6 abstract at the next integrated manuscript freeze.

## Editorial comparison: recent Science versus our paper

| Dimension | Recent Science pattern | Current project response |
|---|---|---|
| Headline | Result/capability, not tool | **Benchmark design redirects biomedical discovery** |
| Opening | Broad scientific bottleneck | AI-selected hypotheses determine what gets tested |
| Technical mechanism | Introduced after significance | Structural observability / benchmark construction |
| Validation | Blind, orthogonal, temporal or real-world evidence | Frozen candidate universes + later BioGRID + independent ChEMBL |
| Main figures | Small number, each advances one step | 4-figure arc locked in v0.6 |
| Negative results | Bound the claim rather than derail it | CtD non-reversal; disease-gene instability; ChEMBL wide-K crossover; Anti-DDI balance limit |
| AI framing | Scientific consequence over model novelty | Benchmark affects scientific allocation of experiments |
| Closing claim | Broad but evidence-bounded | Benchmark design can determine which biology is tested first |

## What must not creep back into the main narrative

- long catalogs of AUCs;
- a separate main figure for every relation family;
- GraphBAN as a new headline simply because it is contemporary;
- Anti-DDI as a second paper embedded inside this paper;
- disease–gene as a clean replication without its instability boundary;
- claims that SVD or degree matching are universally superior;
- terminology such as RIDI in the title/abstract;
- novelty claims around degree bias itself.

## Current fit judgment

**Editorial architecture: strong and plausibly Science-shaped after v0.6.**

The strongest feature is no longer the magnitude of AUC inflation. It is the sequential consequence:

**a seemingly technical benchmark choice changes the model → changes almost every top hypothesis → changes the later/independent evidence encountered first.**

The largest remaining scientific gate is the frozen leakage-free contemporary-model challenge. The largest remaining editorial work is to convert the four figure concepts into publication-quality visual evidence and integrate the ~134-word abstract into the final manuscript.