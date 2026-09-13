# Science editorial-screen playbook — 13 September 2026

## Why this file exists

For *Science*, the first failure mode is editorial screening, not specialist peer review. A 2026 analysis of deidentified AAAS records described 110,303 submissions to *Science* and *Science Advances* from 2015–2020. In the *Science* subset, editors summarily rejected 28.1% of submissions and sent 69.5% to the Board of Reviewing Editors (BoRE); after BoRE advice, editors declined 75% of those manuscripts and sent the remainder to formal peer review. The overall acceptance rate in that historical dataset was 6.1%.

These are historical process data, not a claim about a current 2026 acceptance probability. Their strategic implication is still clear: **the manuscript must establish broad importance before a specialist review ever occurs.**

## The 90-second editorial screen

The title, abstract, first two paragraphs and Fig. 1 must answer four questions without requiring knowledge of link prediction:

1. **What scientific decision is at stake?**  
   Biomedical AI increasingly determines which molecular hypotheses receive experiments.

2. **What is wrong with the current assumption?**  
   Model benchmarks can reward structural observability that is partly independent of pair-specific biology.

3. **Does this change a real scientific decision, not just a score?**  
   Yes. The benchmark changes the selected model and changes 99–100% of the top hypotheses in stability-controlled DTI/PPI analyses.

4. **Is there evidence beyond the benchmark itself?**  
   Yes. Frozen historical BioGRID selection predicts different later-evidence recovery, and independent ChEMBL evidence differs at the early experimental frontier.

If any of these answers requires reading the Supplementary Materials, the first page has failed.

## What the editor should remember after closing the PDF

Not:
- “They propose degree matching.”
- “They benchmark three link-prediction models.”
- “They found negative sampling bias.”

But:

> **A seemingly technical benchmark choice can change which biology an AI-driven discovery program tests first.**

## What Fig. 1 must do

Fig. 1 should be readable as a causal/decision schematic before the reader interprets any AUC:

**same observed biological network**  
→ two plausible benchmark constructions  
→ unequal model sensitivity  
→ different winner  
→ different experimental queue

Cross-domain degree-only AUCs then establish that the first link in the chain is not specific to one dataset.

Avoid using Fig. 1 as a dense benchmarking dashboard.

## What the abstract must not do

Do not:
- list every relation family and every model metric;
- explain degree bins or matching algorithms;
- introduce HT notation;
- mention RIDI;
- introduce Anti-DDI as a second story;
- claim universal superiority of neutralization.

The abstract should contain one decisive number for hypothesis identity and one consequence result. The current preferred ~134-word abstract follows this rule.

## Current best title

**Benchmark design redirects biomedical discovery**

Why it fits recent *Science* style:
- five words;
- states an effect rather than a method;
- keeps “benchmark” because it is the intervention;
- resolves immediately to the scientific consequence, “biomedical discovery.”

Do not revert to titles centered on link prediction, structural shortcuts, degree bias, negative sampling or model evaluation.

## Current relevant Science comparators

### Townley et al., 27 Aug 2026
**De novo design of RNA pseudoknots with deep learning** — *Science* 393, 931–937; DOI 10.1126/science.aeg6829.

Editorial lesson: computational AI is justified by blind challenges and orthogonal biological evidence. The paper's memorable result is what can be designed biologically, not the benchmark score.

**Our analogue:** model selection is frozen before BioGRID future evidence and independent ChEMBL evidence are opened. Treat those as discovery-consequence tests, not routine validation datasets.

### Huang et al., 20 Aug 2026
**Autonomous biomedical research with an artificial intelligence agent** — DOI 10.1126/science.adz4351.

Editorial lesson: heterogeneous benchmarking supports a larger claim about scientific capability; it is not the endpoint.

**Our analogue:** four relation families establish generality of the mechanism; the paper is ultimately about what enters the experimental queue.

### King et al., 6 Aug 2026
**Generative design of bacteriophages with genome language models** — DOI 10.1126/science.aec2657.

Editorial lesson: computational novelty is paired with biological/functional consequence and the title makes that consequence visible.

### Beaglehole et al., 19 Feb 2026
**Toward universal steering and monitoring of AI models** — DOI 10.1126/science.aea6792.

Editorial lesson: wet-lab experiments are not an absolute prerequisite for flagship methodological AI if the question is general, evidence is broad and the consequence is substantive.

## BoRE-level objections we must preempt on page 1

### “Degree bias is already known.”
Response in manuscript logic, not defensive prose:
- concede prior art immediately;
- novelty is the measured downstream decision chain;
- move from structural-only signal to winner reversal to hypothesis identity to opened-later evidence.

### “This is only a benchmark paper.”
Preemption:
- title says discovery;
- Fig. 2 shows changed hypothesis identity;
- Figs. 3–4 show temporal/independent evidence consequence.

### “You cherry-picked a correction that favors SVD.”
Preemption:
- CtD non-reversal is visible;
- disease–gene instability is retained;
- ChEMBL wide-K crossover is visible;
- GraphBAN direction is accepted regardless of outcome;
- no claim that SVD or degree matching is universally superior.

### “The future validation is circular.”
Preemption:
- BioGRID future endpoint ranks the complete fixed historical non-edge universe;
- no degree matching is used in the future-evidence endpoint;
- ChEMBL is independent of model fitting and selection.

### “This is an old-model artifact.”
Preemption:
- leakage-free contemporary GraphBAN challenge is frozen and executing;
- a robust GraphBAN result will be reported as a model-specific boundary, not suppressed.

### “The candidate lists differ just because training is unstable.”
Preemption:
- DTI and PPI use ensemble and leave-one-initialization-out controls;
- cross-selected-model turnover is far larger than within-family turnover.

## Peer-review layer after the editorial screen

Once the paper reaches external review, the detailed defenses move to Supplementary/Extended Data:
- exact matching rules and coverage;
- model hyperparameters and seeds;
- full uncertainty;
- complete candidate-universe definitions;
- temporal release provenance and hashes;
- ChEMBL mapping/evidence rules;
- GraphBAN leakage audit and clean-room protocol;
- Anti-DDI balance diagnostics;
- source data and executable scripts.

This separation is intentional: **broad consequence in the main article; forensic reproducibility in the supplement.**

## Current SEND assessment

Editorial architecture: **PASS** after the v0.6 rewrite.

Scientific packaging: **HOLD** until the leakage-free contemporary-model challenge finishes and is frozen.

Submission package: **HOLD** until four final figures, permanent code/data archive, related-manuscript disclosure and final DOCX/PDF are complete.
