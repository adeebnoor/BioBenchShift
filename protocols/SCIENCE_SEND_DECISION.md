# Science send decision — evidence matrix

**Target:** Science (AAAS flagship)  
**Article type:** Research Article  
**Rule:** this file is a decision gate, not a claim that publication is likely or guaranteed. A SEND decision requires that the scientific story survives the tests below without broadening claims beyond the evidence.

## Editorial thesis

> **Benchmark design can redirect biomedical discovery by changing which AI model wins, which biological hypotheses are prioritized, and which later or independent evidence is concentrated among those priorities.**

The paper is not sold as a new degree-bias method. Its Science-level contribution is the decision-consequence chain:

**structural opportunity → model selection → hypothesis identity → later/independent evidence.**

## Gate matrix

| Gate | Requirement | Current evidence | Status |
|---|---|---|---|
| G1 General mechanism | Structural-only signal must inflate conventional evaluation across independent biomedical relation families | DTI, HuRI PPI, compound–disease and disease–gene all reproduce the effect; 4/4 external relation families | **PASS** |
| G2 Learned-model consequence | Learned model families must respond unequally enough to affect substantive conclusions | NeuralMF, SVD and LightGCN show heterogeneous sensitivity; winner reversals in DTI and disease–gene; non-reversal CtD control | **PASS for consequence; heterogeneous by design** |
| G3 Scientific-decision identity | Winner changes must redirect fixed-universe Top-K hypotheses beyond training noise in at least two relation settings | DTI HT@100=1.000 with much smaller within-family instability; PPI HT@100=0.990 with much smaller within-family instability | **PASS** |
| G4 Non-circular temporal consequence | Model chosen before future evidence must differ in later evidence recovery on the same complete candidate universe | Historical BioGRID: SVD 522 vs NeuralMF 181 later-supported relations in top 50,000; paired bootstrap difference remains positive at every prespecified K | **PASS** |
| G5 Independent non-PPI consequence | External evidence source outside PPI must test the predeclared DTI contrast without participating in model selection | ChEMBL 37: SVD recovers 3/5/6 supported pairs at top 100/500/1000 vs NeuralMF 0/0/2; boundary reverses at wide K and is retained | **PASS with high-priority-frontier boundary** |
| G6 Contemporary architecture | A current stronger architecture must be evaluated under the same frozen, leakage-free dual-regime protocol, regardless of direction | GraphBAN mapping gate passes at 99.7% positive-edge feature coverage; leakage-free TargetDecagon challenge is executing under frozen protocol | **PENDING RESULT** |
| G7 Evidence-state robustness | Anti-DDI should test whether random unknown and curated counter-evidence are interchangeable without overclaiming causal isolation | Overlap weighting leaves residual model-dependent differences, but strict balance gate passes only 5/10 seeds | **PASS as boundary/sensitivity only** |
| G8 Related-paper separation | Science manuscript must be substantively distinct from Anti-DDI and RIDI/Nature work | Explicit related-paper boundary; Science uses a new broad question, external datasets, learned-model experiments, hypothesis identity, temporal and ChEMBL consequences | **PASS subject to final disclosure package** |
| G9 Reproducibility | Headline results must have frozen protocol, provenance/checksums and executable analyses | Core gates, BioGRID snapshots, ChEMBL protocol, GraphBAN protocol amendment and result files are frozen; final release/DOI still needed | **PASS for analysis; release packaging pending** |
| G10 Current Science editorial shape | Title, abstract, first figure and narrative must foreground a broad discovery consequence and fit recent flagship Science papers | v0.6 title is consequence-first; abstract is compressed; main story is four figures; BioGRID and ChEMBL are consequence layers; specialist controls moved to Extended Data | **PASS for editorial architecture; final figures/cover letter pending** |

## Current-Science alignment check — 13 Sep 2026

The manuscript architecture was checked against recent and current *Science* publishing patterns, especially:

- Huang et al., **“Autonomous biomedical research with an artificial intelligence agent”**, *Science* 393, eadz4351 (20 Aug 2026), DOI 10.1126/science.adz4351: broad biomedical bottleneck → general AI mechanism → heterogeneous benchmarking → real-world scientific consequence.
- King et al., **“Generative design of bacteriophages with genome language models”**, *Science* 393, eaec2657 (6 Aug 2026), DOI 10.1126/science.aec2657: broad biological problem → computational advance → direct consequence → general implication.
- Beaglehole et al., **“Toward universal steering and monitoring of AI models”**, *Science* 391, 787–792 (19 Feb 2026), DOI 10.1126/science.aea6792: general model problem → controlled method → broad demonstrations → practical consequence.
- Current issue (10 Sep 2026, vol. 393, issue 6816): result-forward titles continue to dominate, reinforcing the choice to use **“Benchmark design redirects biomedical discovery”** rather than a technical DTI/link-prediction title.

The corresponding manuscript is now `manuscript/SCIENCE_MANUSCRIPT_DRAFT_v0.6.md`.

## Hard SEND rule

### SEND to Science if

1. G6 completes successfully as a valid leakage-free contemporary-model evaluation, **even if GraphBAN is less sensitive than simpler models**;
2. the GraphBAN result is reported without cherry-picking or post-outcome tuning;
3. the final abstract preserves the ChEMBL wide-K boundary and does not convert frontier enrichment into a universal superiority claim;
4. the final related-paper disclosure clearly separates Anti-DDI and RIDI/Nature materials;
5. all headline figure source data and frozen code are bundled in a permanent release/DOI;
6. a final adversarial editorial read concludes that title + abstract + Fig. 1 communicate a broad scientific-decision result to a non-specialist;
7. the main paper remains a **four-figure consequence arc**, rather than expanding back into a catalog of benchmark diagnostics.

### HOLD / redesign before Science if

- the contemporary-model analysis cannot be completed reproducibly under a leakage-free protocol;
- any headline result depends on changing a threshold, split, matching rule or model after inspecting its outcome;
- the manuscript claims that structure-neutralized selection is universally superior rather than showing a decision consequence with boundaries;
- the Science and Anti-DDI/RIDI manuscripts cannot be cleanly distinguished in question, evidence and principal conclusions;
- the external-consequence story collapses after provenance or leakage auditing;
- the first page reads primarily as a negative-sampling or link-prediction methods paper.

## What a favorable GraphBAN result would mean

If GraphBAN also shows a material random-to-structure-neutralized change, the paper can state that benchmark sensitivity extends to a contemporary feature-rich graph architecture. It should remain an Extended Data robustness result unless it changes the central scientific conclusion.

## What an unfavorable GraphBAN result would mean

If GraphBAN is comparatively stable, retain it prominently as a boundary:

> Structural benchmark dependence is model-specific; some contemporary architectures can be more robust, but conventional evaluation can still select substantially different scientific priorities among plausible model families.

A robust GraphBAN result does **not** invalidate the primary chain because that chain is supported by the actual selection reversals and their temporal/independent consequences. It prevents overgeneralization.

## Submission mechanism

Submit through the AAAS **Science Content Tracking System (CTS)** at `https://cts.sciencemag.org/` as a **Research Article**. The final initial-submission file should be optimized for editorial reading, with the four main figures integrated into the narrative and Supplementary Materials supplied separately according to the live author instructions.

Before upload, re-check the live Science pages because technical requirements can change:

- `https://www.science.org/content/page/instructions-preparing-initial-manuscript`
- `https://www.science.org/content/page/science-journals-editorial-policies`

## Current decision

**HOLD only for G6 completion and final submission packaging.**

The scientific case for a *Science* submission is supported by a cross-domain mechanism, two stability-controlled hypothesis-identity replications, a non-circular temporal consequence, and an independent DTI evidence consequence. The presentation has now been recast to match current flagship *Science* style: one broad question, one escalating causal/decision chain, four main figures, explicit negative boundaries, and independent consequence beyond benchmark metrics.