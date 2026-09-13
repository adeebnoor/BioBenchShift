# Science expansion integration plan — 13 September 2026

This file governs how the post-freeze Science expansion is integrated into the manuscript. It is an editorial integration rule, not a new analysis protocol. The original submission baseline remains commit `1298e4c83ff4fa488a784d52219d20e706118a9e` and is not rewritten retrospectively.

## Core claim after expansion

The strongest defensible claim is:

> Benchmark construction is part of the scientific decision process: changing the evaluation/control population can change which model appears best, which experimental hypotheses are prioritized first, and which fine-grained biological programs those queues emphasize.

The manuscript must **not** claim that benchmark design universally changes the best model, that one neutralization rule is uniquely correct, that one queue is biologically correct, or that the study caused later experimental discoveries.

## Main-text evidence to promote

1. **Named temporal BioGRID cases.** Retain all seven preregistered SVD Top-100 later-added interactions, including the four pairs with additional post-freeze XL-MS publication support. Report the cross-model rank displacement explicitly. These are concrete later-supported interactions, not discoveries caused by the study.
2. **Functional biological-program divergence.** Report target-set Jaccard and preregistered Reactome/GO profile divergence at K=100/500/1000. Interpret as different biological prioritization, not correctness.
3. **Matching-rule heterogeneity.** Report that the original joint-log2 degree matching changes the winner in 7/10 seeds, while degree-decile and nearest-neighbor-caliper neutralization do not. Use this to sharpen the conceptual point that benchmark/control-population definition is itself part of model selection.
4. **Contemporary feature-rich challenge.** Incorporate the frozen GraphBAN-style leakage-free adaptation regardless of direction. Its purpose is to bound whether the evaluation effect persists in a contemporary feature-aware graph architecture; it is not a verbatim reproduction of the published GraphBAN framework.

## Evidence to keep secondary / supplementary

1. **BindingDB curated-article validation.** Retain as orthogonal supportive evidence only. It provides sparse evidence (11 mapped primary pairs), no Top-1000 recovery, and publication clustering. Do not describe 9/11 Top-5000 recovery as nine independent discoveries.
2. **Open Targets therapeutic-area analysis.** Retain as a negative scale-boundary result. Coarse therapeutic-area distributions are highly similar despite large target/pathway differences. This directly prohibits claims of broad disease-area redirection.
3. Full per-seed matching results, all enrichment rows, all Open Targets profiles, and pair-level BindingDB audit remain supplementary/source data.

## Required language boundaries

Allowed:
- "redirects experimental prioritization" or "changes the scientific queue" when referring to ranked hypotheses;
- "later-supported interactions" for BioGRID temporal cases;
- "different fine-grained biological programs" for pathway enrichment;
- "benchmark sensitivity is model- and control-definition-dependent" if supported by the contemporary challenge.

Not allowed:
- "caused new discoveries";
- "proved the biologically correct model";
- "redirected therapeutic areas";
- "winner reversal is universal";
- "GraphBAN reproduction" for the current S6 adaptation;
- any implication of wet-lab validation unless real assay data are added by an experimental collaborator.

## Figure architecture target

Aim for four main figures:

- **Fig. 1 — Evaluation changes apparent model quality.** Conventional vs structure-neutralized performance, with honest axes and uncertainty.
- **Fig. 2 — Evaluation changes the selected hypothesis queue.** Model-selection consequence plus HT/rank-displacement summaries.
- **Fig. 3 — Concrete later-supported biological interactions.** Seven fixed BioGRID cases, emphasizing rank displacement and post-freeze evidence without claiming causal discovery.
- **Fig. 4 — Biological consequence and boundary.** Fine-grained pathway divergence contrasted with coarse therapeutic-area similarity; contemporary-model result can occupy a panel or inset if compact.

If the contemporary architecture materially changes model-selection conclusions, Fig. 2/4 should be adjusted mechanically to show that outcome rather than preserving the older narrative.

## Title/abstract rule

Do not finalize title or abstract until the contemporary feature-rich challenge is frozen. The title must be strong but must describe **prioritization/decision consequences**, not experimental discovery itself, unless real prospective assays are added.

## Wet-lab boundary

No experimental assay result exists in this expansion. Do not simulate, imply, or backfill one. If real collaborators later generate assay data with raw measurements and a prespecified protocol, that becomes a new evidence class and requires its own dated amendment.
