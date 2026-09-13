# Submission freeze — 13 September 2026

This addendum records the final submission state without rewriting the original prespecified hypotheses or earlier result checkpoints.

## Closed completion gates

All prespecified submission-completion analyses are now frozen. The final contemporary-model gate used the leakage-free GraphBAN-style protocol defined before outcome inspection, with ChemBERTa drug features, ESM-1b protein features, training-positive-only message passing, fixed seeds 0–2, and the same trained scores evaluated under conventional and structure-neutralized controls.

Mapping retained **18,631 / 18,690 (99.7%)** TargetDecagon positive edges. Across seeds 0, 1 and 2, conventional random-unlabelled evaluation yielded AUROC **0.9972 ± 0.0002** and AUPRC **0.9976**. Structure-neutralized evaluation yielded AUROC **0.9278 ± 0.0141** and AUPRC **0.9425**, for mean differences of **0.0694 AUROC** and **0.0551 AUPRC**. Mean matching coverage was **0.5880**. Held-out positive edges did not enter message passing.

The result supports benchmark sensitivity in a feature-rich contemporary graph architecture while also showing that strong predictive performance remains after structural neutralization. It therefore does not support a claim that the model lacks biological signal.

## Final evidence chain

The submission-level evidence supports the following bounded chain:

**cross-domain structural observability → model-ranking changes → stable selected-model hypothesis-identity changes in DTI and PPI → non-circular later-evidence differences in historical BioGRID → independent DTI evidence concentration at early ChEMBL cutoffs → benchmark sensitivity retained in a leakage-free contemporary feature-rich architecture.**

Important boundaries remain part of the evidence:

- compound–disease prediction did not show a winner reversal;
- single-fit disease–gene rankings were too unstable for clean benchmark-caused hypothesis-turnover attribution;
- ChEMBL evidence favored the neutralized-selected model at early cutoffs but crossed over at broader cutoffs;
- Anti-DDI evidence-state sensitivity did not pass the strict structural-balance gate in all seeds;
- BioGRID later additions reflect research and curation processes rather than unbiased biological truth;
- structural neutralization is not claimed to be uniquely correct or universally beneficial.

## Submission claim ceiling

The final manuscript may claim that benchmark design can materially alter apparent performance, model selection and the identity of biological hypotheses prioritized for follow-up, and that these changes can propagate to later or independent evidence concentration in the tested settings. It may also state that benchmark sensitivity persists in the prespecified leakage-free GraphBAN-style contemporary challenge.

The manuscript should not claim universal prospective superiority of structure-neutralized selection, absence of biological learning, universal structural bias across all architectures or domains, or causal biological truth from database additions alone.