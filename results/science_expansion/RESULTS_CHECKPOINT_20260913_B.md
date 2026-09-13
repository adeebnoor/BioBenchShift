# Science expansion — results checkpoint B (13 September 2026)

Frozen baseline: `1298e4c83ff4fa488a784d52219d20e706118a9e`.

This checkpoint records two preregistered expansion outcomes without rewriting the frozen baseline.

## D. Named later-supported BioGRID cases

Protocol committed before pair identities were inspected: `BIOGRID_NAMED_LATER_CASES_AMENDMENT_20260913.md`.

Immutable workflow artifact:
- artifact id `10313984239`
- digest `sha256:d27592a0012153d322422e805be70401c35ebd9a2e97be9ef94af9e05b7670fe`
- workflow head `157edb984eadb693122a56b3bd6a7da272f78a80`

The reconstruction passed both frozen-count assertions:
- SVD Top-100 later-added closed-world interactions: **7**;
- NeuralMF Top-100 later-added closed-world interactions: **0**.

All seven mechanically selected SVD Top-100 cases are retained:

| SVD rank | Pair | NeuralMF rank | LightGCN rank | Later BioGRID publication evidence | System |
|---:|---|---:|---:|---|---|
| 34 | RPL3–RPS25 | 16,680 | 64,732 | PMID 26344197; PMID 41387687 | Co-fractionation; XL-MS |
| 43 | YAP1–ACTB | 10,815 | 3,890 | PMID 31501420; PMID 38537774 | Affinity Capture-MS |
| 51 | HIST2H3A–HIST1H3A | 21,962 | 20,485 | PMID 35575683; PMID 36479438 | XL-MS |
| 57 | RPL30–RPS19 | 40,920 | 54,806 | PMID 26344197; PMID 41387687 | Co-fractionation; XL-MS |
| 66 | H3F3A–HIST1H2AB | 13,294 | 20,748 | PMID 30021884; PMID 35575683; PMID 36479438; PMID 40593561 | XL-MS |
| 69 | RPS9–RPL23 | 48,996 | 84,007 | PMID 26344197; PMID 41387687 | Co-fractionation; XL-MS |
| 81 | RPL35–RPS4X | 23,035 | 67,876 | PMID 26344197; PMID 41387687 | Co-fractionation; XL-MS |

BioGRID 5.0.250 was compiled on **25 September 2025**. PMID 41387687 (Jiao et al., Nature Communications, DOI 10.1038/s41467-025-66023-0) was published on **12 December 2025**. Thus four members of the mechanically fixed seven-pair case set — RPL3–RPS25, RPL30–RPS19, RPS9–RPL23 and RPL35–RPS4X — have additional experimental XL-MS support from a publication that post-dates the historical freeze. Those four pairs were ranked 34, 57, 69 and 81 by SVD versus 16,680, 40,920, 48,996 and 23,035 by NeuralMF, respectively.

Interpretation boundary: the study did not cause these discoveries. Some selected pairs also have older publication evidence that was not present in the frozen MV-Physical candidate state, and later database additions include curation effects. The defensible statement is that the benchmark-selected queues placed concrete interactions that later entered the BioGRID evidence state — including four with post-freeze XL-MS publication support — at radically different ranks.

## E. Matching-method robustness

Protocol: `MATCHING_ROBUSTNESS_AMENDMENT_20260913.md`.

Immutable workflow artifact:
- artifact id `10313859766`
- digest `sha256:ccde274b500ae7c339222b938b0dd61e60760cafdefa808098e6076f662c7a56`
- workflow head `7fcf9eeb3c732b75e1feb5ccf4e26da2e39db065`

Ten independent split seeds were evaluated under three prespecified neutralization rules.

| Method | Mean matching coverage | Neutralized winner | Winner-change fraction | Mean neutralized AUROC: SVD | NeuralMF | LightGCN |
|---|---:|---|---:|---:|---:|---:|
| Original joint log2 bins | 0.5989 | SVD 7/10; NeuralMF 3/10 | 0.70 | 0.9151 | 0.9093 | 0.8851 |
| Joint degree deciles | 1.0000 | NeuralMF 10/10 | 0.00 | 0.9405 | 0.9790 | 0.9705 |
| Nearest-neighbor caliper | 0.4848 | NeuralMF 10/10 | 0.00 | 0.9046 | 0.9323 | 0.9145 |

Conventional evaluation selected NeuralMF in all 10 seeds for all methods because the conventional arm is shared. The original log2-bin neutralization reproduces the SVD reversal in 7/10 seeds, but the reversal is **not invariant to the neutralization rule**.

Interpretation boundary: this is an important negative/heterogeneity result. It prevents any claim that there is one uniquely correct structural-neutralization procedure or that winner reversal is guaranteed whenever structure is controlled. It strengthens the broader decision-theoretic statement that benchmark construction — including the definition of the control population — is itself part of model selection and can materially change the conclusion.

## Current implication for the Science-only strategy

The expansion now contains:
1. a strong biological-program divergence result (checkpoint A);
2. seven fully mechanical named BioGRID later-supported cases, four with post-freeze XL-MS publication support and very large cross-model rank displacement;
3. a transparent robustness boundary showing that the exact winner reversal depends on how structural neutralization is implemented;
4. mixed/supportive rather than decisive BindingDB evidence.

The principal remaining high-value gates are Open Targets disease-area divergence and the preregistered contemporary architecture ladder. Wet-lab validation remains outside the evidence unless real assays and raw measurements are obtained.
