# H5 pilot — benchmark-induced hypothesis-selection instability

The candidate universe is identical within each seed. Known positives and endpoints absent from training are excluded before either winner ranks candidates.

| Family | Conventional winner | Neutralized winner | Candidate pairs | Spearman | HT@100 | HT@500 | HT@1000 |
|---|---|---|---:|---:|---:|---:|---:|
| DTI | NeuralMF | SVD | 817292 | 0.091 | 1.000 | 1.000 | 0.941 |
| Disease-Gene | NeuralMF | LightGCN | 614361 | 0.258 | 0.942 | 0.899 | 0.876 |

HT@k = 1 - overlap/k. Values are means across frozen split seeds. The JSON records percentile bootstrap 95% intervals over split-level estimates.

Interpretation boundary: turnover is scientific-decision instability, not evidence that either candidate list is biologically correct. Temporal/independent validation is required to assign external validity.
