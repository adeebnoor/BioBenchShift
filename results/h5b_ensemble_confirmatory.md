# H5b confirmatory ensemble analysis

Each model family is represented by the mean score of five independently initialized fits on one fixed training split. Ensemble stability is estimated by leave-one-initialization-out (LOO) ensembles on the same candidate universe.

| Family | K | Cross-selected ensemble HT | Conventional LOO HT mean | Neutralized LOO HT mean | Cross − worst LOO mean |
|---|---:|---:|---:|---:|---:|
| DTI | 100 | 1.000 | 0.085 | 0.042 | 0.915 |
| DTI | 500 | 1.000 | 0.005 | 0.033 | 0.967 |
| DTI | 1,000 | 0.939 | 0.066 | 0.038 | 0.873 |
| Disease-Gene | 100 | 0.800 | 0.378 | 0.219 | 0.422 |
| Disease-Gene | 500 | 0.740 | 0.353 | 0.177 | 0.387 |
| Disease-Gene | 1,000 | 0.717 | 0.325 | 0.191 | 0.392 |

Interpretation rule: a clean H5 replication requires cross-selected-model turnover to remain materially larger than ensemble instability. Ensembling is confirmatory; the original single-fit H5 and its negative/boundary findings remain reported.
