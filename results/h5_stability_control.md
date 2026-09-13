# H5 stability control — model-selection turnover versus initialization instability

All comparisons use one fixed training split and one fixed candidate universe per relation family. Only model initialization/training randomness varies.

| Family | Comparison | HT@100 | HT@500 | HT@1000 |
|---|---|---:|---:|---:|
| DTI | within NeuralMF (mean pairwise) | 0.307 | 0.036 | 0.212 |
| DTI | within SVD (mean pairwise) | 0.138 | 0.136 | 0.155 |
| DTI | conventional winner vs neutralized winner (mean paired-init) | 1.000 | 0.999 | 0.945 |
| Disease-Gene | within NeuralMF (mean pairwise) | 0.934 | 0.877 | 0.849 |
| Disease-Gene | within LightGCN (mean pairwise) | 0.684 | 0.658 | 0.642 |
| Disease-Gene | conventional winner vs neutralized winner (mean paired-init) | 0.944 | 0.898 | 0.872 |

A benchmark-selection effect is strongest when cross-selected-model turnover materially exceeds ordinary within-model initialization turnover on the same candidate universe. No biological correctness is inferred from stability alone.
