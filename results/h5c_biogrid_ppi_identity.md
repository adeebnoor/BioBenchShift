# PPI H5 — BioGRID hypothesis identity with ensemble stability

Both selected model families are trained on the complete frozen BioGRID 5.0.250 human MV-Physical graph. Each family is a three-initialization score ensemble. Leave-one-initialization-out ensembles quantify within-family ranking instability on the identical full historical non-edge universe.

| K | NeuralMF vs SVD ensemble HT | NeuralMF LOO HT mean | SVD LOO HT mean |
|---:|---:|---:|---:|
| 100 | 0.990 | 0.137 | 0.050 |
| 500 | 0.986 | 0.139 | 0.035 |
| 1,000 | 0.976 | 0.130 | 0.038 |

Candidate universe: **70,041,100** historical non-edges among **11,844** historical proteins.

This experiment measures scientific-decision identity only. Later-evidence correctness is evaluated separately in H4/H4b.
