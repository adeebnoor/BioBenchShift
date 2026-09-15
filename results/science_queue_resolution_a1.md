# Science A1/A2 same-fit queue-resolution replay

Protocol: `BioBenchShift-Science-A1-A2-secondary-20260915-v2-replay`
Candidate universe: **1,005,757** mapped unknown drug-target pairs.
Frozen-vs-replay status: **distinct_replay_execution**. Frozen metrics are retained separately and are not substituted into this analysis.

## GraphBAN-NeuralMF nominated comparison

| Seed | HT@100 | HT@500 | HT@1000 | |ΔAUROC| random | |ΔAUROC| matched |
|---:|---:|---:|---:|---:|---:|
| 0 | 1.000 | 0.996 | 0.903 | 0.001293 | 0.040556 |
| 1 | 1.000 | 0.844 | 0.815 | 0.000616 | 0.052321 |
| 2 | 1.000 | 0.994 | 0.881 | 0.000518 | 0.026038 |
| 3 | 1.000 | 0.984 | 0.844 | 0.000278 | 0.035211 |
| 4 | 1.000 | 0.998 | 0.876 | 0.001987 | 0.006938 |

Mean HT@100 = **1.000**; mean HT@500 = **0.963**; mean HT@1000 = **0.864**.
Mean absolute AUROC gap = **0.000939** under random controls and **0.032213** under matched controls.

All six model pairs are retained in the pairwise CSV. These are descriptive same-fit results. They do not establish statistical equivalence or a universal benchmark-resolution threshold.
