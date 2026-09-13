# H3 pilot — model-rank instability under structure-neutralized evaluation

## Result

Across all three bipartite biomedical relation tasks evaluated with the same train/test and degree-matched protocol, changing only the negative-evaluation design changed at least one pairwise model ordering.

| Task | Conventional random-negative ranking | Degree-matched ranking | Explicit reversal |
|---|---|---|---|
| BioSNAP drug–target | NeuralMF (0.997) > LightGCN (0.989) > SVD (0.950) | SVD (0.914) > NeuralMF (0.908) > LightGCN (0.883) | SVD moves 3rd → 1st |
| Hetionet drug–disease treatment | LightGCN (0.902) > NeuralMF (0.879) > SVD (0.717) | LightGCN (0.882) > SVD (0.705) > NeuralMF (0.560) | NeuralMF moves 2nd → 3rd |
| Hetionet disease–gene | NeuralMF (0.863) > LightGCN (0.794) > SVD (0.620) | LightGCN (0.721) > NeuralMF (0.640) > SVD (0.604) | NeuralMF and LightGCN swap 1st/2nd |

## Rank-correlation summary

With ranks encoded from best=1 to worst=3:

- DTI: Spearman rank correlation = **−0.50**; Kendall tau = **−0.33**.
- Drug–disease: Spearman rank correlation = **0.50**; Kendall tau = **0.33**.
- Disease–gene: Spearman rank correlation = **0.50**; Kendall tau = **0.33**.

The purpose of this pilot is not to crown any of these models. It establishes that a benchmark-design intervention can alter comparative conclusions even when all models are trained on exactly the same training edges.

## Boundary

This is a **pilot H3 result**, not yet the final Science claim. The model panel is intentionally small (SVD, neural matrix factorization, LightGCN-style message passing). The result must be challenged with stronger KG/GNN and task-specific models and then related to temporal/independent generalization before the manuscript can claim a broadly applicable evaluation standard.
