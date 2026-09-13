# Gate 2 — LightGCN on Hetionet Disease-associates-Gene (DaG)

- Random-negative AUC: **0.794 ± 0.012**
- Degree-matched AUC: **0.721 ± 0.016**
- Mean AUC drop: **0.073**
- Matching coverage: **1.000**
- Seen-endpoint test fraction: **0.716**

Adjacency and degree are constructed from training edges only; held-out positive edges are excluded from message passing.
