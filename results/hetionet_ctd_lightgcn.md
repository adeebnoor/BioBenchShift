# Gate 2 — LightGCN on Hetionet Compound-treats-Disease (CtD)

- Random-negative AUC: **0.902 ± 0.005**
- Degree-matched AUC: **0.882 ± 0.015**
- Mean AUC drop: **0.020**
- Matching coverage: **1.000**
- Seen-endpoint test fraction: **0.636**

Adjacency and degree are constructed from training edges only; held-out positive edges are excluded from message passing.
