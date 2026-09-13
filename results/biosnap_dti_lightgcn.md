# Gate 2 — LightGCN on BioSNAP TargetDecagon DTI

- Random-negative AUC: **0.989 ± 0.001**
- Degree-matched AUC: **0.883 ± 0.005**
- Mean AUC drop: **0.105**
- Matching coverage: **0.599**
- Seen-endpoint test fraction: **0.838**

Adjacency and degree are constructed from training edges only; held-out positive edges are excluded from message passing.
