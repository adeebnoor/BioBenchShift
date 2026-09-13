# Gate 1 — BioSNAP DTI result

- Edges: **18,690**
- Drugs: **284**
- Targets: **3,648**
- Degree-only AUC, conventional random negatives: **0.983 ± 0.001**
- Degree-only AUC, targeted degree-matched negatives: **0.620 ± 0.012**
- Mean AUC inflation: **0.364**
- Mean matched fraction: **0.599**
- Mean evaluable seen-endpoint fraction: **0.837**
- Gate-1 DTI structural-inflation signal: **PASS**

Matching uses nonedges drawn from the same prespecified log2 degree-bin pair as each positive test edge, with endpoint degrees computed from training edges only. The internal pass rule is conventional degree-only AUC >= 0.60, AUC inflation >= 0.08, and >= 50% matching coverage.
