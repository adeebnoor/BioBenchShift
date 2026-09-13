# Anti-DDI evidence-state overlap-weighted sensitivity

- Frozen positives: **8,094** pairs / **437** ATC5 nodes; SHA verified.
- Eligible curated T1/T2 pairs in positive vocabulary: **621**.
- Seeds meeting prespecified balance + ESS gates: **5/10**.

| Model | Ordinary AUC random | Ordinary AUC curated | Overlap AUC random | Overlap AUC curated | Overlap Δ random−curated (95% seed-bootstrap CI) |
|---|---:|---:|---:|---:|---:|
| SVD | 0.771 | 0.730 | 0.697 | 0.689 | +0.008 [+0.003, +0.012] |
| NeuralMF | 0.919 | 0.830 | 0.798 | 0.750 | +0.048 [+0.040, +0.057] |
| LightGCN | 0.864 | 0.789 | 0.737 | 0.706 | +0.031 [+0.027, +0.036] |

Strict balance-pass subset (prespecified max |SMD| < 0.10 and both ESS >=100):

| Model | Overlap Δ random−curated among balance-pass seeds (95% seed-bootstrap CI) |
|---|---:|
| SVD | +0.011 [+0.004, +0.016] |
| NeuralMF | +0.049 [+0.037, +0.061] |
| LightGCN | +0.034 [+0.027, +0.041] |

**Boundary:** only 5/10 frozen seeds met the full prespecified balance gate. The direction among those strict-pass seeds is consistent with a residual evidence-state contribution for the learned neural models, but the analysis does not satisfy the protocol requirement for a clean isolated evidence-state claim across all seeds. It remains a sensitivity/boundary result.

Overlap weighting standardizes measured endpoint-degree exposure only; it does not prove that all structural confounding is removed or that curated pairs are clinically safe.
