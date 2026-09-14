# Fair DTI model-selection panel on the GraphBAN-mapped universe

All four models are compared on the same mapped TargetDecagon relation universe, identical frozen splits, identical conventional controls, and identical degree-matched controls used by the leakage-free GraphBAN S6 challenge.

| Model | Conventional AUROC | Neutralized AUROC | Conventional AUPRC | Neutralized AUPRC |
|---|---:|---:|---:|---:|
| svd | 0.951 ± 0.006 | 0.915 ± 0.006 | 0.973 | 0.936 |
| neuralmf | 0.996 ± 0.001 | 0.907 ± 0.005 | 0.996 | 0.887 |
| lightgcn | 0.989 ± 0.001 | 0.881 ± 0.006 | 0.991 | 0.869 |
| graphban | 0.997 ± 0.000 | 0.927 ± 0.012 | 0.997 | 0.942 |

Mean matching coverage: **0.592**.
Conventional winner by mean AUROC: **graphban**.
Structure-neutralized winner by mean AUROC: **graphban**.
Full-panel winner reversal: **NO**.
Seed-wise conventional winners: ['graphban', 'neuralmf', 'graphban', 'graphban', 'graphban'].
Seed-wise neutralized winners: ['graphban', 'svd', 'graphban', 'graphban', 'graphban'].

This panel supersedes comparisons that mixed different relation universes or evaluation-control draws. It does not erase the historical baseline-ladder result; it determines the contemporary full-panel model-selection claim.
