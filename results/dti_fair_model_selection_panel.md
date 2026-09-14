# Fair DTI model-selection panel on the GraphBAN-mapped universe

All four models are compared on the same mapped TargetDecagon relation universe, identical frozen splits, identical conventional controls, and identical degree-matched controls used by the leakage-free GraphBAN S6 challenge.

| Model | Conventional AUROC | Neutralized AUROC | Conventional AUPRC | Neutralized AUPRC |
|---|---:|---:|---:|---:|
| svd | 0.954 ± 0.004 | 0.917 ± 0.007 | 0.975 | 0.936 |
| neuralmf | 0.997 ± 0.001 | 0.908 ± 0.004 | 0.996 | 0.889 |
| lightgcn | 0.989 ± 0.001 | 0.883 ± 0.007 | 0.991 | 0.872 |
| graphban | 0.997 ± 0.000 | 0.928 ± 0.014 | 0.998 | 0.943 |

Mean matching coverage: **0.588**.
Conventional winner by mean AUROC: **graphban**.
Structure-neutralized winner by mean AUROC: **graphban**.
Full-panel winner reversal: **NO**.
Seed-wise conventional winners: ['graphban', 'neuralmf', 'graphban'].
Seed-wise neutralized winners: ['graphban', 'svd', 'graphban'].

This panel supersedes comparisons that mixed different relation universes or evaluation-control draws. It does not erase the historical baseline-ladder result; it determines the contemporary full-panel model-selection claim.
