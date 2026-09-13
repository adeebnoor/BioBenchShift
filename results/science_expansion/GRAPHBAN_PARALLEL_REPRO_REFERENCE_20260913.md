# GraphBAN-style S6 parallel reproducibility reference — 13 September 2026

This reference was recorded on the Science expansion branch before the current parallel verification run completed. It copies the already-frozen S6 target values from baseline submission commit `1298e4c83ff4fa488a784d52219d20e706118a9e`. The current run is therefore a reproducibility check of an existing result, not a new directional experiment.

## Frozen baseline identity

- Benchmark: BioSNAP TargetDecagon DTI
- Architecture label: leakage-free GraphBAN-style transductive GraphSAGE + ChemBERTa + ESM-1b
- Positive edges: 18,690
- Mapped positive edges: 18,631
- Mapping edge fraction: 0.9968432316746924
- Seeds: 0, 1, 2
- Held-out positive edges in message passing: false
- Epochs: 20
- Hidden dimension: 256

## Frozen aggregate target values

- AUROC conventional mean: 0.9972193162808759
- AUROC conventional SD: 0.0001582502911446738
- AUROC degree-matched mean: 0.9278350215958379
- AUROC degree-matched SD: 0.014142067145034874
- AUROC conventional-minus-matched mean: 0.069384294685038
- AUPRC conventional mean: 0.9975797767947417
- AUPRC degree-matched mean: 0.9425010059661556
- AUPRC conventional-minus-matched mean: 0.05507877082858607
- Mean matching coverage: 0.5880116640225057

## Frozen seed-level target values

| seed | n_train | n_test_mapped | n_seen | n_matched | match_fraction | AUROC random | AUPRC random | AUROC matched | AUPRC matched | AUROC drop | AUPRC drop |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 14904 | 3727 | 3129 | 1929 | 0.6164908916586769 | 0.9970762894891151 | 0.9975422173389153 | 0.9440215789076228 | 0.9563752026085196 | 0.05305471058149225 | 0.041167014730395746 |
| 1 | 14912 | 3719 | 3115 | 1813 | 0.5820224719101124 | 0.9973893215915327 | 0.997729156683638 | 0.917872362045398 | 0.9338487184831499 | 0.07951695954613469 | 0.06388043820048817 |
| 2 | 14904 | 3727 | 3144 | 1778 | 0.5655216284987278 | 0.9971923377619796 | 0.9974679563616715 | 0.9216111238344925 | 0.9372790968067972 | 0.07558121392748707 | 0.0601888595548743 |

## Comparison rule for the current parallel run

The parallel run must be compared against these values without choosing whichever execution is more favorable. Exact agreement is preferred. Small numerical differences attributable to deterministic feature materialization / library execution are reported quantitatively. Material disagreement triggers investigation of execution fidelity; it cannot be resolved by selecting the result that better supports the manuscript.

The architecture-fidelity label remains **GraphBAN-style leakage-free adaptation**, not a verbatim full GraphBAN reproduction.
