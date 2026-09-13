# Science expansion — results checkpoint D (13 September 2026)

Frozen baseline: `1298e4c83ff4fa488a784d52219d20e706118a9e`.

## G. Contemporary feature-rich challenge changes the claim boundary

The primary contemporary challenge was already frozen in the baseline repository as `results/graphban_targetdecagon_clean_summary.json`. It is a leakage-free **GraphBAN-style** transductive graph model using ChemBERTa compound features, ESM-1b protein features and graph message passing. It is not claimed to be a verbatim reproduction of the published GraphBAN implementation.

Frozen result across seeds 0,1,2:

| Metric | Conventional random controls | Structure-neutralized controls | Difference |
|---|---:|---:|---:|
| AUROC | 0.99722 ± 0.00016 | 0.92784 ± 0.01414 | -0.06938 |
| AUPRC | 0.99758 | 0.94250 | -0.05508 |

Mapping coverage was 0.99684 of BioSNAP TargetDecagon edges; mean matching coverage was 0.58801; held-out positive edges were excluded from message passing.

For the controlled structural DTI ladder, the frozen simple-model means were:
- SVD: AUROC 0.950 conventional, 0.914 matched;
- NeuralMF: 0.997 conventional, 0.908 matched;
- LightGCN: 0.989 conventional, 0.883 matched.

The feature-rich GraphBAN-style model therefore remains materially benchmark-sensitive, but its matched AUROC (0.92784) exceeds the matched AUROC of the simple structural models. Consequently, the older NeuralMF -> SVD winner reversal does **not** transfer straightforwardly to the expanded feature-rich panel.

## Interpretation boundary

This is a strengthening result for mechanism, not for universality. It shows that benchmark/control-population design can materially change apparent performance even for a contemporary feature-aware graph architecture, while richer molecular/protein features can reduce susceptibility enough to alter the model-selection consequence.

The manuscript must therefore distinguish two claims:

1. **Controlled structural ladder:** evaluation design can reverse the selected model and almost completely change the resulting hypothesis queue.
2. **Contemporary feature-rich challenge:** evaluation design still changes apparent performance, but the feature-rich model remains comparatively robust and the simple-model winner reversal should not be generalized to all modern architectures.

No GraphBAN hypothesis-turnover analysis is promoted unless the frozen model-selection rule identifies different contemporary winners under the two evaluation regimes. The present evidence does not justify such a claim.

## Science implication

The central paper should no longer be framed as a universal 'benchmark changes the winner' result. The stronger and more defensible general statement is that **benchmark construction is part of the scientific allocation decision**: its effect on model ranking and downstream hypotheses is substantial but depends on both the control definition and the model class.
