# GraphBAN-style model-selection gate — 13 September 2026

Baseline frozen submission commit: `1298e4c83ff4fa488a784d52219d20e706118a9e`.

This gate applies the conditional rule already frozen in `protocols/CONTEMPORARY_MODEL_CHALLENGE.md`: after the contemporary challenge is frozen, add the GraphBAN-style model to the SVD/NeuralMF/LightGCN panel under the same evaluation contract; if it becomes the winner under at least one frozen regime, promote it into the conditional H5 scientific-decision comparison and the unchanged ChEMBL external-validation protocol.

No model was retrained and no outcome-dependent threshold was introduced to make this decision. Values below combine the already-frozen GraphBAN-style S6 seed results with the same-seed (`0,1,2`) log2-degree-bin results from the preregistered matching-robustness analysis.

## Same-seed AUROC comparison

| Model | Conventional mean AUROC, seeds 0–2 | Log2 degree-matched mean AUROC, seeds 0–2 |
|---|---:|---:|
| SVD | 0.953061 | 0.917098 |
| NeuralMF | 0.997213 | 0.910563 |
| LightGCN | 0.988671 | 0.884024 |
| GraphBAN-style adaptation | **0.997219** | **0.927835** |

### Seed-level caution

GraphBAN-style does not dominate every seed. Under conventional evaluation it is slightly below NeuralMF at seed 1; under matched evaluation it is slightly below SVD at seed 1. The conventional mean advantage over NeuralMF is only approximately `+0.0000064`, effectively a near-tie for scientific interpretation. The matched mean advantage over SVD is approximately `+0.01074`.

## Gate decision

By the frozen mean-AUROC model-selection convention, the GraphBAN-style adaptation is the numerical winner in **both** evaluation regimes for seeds 0–2. Therefore the conditional branch in `CONTEMPORARY_MODEL_CHALLENGE.md` is activated:

1. add the GraphBAN-style model to the H5 hypothesis-identity / Top-K decision comparison on the identical frozen unknown TargetDecagon universe;
2. quantify within-GraphBAN initialization/ensemble instability and compare it with cross-model queue turnover;
3. add the GraphBAN-style full-universe ranking to the already-frozen ChEMBL external-validation protocol without changing ChEMBL mappings, thresholds, evidence definitions or K values.

## Interpretation boundary

This gate is **not** evidence that GraphBAN-style is materially superior to NeuralMF under the conventional benchmark; the conventional mean difference is negligible. Its purpose is mechanical compliance with the prespecified conditional analysis rule. The architecture remains labeled a `GraphBAN-style leakage-free adaptation`, not a full reproduction of the published GraphBAN framework.
