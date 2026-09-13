# H4 pilot — BioGRID temporal external validity

Historical release: **5.0.250**; later release: **5.0.261**.
Historical human multi-validated physical network: **93,146 edges / 11,844 nodes**.
Closed-world later-added edges: **5,635**.

## Model selection and future performance

| Model | Internal random AUC | Internal degree-matched AUC | Future AUC vs degree-matched persistent controls |
|---|---:|---:|---:|
| LightGCN | 0.908 | 0.595 | 0.589 |
| NeuralMF | 0.931 | 0.552 | 0.549 |
| SVD | 0.877 | 0.766 | 0.703 |

Conventional benchmark winner: **NeuralMF**.
Structure-neutralized winner: **SVD**.
Future matched-control AUC difference (neutralized-selected minus conventional-selected): **+0.154**.
Pilot verdict: **NEUTRALIZED_WINS_FUTURE**.

Future-added BioGRID edges provide temporal external evidence, not an exhaustive truth set; persistent unobserved pairs are not asserted to be true negatives.
