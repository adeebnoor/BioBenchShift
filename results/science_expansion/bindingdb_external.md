# BindingDB curated-article external validation

BindingDB input SHA-256: `3a5153d2d645eaba4e30aa7c9bd67260ec39eb8ee18ec681b46b4d6091f4a813`
BioSNAP input SHA-256: `fb4867f90dd8b26383689c9ebce77cd2746475c585a3dff87de83bf0954c69db`
Primary externally supported candidate pairs (<=1,000 nM): **11**.
Strict <=100 nM supported candidate pairs: **2**.

## Primary ranking endpoint

| Model | K | BindingDB-supported hits | Precision lower bound | Recall supported | Enrichment vs uniform |
|---|---:|---:|---:|---:|---:|
| NeuralMF | 100 | 0 | 0.00000 | 0.00000 | 0.00x |
| NeuralMF | 500 | 0 | 0.00000 | 0.00000 | 0.00x |
| NeuralMF | 1,000 | 0 | 0.00000 | 0.00000 | 0.00x |
| NeuralMF | 5,000 | 0 | 0.00000 | 0.00000 | 0.00x |
| NeuralMF | 10,000 | 1 | 0.00010 | 0.09091 | 9.25x |
| NeuralMF | 50,000 | 10 | 0.00020 | 0.90909 | 18.50x |
| SVD | 100 | 0 | 0.00000 | 0.00000 | 0.00x |
| SVD | 500 | 0 | 0.00000 | 0.00000 | 0.00x |
| SVD | 1,000 | 0 | 0.00000 | 0.00000 | 0.00x |
| SVD | 5,000 | 9 | 0.00180 | 0.81818 | 166.47x |
| SVD | 10,000 | 9 | 0.00090 | 0.81818 | 83.24x |
| SVD | 50,000 | 9 | 0.00018 | 0.81818 | 16.65x |
| LightGCN | 100 | 0 | 0.00000 | 0.00000 | 0.00x |
| LightGCN | 500 | 0 | 0.00000 | 0.00000 | 0.00x |
| LightGCN | 1,000 | 0 | 0.00000 | 0.00000 | 0.00x |
| LightGCN | 5,000 | 0 | 0.00000 | 0.00000 | 0.00x |
| LightGCN | 10,000 | 1 | 0.00010 | 0.09091 | 9.25x |
| LightGCN | 50,000 | 9 | 0.00018 | 0.81818 | 16.65x |

## Degree/popularity-matched null

| Model | K | Observed | Matched-null mean | z | empirical two-sided p |
|---|---:|---:|---:|---:|---:|
| SVD | 100 | 0 | 0.00 | -0.01 | 1.00000 |
| SVD | 500 | 0 | 0.01 | -0.11 | 1.00000 |
| SVD | 1,000 | 0 | 0.03 | -0.16 | 1.00000 |
| NeuralMF | 100 | 0 | 0.00 | nan | 1.00000 |
| NeuralMF | 500 | 0 | 0.00 | nan | 1.00000 |
| NeuralMF | 1,000 | 0 | 0.01 | -0.08 | 1.00000 |

Mechanically selected later-published discordant case-study rows: **0**.

No absent BindingDB pair is interpreted as a negative interaction.
