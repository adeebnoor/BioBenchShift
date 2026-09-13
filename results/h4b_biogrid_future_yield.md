# H4b — non-circular prospective discovery yield

Frozen 2025 candidate universe: **70,041,100** unknown human protein pairs among historical nodes; **5,635** became BioGRID MV-Physical edges by 5.0.261.

| Model | K | Later-added hits | Precision lower bound | Recall of later-added edges | Enrichment vs uniform |
|---|---:|---:|---:|---:|---:|
| SVD | 100 | 7 | 0.0700 | 0.0012 | 870.1× |
| SVD | 500 | 18 | 0.0360 | 0.0032 | 447.5× |
| SVD | 1,000 | 26 | 0.0260 | 0.0046 | 323.2× |
| SVD | 5,000 | 87 | 0.0174 | 0.0154 | 216.3× |
| SVD | 10,000 | 159 | 0.0159 | 0.0282 | 197.6× |
| SVD | 50,000 | 522 | 0.0104 | 0.0926 | 129.8× |
| NeuralMF | 100 | 0 | 0.0000 | 0.0000 | 0.0× |
| NeuralMF | 500 | 1 | 0.0020 | 0.0002 | 24.9× |
| NeuralMF | 1,000 | 5 | 0.0050 | 0.0009 | 62.1× |
| NeuralMF | 5,000 | 22 | 0.0044 | 0.0039 | 54.7× |
| NeuralMF | 10,000 | 42 | 0.0042 | 0.0075 | 52.2× |
| NeuralMF | 50,000 | 181 | 0.0036 | 0.0321 | 45.0× |
| LightGCN | 100 | 0 | 0.0000 | 0.0000 | 0.0× |
| LightGCN | 500 | 4 | 0.0080 | 0.0007 | 99.4× |
| LightGCN | 1,000 | 5 | 0.0050 | 0.0009 | 62.1× |
| LightGCN | 5,000 | 30 | 0.0060 | 0.0053 | 74.6× |
| LightGCN | 10,000 | 66 | 0.0066 | 0.0117 | 82.0× |
| LightGCN | 50,000 | 241 | 0.0048 | 0.0428 | 59.9× |

No degree matching or future-control sampling is used in this ranking analysis. Every model ranks the same full historical non-edge universe. Later-added BioGRID edges are treated as future evidence; persistent-unobserved pairs are not called biological negatives.
