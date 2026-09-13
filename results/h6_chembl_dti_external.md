# Independent DTI external validation — ChEMBL 37

ChEMBL release: **ChEMBL_37**. Models were fitted on BioSNAP only; ChEMBL was opened only after rankings were defined according to the frozen protocol.

Mapping coverage: **264/284 drugs (93.0%)** and **1897/3648 genes (52.0%)**.
Primary externally supported unknown BioSNAP candidate pairs (binding, confidence 9, human single protein, pChEMBL >= 6): **66**.

| Model | K | ChEMBL-supported hits | Precision lower bound | Recall of supported pairs | Enrichment vs uniform |
|---|---:|---:|---:|---:|---:|
| NeuralMF | 100 | 0 | 0.0000 | 0.0000 | 0.0× |
| NeuralMF | 500 | 0 | 0.0000 | 0.0000 | 0.0× |
| NeuralMF | 1,000 | 2 | 0.0020 | 0.0303 | 30.8× |
| NeuralMF | 5,000 | 12 | 0.0024 | 0.1818 | 37.0× |
| NeuralMF | 10,000 | 19 | 0.0019 | 0.2879 | 29.3× |
| NeuralMF | 50,000 | 35 | 0.0007 | 0.5303 | 10.8× |
| SVD | 100 | 3 | 0.0300 | 0.0455 | 462.4× |
| SVD | 500 | 5 | 0.0100 | 0.0758 | 154.1× |
| SVD | 1,000 | 6 | 0.0060 | 0.0909 | 92.5× |
| SVD | 5,000 | 13 | 0.0026 | 0.1970 | 40.1× |
| SVD | 10,000 | 15 | 0.0015 | 0.2273 | 23.1× |
| SVD | 50,000 | 23 | 0.0005 | 0.3485 | 7.1× |
| LightGCN | 100 | 0 | 0.0000 | 0.0000 | 0.0× |
| LightGCN | 500 | 2 | 0.0040 | 0.0303 | 61.7× |
| LightGCN | 1,000 | 5 | 0.0050 | 0.0758 | 77.1× |
| LightGCN | 5,000 | 11 | 0.0022 | 0.1667 | 33.9× |
| LightGCN | 10,000 | 15 | 0.0015 | 0.2273 | 23.1× |
| LightGCN | 50,000 | 33 | 0.0007 | 0.5000 | 10.2× |

Primary confirmatory contrast was frozen before outcome inspection: SVD (structure-neutralized-selected) minus NeuralMF (conventional-selected) supported-hit yield at prespecified K. Unlabeled ChEMBL pairs are not called negatives.
