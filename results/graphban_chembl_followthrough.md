# Fair-universe GraphBAN -> ChEMBL 37 follow-through

Mapped TargetDecagon coverage: **99.7%**. Candidate universe: **1,005,757** unknown mapped drug-target pairs.
ChEMBL primary supported candidates in this universe: **61**; strict pChEMBL >= 7: **22**.

GraphBAN is the mean-AUROC winner under both conventional and structure-neutralized evaluation in the frozen fair four-model panel. Therefore the contemporary DTI panel has **no benchmark-induced model-identity switch**; the same GraphBAN ranking is the relevant external-evidence ranking for both selection regimes.

## Primary ChEMBL support (pChEMBL >= 6)

| Model | @100 | @500 | @1,000 | @5,000 | @10,000 | @50,000 |
|---|---:|---:|---:|---:|---:|---:|
| SVD | 2 | 3 | 4 | 10 | 12 | 20 |
| NeuralMF | 0 | 0 | 1 | 9 | 13 | 28 |
| LightGCN | 0 | 1 | 4 | 7 | 9 | 28 |
| GraphBAN | 0 | 2 | 4 | 11 | 12 | 26 |

## GraphBAN ranking stability

Leave-one-fit-out turnover relative to the five-fit ensemble: HT@100 mean **0.136**, HT@500 **0.105**, HT@1000 **0.095**.

**Interpretation boundary.** The full contemporary DTI panel does not support a benchmark-induced winner reversal: GraphBAN remains selected under both evaluation regimes despite a material performance drop after structural neutralization. ChEMBL therefore characterizes the independently supported content of the common GraphBAN discovery frontier rather than a benchmark-induced DTI switch. The PPI experiment remains the clean full chain from benchmark rule to winner reversal, hypothesis turnover, and later-evidence consequence.
