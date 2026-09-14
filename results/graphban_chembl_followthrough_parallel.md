# Fair-universe GraphBAN -> ChEMBL 37 follow-through

Mapped TargetDecagon coverage: **99.7%**. Candidate universe: **1,005,757** unknown mapped pairs.

GraphBAN is selected under both frozen DTI evaluation regimes; therefore the contemporary DTI panel has **no benchmark-induced model-identity switch**.

| Model | @100 | @500 | @1,000 | @5,000 | @10,000 | @50,000 |
|---|---:|---:|---:|---:|---:|---:|
| SVD | 2 | 3 | 4 | 10 | 12 | 20 |
| NeuralMF | 0 | 0 | 1 | 9 | 13 | 28 |
| LightGCN | 0 | 1 | 4 | 7 | 9 | 28 |
| GraphBAN | 2 | 2 | 4 | 12 | 13 | 24 |

GraphBAN leave-one-fit-out turnover: HT@100 **0.170**, HT@500 **0.115**, HT@1000 **0.105**.

**Boundary:** ChEMBL characterizes the independently supported content of the common GraphBAN frontier; it does not create a DTI benchmark-induced switch because GraphBAN is the selected model under both evaluation rules. The historical PPI experiment remains the full benchmark -> winner -> hypotheses -> later-evidence chain.
