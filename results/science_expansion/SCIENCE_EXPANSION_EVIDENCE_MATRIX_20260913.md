# Science expansion evidence matrix — 13 September 2026

Baseline frozen submission commit: `1298e4c83ff4fa488a784d52219d20e706118a9e`

This file synthesizes expansion results without changing the frozen Science submission. It records supportive evidence, boundary conditions, and prohibited over-interpretations. The contemporary GraphBAN challenge remains pending and must be incorporated regardless of direction.

| Evidence stream | Prespecified result | What it supports | Boundary / what it does not support |
|---|---|---|---|
| BioGRID temporal future evidence | Seven later-added human physical interactions fall in SVD ranks 34–81 but NeuralMF ranks 10,815–48,996; all seven were retained by a mechanical, all-cases extraction rule. | Benchmark choice can radically change which subsequently supported relations enter an early experimental queue. | Later-added BioGRID records are temporal support, not exhaustive truth and not proof that persistent-unobserved pairs are false. Do not call these interactions discoveries by this study. |
| Functional consequence (g:Profiler) | Target Jaccard is 0.127, 0.193, 0.297 at K=100,500,1000; enrichment-profile JS divergence is 0.675, 0.701, 0.594. | Different benchmark-selected queues emphasize substantially different target sets and biological programs. | Enrichment does not establish causal correctness or validate individual predicted pairs. |
| Open Targets therapeutic-area context | Disease-area JS divergence is only 0.00432, 0.00186, 0.00155 at K=100,500,1000, with top-5 therapeutic-area overlap coefficient 0.6 throughout. | The strongest redirection is at target/program resolution rather than broad therapeutic-area allocation. | Do not claim that benchmark choice broadly redirects disease domains in this analysis. |
| BindingDB curated-article external support | 11 candidate pairs supported at <=1,000 nM (2 at <=100 nM). No model retrieves a supported pair in Top-1,000. At Top-5,000 SVD retrieves 9/11 whereas NeuralMF and LightGCN retrieve 0/11; by Top-50,000 NeuralMF retrieves 10/11 and SVD/LightGCN 9/11. | Independent orthogonal evidence that supported pairs can occupy very different rank regions across models. | Sparse and clustered support; pair audit shows non-independence/concentration. Not evidence for a general early-frontier advantage and not a basis for claiming absent BindingDB pairs are negatives. |
| Matching-definition robustness | With log2 degree bins, the neutralized winner changes in 7/10 seeds (SVD 7, NeuralMF 3; mean matching coverage 0.599). With degree deciles and nearest-caliper matching, winner change is 0/10; NeuralMF remains winner in all seeds (coverage 1.000 and 0.485 respectively). | The benchmark-framing effect has a measurable dependence on the exact structural-neutralization design; sensitivity itself is an important benchmark-design boundary condition. | Do not claim invariant winner reversal across all reasonable degree-matching schemes. |
| Contemporary GraphBAN challenge | **PENDING — retain outcome regardless of direction.** Frozen ChemBERTa drug representations + ESM-1b protein representations, bilinear interaction head, three seeds, random and cold-both evaluation. | If cold-both degradation is large, stronger representations do not remove framing sensitivity; if small, GraphBAN is a genuine boundary case showing reduced sensitivity. | No architecture, split, representation, or interpretation rule may be changed after outcome inspection to obtain a preferred direction. |

## Claim hierarchy after completed expansion streams

**Strongly supported now:** benchmark design can redirect the identity of high-priority relations, targets, and functional programs; retrospective model ranking need not identify the model whose earliest queue contains later-supported relations.

**Supported with qualification:** independent BindingDB evidence shows rank-location differences, but the validation set is sparse and clustered. Structural-neutralization conclusions depend on the matching definition.

**Not supported:** a universal claim that one model is biologically correct; a claim that every neutralization scheme reverses model ranking; a claim that benchmark choice substantially changes broad therapeutic-area allocation; or a claim that absent database records are biological negatives.

## Locked synthesis rule

The final expansion narrative must include all completed prespecified analyses, including null, mixed, or contradictory outcomes. GraphBAN will be added to this matrix when its frozen S6 run finishes; its direction cannot trigger replacement by another contemporary model.
