# Science expansion — results checkpoint C (13 September 2026)

Frozen baseline: `1298e4c83ff4fa488a784d52219d20e706118a9e`.

## F. Open Targets therapeutic-area context

Protocol committed before outcome inspection: `OPEN_TARGETS_DISEASE_AREA_AMENDMENT_20260913.md`.

Immutable workflow artifact:
- workflow run `34747352215`
- artifact id `10313809941`
- digest `sha256:d162241126044b6337248726944082bddb2e54271864ed23d3fb05fead566a11`
- workflow head `5be39822517fb18eca1df19aa8b976de56836565`

Release validation performed inside the run:
- Open Targets API version **26.6.3**
- Open Targets data version **26.06**
- 780/780 unique queue targets mapped to Ensembl
- 25 high-level therapeutic areas
- 2,000 target-bootstrap replicates per K

The preregistered equal-target therapeutic-area profiles were highly similar despite the much larger Reactome/GO pathway differences reported in checkpoint A.

| K | SVD annotated targets | NeuralMF annotated targets | JS divergence | Bootstrap median JS | 95% bootstrap interval | Top-5 area overlap coefficient |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 78/78 | 100/100 | 0.00432 | 0.00652 | 0.00324–0.01132 | 0.60 |
| 500 | 250/250 | 497/498 | 0.00186 | 0.00263 | 0.00155–0.00418 | 0.60 |
| 1,000 | 342/342 | 669/670 | 0.00155 | 0.00215 | 0.00127–0.00336 | 0.60 |

## Interpretation

This is a **negative/scale-dependent result**. The benchmark-selected queues differ strongly in target identity and fine-grained Reactome/GO functional enrichment, but those differences do not propagate to large changes in the coarse Open Targets therapeutic-area distribution under the preregistered equal-target aggregation.

The manuscript therefore must not claim that the two benchmark designs redirect attention between broad disease areas. The defensible biological-consequence statement is narrower and stronger: they redirect the concrete targets and functional biological programs encountered first, while much broader therapeutic-area allocation remains similar.

This negative result is retained because it defines the scale at which benchmark-induced redirection is observed and prevents overgeneralization.
