# Science expansion — results checkpoint A (13 September 2026)

Frozen baseline: `1298e4c83ff4fa488a784d52219d20e706118a9e`.

This checkpoint records expansion outcomes without changing the original frozen submission evidence.

## A. BindingDB curated-article direct DTI evidence

Immutable workflow artifact: GitHub Actions artifact `10313888519`, digest `sha256:a6486bd40bcbae7cdebbec8d066e6bf245b1fa0e03213283144483730df01740`.

Inputs:
- BindingDB September 2026 article-curated subset, excluding the dedicated BindingDB ChEMBL-import subset from the primary source lineage.
- extracted TSV SHA-256: `3a5153d2d645eaba4e30aa7c9bd67260ec39eb8ee18ec681b46b4d6091f4a813`;
- BioSNAP input SHA-256: `fb4867f90dd8b26383689c9ebce77cd2746475c585a3dff87de83bf0954c69db`.

Mapped quantitative human single-protein evidence produced **11** previously unlabelled BioSNAP candidate pairs at <=1,000 nM and **2** at <=100 nM.

Frozen-ranking recovery among the 11 primary pairs:

| K | SVD | NeuralMF | LightGCN |
|---:|---:|---:|---:|
| 100 | 0 | 0 | 0 |
| 500 | 0 | 0 | 0 |
| 1,000 | 0 | 0 | 0 |
| 5,000 | 9 | 0 | 0 |
| 10,000 | 9 | 1 | 1 |
| 50,000 | 9 | 10 | 9 |

Interpretation: this source does **not** reproduce the ChEMBL early-frontier result. It instead shows a pronounced mid-rank concentration difference that disappears at broad cutoffs. No post-2018 mechanically eligible Top-1,000 named case was found. Absence from BindingDB is not a negative label.

A transparent pair audit (artifact `10313918696`, digest `sha256:11469cacc365a8bbcc49239fbda83512da9ced60e357908f2139b5d2d33f88a7`) further shows that much of the SVD Top-5,000 recovery is clustered in carbonic-anhydrase interactions from a shared older publication. The 9/11 count must therefore **not** be described as nine independent discoveries.

## B. Biological-program consequence of the competing queues

Immutable workflow artifact: GitHub Actions artifact `10314526509`, digest `sha256:b16ea6d388d86073fc4c7ea4ce53a3ab5206bad8026110fe239046e68c59c393`.

g:Profiler version recorded by the run: `e114_eg62_p19_27110d83`; custom background = 3,648 BioSNAP targets; sources = Reactome + GO Biological Process; FDR < 0.05; no manual term selection.

| K | SVD unique targets | NeuralMF unique targets | Target Jaccard | SVD sig. terms | NeuralMF sig. terms | Shared terms | Term overlap coefficient | JS divergence |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 78 | 100 | 0.1266 | 439 | 186 | 73 | 0.3925 | 0.6751 |
| 500 | 250 | 498 | 0.1930 | 734 | 846 | 223 | 0.3038 | 0.7011 |
| 1,000 | 342 | 670 | 0.2974 | 791 | 1,362 | 392 | 0.4956 | 0.5938 |

Interpretation: the benchmark-selected queues emphasize materially different target sets and different functional-enrichment profiles. This is evidence of **different biological prioritization**, not evidence that one queue is biologically correct or that an enriched pathway validates an individual pair.

## C. Pending preregistered gates

- 10-seed matching-method robustness: running under `MATCHING_ROBUSTNESS_AMENDMENT_20260913.md`.
- Named BioGRID later-supported case extraction: protocol committed before identities were inspected as `BIOGRID_NAMED_LATER_CASES_AMENDMENT_20260913.md`.
- Open Targets disease-area context: not yet opened; requires a separate aggregation amendment before outcome inspection.
- Contemporary architecture ladder: not yet opened. Modern hypothesis turnover will only be computed if the modern benchmark winner changes under the common evaluation contract.
- Prospective wet-lab validation: unavailable unless a real experimental collaborator, assay protocol and raw measurements are obtained; no such evidence will be simulated or implied.

## Current Science interpretation

Expansion B is a clear positive evidence class: the benchmark-induced model change produces substantially different biological-program prioritization. Expansion A is mixed/supportive rather than decisive because the direct BindingDB evidence is sparse, absent from the first 1,000 ranks, and publication-clustered. The strongest remaining opportunity for the Science-only case is therefore the preregistered named temporal BioGRID case set plus the contemporary-model and matching-robustness gates.
