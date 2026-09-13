# Science Expansion Results Checkpoint — 13 September 2026

This checkpoint records observed expansion results after the preregistered protocols were committed. It does not modify the original frozen submission commit `1298e4c83ff4fa488a784d52219d20e706118a9e` and does not retroactively alter any gate or cutoff.

## A. BindingDB curated-article external validation

GitHub Actions run: `34746549599`  
Immutable artifact: `science-expansion-bindingdb-1700028a4dee4757c4cbda652b41db98a0157f89`  
Artifact digest: `sha256:a6486bd40bcbae7cdebbec8d066e6bf245b1fa0e03213283144483730df01740`

The analysis, checksum capture and artifact upload steps all completed successfully. The workflow was marked failed only because the final `git push` occurred after the branch had advanced with a separately committed protocol amendment. The principal result files were subsequently copied to the expansion branch without rerunning or changing the analysis.

### Frozen inputs and evidence set

- BioSNAP input SHA-256: `fb4867f90dd8b26383689c9ebce77cd2746475c585a3dff87de83bf0954c69db`
- BindingDB curated-articles TSV SHA-256: `3a5153d2d645eaba4e30aa7c9bd67260ec39eb8ee18ec681b46b4d6091f4a813`
- BindingDB rows inspected: 93,712
- mapped human single-protein quantitative rows: 607
- externally supported BioSNAP candidate pairs at <=1,000 nM: 11
- strict <=100 nM supported pairs: 2

### Primary <=1,000 nM support counts

| K | SVD | NeuralMF | LightGCN | SVD − NeuralMF |
|---:|---:|---:|---:|---:|
| 100 | 0 | 0 | 0 | 0 |
| 500 | 0 | 0 | 0 | 0 |
| 1,000 | 0 | 0 | 0 | 0 |
| 5,000 | 9 | 0 | 0 | +9 |
| 10,000 | 9 | 1 | 1 | +8 |
| 50,000 | 9 | 10 | 9 | −1 |

At the preregistered degree/popularity-control cutoffs (100, 500, 1,000), all observed support counts were zero and empirical two-sided p values were 1.0. Therefore the primary early-frontier BindingDB confirmation **did not pass**. The large SVD advantage observed at K=5,000 and 10,000 is retained as an observed secondary result but is not retroactively promoted to the preregistered early-frontier confirmatory gate.

### Pair-level audit

GitHub Actions run: `34746792813`  
Immutable artifact digest: `sha256:11469cacc365a8bbcc49239fbda83512da9ced60e357908f2139b5d2d33f88a7`

A transparent audit enumerated all 11 supported pairs and their ranks under all models. The mechanically prespecified post-2018 named-case gate returned **0 cases**. All 11 supporting BindingDB publications predated 2018 (2004–2012), so these data cannot be presented as later discoveries.

Nine of the 11 <=1,000 nM pairs were within the SVD top 5,000. Most of these were thiazide/thiazide-like diuretic interactions with carbonic-anhydrase isoforms (including CA4, CA9, CA12, CA5A and CA5B). This is biologically coherent external evidence, but it is historical biochemical evidence rather than a prospective or temporally later discovery result.

## B. Biological consequence of different hypothesis queues

GitHub Actions run: `34746905506`  
Immutable artifact: `science-expansion-biological-consequence-c80348e702ac687203bf34be2d7a6b24c966aa95`  
Artifact digest: `sha256:b16ea6d388d86073fc4c7ea4ce53a3ab5206bad8026110fe239046e68c59c393`

The preregistered Reactome + GO Biological Process analysis used the full BioSNAP target set as the custom background and FDR correction. No pathway term was manually selected or deleted.

| K | SVD unique targets | NeuralMF unique targets | Target Jaccard | SVD sig. terms | NeuralMF sig. terms | Shared terms | Term overlap coefficient | J–S divergence |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 78 | 100 | 0.1266 | 439 | 186 | 73 | 0.3925 | 0.6751 |
| 500 | 250 | 498 | 0.1930 | 734 | 846 | 223 | 0.3038 | 0.7011 |
| 1,000 | 342 | 670 | 0.2974 | 791 | 1,362 | 392 | 0.4956 | 0.5938 |

The prespecified functional-profile divergence is large at all three decision-relevant cutoffs. A source-stratified descriptive check (not a new confirmatory gate) also shows that the pattern is not generated solely by GO hierarchy: Reactome-only Jensen–Shannon divergence is approximately 0.695, 0.882 and 0.799 at K=100, 500 and 1,000, respectively.

Interpretation is deliberately bounded: the result shows that the benchmark-selected queues emphasize materially different biological programs. It does **not** prove that either queue is biologically correct or validate any individual predicted interaction.

## C. Science-only stop/go status at this checkpoint

The preregistered Science-expansion decision rule requires at least two of three strengthened evidence classes to survive without claim inflation.

1. **Independent direct BindingDB early-frontier evidence after popularity control:** NOT PASSED. The prespecified early cutoffs contained no supported pairs. A secondary mid-list difference exists but is not promoted post hoc.
2. **Biological-program prioritization difference:** PASSED. Target composition and functional enrichment profiles differ strongly at all three prespecified cutoffs.
3. **Contemporary multi-architecture winner-change / decision consequence:** NOT YET ESTABLISHED at this checkpoint.

Therefore the expanded evidence has improved the manuscript materially, but it has **not yet crossed the preregistered Science-only go threshold**. The original CTS-ready package should remain frozen and should not be rewritten to imply that this threshold has already been met.
