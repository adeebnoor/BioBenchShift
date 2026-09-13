# Science project: structural shortcuts in biomedical AI

## Working premise

This project does **not** claim that degree bias is newly discovered. Recent work has already documented rich-node / prior / degree effects in network biology, drug-target interaction prediction, and biomedical knowledge-graph evaluation. The target contribution is broader and more falsifiable:

> **Apparent predictive intelligence in biomedical relation models can be decomposed into biological signal and structural shortcut signal; standard benchmark construction can systematically reward the latter across distinct biomedical tasks.**

The Science-level question is therefore not "does degree bias exist?" but:

1. How widespread is structural shortcut learning across biomedical relation types and model classes?
2. How much reported predictive performance survives when structural opportunity is neutralized?
3. Does structural correction change which models appear best?
4. Which evaluation design best predicts performance on independent, temporal, or genuinely novel relationships?

Anti-DDI is the **seed observation**, not the whole paper.

## Existing seed result

Anti-DDI v3.0.1 contains a reproducible ATC5 diagnostic in which a popularity-only score with no pair-specific pharmacology achieves approximately 0.90 AUC under conventional/random or curated-negative evaluation and approximately 0.50 after degree matching. This is a model-free proof that benchmark construction can make a non-biological signal look highly predictive.

The seed result must remain clearly separated from clinical non-interaction claims.

## Primary hypotheses (pre-specified)

### H1 — Cross-domain structural inflation
Across at least four independent biomedical relation families, a degree/popularity-only null will perform materially above chance under conventional random-negative evaluation and lose a substantial fraction of that performance after structure-matched evaluation.

### H2 — Model performance inflation
For learned models, conventional random-negative evaluation will overestimate performance relative to structure-matched evaluation. The magnitude of inflation will vary by model family and task.

### H3 — Model-rank instability
Correcting structural opportunity will change the rank ordering of at least some model families. A model that wins under random evaluation need not win under structure-neutralized evaluation.

### H4 — External-validity test
Performance under structure-neutralized evaluation will be more predictive than conventional random-split performance of performance on independent, temporal, cold-start, or otherwise distribution-shifted validation sets.

H4 is the key translational/generalization claim. If H4 fails, the paper should be reframed as an evaluation-bias atlas rather than a universal benchmark-correction claim.

## Benchmark families

The first-pass benchmark matrix is intentionally heterogeneous so that no single database lineage can drive the conclusion.

| Family | Relation | Preferred public source | Role |
|---|---|---|---|
| DDI | drug–drug interaction | Anti-DDI + frozen GoldD2 reference; later DDInter if licensing permits redistribution-safe processing | Seed / mechanism-aware case |
| DTI | drug–target interaction | BindingDB and/or a redistribution-safe BioSNAP-derived split | Bipartite molecular interaction |
| PPI | protein–protein interaction | HuRI; optional OGB protein-association benchmark | Network-biology replication |
| Drug–disease | indication / treatment | PrimeKG/OpenBioLink or another open relation subset | Biomedical KG task |
| Multi-relational KG | heterogeneous biomedical relations | OGB BioKG/OpenBioLink | Stress test across relation types |

A dataset enters the primary paper only if its license, version, split construction, identifiers, and full provenance can be frozen and redistributed or deterministically rebuilt.

## Evaluation designs

Every benchmark should expose the same ladder whenever technically meaningful:

1. **Random unlabelled negatives** — conventional optimistic baseline.
2. **Curated negatives** — where a defensible negative/counter-evidence set exists.
3. **Degree-matched negatives** — match endpoint degree bins without using pair-specific biology.
4. **Exact/near-exact degree balancing** — sensitivity analysis using tighter matching or propensity weighting.
5. **Endpoint-disjoint / cold-start split** — test unseen nodes where applicable.
6. **Temporal or independent holdout** — preferred external-validity endpoint.
7. **Degree-preserving rewiring/permutation** — causal diagnostic: preserve structural opportunity while disrupting biological pair identity.

No one correction is treated as perfect. Agreement across multiple interventions is the desired evidence.

## Model ladder

The paper should compare increasingly expressive models rather than only one SOTA system.

### Structure-only nulls
- endpoint degree product / log-degree product
- preferential attachment
- common-neighbor / Adamic–Adar where graph type permits

### Classical learned baselines
- matrix factorization / logistic edge model
- DistMult / ComplEx for multi-relational settings

### Modern graph models
- GraphSAGE/GCN-style link predictor for homogeneous graphs
- NBFNet or an equivalently strong relational GNN for heterogeneous KGs

### Task-specific strong model
At least one competitive task-specific architecture on DTI/DDI or drug–disease prediction, selected only after the benchmark pipeline is frozen.

The central result is **not** absolute SOTA. It is how model conclusions change under evaluation interventions.

## Primary metrics

For each task/model/split:

- AUROC
- AUPRC
- Hits@K / MRR where appropriate
- bootstrap confidence intervals over test edges/nodes
- repeated split/seeding uncertainty

### Structural inflation

`Inflation_AUC = AUC_random - AUC_structure_matched`

### Normalized collapse fraction

For models with `AUC_random > 0.55`:

`CollapseFraction = (AUC_random - AUC_structure_matched) / (AUC_random - 0.5)`

This is descriptive, not a new universal metric. Report raw performance beside it.

### Model-rank instability

Use Kendall/Spearman rank correlation and explicit pairwise rank reversals between conventional and structure-neutralized evaluations.

## Falsification / go-no-go rules

This project should be killable by data.

### Gate 1 — cross-domain signal
Proceed to the full Science build only if at least **3 of 4 independent relation families** show a reproducible structural-inflation signal in a model-free degree/popularity baseline.

### Gate 2 — learned models
Proceed only if the bias is not confined to a trivial null model; at least two learned model families must show meaningful evaluation sensitivity in at least three relation families.

### Gate 3 — scientific consequence
For a Science-level claim, we need one of the following:

- material model-rank reversals after correction, or
- strong evidence that conventional evaluation is poorly related to independent/temporal generalization while corrected evaluation is better, or
- a large and reproducible cross-domain performance collapse that changes substantive conclusions about claimed biological prediction.

If none occurs, target a strong methods/bioinformatics journal instead of forcing a Science claim.

## Figure architecture

### Figure 1 — The illusion
Anti-DDI seed result: a no-pharmacology popularity score looks excellent under conventional evaluation and collapses after degree matching. Include the causal logic, not clinical claims.

### Figure 2 — Cross-domain atlas
Same structural-null test across DDI, DTI, PPI, drug–disease, and multi-relational KG tasks.

### Figure 3 — Modern models
Performance of classical and modern learned models under random vs structure-neutralized evaluation.

### Figure 4 — Rank reversals
Show which models change rank and which task families are most vulnerable.

### Figure 5 — External validity
Compare random-evaluation and corrected-evaluation scores with temporal/independent/cold-start generalization.

### Figure 6 — Practical standard
A compact evaluation protocol: report structural null, structure-matched test, cold-start/temporal test, and independent validation before claiming biological relation prediction.

## Critical threats to validity

1. **Prior art:** degree/rich-node/prior bias is already known in several domains. Novelty must come from a general cross-domain causal/evaluation result, not rediscovery.
2. **Negative labels:** unknown pairs are not true negatives. The paper must distinguish unobserved, curated negative, counter-evidence, and contradicted states.
3. **Licensing:** DrugBank-derived data may not be redistributable. Prefer open data or deterministic retrieval scripts plus manifests.
4. **Leakage:** graph construction, ontological duplicates, relation inverses, homologs, and temporal contamination can inflate results independently of degree.
5. **Hyperparameter leakage:** tune only on training/validation constructed under the same prespecified protocol; never tune on independent holdout.
6. **Anti-DDI lineage:** the retracted predecessor and the withdrawn leaked classifier demonstration must remain transparently disclosed. The Science project must use only audited v3 assets and new frozen analyses.
7. **Clinical overreach:** none of the structural tests establishes drug safety or authorizes alert suppression.

## Immediate work sequence

### Phase A — model-free audit
Build a generic structural-shortcut audit that accepts a positive edge list and a candidate-negative/unlabelled pool, then reports random, curated, and degree-matched performance for degree-only scores with repeated seeds.

### Phase B — four public relation families
Freeze four redistribution-safe datasets and run the model-free audit. Do not invest in GPUs before Gate 1 is met.

### Phase C — learned-model benchmark
Only after Gate 1, add a controlled model ladder and pre-register split construction and metrics in the repository.

### Phase D — external validity
Add temporal/cold-start/independent validation. This is required for the strongest claim.

### Phase E — manuscript
Write the paper around the cross-domain result. Anti-DDI becomes Figure 1 / motivating case, not the title or central scope.

## Working title candidates

1. **Structural shortcuts inflate predictive performance across biomedical networks**
2. **When biomedical AI learns popularity instead of biology**
3. **Structural opportunity masquerades as biological prediction in biomedical networks**

The title will be frozen only after the cross-domain results are known.
