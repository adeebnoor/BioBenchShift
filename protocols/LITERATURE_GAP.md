# Literature gap and novelty boundary

## Why this file exists

The Science project must not claim that degree bias, rich-node bias, prior bias, negative-sampling bias, temporal benchmarking, or biomedical hypothesis-generation evaluation are new. Several strong papers already establish parts of that story. The project is viable only if it connects these pieces into a new, falsifiable statement about **benchmark-induced model selection and scientific decision identity** across biomedical relation families.

## Closest prior art identified and rechecked (13 September 2026)

### 1. Generic graph link prediction

**Aiyappa et al., ICML 2025 — “Implicit degree bias in the link prediction task.”**

Core result: conventional edge-vs-random-nonedge sampling favors high-degree nodes; a degree-only null can approach optimal performance; the authors propose degree-corrected link prediction and show improved alignment with recommendation-task performance.

Implication for us: we cannot claim invention of degree-corrected link prediction. Our contribution must be biomedical, cross-domain, and tied to substantive biological/external-validation consequences.

Source: https://proceedings.mlr.press/v267/aiyappa25a.html

### 2. Network biology

**Yılmaz, Yorgancioglu & Koyutürk, PNAS 2025 — “Bias-aware training and evaluation of link prediction algorithms in network biology.”**

Core result: uniform random edge-based evaluation favors rich/high-degree nodes, including across network snapshots; the paper proposes AWARE strategies emphasizing low-degree nodes and evaluates BioGRID over time as well as STRING evidence sources.

Implication for us: PPI/network-biology degree bias **and temporal BioGRID analysis are already established**. A BioGRID time-slice result by itself is not novel. Our PPI temporal analysis matters only as one external-validity link in a broader chain: cross-domain structural inflation → model-rank reversal → hypothesis-selection turnover → later-evidence consequence.

Source: https://doi.org/10.1073/pnas.2416646122

### 3. Drug–target interaction prediction

**Lin et al., Nature Communications 2025 — TAPB.**

Core result: target prior tendency can act as a confounder in DTI datasets; models can exploit target-specific label priors rather than genuine interaction mechanisms; the authors introduce an interventional debiasing framework.

Implication for us: DTI prior bias is established. A DTI-only paper will not be novel enough. DTI should be one relation family in a broader test.

Source: https://www.nature.com/articles/s41467-025-66915-1

### 4. Biomedical knowledge-graph link prediction

**Baudot and colleagues, Bioinformatics 2026 — “Benchmarking the impact of data leakage on the performance of knowledge graph embedding models for biomedical link prediction.”**

Core result: evaluates train-test redundancy, degree-related illegitimate features, and test-set sampling/distribution shift across biomedical KGs. Degree-preserving permutation did not show that decoder-only models were driven by degree alone. Random and cold-start evaluations were also compared with an independent Orphanet indication set, with substantial external-performance loss.

Implication for us: “biomedical KG models may use degree” and “random splits overestimate external generalization” are both insufficient novelty claims. We need the model-selection and downstream-hypothesis consequence, not only another leakage audit.

Source: https://academic.oup.com/bioinformatics/article/42/8/btag608/8767164

### 5. Modern biomedical relational GNNs

**BioPathNet, Nature Biomedical Engineering 2026.**

The paper explicitly analyzes high-degree bias and uses stringent negative sampling while reporting strong performance across several biomedical relation tasks.

Implication for us: strong modern models are already bias-aware. The Science project must show either (a) residual benchmark shortcut effects even in strong contemporary models, (b) model-rank reversals under standardized structural controls, or (c) that corrected evaluation better predicts genuinely independent/temporal performance.

Source: https://www.nature.com/articles/s41551-025-01598-z

### 6. Temporal biomedical hypothesis-generation benchmarks

**Tyagin & Safro, BMC Bioinformatics 2024 — Dyport.**

Core result: proposes a dynamic biomedical hypothesis-generation benchmark using time-sliced curated knowledge, evaluates link-prediction systems under realistic temporal conditions, and introduces discovery-importance weighting to move beyond static link-prediction metrics.

Implication for us: temporal validation of generated biomedical hypotheses is not new. Our RIDI-inspired H5 should therefore **not** be sold as the invention of temporal hypothesis benchmarking. Its defensible target is narrower and different: a benchmark intervention changes which model wins; holding the candidate universe fixed, that model-selection reversal changes the identity of prioritized hypotheses; external evidence then tests the consequence of the selection decision.

Source: https://doi.org/10.1186/s12859-024-05812-8

### 7. Earlier systematic bias auditing in biological ML

**Communications Biology 2021 — “Systematic auditing is essential to debiasing machine learning in biology.”**

Core result: paired-input biological predictors can attain high in-network performance without learning generalizable biology; the work proposes systematic auditing of predictor bias.

Implication for us: the broad warning that biological ML can exploit dataset biases is established. We need quantitative, cross-domain decision consequences rather than another generic auditing manifesto.

Source: https://www.nature.com/articles/s42003-021-01674-5

## The defensible novelty target

The project should aim to establish all four levels below.

### Level A — cross-domain prevalence

Use the same prespecified model-free structural audit across independent relation families. This yields a comparable atlas rather than a collection of anecdotes.

### Level B — model-selection consequence

Demonstrate that structural correction changes conclusions about learned models, including rank ordering. Merely showing a degree-only null performs well is not enough after ICML 2025 and PNAS 2025.

### Level C — scientific decision identity

When the benchmark chooses a different winner, rank the **same frozen biological candidate universe** and quantify how much the prioritized hypothesis set changes. This is the RIDI-inspired layer. Its novelty is not ranking hypotheses per se; it is tracing a benchmark intervention through model selection into hypothesis identity.

### Level D — external consequence

Test whether the evaluation regime that selects a model is associated with better recovery of genuinely later or independent evidence. The strongest version is non-circular: select models historically, rank the full identical historical unknown-pair universe, then open a later snapshot and count later-supported relations without constructing the future test by the same structural matching rule.

If Levels A–D survive independent domains and stronger models, the result is substantially stronger than the current prior art because it connects benchmark design to the identity and eventual support of biomedical scientific decisions.

## Anti-DDI’s role

Anti-DDI remains valuable because it contains a particularly clean motivating contradiction:

- a popularity-only score with no pair-specific pharmacology gives ~0.90 AUC under conventional/curated ATC5 evaluation;
- after degree matching, the same score is ~0.50;
- the resource explicitly distinguishes evidence states rather than equating missing edges with negative biology.

But Anti-DDI should be a motivating/evidence-state case, not the scope of the paper.

## Claims we must avoid

Do **not** claim:

- “we discovered degree bias in link prediction”;
- “we invented temporal biomedical hypothesis evaluation”;
- “biomedical AI generally learns no biology”;
- “all reported biomedical link-prediction performance is invalid”;
- “degree matching is the uniquely correct evaluation design”;
- “Anti-DDI candidates are safe drug combinations.”

Any final claim must be proportional to the cross-domain, model-selection, hypothesis-identity and external-validation evidence actually obtained.
