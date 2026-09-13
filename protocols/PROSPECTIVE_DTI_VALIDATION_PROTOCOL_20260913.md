# Prospective DTI Validation Protocol — Science Expansion

**Frozen:** 13 September 2026  
**Branch:** `science-expansion-20260913`  
**Immutable baseline:** `1298e4c83ff4fa488a784d52219d20e706118a9e`

## Scientific objective

Test prospectively whether two defensible benchmark designs, by selecting different models, direct a fixed experimental budget toward different drug–target hypotheses and therefore toward different experimentally validated biology.

This protocol is committed **before any prospective wet-lab result is available**. Null, mixed, and adverse results will be retained.

## Primary causal chain

`benchmark definition -> selected model -> top-B hypothesis queue -> blinded experiment -> validated discoveries`

The intervention is the benchmark/evaluation rule used upstream of model selection. The experiment does not assume that structural neutralization is universally superior.

## Primary decision policies

The primary prospective comparison preserves the already frozen DTI model-selection result:

- **Conventional policy:** model selected under random-unlabelled evaluation = NeuralMF.
- **Structure-neutralized policy:** model selected under the original joint endpoint-degree-matched evaluation = SVD.

The model identities above are frozen prior results and will not be reselected using prospective assay outcomes.

A separate contemporary-model panel is treated as a robustness/moderation experiment. It does not replace the primary prospective contrast after results are known.

## Experimental budget

The primary budget is **B = 100 assay slots per policy**.

- Conventional policy: top 100 assay-eligible novel candidate pairs.
- Structure-neutralized policy: top 100 assay-eligible novel candidate pairs.
- If the two lists overlap, shared pairs are assayed once but count toward both policy budgets.
- Each policy also has a mechanically ordered reserve list of 25 candidates for pre-assay failures such as compound unavailability, failed protein production, or platform incompatibility.
- Reserve substitution follows rank order only and must occur before arm unblinding by the laboratory.

The top-100 budget is chosen because HT@100 is a frozen primary decision-scale endpoint in the computational study. With 100 candidates per policy, the design has approximately 81% asymptotic two-sided power at alpha=0.05 to distinguish a 10% from a 25% validation rate; the study is not powered only around that one effect size and all confidence intervals will be reported.

## Prospective novelty and assay-eligibility universe

The assay-eligible universe is defined **before score-based selection**. A candidate pair is eligible only if all of the following hold:

1. human drug–target pair maps unambiguously to the frozen BioSNAP identifiers;
2. the pair is absent from the frozen BioSNAP positive edge set;
3. no qualifying direct interaction is present in the frozen ChEMBL 37 or BindingDB 2026-09 evidence states used by this project before panel lock;
4. compound can be sourced at >=95% stated purity or equivalent CRO quality standard;
5. target has a validated direct biochemical/binding assay available through the selected independent laboratory/CRO;
6. assay conditions permit a quantitative concentration-response measurement or direct binding estimate;
7. the pair is not excluded by a prespecified safety, handling, or solubility feasibility criterion supplied by the laboratory before it sees model identity.

The laboratory may define an assay catalog before scoring. Restricting to that catalog is allowed if the catalog is frozen before model ranks are opened.

## Candidate selection

Within the assay-eligible universe:

1. rank all candidates by the frozen conventional-policy model score and take the first 100;
2. independently rank all candidates by the frozen structure-neutralized-policy model score and take the first 100;
3. do not hand-pick candidates for biological attractiveness;
4. do not remove a candidate because its prediction seems implausible;
5. all exclusions and reserve substitutions must carry a machine-readable reason code.

For each selected pair, retain both policy ranks and both model scores so that cross-policy counterfactual rank can be reported after unblinding.

## Blinding

The wet-lab/CRO receives only:

- coded sample identifier (`BBX-P####`),
- compound identity needed to execute the assay,
- target identity needed to execute the assay,
- assay instructions.

It does **not** receive:

- benchmark arm,
- model name,
- prediction score,
- policy rank,
- expected direction of the result.

The computational arm key is stored separately from the laboratory result table. Raw laboratory results are frozen and checksummed before the arm key is joined back to the data.

## Assay hierarchy

### Primary assay

Use a direct biochemical or biophysical interaction assay appropriate to the target class, such as SPR, BLI, MST, competition binding, radioligand binding, biochemical enzyme binding/activity with direct target attribution, or an equivalently validated CRO platform.

The same target should use the same primary assay format across both policy arms whenever technically possible.

### Confirmation

Every primary positive must undergo an orthogonal confirmation assay when a technically valid orthogonal method exists. Cell-based functional confirmation is strongly preferred for biologically actionable positives but is a secondary endpoint unless prespecified for a target class before unblinding.

## Primary validation rule

A pair is a **primary validated interaction** when all prespecified platform QC criteria pass and a reproducible concentration-dependent interaction is observed at a potency/affinity equivalent to <=1,000 nM in the primary assay.

A stricter <=100 nM threshold is a frozen sensitivity analysis.

If a platform does not report nM-equivalent affinity/potency, that target class must have its binary/quantitative rule frozen in a laboratory appendix before arm unblinding and must not be pooled into the primary <=1,000 nM endpoint unless scientifically commensurate.

## Primary endpoints

For each policy at B=100 report:

1. number and fraction of primary validated interactions;
2. risk difference and risk ratio with 95% confidence intervals;
3. Fisher exact test for validation-rate difference, treated as one inferential summary rather than the sole success criterion;
4. number of validated interactions passing orthogonal confirmation;
5. cross-policy rank of every validated interaction;
6. **counterfactual experimental miss:** for each validated pair selected by one policy, whether the competing policy would have placed it outside its top-100 budget;
7. number of distinct targets and compounds represented among validated interactions.

## Scientific-program endpoints

Among validated interactions, compare the two policy portfolios using prespecified mappings:

- target-set Jaccard similarity;
- Reactome and GO Biological Process profiles using the same frozen enrichment pipeline as the computational study;
- Jensen–Shannon divergence of normalized enrichment profiles;
- high-level Open Targets disease-area profiles as interpretive context only.

The claim of redirected biological attention is made at the finest level directly supported by the data. Similarity at a broad therapeutic-area scale will be reported if present.

## Statistical analysis

- Analyze all successfully completed assay slots under the policy that selected them.
- Shared candidates count once experimentally but contribute to both policy portfolios.
- Missing assays are reported with reason codes; no outcome-dependent replacement is permitted.
- Primary binary validation uses exact binomial confidence intervals within each policy and Fisher exact comparison between policies.
- Rank/yield curves are reported over available prefix budgets where feasible.
- Enrichment/profile comparisons retain multiple-testing correction exactly as in the frozen biological-consequence analysis.
- No candidate, target class, or assay platform may be removed post hoc to improve the headline contrast.

## Independent execution requirement

The preferred design uses an independent laboratory or CRO that did not develop the computational models and is blinded to policy identity. Raw instrument exports, plate maps, QC flags, concentration-response tables, and analysis scripts must be retained.

No collaborator will be listed as an author solely for credibility; authorship requires a qualifying contribution and manuscript approval.

## Contemporary-model extension

If the exact published contemporary-model ladder changes winner identity between conventional and neutralized evaluation, an additional **secondary** prospective panel of 25 candidates per modern selected policy is opened under the same blinding and assay rules.

If modern winner identity does not change, that null result is retained and no post hoc modern experimental contrast is manufactured.

## Stop/go interpretation

The wet-lab experiment is considered scientifically informative regardless of direction.

- Different validated portfolios at equal budget support the claim that benchmark design changes experimentally realized scientific attention.
- Similar validation yield but biologically distinct validated portfolios support a scale-dependent prioritization effect.
- Similar yield and similar validated portfolios weaken the discovery-redirection claim and must narrow the manuscript.
- Better conventional-policy yield is not suppressed; it would show that structural neutralization can trade conventional predictive utility for a different search frontier.

## Data release

After laboratory unblinding and collaborator approval, release:

- frozen candidate manifest;
- arm-blind laboratory manifest;
- raw and processed assay results where licensing permits;
- arm key;
- analysis code;
- checksums and commit identifiers;
- all exclusion and substitution reason codes.

No wet-lab conclusion will be added to the manuscript before the raw result package is frozen and auditable.
