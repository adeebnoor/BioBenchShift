# Frozen external DTI validation protocol — ChEMBL 37

**Protocol freeze:** 2026-09-13  
**Branch:** `main`  
**Purpose:** independent-source validation of the benchmark-selection consequence in drug–target interaction prediction.  
**Critical rule:** this protocol is frozen before inspecting external ChEMBL yield by model.

## Scientific question

BioSNAP TargetDecagon is used for model development and for choosing a winner under conventional versus structure-neutralized evaluation. ChEMBL is then opened only as an independent evidence source.

Primary question:

> When the conventionally selected and structure-neutralized-selected BioSNAP models rank the same frozen unknown drug–target universe, which ranking places more independently supported ChEMBL relations near the top?

ChEMBL information MUST NOT be used for model fitting, hyperparameter selection, benchmark winner selection, threshold selection, or candidate-universe construction except for post hoc evidence labeling after rankings are frozen.

## Frozen source

- ChEMBL release: **ChEMBL 37**.
- Query/retrieval date is recorded by the workflow.
- BioSNAP source remains the already frozen TargetDecagon edge list and checksum.
- All raw mapping/query outputs used in the final analysis are checksumed or summarized in a deterministic manifest.

## Identifier mapping

### Compounds

BioSNAP drugs are PubChem CIDs encoded as `CID#########`.

1. Strip the `CID` prefix and leading zeros to recover the integer PubChem CID.
2. Map PubChem -> ChEMBL by **UniChem**, source ID 22 (PubChem) to source ID 1 (ChEMBL).
3. Preserve all current ChEMBL mappings returned by UniChem.
4. If one PubChem CID maps to multiple ChEMBL molecule IDs, all mapped molecule IDs are eligible evidence aliases; the BioSNAP drug is supported if any eligible alias supports the target.
5. Unmapped compounds remain in the ranking universe but cannot receive a ChEMBL support label; mapping coverage is reported explicitly.

### Targets

BioSNAP targets are NCBI GeneID/Entrez identifiers.

1. Map GeneID -> UniProtKB using the official UniProt ID-mapping service.
2. Restrict ChEMBL target evidence to **Homo sapiens** (`tax_id = 9606`).
3. Map eligible human UniProt accessions to ChEMBL target components/targets.
4. Only ChEMBL targets whose target type is **SINGLE PROTEIN** are eligible for the primary analysis.
5. If a GeneID maps to multiple human UniProt entries/ChEMBL single-protein targets, all eligible mappings are retained and the BioSNAP gene is supported if any eligible mapping is supported.
6. Mapping coverage and multiplicity are reported; unmapped BioSNAP genes remain in the ranking universe but cannot receive a ChEMBL support label.

## Primary ChEMBL evidence definition

A BioSNAP drug–target candidate is **independently supported** if at least one ChEMBL 37 activity satisfies ALL of the following:

- target organism: **Homo sapiens**;
- target type: **SINGLE PROTEIN**;
- assay type: **B** (binding);
- ChEMBL target-assignment confidence score: **9** (direct single-protein target);
- `pchembl_value` is present and numeric;
- **primary potency threshold: pChEMBL >= 6.0** (approximately <= 1 micromolar for standardized molar potency values);
- activity has no explicit invalid-data flag that indicates the numeric record should be excluded;
- the activity maps to an eligible BioSNAP PubChem drug and GeneID target through the frozen mappings above.

No requirement on mechanism-of-action direction (agonist/inhibitor/etc.) is imposed because BioSNAP TargetDecagon represents drug–target interaction rather than a signed pharmacological mechanism.

## Prespecified sensitivity analyses

These are secondary and will NOT replace the primary result based on outcome favorability.

1. **Stringent potency:** pChEMBL >= 7.0 (approximately <= 100 nM).
2. **Exact-relation sensitivity:** where relation metadata permit, restrict to exact/equivalent quantitative relations rather than ambiguous bounds.
3. **Reviewed-target sensitivity:** restrict GeneID mappings to reviewed UniProtKB/Swiss-Prot where mapping information permits.
4. **Mapped-universe sensitivity:** compute metrics both on the full BioSNAP candidate universe and on the subset for which both endpoints have valid external mappings. The mapped-universe analysis is a sensitivity check; the primary ranking itself is not rebuilt after seeing ChEMBL.

## Frozen candidate universe

- Models are trained on BioSNAP only.
- The unknown candidate universe is the Cartesian product of BioSNAP drug nodes and target nodes excluding all BioSNAP-observed positive edges.
- Every compared model ranks the identical candidate universe.
- Rankings are generated before any ChEMBL evidence labels are used.

## Models and selection regimes

Primary comparison follows the already frozen BioSNAP benchmark-selection result:

- conventional-evaluation winner: **NeuralMF**;
- structure-neutralized-evaluation winner: **SVD**.

LightGCN is retained as a prespecified comparison model.

For stability, primary ranking scores should use the frozen ensemble convention already established in H5b where computationally feasible; single-fit results, if shown, are secondary and cannot replace ensemble results based on favorability.

## Primary endpoints

For each model, compute ChEMBL-supported relations among:

- Top 100;
- Top 500;
- Top 1,000;
- Top 5,000;
- Top 10,000;
- Top 50,000 or the maximum available universe if smaller.

Report:

- supported hits at K;
- precision lower bound at K (supported hits / K; called a lower bound because ChEMBL is incomplete);
- recall of externally supported candidate relations;
- enrichment versus uniform ranking over the same evaluable universe;
- cumulative hit curve;
- rank distribution / average precision using ChEMBL support only as a positive evidence label, never interpreting unlabeled pairs as biological negatives.

## Primary confirmatory contrast

The predeclared confirmatory comparison is:

> **ChEMBL-supported hit yield of the structure-neutralized-selected SVD ensemble minus the conventional-selected NeuralMF ensemble on the same frozen candidate universe.**

The direction is hypothesized to be positive, but a zero or negative result is retained and reported.

## Uncertainty

- Bootstrap externally supported relations to obtain 95% confidence intervals for hit-rate/recall contrasts at prespecified K values.
- Report exact counts alongside intervals.
- Where appropriate, use a hypergeometric/binomial enrichment calculation only as a descriptive null against uniform ranking; do not treat ChEMBL unlabeled pairs as confirmed negatives.

## Leakage / overlap safeguards

- ChEMBL is not used to train BioSNAP models.
- The analysis concerns independent database support, not temporal novelty unless a ChEMBL release-date-specific historical comparison is separately performed.
- Because BioSNAP and ChEMBL may share upstream literature or database ancestry, quantify overlap/provenance where feasible and call this **independent-source evidence**, not independent experimental validation.
- A stronger orthogonal validation layer (e.g., manually verified primary experiments or another non-ChEMBL source) remains desirable for a flagship Science claim.

## Failure and interpretation rules

A negative ChEMBL result does not invalidate H1–H5 or the BioGRID temporal result; it caps the generality claim.

If SVD does not outperform NeuralMF on the frozen primary ChEMBL endpoint, we will NOT alter the potency threshold, mappings, Top-K cutoffs, or evidence definition to manufacture a positive replication.

If mapping coverage is insufficient, the result is labeled inconclusive and coverage itself is reported.

## Claim boundary if primary replication succeeds

Allowed:

> In two distinct biomedical relation settings, benchmark selection changed the model chosen for discovery, and the structure-neutralized-selected model showed greater later/independent evidence yield on a fixed candidate universe.

Still not allowed:

- all biomedical AI benchmarks are invalid;
- degree matching is uniquely correct;
- every structure-neutralized winner generalizes better;
- ChEMBL-supported pairs are necessarily clinically relevant interactions;
- ChEMBL-unlabeled pairs are negative interactions.
