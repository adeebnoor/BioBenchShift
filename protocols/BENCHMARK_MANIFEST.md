# Gate 1 benchmark manifest

Status as of 13 September 2026. This is a working acquisition/provenance manifest, not a claim that every candidate will enter the final paper.

| ID | Relation family | Candidate source | Scale / rationale | Status | Primary role |
|---|---|---|---|---|---|
| DDI-ANTI | drug–drug | Anti-DDI v3.0.1 + frozen GoldD2 ATC5 positive reference | Existing seed diagnostic; audited evidence states and reproduced ~0.90→~0.50 degree-bias collapse | **READY** | Figure 1 seed observation |
| DTI-BIOSNAP | drug–target | BioSNAP TargetDecagon (STITCH-derived experimentally supported targets) | 18,690 edges; 284 drugs; 3,648 genes | **SOURCE IDENTIFIED**; freeze checksum + license/provenance before analysis | Independent bipartite relation family |
| PPI-HURI | protein–protein | HuRI / Human Reference Interactome | Independent physical-interaction network; avoids DDI/DTI lineage | **ACQUISITION + LICENSE CHECK** | Independent network-biology family |
| DRUGDIS-HET | compound–treats–disease | Hetionet v1.0 CtD / underlying PharmacotherapyDB relation | 755 CtD edges in published Hetionet; mechanistically different endpoint types | **SOURCE IDENTIFIED**; use edge-level provenance/license metadata | Drug–disease family |
| GENEDIS-HET | disease–associates–gene | Hetionet v1.0 DaG | 12,623 edges; useful sensitivity relation but shares Hetionet integration pipeline | **SECONDARY** | Within-KG relation-type sensitivity |
| PPI-HET | gene–interacts–gene | Hetionet v1.0 GiG | 147,164 edges; large PPI relation but not independent of Hetionet | **SECONDARY** | Within-KG relation-type sensitivity |
| OBL-HQ | heterogeneous biomedical KG | OpenBioLink2020 directed high-quality | >8.5M train edges; 28 relations; includes positive and true-negative resources and standardized splits | **SOURCE IDENTIFIED**; large-scale stress test after Gate 1 | Multi-relational stress test |

## Gate 1 admission rules

A benchmark may contribute to the **primary cross-domain claim** only if all are true:

1. exact version/release is frozen;
2. source and endpoint identifiers are documented;
3. positive-edge definition is explicit;
4. negative/unlabelled construction is explicit;
5. license permits the intended analysis and reproducible release strategy;
6. no test information is used to tune the structural matching procedure;
7. relation family is not counted as “independent” merely because it is a different edge type from the same integrated graph.

## Immediate acquisition order

1. BioSNAP DTI — first independent bipartite test.
2. HuRI PPI — first independent homogeneous biological-network test.
3. Hetionet CtD — first drug–disease test, with license/provenance retained at edge level.
4. OpenBioLink HQ — large-scale multi-relational stress test only after the low-cost audits establish effect sizes.

## Current environment note

The canonical BioSNAP file endpoint has been identified, but binary dataset retrieval is not available through the current execution runtime. This is an execution-environment constraint, not a scientific blocker. The analysis code is already generic and ready for the frozen edge list once bytes are available through a mounted/uploaded artifact or a text-accessible source.
