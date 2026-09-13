# Frozen contemporary-model challenge: GraphBAN on BioSNAP

Status: preregistered before inspecting any GraphBAN structure-matched result.

## Purpose

Challenge the structural-benchmark conclusion with a contemporary compound–protein interaction model published with BioSNAP experiments. This is a robustness test, not a new leaderboard.

## External model

GraphBAN repository: `HamidHadipour/GraphBAN`
Pinned upstream commit: `2dded31a46edfb18d685fb1bb361a95e99eaae49`
Primary setting: the authors' published BioSNAP **transductive** splits and pretrained checkpoints.

The upstream implementation supplies BioSNAP train/validation/test CSVs with `SMILES`, `Protein`, `Y`, plus pretrained transductive checkpoints for seeds 12, 14, 16, 18 and 20.

## Primary audit

For each upstream seed:

1. Use the upstream train/test split without relabeling or resplitting.
2. Run the corresponding published pretrained GraphBAN checkpoint on its official test set.
3. Verify that exported prediction targets exactly match test `Y` row order.
4. Define training endpoint degree from **positive (`Y=1`) training relations only**.
5. Assign each drug and protein to the same log2 degree-bin rule used elsewhere in this project: bin 0 for unseen/zero degree; otherwise `floor(log2(degree))+1`.
6. Conventional AUROC is computed on the complete official test set.
7. Structure-matched AUROC is computed by matching positive and negative test examples within the identical `(drug-degree-bin, protein-degree-bin)` stratum.
8. Because some strata can contain unequal positive/negative counts, repeat deterministic within-stratum downsampling for **100 matching replicates** and report mean, SD and range.
9. Report matching coverage as matched positives / all test positives.
10. Compute a degree-only score on the same conventional and matched examples as a manipulation check. A valid neutralization should substantially reduce degree-only discrimination.

## Interpretation locked before outcome inspection

- A GraphBAN conventional-to-matched AUROC change of >=0.03 with adequate matching coverage is treated as material evaluation sensitivity.
- A smaller change is a boundary result and will be reported, not hidden.
- This challenge does **not** require GraphBAN to lose to SVD or any other model. The scientific question is whether a strong contemporary CPI architecture is itself sensitive to benchmark structure.
- No GraphBAN hyperparameter tuning is allowed for this audit because published pretrained checkpoints are used.
- The five upstream seeds are the complete planned replication set; no seed will be discarded based on outcome.

## Portability-only execution patches

The GitHub Actions runner is CPU-only. The upstream prediction script may be patched only to:

- load the checkpoint with `map_location=cpu`;
- release ESM memory before loading ChemBERTa;
- release ChemBERTa memory after feature extraction.

These patches change device/memory handling only, not model weights, architecture, data, scores or thresholds.

## Reporting

The primary manuscript text will report all five seeds if execution is successful. A one-seed smoke test is operational only and cannot by itself satisfy the contemporary-model send gate.
