# Gate 2 pilot — learned models on BioSNAP DTI

| Model | Random-negative AUC | Degree-matched AUC | Mean drop |
|---|---:|---:|---:|
| Truncated-SVD latent factors | 0.950 ± 0.005 | 0.914 ± 0.007 | 0.036 |
| Logistic neural matrix factorization | 0.997 ± 0.001 | 0.908 ± 0.008 | 0.089 |

Mean degree-matching coverage: **0.599**. These are structural latent models, a Gate-2 pilot rather than the final SOTA model panel.

These values were reproduced on multiple GitHub Actions runs from `analysis/run_biosnap_dti_gate2_models.py`; this summary is frozen separately because concurrent result-writing workflows caused the original result commit to race even though the analysis step itself completed successfully.
