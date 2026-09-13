# H4c — uncertainty for non-circular BioGRID future yield

Paired bootstrap over **5,635** later-added relations; **20,000** resamples. All rankings are generated from the frozen historical network before later-edge labels are used.

| K | SVD hits | NeuralMF hits | Recall difference | Paired bootstrap 95% CI | P(bootstrap difference <= 0) |
|---:|---:|---:|---:|---:|---:|
| 100 | 7 | 0 | 0.0012 | [0.0004, 0.0023] | 0.0008 |
| 500 | 18 | 1 | 0.0030 | [0.0016, 0.0046] | 5e-05 |
| 1,000 | 26 | 5 | 0.0037 | [0.0020, 0.0057] | 5e-05 |
| 5,000 | 87 | 22 | 0.0115 | [0.0082, 0.0151] | 5e-05 |
| 10,000 | 159 | 42 | 0.0208 | [0.0163, 0.0252] | 5e-05 |
| 50,000 | 522 | 181 | 0.0605 | [0.0531, 0.0680] | 5e-05 |

Hypergeometric enrichment probabilities against a uniformly random top-K set are stored in the model table. The paired bootstrap is the primary uncertainty summary for the prespecified SVD-versus-NeuralMF future-evidence contrast.
