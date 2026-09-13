Dear Editors,

We submit the Research Article, **“Benchmark design redirects biomedical discovery,”** for consideration in *Science*. Biomedical AI increasingly determines which molecular relations are prioritized for experimental follow-up, yet model selection is commonly treated as though the benchmark were a passive measuring instrument. We test a broader possibility: that benchmark construction itself can change which model is selected and therefore which biology is tested first.

Across drug–target, protein–protein, compound–disease and disease–gene relations, structural observability alone produced strong apparent discrimination under conventional sampled-unknown evaluation and largely disappeared after structural matching. Learned models depended on this signal unequally. In drug–target and historical protein-interaction analyses, the benchmark change reversed the selected model; stability-controlled ensembles then disagreed on 100% and 99% of top-100 hypotheses, respectively, far beyond within-model variability. We next froze a historical BioGRID candidate universe before opening later evidence: the structure-neutralized-selected model recovered 522 of 5,635 later-added interactions within its top 50,000 predictions compared with 181 for the conventional winner, with paired bootstrap differences positive at every prespecified cutoff. In a separately curated ChEMBL 37 evidence source that was not used for model fitting, selection or threshold tuning, the same evaluation principle concentrated more supported drug–target relations at the earliest experimental cutoffs, while a broad-cutoff crossover was retained as an explicit boundary.

The advance is not a new degree-bias correction. Degree/rich-node bias, prior bias and benchmark leakage have important precedents. The contribution is the experimentally linked scientific-decision chain: **benchmark construction → model identity → hypothesis identity → later or independent evidence**. The work spans multiple biomedical relation families, includes non-reversal and instability controls, preserves negative/boundary results, and separates the central claim from a related Anti-DDI resource and a separate allocation-identity manuscript. Code, frozen protocols, seed-level outputs, source data and provenance are maintained in a public reproducibility branch and will be archived to a persistent release before submission.

We believe the manuscript is appropriate for *Science* because it concerns how AI evaluation choices can redirect the allocation of experimental attention across biomedical discovery, rather than the performance of a single algorithm or benchmark. The manuscript is not under consideration elsewhere. Related work and resource lineage will be disclosed transparently in the submission system and cover materials.

Sincerely,

**Adeeb Noor**  
Department of Information Technology  
Faculty of Computing and Information Technology  
King Abdulaziz University  
Jeddah, Saudi Arabia  
arnoor@kau.edu.sa
