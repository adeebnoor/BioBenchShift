# Matching robustness amendment — 13 September 2026

This amendment refines Expansion D in `SCIENCE_EXPANSION_PROTOCOL_20260913.md` before any matching-robustness outcome has been generated or inspected.

The three neutralization variants are fixed as follows:

1. **Original joint log2 degree-bin matching**: exactly the frozen project rule (`0` for degree 0; otherwise `floor(log2(d))+1`) applied separately to the two endpoints, with joint-bin matching.
2. **Joint degree-decile matching**: compute `log1p` endpoint degree from the training graph only. Learn decile cut points separately for left and right endpoint distributions from the training nodes; collapse duplicate cut points deterministically; match positives to nonedges in the same pair of resulting endpoint bins.
3. **Nearest-neighbour degree matching**: standardize the two-dimensional vector `(log1p left degree, log1p right degree)` using means and standard deviations from training-graph endpoint degrees. For each positive, select an unused nonedge minimizing Euclidean distance subject to a fixed **0.25-standard-deviation caliper on each coordinate**. Positives without an eligible match are excluded and matching coverage is reported.

For all methods:

- endpoint degrees are computed from training edges only;
- held-out positives and all known positives are excluded from negative pools;
- one negative is used per matched positive;
- the same train/test splits and model fits are reused across neutralization methods within a seed so only the evaluation control changes;
- model comparison uses SVD, NeuralMF and LightGCN on **10 fixed seeds (0–9)**;
- conventional random-unlabelled evaluation remains identical across the three neutralization variants;
- report per-seed AUROC, neutralized AUROC, matching coverage and model winner; summarize the frequency of winner identity under each regime;
- the original log2-bin method remains the frozen primary method. The two added methods are robustness analyses, not replacements selected after outcome inspection.
