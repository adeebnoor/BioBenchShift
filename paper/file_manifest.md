# Science submission file manifest

## Source files

- `paper/manuscript.md` — current main-text source.
- `paper/supplementary_materials.md` — current Supplementary source.
- `paper/cover_letter.md` — current cover-letter source.
- `submission/SCIENCE_CTS_METADATA.md` — CTS entry sheet.
- `submission/RELATED_WORK_DISCLOSURE.md` — related-work firewall/disclosure wording.
- `figures/source_data_main.csv` — source data for the four main figures.
- `analysis/build_science_figures.py` — reproducible main-figure builder.

## Generated files

The workflow `.github/workflows/build_science_submission_package.yml` generates:

- `figures/generated/Figure1_structural_observability.png`
- `figures/generated/Figure1_structural_observability.pdf`
- `figures/generated/Figure2_hypothesis_identity.png`
- `figures/generated/Figure2_hypothesis_identity.pdf`
- `figures/generated/Figure3_biogrid_future_evidence.png`
- `figures/generated/Figure3_biogrid_future_evidence.pdf`
- `figures/generated/Figure4_chembl_frontier.png`
- `figures/generated/Figure4_chembl_frontier.pdf`
- `paper/generated/Science_Research_Article_PRE_SUBMISSION.docx`
- `paper/generated/Science_Supplementary_Materials_PRE_SUBMISSION.docx`
- `paper/generated/Science_Cover_Letter.docx`

## Final-freeze actions still required

1. Insert the frozen leakage-free GraphBAN result without changing its predeclared criterion.
2. Run the Science editorial red-team screen once more.
3. Verify all dataset/database citations and Science house reference formatting.
4. Replace `PRE_SUBMISSION` labels only after the send gate passes.
5. Freeze the exact Git commit/tag used for submission.
6. Archive that release under a persistent DOI and update data/code availability wording.
7. Verify current CTS technical instructions on the day of upload.
8. Upload only files generated from the final frozen source commit.
