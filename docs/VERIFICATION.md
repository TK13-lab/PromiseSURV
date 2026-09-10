# Verification record

Verified on 2026-09-03 against the frozen analysis tree used for the revision.

- All eight R files parsed successfully and contain only two top-level expressions: a local `main` function and its invocation.
- No R figure script contains `source()` or relies on a shared plotting configuration.
- Figures 1–6 all ran successfully and generated 400-dpi PNG, PDF, and SVG outputs.
- PNG dimensions were 4600 x 3360 (Figure 1), 4480 x 2480 (Figure 2), 3600 x 3200 (Figure 3), 3600 x 3360 (Figure 4), 3600 x 4640 (Figure 5), and 3600 x 3840 (Figure 6).
- Visual inspection covered all six figures. The final checks specifically confirmed enlarged text in Figures 3–6, an unclipped Figure 5 Panel F header, protected Figure 6 Panel A labels, the original Figure 6 Panel C height, reduced Panel C bubble sizes, and a visible external size legend.
- Every Python file passed byte-code compilation.
- An AST audit found no module-level Python configuration assignments or executable setup calls.
- The Figure 1 Python builder generated a CSV byte-for-byte identical to the frozen panel-data CSV.
- Runtime tests completed for the external TCGA benchmark, GEO survival/marker tests, Figure 3 biology tests, and the internal benchmark `--quick-check` path.
- The canonical R reviewer audit ran successfully with 5,000 internal split-level bootstrap resamples and 2,000 external patient-level bootstrap resamples, reproducing the revised manuscript Table 2 values from the locked student-strengthening metrics.
- The full Python internal benchmark remains the default (10 splits, 5,000 split-level bootstrap samples, and 2,000 patient-level bootstrap samples); it is computationally intensive and therefore has an explicitly labeled quick-check mode for installation verification.

The optional Python internal refit uses the current 314-patient complete-case input and is not claimed to be byte-identical to a superseded 316-patient historical Task 3 snapshot. The canonical manuscript Table 2 route is `R/analysis/reviewer_statistical_analysis.R` with the locked student-strengthening metrics.

Expected benign messages are limited to scale replacement from `survminer`/`ggsurvplot` and a ggpubr line-size deprecation warning. They do not alter the exported panels.
