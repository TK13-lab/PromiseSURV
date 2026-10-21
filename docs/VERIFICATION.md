# Verification scope

Checks performed on 21 October 2026:

- All eight R and five Python scripts passed syntax parsing.
- The Figure 1 input builder reproduced the frozen panel CSV byte for byte (SHA-256 `1cff6e3199192b75e6eff44d36461290c0e64b5c32748d63b52a3c88f97a7366`).
- Figures 1 and 2 ran successfully and exported PNG, PDF and SVG files.
- The unchanged Table 2–3 calculation blocks from the R statistical audit were run separately against the frozen inputs, using 5,000 split-level and 2,000 patient-level bootstrap resamples. All four Table 2 rows and nine Table 3 rows matched the stored outputs within a numerical tolerance of 1e-12.
- Data retrieval dates were checked against original acquisition records.

The complete R audit could not be rerun in this checkout because the TCGA-CDR workbook was unavailable. Supply that workbook with `--tcga-cdr-file=PATH`. The full internal Python model refit and complete raw-data pipeline were not rerun for this update. These checks do not establish end-to-end reproduction from a fresh public clone: the frozen input tree remains a separate requirement.

Historical project records document successful rendering of Figures 1–6 and benchmark checks on 3 September 2026; those broader runs were not repeated for this documentation and missing-script update. Existing analysis implementations are unchanged.
