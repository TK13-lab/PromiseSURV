# Statistical methods and reproducibility

All inferential p-values produced by the statistical audit and Python benchmark scripts are two-sided.

## Manuscript statistical audit

`R/analysis/reviewer_statistical_analysis.R` assembles Tables 1–3 from frozen analysis artifacts. Its `table4_external_cox_associations.csv` exports stored continuous TCGA Cox associations; it is not manuscript Table 4, which combines TCGA and GEO survival contrasts. The external TCGA and GEO Python scripts expose the corresponding Cox and log-rank calculations.

1. Internal performance and uncertainty
   - Ten paired repeated-split C-index estimates for the blended RNA+clinical student and three comparators.
   - Mean C-index 95% percentile intervals from 5,000 nonparametric resamples of the ten split-level estimates.
   - Blended-minus-comparator mean differences with 95% percentile intervals from 5,000 paired split resamples.
   - Two-sided paired Wilcoxon signed-rank tests for model comparisons.
   - Random seed 20260903.
2. External TCGA performance
   - Cancer-specific Harrell C-index for HNSC, KIRC, and LUAD.
   - C-index 95% percentile intervals from 2,000 patient-level nonparametric bootstrap resamples.
   - Locked 1-, 3-, and 5-year AUC and integrated Brier summaries are joined from the original benchmark artifact.
3. External Cox associations
   - Stored continuous PROMISE risk associations per standard deviation, assessed with two-sided Cox tests; the R audit exports these estimates without refitting.
   - Both unadjusted estimates and estimates adjusted for age and stage are reported with hazard-ratio 95% confidence intervals.

## Transparent Python recomputations

1. Internal CPTAC model benchmark
   - Ten stratified 70/30 train/test splits, seeds 20260623–20260632.
   - Penalized Cox models: penalizer 0.50 and L1 ratio 0.00.
   - RNA+clinical distilled student: ridge alpha 10.0.
   - Fixed-split patient-level C-index uncertainty: 2,000 nonparametric bootstrap samples and 95% percentile intervals.
   - Mean repeated-split C-index and paired mean-difference uncertainty: 5,000 split-level bootstrap samples and 95% percentile intervals.
   - Across-split model comparison: paired two-sided Wilcoxon signed-rank tests with Holm correction.
2. External TCGA benchmark
   - Bootstrap base seed 20260623 plus the zero-based cancer index in the script's declared order.
   - Cancer-specific and pooled Harrell C-index.
   - Two-sided univariable Cox test per one-standard-deviation increase in score, with 95% confidence interval.
   - Median-split Kaplan–Meier comparison using a two-sided log-rank test.
   - Paired patient-level C-index differences between PROMISE and the clinical scaffold: 2,000 bootstrap samples, 95% percentile intervals, two-sided bootstrap p-values, and Holm correction across cancers.
3. GEO validation
   - Two-sided univariable Cox and log-rank tests for each locked cohort/analysis set.
   - Two-sided Welch tests for high-versus-low standardized marker expression, followed by Benjamini–Hochberg adjustment within each dataset.
4. Figure 3 biology
   - Two-sided Spearman rank correlations for teacher-to-student transfer and RNA/RPPA program associations.
   - Benjamini–Hochberg adjustment within the declared test families.

## Original analysis and reproducibility limits

The training-selected blended model underlying Table 2 used outer seeds 20260623–20260632 and inner seeds equal to the outer seed plus 1000 plus the zero-based inner-repeat index. Its frozen `student_strengthening` metrics are the canonical Table 2 input. The fixed-penalty Python refit is a separate analysis, not a recreation of that tuning and blending workflow. The reported marker negative controls used seed 20260702; their original generation pipeline is not included in this package.

`R/analysis/shared_marker_enrichment.R` uses seed 20260225, serial `fgsea` execution, gene-set sizes 3–200 and `eps = 0`. Enrichment also depends on the MSigDB snapshot supplied through `msigdbr`.

Initial cohort-level missingness filtering, normalization, feature-space construction and imputation preceded model splitting. Subsequent model-stage scaling, low-variance filtering and outcome-guided compact-teacher selection used training partitions. The complete pipeline should therefore not be described as nested preprocessing from raw data.

Reproduction requires identical ordered inputs and software/gene-set versions. Seeds control the specified random operations; solver convergence, BLAS libraries, hardware and parallel numerical execution can still cause small differences. Deterministic table assembly and plotting do not require random seeds.

