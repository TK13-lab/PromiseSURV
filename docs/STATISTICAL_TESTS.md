# Statistical tests exposed for review

All inferential p-values produced by the reviewer-facing scripts are two-sided.

## Canonical revised-manuscript audit

`R/analysis/reviewer_statistical_analysis.R` reconstructs the revised cohort table and Tables 2–4 from locked analysis artifacts.

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
   - Continuous PROMISE risk per standard deviation, assessed with two-sided Cox tests.
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

The exact values, test families, grouping variables, and random seeds are stated locally inside each file rather than imported from a global configuration file.
