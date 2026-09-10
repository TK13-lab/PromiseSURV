# Input data contract

This repository is intentionally code-only. Raw CPTAC, TCGA, and GEO data are not duplicated. 

## Figure scripts

- Figure 3:
  - `outputs/task6/task6_biology_audit.csv`
  - `outputs/task6/cptac/task6_cptac_teacher_student_risk_scores.csv`
  - `outputs/task6/tcga/task6_tcga_rna_pathway_risk_associations.csv`
  - `outputs/task6/tcga/task6_tcga_rna_pathway_scores.csv`
  - `outputs/task6b/rppa/task6b_rppa_module_risk_associations.csv`
- Figure 4:
  - `outputs/task5/task5_tcga_clinical_utility_input.csv`
  - `outputs/task4/geo_external_validation/plan2_geo_marker_score_validation.csv`
  - `outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_metrics.csv`
  - `outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_predictions.csv`
- Figure 5: the following files under `outputs/task6/shared_markers/`:
  - `task6_shared_marker_core.csv`
  - `task6_shared_marker_model_contrast.csv`
  - `task6_shared_marker_summary.csv`
  - `task6_shared_marker_overlap_membership.csv`
  - `task6_shared_marker_selection_steps.csv`
  - `task6_shared_marker_enrichment_panel_b.csv`
  - `task6_shared_marker_enrichment_ranked_gene_scores.csv`
  - `task6_shared_marker_enrichment_gsea.csv`
  - `task6_shared_marker_enrichment_ora.csv`
  - `task6_shared_marker_enrichment_curve_data.csv`
  - `task6_shared_program_summary.csv`
  - `task6_shared_program_evidence.csv`
  - `task6_shared_marker_application_summary.csv`
  - `task6_shared_marker_selection_config.json`
- Figure 6:
  - `outputs/task4/geo_external_validation/plan2_geo_marker_score_validation.csv`
  - `outputs/task4/geo_external_validation/plan2_geo_marker_score_predictions.csv`
  - `outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_metrics.csv`
  - `outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_predictions.csv`
  - `outputs/task6/shared_markers/task6_shared_marker_enrichment_panel_b.csv`
  - `outputs/task6/shared_markers/task6_shared_marker_enrichment_curve_data.csv`

## Python scripts

- `Figure_1_build_panel_data.py` reads the locked Task 1–4 cohort, missingness, complete-case, and TCGA deployment tables named directly in that file.
- `01_internal_survival_benchmarks.py` reads `outputs/task3/model_inputs/cptac_mvp_complete_case_input.parquet`.
- `02_external_TCGA_benchmarks.py` reads the two locked prediction tables under `outputs/task4/external_tcga/` named directly in that file.
- `03_GEO_survival_and_marker_tests.py` reads the two locked GEO prediction tables named directly in that file.
- `04_Figure_3_biology_tests.py` reads the CPTAC risk-transfer, TCGA RNA pathway-score, and TCGA RPPA risk-input tables named directly in that file.

## Rebuilding an intermediate

`R/analysis/shared_marker_enrichment.R` regenerates the GSEA/ORA and curve tables for Figures 5 and 6 from `outputs/task6/shared_markers/task6_shared_marker_overlap_membership.csv`. MSigDB retrieval performed by `msigdbr` requires an available package database/cache.

## Reviewer statistical audit

`R/analysis/reviewer_statistical_analysis.R` reads these locked inputs:

- `outputs/task3/model_inputs/cptac_mvp_complete_case_input.csv`
- `outputs/task2/processed_modalities/clinical/clinical_features_cptac.csv`
- `outputs/task3/student_strengthening/task3_student_strengthening_metrics.csv`
- `outputs/task4/external_tcga/task4_tcga_model_predictions.csv`
- `outputs/task4/external_tcga/task4_tcga_model_comparison_metrics.csv`
- `outputs/task5/task5_tcga_clinical_utility_input.csv`
- `outputs/task5/task5_tcga_multivariable_cox.csv`
- TCGA-CDR Supplemental Table S1, supplied explicitly with `--tcga-cdr-file=PATH`

The canonical revised-manuscript Table 2 is reconstructed from the locked `student_strengthening` metrics. The optional Python internal refit instead uses the current 314-patient complete-case matrix. It is provided as a transparent model recomputation and is not claimed to be byte-identical to the superseded historical Task 3 repeated-split files, which reflect an earlier 316-patient snapshot.
