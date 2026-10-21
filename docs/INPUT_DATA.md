# Input data contract

This repository is intentionally code-only. Raw CPTAC, TCGA, and GEO data are not duplicated. Pass the root of the frozen analysis tree with `--input-root=PATH` (R) or `--input-root PATH` (Python). The tree must retain the relative paths below.

The frozen processed inputs are not distributed in this repository. The listed public sources identify the study data, but this package does not provide the complete acquisition, preprocessing and original model-training pipeline needed to recreate every frozen artifact.

## Source data

Retrieval dates below were checked against original timestamped acquisition records. They refer to source retrieval, not file modification times or dataset release dates.

| Source | Accession or resource | Retrieved |
| --- | --- | --- |
| [PDC](https://pdc.cancer.gov/pdc/) | LUAD: PDC000153/PDC000149; KIRC: PDC000127/PDC000128 (proteome/phosphoproteome, clinical records and sample maps) | 23 June 2026 |
| [PDC](https://pdc.cancer.gov/pdc/) | HNSC: PDC000221/PDC000222 | 24 June 2026 |
| [cBioPortal API](https://www.cbioportal.org/api) | Matched CPTAC RNA profiles: `luad_cptac_gdc_mrna_seq_tpm`, `rcc_cptac_gdc_mrna_seq_tpm`, `ohnca_cptac_gdc_mrna_seq_tpm` | 24 June 2026 |
| [GDC PanCancer Atlas](https://gdc.cancer.gov/about-data/publications/pancanatlas) and [TCGA-CDR](https://doi.org/10.1016/j.cell.2018.02.052) | TCGA clinical/CDR; RNA file `3586c0da-64d0-4b74-a449-5ff4d9136611`; RPPA file `fcbb373e-28d4-4818-92f3-601ede3da5e1` | 23 June 2026 |
| [UCSC Xena](https://xenabrowser.net/datapages/) | `TCGA_phenotype_denseDataOnlyDownload.tsv.gz` | 23 June 2026 |
| [cBioPortal API](https://www.cbioportal.org/api) | Subtypes: `luad_tcga_pan_can_atlas_2018`, `kirc_tcga_pan_can_atlas_2018` | 23 June 2026 |
| [cBioPortal API](https://www.cbioportal.org/api) | Subtypes: `hnsc_tcga_pan_can_atlas_2018` | 24 June 2026 |
| [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) | Series GSE31210, GSE68465, GSE41613, GSE65858, GSE29609, GSE22541; platforms GPL96, GPL570, GPL10558, GPL1708 | 2 July 2026 |

These analyses used public processed resources. PDC reuse follows its [CC BY 4.0 and attribution requirements](https://pdc.cancer.gov/pdc/data-use-guidelines); GDC distinguishes [open and controlled access](https://gdc.cancer.gov/access-data/data-access-processes-and-tools); GEO reuse is subject to its [disclaimer and submitter rights](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html). Source-specific conditions also apply to cBioPortal-derived resources.

## Figure scripts

- Figure 1: `outputs/figures/figure2/figure2_cohort_construction_panel_data.csv`
- Figure 2: no data file; all workflow nodes, connectors, coordinates, fonts, and colors are explicit in the script.
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

## Manuscript statistical audit

`R/analysis/reviewer_statistical_analysis.R` reads these locked inputs:

- `outputs/task3/model_inputs/cptac_mvp_complete_case_input.csv`
- `outputs/task2/processed_modalities/clinical/clinical_features_cptac.csv`
- `outputs/task3/student_strengthening/task3_student_strengthening_metrics.csv`
- `outputs/task4/external_tcga/task4_tcga_model_predictions.csv`
- `outputs/task4/external_tcga/task4_tcga_model_comparison_metrics.csv`
- `outputs/task5/task5_tcga_clinical_utility_input.csv`
- `outputs/task5/task5_tcga_multivariable_cox.csv`
- TCGA-CDR Supplemental Table S1, supplied explicitly with `--tcga-cdr-file=PATH`

The canonical Table 2 is reconstructed from the frozen `student_strengthening` metrics for the 314-patient analysis. The optional Python internal refit uses the same current complete-case input with fixed penalties; it does not reproduce the original training-selected blending procedure. Superseded historical Task 3 repeated-split files reflect an earlier 316-patient snapshot and should not be substituted for the canonical metrics.
