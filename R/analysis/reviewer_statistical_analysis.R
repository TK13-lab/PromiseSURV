main <- function() {
  # Purpose: Generate reviewer-requested cohort, uncertainty, baseline-comparison,
  # and inferential-statistics tables from locked PROMISE-Surv artifacts.
  # Usage:
  #   Rscript R/analysis/reviewer_statistical_analysis.R --input-root=PATH --tcga-cdr-file=FILE

  args_file <- commandArgs(FALSE)
  file_arg <- args_file[grepl("^--file=", args_file)]
  script_path <- if (length(file_arg) > 0) {
    sub("^--file=", "", file_arg[[1]])
  } else {
    "R/analysis/reviewer_statistical_analysis.R"
  }
  package_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  trailing_args <- commandArgs(trailingOnly = TRUE)
  argument_value <- function(prefix, default) {
    hit <- trailing_args[startsWith(trailing_args, prefix)]
    if (length(hit) == 0) default else sub(prefix, "", hit[[1]], fixed = TRUE)
  }
  project_dir <- normalizePath(
    argument_value("--input-root=", package_dir),
    winslash = "/", mustWork = TRUE
  )
  out_dir <- argument_value(
    "--output-dir=", file.path(package_dir, "outputs", "reviewer_statistics")
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  suppressPackageStartupMessages({
    library(survival)
    library(readxl)
  })

  set.seed(20260903)

  fmt_n_pct <- function(n, denominator) {
    if (is.na(n) || denominator == 0) {
      return("NA")
    }
    sprintf("%d (%.1f%%)", as.integer(n), 100 * n / denominator)
  }

  fmt_median_iqr <- function(x, digits = 1) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) {
      return("NA")
    }
    q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    sprintf(paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]"), q[2], q[1], q[3])
  }

  harrell_c <- function(time, status, risk) {
    ok <- is.finite(time) & is.finite(status) & is.finite(risk)
    if (sum(ok) < 5 || sum(status[ok] == 1) < 2) {
      return(NA_real_)
    }
    as.numeric(concordance(Surv(time[ok], status[ok]) ~ risk[ok], reverse = TRUE)$concordance)
  }

  bootstrap_stat <- function(n, statistic, iterations = 2000L) {
    vals <- replicate(iterations, statistic(sample.int(n, n, replace = TRUE)))
    vals <- vals[is.finite(vals)]
    c(mean = mean(vals), lower95 = unname(quantile(vals, 0.025)), upper95 = unname(quantile(vals, 0.975)), used = length(vals))
  }

  # Cohort characteristics: exact complete-case discovery IDs and deployed TCGA IDs.
  model_input <- read.csv(file.path(project_dir, "outputs/task3/model_inputs/cptac_mvp_complete_case_input.csv"), check.names = FALSE)
  cptac_clin <- read.csv(file.path(project_dir, "outputs/task2/processed_modalities/clinical/clinical_features_cptac.csv"), check.names = FALSE)
  cptac <- merge(
    unique(model_input[, c("patient_id", "cancer_label", "OS", "OS_time")]),
    unique(cptac_clin[, c("patient_id", "cancer_label", "age_raw", "gender_clean", "stage_clean")]),
    by = c("patient_id", "cancer_label"), all.x = TRUE, suffixes = c("", ".clinical")
  )

  tcga_input <- read.csv(file.path(project_dir, "outputs/task5/task5_tcga_clinical_utility_input.csv"), check.names = FALSE)
  tcga_input <- tcga_input[tcga_input$model == "promise_rna_clinical_student", ]
  tcga_cdr_path <- argument_value(
    "--tcga-cdr-file=",
    file.path(project_dir, "data_links", "raw", "TCGA-CDR", "TCGA-CDR-SupplementalTableS1.xlsx")
  )
  if (!file.exists(tcga_cdr_path)) {
    stop("TCGA-CDR workbook not found; pass --tcga-cdr-file=PATH: ", tcga_cdr_path)
  }
  tcga_cdr <- suppressWarnings(as.data.frame(read_excel(tcga_cdr_path, sheet = "TCGA-CDR")))
  tcga <- merge(
    tcga_input,
    tcga_cdr[, c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage")],
    by.x = "patient_id", by.y = "bcr_patient_barcode", all.x = TRUE
  )

  summarize_cohort <- function(df, source_label, raw_age, male, early, late, missing_stage) {
    do.call(rbind, lapply(c("LUAD", "KIRC", "HNSC"), function(cancer) {
      x <- df[df$cancer_label == cancer, ]
      n <- nrow(x)
      data.frame(
        Cancer = cancer,
        Cohort = source_label,
        N = n,
        OS_events = sum(x$OS == 1, na.rm = TRUE),
        Median_follow_up_years = sprintf("%.2f", median(x$OS_time, na.rm = TRUE) / 365.25),
        Age_years_median_IQR = fmt_median_iqr(raw_age(x), 1),
        Male_n_percent = fmt_n_pct(sum(male(x), na.rm = TRUE), n),
        Stage_I_II_n_percent = fmt_n_pct(sum(early(x), na.rm = TRUE), n),
        Stage_III_IV_n_percent = fmt_n_pct(sum(late(x), na.rm = TRUE), n),
        Stage_missing_n_percent = fmt_n_pct(sum(missing_stage(x), na.rm = TRUE), n),
        Modalities = if (source_label == "CPTAC/PDC discovery") "Clinical; RNA; proteomics; phosphoproteomics" else "Clinical; RNA",
        stringsAsFactors = FALSE
      )
    }))
  }

  cptac_table <- summarize_cohort(
    cptac, "CPTAC/PDC discovery",
    function(x) x$age_raw,
    function(x) toupper(x$gender_clean) == "MALE",
    function(x) grepl("STAGE_I$|STAGE_II$", x$stage_clean),
    function(x) grepl("STAGE_III$|STAGE_IV$", x$stage_clean),
    function(x) is.na(x$stage_clean) | x$stage_clean == "MISSING"
  )

  tcga_stage <- toupper(as.character(tcga$ajcc_pathologic_tumor_stage))
  tcga_table <- summarize_cohort(
    tcga, "TCGA deployment",
    function(x) x$age_at_initial_pathologic_diagnosis,
    function(x) toupper(x$gender) == "MALE",
    function(x) grepl("STAGE I($|[ABC])|STAGE II($|[ABC])", toupper(x$ajcc_pathologic_tumor_stage)),
    function(x) grepl("STAGE III($|[ABC])|STAGE IV($|[ABC])", toupper(x$ajcc_pathologic_tumor_stage)),
    function(x) is.na(x$ajcc_pathologic_tumor_stage) | !grepl("STAGE (I|II|III|IV)", toupper(x$ajcc_pathologic_tumor_stage))
  )
  cohort_table <- rbind(cptac_table, tcga_table)
  write.csv(cohort_table, file.path(out_dir, "table1_cohort_characteristics.csv"), row.names = FALSE)

  # Internal repeated-split summaries and paired tests against the blended student.
  strength <- read.csv(file.path(project_dir, "outputs/task3/student_strengthening/task3_student_strengthening_metrics.csv"), check.names = FALSE)
  strength <- strength[strength$split_label == "repeated", ]
  model_labels <- c(
    cv_blended_student = "Blended RNA+clinical student",
    clinical_cox = "Clinical Cox scaffold",
    compact_multiomics_teacher = "Compact multi-omics teacher",
    distilled_student = "Direct distilled RNA+clinical student"
  )
  models <- names(model_labels)
  baseline <- strength[strength$model == "cv_blended_student", c("repeat", "test_c_index")]
  names(baseline)[2] <- "blended_c"

  internal_rows <- lapply(models, function(model_name) {
    x <- strength[strength$model == model_name, c("repeat", "test_c_index")]
    x <- x[order(x[["repeat"]]), ]
    boot <- bootstrap_stat(nrow(x), function(ii) mean(x$test_c_index[ii]), 5000L)
    if (model_name == "cv_blended_student") {
      delta <- rep(0, nrow(x))
      p <- NA_real_
    } else {
      paired <- merge(baseline, x, by = "repeat")
      delta <- paired$blended_c - paired$test_c_index
      p <- wilcox.test(paired$blended_c, paired$test_c_index, paired = TRUE, exact = FALSE)$p.value
    }
    delta_boot <- if (model_name == "cv_blended_student") {
      c(mean = 0, lower95 = 0, upper95 = 0, used = nrow(x))
    } else {
      bootstrap_stat(length(delta), function(ii) mean(delta[ii]), 5000L)
    }
    data.frame(
      Model = unname(model_labels[model_name]),
      Repeated_splits = nrow(x),
      Mean_C_index = mean(x$test_c_index),
      SD_C_index = sd(x$test_c_index),
      Mean_C_index_lower95 = boot["lower95"],
      Mean_C_index_upper95 = boot["upper95"],
      Blended_minus_comparator = ifelse(model_name == "cv_blended_student", NA, mean(delta)),
      Difference_lower95 = ifelse(model_name == "cv_blended_student", NA, delta_boot["lower95"]),
      Difference_upper95 = ifelse(model_name == "cv_blended_student", NA, delta_boot["upper95"]),
      Paired_Wilcoxon_P = p,
      stringsAsFactors = FALSE
    )
  })
  internal_table <- do.call(rbind, internal_rows)
  write.csv(internal_table, file.path(out_dir, "table2_internal_performance_uncertainty.csv"), row.names = FALSE)

  # External TCGA bootstrap confidence intervals for C-index, plus AUC and Brier summaries.
  pred <- read.csv(file.path(project_dir, "outputs/task4/external_tcga/task4_tcga_model_predictions.csv"))
  metrics <- read.csv(file.path(project_dir, "outputs/task4/external_tcga/task4_tcga_model_comparison_metrics.csv"))
  external_rows <- list()
  counter <- 1L
  for (cancer in c("HNSC", "KIRC", "LUAD")) {
    for (model_name in c("promise_rna_clinical_student", "rna_clinical_cox", "common_clinical_cox")) {
      x <- pred[pred$cancer_label == cancer & pred$model == model_name, ]
      boot <- bootstrap_stat(nrow(x), function(ii) harrell_c(x$OS_time[ii], x$OS[ii], x$risk_score[ii]), 2000L)
      m <- metrics[metrics$cancer_label == cancer & metrics$model == model_name, ][1, ]
      external_rows[[counter]] <- data.frame(
        Cancer = cancer,
        Model = c(
          promise_rna_clinical_student = "PROMISE RNA+clinical student",
          rna_clinical_cox = "RNA+clinical Cox",
          common_clinical_cox = "Clinical Cox scaffold"
        )[[model_name]],
        N = nrow(x),
        OS_events = sum(x$OS == 1),
        C_index = harrell_c(x$OS_time, x$OS, x$risk_score),
        C_index_lower95 = boot["lower95"],
        C_index_upper95 = boot["upper95"],
        AUC_1_year = m$auc_365d,
        AUC_3_year = m$auc_1095d,
        AUC_5_year = m$auc_1825d,
        Integrated_Brier_mean = m$integrated_brier_mean,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  external_table <- do.call(rbind, external_rows)
  write.csv(external_table, file.path(out_dir, "table3_external_performance_with_uncertainty.csv"), row.names = FALSE)

  # Continuous PROMISE risk associations, unadjusted and age/stage-adjusted.
  cox_results <- read.csv(file.path(project_dir, "outputs/task5/task5_tcga_multivariable_cox.csv"))
  cox_results <- cox_results[cox_results$model %in% c("risk_only", "risk_plus_age_stage"), ]
  cox_results$Analysis <- ifelse(cox_results$model == "risk_only", "Unadjusted PROMISE risk (per SD)", "PROMISE risk adjusted for age and stage (per SD)")
  cox_table <- cox_results[, c("cancer_label", "Analysis", "n", "events", "primary_hr", "primary_hr_lower95", "primary_hr_upper95", "primary_p", "concordance", "ph_global_min_p")]
  names(cox_table) <- c("Cancer", "Analysis", "N", "OS_events", "HR", "HR_lower95", "HR_upper95", "P_value", "Model_C_index", "PH_global_min_P")
  write.csv(cox_table, file.path(out_dir, "table4_external_cox_associations.csv"), row.names = FALSE)

  # Machine-readable audit summary.
  audit <- data.frame(
    Item = c("Analysis seed", "Internal mean bootstrap iterations", "External C-index bootstrap iterations", "Internal paired test", "External concordance", "TCGA clinical source"),
    Value = c("20260903", "5000", "2000", "two-sided paired Wilcoxon signed-rank", "Harrell C via survival::concordance", tcga_cdr_path),
    stringsAsFactors = FALSE
  )
  write.csv(audit, file.path(out_dir, "revision_statistical_audit.csv"), row.names = FALSE)

  message("Wrote reviewer statistical outputs to: ", out_dir)
}

main()
