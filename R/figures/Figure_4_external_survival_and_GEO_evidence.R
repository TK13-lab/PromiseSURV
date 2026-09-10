# Purpose: Reproduce manuscript Figure 4 (TCGA survival and GEO evidence).


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  argument_value <- function(prefix, default) {
    hit <- args[startsWith(args, prefix)]
    if (length(hit) == 0) default else sub(prefix, "", hit[[1]], fixed = TRUE)
  }
  script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][[1]])
  package_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
  input_root <- normalizePath(argument_value("--input-root=", package_root), winslash = "/", mustWork = TRUE)
  output_dir <- argument_value("--output-dir=", file.path(package_root, "outputs", "figures"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(grid)
    library(survival)
    library(survminer)
  })

  utility_file <- file.path(input_root, "outputs/task5/task5_tcga_clinical_utility_input.csv")
  geo_summary_file <- file.path(input_root, "outputs/task4/geo_external_validation/plan2_geo_marker_score_validation.csv")
  kirc_summary_file <- file.path(input_root, "outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_metrics.csv")
  kirc_predictions_file <- file.path(input_root, "outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_predictions.csv")
  required_files <- c(utility_file, geo_summary_file, kirc_summary_file, kirc_predictions_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) stop("Missing Figure 4 inputs: ", paste(missing_files, collapse = ", "))

  utility <- read.csv(utility_file, stringsAsFactors = FALSE)
  geo_summary <- read.csv(geo_summary_file, stringsAsFactors = FALSE)
  kirc_summary <- read.csv(kirc_summary_file, stringsAsFactors = FALSE)
  kirc_predictions <- read.csv(kirc_predictions_file, stringsAsFactors = FALSE)

  fmt2 <- function(x) sprintf("%.2f", as.numeric(x))
  fmt_p <- function(x) {
    x <- as.numeric(x)
    ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
  }

  # Panel A: cancer-specific TCGA Kaplan-Meier curves.
  utility <- utility[
    utility$cancer_label %in% c("HNSC", "KIRC", "LUAD") &
      complete.cases(utility[, c("OS_time", "OS", "risk_score")]) &
      utility$OS_time > 0,
  ]
  utility$OS_years <- utility$OS_time / 365.25
  utility$cancer_label <- factor(utility$cancer_label, levels = c("HNSC", "KIRC", "LUAD"))

  make_median_result <- function(data_for_one_cancer) {
    cutoff <- median(data_for_one_cancer$risk_score, na.rm = TRUE)
    data_for_one_cancer$risk_group <- ifelse(
      data_for_one_cancer$risk_score >= cutoff,
      "High risk",
      "Low risk"
    )
    data_for_one_cancer$risk_group <- factor(
      data_for_one_cancer$risk_group,
      levels = c("Low risk", "High risk")
    )
    survival_fit <- survfit(Surv(OS_years, OS) ~ risk_group, data = data_for_one_cancer)
    cox_fit <- coxph(Surv(OS_years, OS) ~ risk_group, data = data_for_one_cancer)
    cox_summary <- summary(cox_fit)
    logrank_fit <- survdiff(Surv(OS_years, OS) ~ risk_group, data = data_for_one_cancer)
    list(
      data = data_for_one_cancer,
      fit = survival_fit,
      hr = cox_summary$conf.int[1, "exp(coef)"],
      hr_low = cox_summary$conf.int[1, "lower .95"],
      hr_high = cox_summary$conf.int[1, "upper .95"],
      logrank_p = pchisq(logrank_fit$chisq, length(logrank_fit$n) - 1, lower.tail = FALSE)
    )
  }

  median_results <- lapply(c("HNSC", "KIRC", "LUAD"), function(cancer) {
    make_median_result(utility[utility$cancer_label == cancer, ])
  })
  names(median_results) <- c("HNSC", "KIRC", "LUAD")

  make_km_panel <- function(cancer, result, show_legend, show_y_axis, panel_tag) {
    annotation <- paste0(
      "HR ", fmt2(result$hr), " (", fmt2(result$hr_low), "-", fmt2(result$hr_high), ")",
      "\nlog-rank p=", fmt_p(result$logrank_p)
    )
    survival_plot <- ggsurvplot(
      result$fit,
      data = result$data,
      conf.int = FALSE,
      risk.table = FALSE,
      pval = FALSE,
      censor = FALSE,
      size = 0.86,
      break.time.by = 2,
      palette = c("#2E78B7", "#C44E52"),
      legend.title = "",
      legend.labs = c("Low risk", "High risk"),
      ggtheme = theme_classic(base_family = "Arial", base_size = 11.5)
    )
    survival_plot$plot +
      annotate(
        "label",
        x = 0.10, y = 0.19, label = annotation,
        hjust = 0, vjust = 0, family = "Arial", size = 4.00,
        lineheight = 0.98, fill = scales::alpha("white", 0.82), color = "#222222"
      ) +
      scale_x_continuous(breaks = c(0, 1, 3, 5), expand = c(0.01, 0.01)) +
      scale_y_continuous(breaks = c(0, 0.5, 1.0), limits = c(0, 1), expand = c(0.01, 0.01)) +
      coord_cartesian(xlim = c(0, 5), ylim = c(0, 1), clip = "on") +
      labs(
        title = cancer,
        x = "Years from diagnosis",
        y = if (show_y_axis) "Overall survival" else NULL,
        tag = panel_tag
      ) +
      theme(
        text = element_text(family = "Arial", size = 11.5, color = "#222222"),
        axis.text = element_text(size = 11.5, color = "#222222"),
        axis.title = element_text(size = 12.7, color = "#222222"),
        plot.title = element_text(
          family = "Arial", face = "bold", size = 13.3, hjust = 0.5,
          color = c(HNSC = "#C45A6A", KIRC = "#2F9C8B", LUAD = "#6A73B7")[[cancer]]
        ),
        legend.position = if (show_legend) "bottom" else "none",
        legend.direction = "horizontal",
        legend.text = element_text(size = 10.5),
        legend.key.width = unit(0.22, "inches"),
        plot.tag = element_text(family = "Arial", face = "bold", size = 24, color = "#111111"),
        plot.tag.position = c(0.002, 0.998),
        plot.margin = margin(7, 8, 7, 8)
      )
  }

  panel_a <- wrap_plots(
    list(
      make_km_panel("HNSC", median_results[["HNSC"]], FALSE, TRUE, "A"),
      make_km_panel("KIRC", median_results[["KIRC"]], TRUE, FALSE, NULL),
      make_km_panel("LUAD", median_results[["LUAD"]], FALSE, FALSE, NULL)
    ),
    nrow = 1
  )

  # Panel B: locked external-marker Cox associations, including endpoint labels.
  geo_rows <- geo_summary[, c("accession", "cancer_label", "n", "events", "cox_hr_per_sd", "cox_p")]
  geo_rows$endpoint <- "OS"
  geo_rows$evidence <- ifelse(
    geo_rows$cox_hr_per_sd > 1 & geo_rows$cox_p < 0.05,
    "supports",
    ifelse(geo_rows$cox_hr_per_sd > 1, "suggestive", "opposite/weak")
  )
  kirc_rows_source <- kirc_summary[
    kirc_summary$analysis_set %in% c("primary_all", "primary_m0_sensitivity"),
  ]
  kirc_rows <- data.frame(
    accession = ifelse(
      kirc_rows_source$analysis_set == "primary_all",
      "GSE22541 primary",
      "GSE22541 M0-only"
    ),
    cancer_label = "KIRC",
    n = kirc_rows_source$n,
    events = kirc_rows_source$events,
    cox_hr_per_sd = kirc_rows_source$cox_hr_per_sd,
    cox_p = kirc_rows_source$cox_p,
    endpoint = "DFS/progression",
    evidence = "supports",
    stringsAsFactors = FALSE
  )
  geo_plot <- rbind(geo_rows, kirc_rows)
  geo_plot$facet_label <- ifelse(
    grepl("^GSE22541 ", geo_plot$accession),
    paste0(geo_plot$cancer_label, "\nGSE22541\n", sub("^GSE22541 ", "", geo_plot$accession)),
    paste0(geo_plot$cancer_label, "\n", geo_plot$accession)
  )
  geo_plot$endpoint_short <- ifelse(geo_plot$endpoint == "OS", "OS", "DFS/prog.")
  geo_plot$mini_label <- paste0(
    geo_plot$endpoint_short,
    "\nn=", geo_plot$n, "/", geo_plot$events,
    "\nHR=", fmt2(geo_plot$cox_hr_per_sd),
    "\np=", fmt_p(geo_plot$cox_p)
  )

  panel_b <- ggplot(geo_plot, aes(x = cox_hr_per_sd, y = 1, color = evidence, shape = endpoint)) +
    geom_vline(xintercept = 1, color = "#777777", linewidth = 0.45, linetype = "dashed") +
    geom_point(size = 5.2) +
    geom_label(
      aes(x = 1, y = 0.72, label = mini_label),
      family = "Arial", size = 3.20, lineheight = 0.90, color = "#222222",
      fill = scales::alpha("white", 0.82), linewidth = 0,
      label.padding = unit(0.03, "lines")
    ) +
    facet_wrap(~facet_label, ncol = 4) +
    scale_x_log10(limits = c(0.45, 3.8), breaks = c(0.5, 1, 2.5), minor_breaks = NULL) +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_color_manual(values = c(
      supports = "#4E9F8E", suggestive = "#C9A227", `opposite/weak` = "#B65C6A"
    )) +
    coord_cartesian(ylim = c(0.51, 1.25), clip = "off") +
    labs(x = "Locked marker-score HR per SD", tag = "B", color = NULL, shape = NULL) +
    theme_classic(base_family = "Arial", base_size = 10.5) +
    theme(
      axis.text = element_text(size = 9.7, color = "#222222"),
      axis.title = element_text(size = 11.0, color = "#222222"),
      strip.text = element_text(size = 9.8, face = "bold", lineheight = 0.88),
      legend.position = "bottom",
      legend.text = element_text(size = 10.0),
      panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.28),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(8, 9, 8, 8)
    )

  # Panel C: GSE22541 primary and M0-only score distributions.
  box_data <- kirc_predictions[
    kirc_predictions$analysis_set %in% c("primary_all", "primary_m0_sensitivity") &
      kirc_predictions$dfi_group %in% c("long_dfs", "short_dfs"),
  ]
  box_data$dfi_group <- factor(
    box_data$dfi_group,
    levels = c("long_dfs", "short_dfs"),
    labels = c("Long\nDFS", "Short\nDFS")
  )
  box_data$analysis_set <- factor(
    box_data$analysis_set,
    levels = c("primary_all", "primary_m0_sensitivity"),
    labels = c("Primary ccRCC", "M0-only sensitivity")
  )
  panel_c <- ggplot(box_data, aes(x = dfi_group, y = promise_marker_score, fill = dfi_group)) +
    geom_boxplot(
      width = 0.55, outlier.shape = NA, alpha = 0.25,
      color = "#333333", linewidth = 0.45
    ) +
    geom_point(
      aes(color = dfi_group),
      position = position_jitter(width = 0.08, height = 0, seed = 20260903),
      size = 2.8, alpha = 0.9, show.legend = FALSE
    ) +
    facet_wrap(~analysis_set, nrow = 1) +
    scale_fill_manual(values = c("Long\nDFS" = "#4E79A7", "Short\nDFS" = "#B65C6A")) +
    scale_color_manual(values = c("Long\nDFS" = "#4E79A7", "Short\nDFS" = "#B65C6A")) +
    labs(x = NULL, y = "Locked PROMISE marker score", tag = "C") +
    theme_classic(base_family = "Arial", base_size = 10.8) +
    theme(
      axis.text.x = element_text(size = 9.3, color = "#222222"),
      axis.text.y = element_text(size = 10.5, color = "#222222"),
      axis.title = element_text(size = 10.2, color = "#222222"),
      strip.text = element_text(size = 9.8, face = "bold", lineheight = 0.88),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.28),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(7, 7, 7, 22)
    )

  figure4 <- panel_a / (panel_b + panel_c + plot_layout(widths = c(1.42, 0.78))) +
    plot_layout(heights = c(0.95, 1.05)) &
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998)
    )

  png_file <- file.path(output_dir, "Figure_4_external_survival_and_geo_evidence_large_text_400dpi.png")
  svg_file <- file.path(output_dir, "Figure_4_external_survival_and_geo_evidence_large_text.svg")
  pdf_file <- file.path(output_dir, "Figure_4_external_survival_and_geo_evidence_large_text.pdf")
  ggsave(png_file, figure4,
    width = 9.0, height = 8.4, units = "in", dpi = 400,
    bg = "white", device = ragg::agg_png, limitsize = FALSE
  )
  ggsave(svg_file, figure4,
    width = 9.0, height = 8.4, units = "in",
    bg = "white", device = svglite::svglite, limitsize = FALSE
  )
  ggsave(pdf_file, figure4,
    width = 9.0, height = 8.4, units = "in",
    bg = "white", device = grDevices::cairo_pdf, limitsize = FALSE
  )
  message("Wrote ", png_file)
}

main()
