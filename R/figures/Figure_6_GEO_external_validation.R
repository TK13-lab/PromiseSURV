# Purpose: Reproduce manuscript Figure 6 (locked GEO validation).


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

  geo_summary_file <- file.path(input_root, "outputs/task4/geo_external_validation/plan2_geo_marker_score_validation.csv")
  geo_predictions_file <- file.path(input_root, "outputs/task4/geo_external_validation/plan2_geo_marker_score_predictions.csv")
  kirc_summary_file <- file.path(input_root, "outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_metrics.csv")
  kirc_predictions_file <- file.path(input_root, "outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_predictions.csv")
  enrichment_panel_file <- file.path(input_root, "outputs/task6/shared_markers/task6_shared_marker_enrichment_panel_b.csv")
  enrichment_curve_file <- file.path(input_root, "outputs/task6/shared_markers/task6_shared_marker_enrichment_curve_data.csv")
  required_files <- c(
    geo_summary_file, geo_predictions_file, kirc_summary_file,
    kirc_predictions_file, enrichment_panel_file, enrichment_curve_file
  )
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) stop("Missing Figure 6 inputs: ", paste(missing_files, collapse = ", "))

  geo_summary <- read.csv(geo_summary_file, stringsAsFactors = FALSE)
  geo_predictions <- read.csv(geo_predictions_file, stringsAsFactors = FALSE)
  kirc_summary <- read.csv(kirc_summary_file, stringsAsFactors = FALSE)
  kirc_predictions <- read.csv(kirc_predictions_file, stringsAsFactors = FALSE)
  enrichment_panel <- read.csv(enrichment_panel_file, stringsAsFactors = FALSE)
  enrichment_curve <- read.csv(enrichment_curve_file, stringsAsFactors = FALSE)

  fmt2 <- function(x) sprintf("%.2f", as.numeric(x))
  fmt_p <- function(x) {
    x <- as.numeric(x)
    ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
  }

  # Panel A: continuous-score Cox associations for every external cohort.
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
  forest_data <- rbind(geo_rows, kirc_rows)
  forest_data$endpoint_short <- ifelse(forest_data$endpoint == "OS", "OS", "DFS/prog.")
  dataset_labels <- paste(forest_data$cancer_label, forest_data$accession)
  forest_data$dataset_short <- factor(
    dataset_labels,
    levels = rev(dataset_labels[order(forest_data$cancer_label, forest_data$accession)])
  )

  panel_a <- ggplot(
    forest_data,
    aes(x = cox_hr_per_sd, y = dataset_short, color = evidence, shape = endpoint)
  ) +
    geom_vline(xintercept = 1, color = "#777777", linewidth = 0.45, linetype = "dashed") +
    geom_segment(
      aes(x = 1, xend = cox_hr_per_sd, yend = dataset_short),
      color = "#B8B8B8", linewidth = 0.45
    ) +
    geom_point(size = 5.0) +
    geom_label(
      aes(
        label = paste0(endpoint_short, "; HR=", fmt2(cox_hr_per_sd)),
        hjust = ifelse(cox_hr_per_sd > 1.25, 1.08, -0.08)
      ),
      family = "Arial", size = 3.25, color = "#222222",
      fill = scales::alpha("white", 0.94), linewidth = 0,
      label.padding = unit(0.035, "lines"), label.r = unit(0.02, "lines")
    ) +
    scale_x_log10(limits = c(0.45, 4.3), breaks = c(0.5, 1, 2, 3.5), minor_breaks = NULL) +
    scale_color_manual(values = c(
      supports = "#4E9F8E", suggestive = "#C9A227", `opposite/weak` = "#B65C6A"
    )) +
    labs(x = "Locked marker-score HR per SD", y = NULL, tag = "A", color = NULL, shape = NULL) +
    theme_classic(base_family = "Arial", base_size = 10.5) +
    guides(
      color = guide_legend(order = 1, nrow = 1),
      shape = guide_legend(order = 2, nrow = 1)
    ) +
    theme(
      axis.text = element_text(size = 9.6, color = "#222222"),
      axis.title = element_text(size = 10.8, color = "#222222"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 9.4),
      panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.28),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(8, 44, 8, 8)
    )

  # Panel B: four endpoint-labeled Kaplan-Meier stress tests.
  make_survival_panel <- function(data_for_panel, title, endpoint, hr, p_value, x_limit, show_y_axis) {
    data_for_panel <- data_for_panel[
      !is.na(data_for_panel$time_plot) &
        !is.na(data_for_panel$event_plot) &
        data_for_panel$risk_group %in% c("low", "high"),
    ]
    data_for_panel$risk_group <- factor(
      data_for_panel$risk_group,
      levels = c("low", "high"),
      labels = c("Low score", "High score")
    )
    survival_fit <- survfit(Surv(time_plot, event_plot) ~ risk_group, data = data_for_panel)
    x_upper <- if (is.null(x_limit)) max(data_for_panel$time_plot, na.rm = TRUE) else x_limit[[2]]
    label_x <- if (is.null(x_limit)) 0.06 * x_upper else x_limit[[1]] + 0.06 * diff(x_limit)
    survival_plot <- ggsurvplot(
      survival_fit,
      data = data_for_panel,
      conf.int = FALSE,
      risk.table = FALSE,
      pval = FALSE,
      censor = FALSE,
      size = 1.10,
      break.time.by = ifelse(x_upper <= 5.5, 1, 2),
      palette = c("#2F78B7", "#C44E52"),
      legend.title = "",
      legend.labs = c("Low score", "High score"),
      ggtheme = theme_classic(base_family = "Arial", base_size = 10.2)
    )
    plot <- survival_plot$plot +
      annotate(
        "label",
        x = label_x, y = 0.25,
        label = paste0(endpoint, "\nHR=", fmt2(hr), "; p=", fmt_p(p_value)),
        hjust = 0, family = "Arial", size = 3.25, lineheight = 0.94,
        fill = scales::alpha("white", 0.88), color = "#222222"
      ) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0.01, 0.01)) +
      labs(x = "Years", y = if (show_y_axis) "Survival probability" else NULL, title = title) +
      theme(
        plot.title = element_text(family = "Arial", face = "bold", hjust = 0.5, size = 10.8, lineheight = 0.90),
        axis.title = element_text(family = "Arial", size = 10.5),
        axis.text = element_text(family = "Arial", size = 9.5),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5)
      )
    if (!is.null(x_limit)) {
      plot <- plot + coord_cartesian(xlim = x_limit, ylim = c(0, 1), clip = "on")
    }
    plot
  }

  km_plots <- list()
  for (accession in c("GSE31210", "GSE68465")) {
    cohort_predictions <- geo_predictions[geo_predictions$accession == accession, ]
    cohort_predictions$time_plot <- cohort_predictions$OS_time / 365.25
    cohort_predictions$event_plot <- cohort_predictions$OS
    cohort_summary <- geo_summary[geo_summary$accession == accession, ][1, ]
    km_plots[[accession]] <- make_survival_panel(
      cohort_predictions,
      paste(cohort_summary$cancer_label, accession),
      "OS",
      cohort_summary$median_high_vs_low_hr,
      cohort_summary$logrank_p,
      c(0, 5),
      accession == "GSE31210"
    )
  }
  for (analysis_set in c("primary_all", "primary_m0_sensitivity")) {
    cohort_predictions <- kirc_predictions[kirc_predictions$analysis_set == analysis_set, ]
    cohort_predictions$time_plot <- cohort_predictions$dfs_time_months / 12
    cohort_predictions$event_plot <- cohort_predictions$dfs_event
    cohort_summary <- kirc_summary[kirc_summary$analysis_set == analysis_set, ][1, ]
    km_title <- ifelse(
      analysis_set == "primary_all",
      "KIRC GSE22541\nprimary",
      "KIRC GSE22541\nM0-only"
    )
    km_plots[[analysis_set]] <- make_survival_panel(
      cohort_predictions,
      km_title,
      "DFS/progression",
      cohort_summary$cox_hr_per_sd,
      cohort_summary$cox_p,
      NULL,
      analysis_set == "primary_all"
    )
  }
  km_grid <- wrap_plots(km_plots, ncol = 2) +
    plot_annotation(caption = "Blue = low marker score; red = high marker score") &
    theme(plot.caption = element_text(family = "Arial", size = 10.2, hjust = 0.5))
  panel_b <- wrap_elements(full = km_grid) +
    labs(tag = "B") +
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998)
    )

  # Panel C: direction-consistent shifts for all eight locked marker genes.
  marker_genes <- c("NAA38", "RAN", "MRPL42", "CDC45", "MAF", "FMO1", "ZBTB4", "STAT5B")
  positive_genes <- c("NAA38", "RAN", "MRPL42", "CDC45")
  zscore <- function(x) {
    numeric_x <- as.numeric(x)
    standard_deviation <- sd(numeric_x, na.rm = TRUE)
    if (is.na(standard_deviation) || standard_deviation < 1e-8) {
      return(rep(NA_real_, length(numeric_x)))
    }
    (numeric_x - mean(numeric_x, na.rm = TRUE)) / standard_deviation
  }
  gene_effect_one <- function(data_for_cohort, dataset_label) {
    data_for_cohort <- data_for_cohort[data_for_cohort$risk_group %in% c("low", "high"), ]
    do.call(rbind, lapply(marker_genes, function(gene) {
      if (!gene %in% names(data_for_cohort)) {
        return(NULL)
      }
      standardized_expression <- zscore(data_for_cohort[[gene]])
      raw_shift <- mean(
        standardized_expression[data_for_cohort$risk_group == "high"],
        na.rm = TRUE
      ) - mean(
        standardized_expression[data_for_cohort$risk_group == "low"],
        na.rm = TRUE
      )
      data.frame(
        dataset = dataset_label,
        gene_symbol = gene,
        direction_effect = ifelse(gene %in% positive_genes, raw_shift, -raw_shift),
        stringsAsFactors = FALSE
      )
    }))
  }
  gene_effects <- do.call(rbind, lapply(
    c("GSE31210", "GSE68465", "GSE41613", "GSE65858", "GSE29609"),
    function(accession) {
      cohort_predictions <- geo_predictions[geo_predictions$accession == accession, ]
      gene_effect_one(
        cohort_predictions,
        paste(unique(cohort_predictions$cancer_label)[1], accession)
      )
    }
  ))
  gene_effects <- rbind(
    gene_effects,
    gene_effect_one(
      kirc_predictions[kirc_predictions$analysis_set == "primary_all", ],
      "KIRC GSE22541 primary"
    ),
    gene_effect_one(
      kirc_predictions[kirc_predictions$analysis_set == "primary_m0_sensitivity", ],
      "KIRC GSE22541 M0-only"
    )
  )
  gene_effects$gene_symbol <- factor(gene_effects$gene_symbol, levels = marker_genes)
  gene_effects$dataset <- factor(gene_effects$dataset, levels = rev(unique(gene_effects$dataset)))
  gene_effects$abs_effect <- pmin(abs(gene_effects$direction_effect), 1.6)

  panel_c <- ggplot(
    gene_effects,
    aes(x = gene_symbol, y = dataset, fill = direction_effect, size = abs_effect)
  ) +
    geom_point(shape = 21, color = "white", stroke = 0.65, alpha = 0.95, na.rm = TRUE) +
    scale_fill_gradient2(
      low = "#4E79A7", mid = "#F5F2EA", high = "#B65C6A", midpoint = 0,
      limits = c(-1.6, 1.6), oob = scales::squish,
      name = "Direction-\nconsistent\nshift"
    ) +
    scale_size_continuous(range = c(1.2, 6.1), breaks = c(0.5, 1.0, 1.5), name = "|shift|") +
    labs(x = "PROMISE marker genes", y = NULL, tag = "C") +
    theme_classic(base_family = "Arial", base_size = 10.2) +
    guides(
      fill = guide_colorbar(
        order = 1, title.position = "top", title.hjust = 0.5,
        barheight = unit(0.72, "inches"), barwidth = unit(0.10, "inches")
      ),
      size = guide_legend(
        order = 2, title.position = "top", title.hjust = 0.5,
        override.aes = list(
          shape = 21, fill = "#C98C96", color = "#555555",
          stroke = 0.35, alpha = 1
        )
      )
    ) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 9.7),
      axis.text.y = element_text(size = 9.5),
      axis.title = element_text(size = 10.5),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "#222222", fill = NA, linewidth = 0.6),
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 7.6),
      legend.title = element_text(size = 8.0, lineheight = 0.86),
      legend.key.height = unit(0.20, "inches"),
      legend.key.width = unit(0.15, "inches"),
      legend.spacing.y = unit(0.03, "inches"),
      legend.box.spacing = unit(0.03, "inches"),
      legend.margin = margin(1, 1, 1, 1),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(7, 3, 7, 6)
    )

  # Panel D: four top GSEA running-score curves from the locked CPTAC RNA ranking.
  wrap_label <- function(x, width = 24) {
    vapply(x, function(label) paste(strwrap(label, width = width), collapse = "\n"), character(1))
  }
  fmt_q <- function(x) {
    x <- as.numeric(x)
    ifelse(is.na(x), "=NA", ifelse(x < 0.001, "<0.001", paste0("=", sprintf("%.3f", x))))
  }
  enrichment_panel <- enrichment_panel[order(enrichment_panel$panel_rank), , drop = FALSE]
  enrichment_panel$term_label_wrapped <- wrap_label(enrichment_panel$term_label, width = 24)
  enrichment_panel$term_label_wrapped <- factor(
    enrichment_panel$term_label_wrapped,
    levels = rev(enrichment_panel$term_label_wrapped)
  )
  enrichment_panel$fdr_label <- paste0("q", fmt_q(enrichment_panel$padj))
  enrichment_colors <- c(
    Reactome = "#4E79A7",
    `GO biological process` = "#59A14F",
    Hallmark = "#C45A6A"
  )
  top_terms <- enrichment_panel[enrichment_panel$panel_rank <= 4, , drop = FALSE]
  curve_top <- enrichment_curve[
    enrichment_curve$panel_rank %in% top_terms$panel_rank, ,
    drop = FALSE
  ]
  curve_top <- merge(
    curve_top,
    top_terms[, c("pathway", "term_label_wrapped", "fdr_label", "neglog10_fdr")],
    by = "pathway",
    all.x = TRUE
  )
  curve_top$term_label_wrapped <- factor(
    curve_top$term_label_wrapped,
    levels = rev(levels(top_terms$term_label_wrapped))
  )
  curve_hits <- curve_top[curve_top$hit, , drop = FALSE]
  top_gene_labels <- do.call(rbind, lapply(split(curve_hits, curve_hits$pathway), function(data_for_pathway) {
    data_for_pathway <- data_for_pathway[
      order(data_for_pathway$rank, -data_for_pathway$rank_score, data_for_pathway$gene_symbol), ,
      drop = FALSE
    ]
    data_for_pathway <- data_for_pathway[seq_len(min(3, nrow(data_for_pathway))), , drop = FALSE]
    data_for_pathway$top_gene_rank <- seq_len(nrow(data_for_pathway))
    data_for_pathway
  }))
  top_gene_labels$term_label_wrapped <- factor(
    top_gene_labels$term_label_wrapped,
    levels = levels(curve_top$term_label_wrapped)
  )
  top_gene_labels$label_x <- ifelse(
    top_gene_labels$rank < 35,
    top_gene_labels$rank + 14 + (top_gene_labels$top_gene_rank - 1) * 28,
    top_gene_labels$rank + 5
  )
  top_gene_labels$label_x <- pmin(
    top_gene_labels$label_x,
    max(curve_top$rank, na.rm = TRUE) - 28
  )
  top_gene_labels$label_y <- c(0.78, 0.58, 0.38)[top_gene_labels$top_gene_rank]

  curve_labels <- unique(curve_top[, c(
    "pathway", "term_label_wrapped", "collection_label", "NES", "fdr_label"
  )])
  curve_y <- aggregate(running_es ~ pathway, data = curve_top, max)
  curve_labels <- merge(curve_labels, curve_y, by = "pathway", all.x = TRUE)
  curve_labels$x <- max(curve_top$rank, na.rm = TRUE) * 0.96
  curve_labels$y <- curve_labels$running_es * 0.88
  curve_labels$label <- paste0(
    "NES ", fmt2(curve_labels$NES), "\n", curve_labels$fdr_label
  )

  rank_score_scale <- max(abs(curve_top$rank_score), na.rm = TRUE)
  rank_track <- unique(curve_top[, c(
    "pathway", "term_label_wrapped", "rank", "rank_score"
  )])
  rank_track$rank_base <- -0.68
  rank_track$rank_y <- rank_track$rank_base + 0.25 * rank_track$rank_score / rank_score_scale
  rank_track$term_label_wrapped <- factor(
    rank_track$term_label_wrapped,
    levels = levels(curve_top$term_label_wrapped)
  )
  track_labels <- unique(curve_top[, c("pathway", "term_label_wrapped")])
  track_labels$x <- 6
  track_labels$hit_y <- -0.18
  track_labels$rank_y <- -0.68
  track_labels$term_label_wrapped <- factor(
    track_labels$term_label_wrapped,
    levels = levels(curve_top$term_label_wrapped)
  )

  panel_d <- ggplot() +
    geom_hline(yintercept = 0, color = "#777777", linewidth = 0.25, linetype = "22") +
    geom_tile(
      data = rank_track,
      aes(x = rank, y = -0.34, fill = rank_score),
      width = 1.02, height = 0.10, inherit.aes = FALSE
    ) +
    geom_segment(
      data = curve_hits,
      aes(x = rank, xend = rank, y = -0.29, yend = -0.08),
      inherit.aes = FALSE, color = "#222222", linewidth = 0.16, alpha = 0.58
    ) +
    geom_segment(
      data = rank_track,
      aes(x = rank, xend = rank, y = rank_base, yend = rank_y),
      inherit.aes = FALSE, color = "#B8B8B8", linewidth = 0.22, alpha = 0.72
    ) +
    geom_hline(yintercept = -0.68, color = "#555555", linewidth = 0.22, linetype = "22") +
    geom_line(
      data = curve_top,
      aes(x = rank, y = running_es, color = collection_label, group = pathway),
      linewidth = 0.76, lineend = "round"
    ) +
    geom_segment(
      data = top_gene_labels,
      aes(x = rank, y = running_es, xend = label_x, yend = label_y),
      inherit.aes = FALSE, color = "#222222", linewidth = 0.22,
      arrow = arrow(length = unit(0.035, "inches"), type = "open")
    ) +
    geom_text(
      data = top_gene_labels,
      aes(x = label_x, y = label_y, label = gene_symbol),
      inherit.aes = FALSE, family = "Arial", color = "#111111",
      size = 3.20, fontface = "italic", hjust = 0
    ) +
    geom_text(
      data = curve_labels,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE, family = "Arial", size = 3.20,
      lineheight = 0.90, hjust = 1, color = "#222222"
    ) +
    geom_text(
      data = track_labels,
      aes(x = x, y = hit_y, label = "hits"),
      inherit.aes = FALSE, family = "Arial", size = 3.20,
      hjust = 0, color = "#444444"
    ) +
    geom_text(
      data = track_labels,
      aes(x = x, y = rank_y, label = "ranked score"),
      inherit.aes = FALSE, family = "Arial", size = 3.20,
      hjust = 0, color = "#444444"
    ) +
    facet_wrap(~term_label_wrapped, ncol = 2, scales = "free_y") +
    scale_color_manual(values = enrichment_colors, name = NULL) +
    scale_fill_gradient2(
      low = "#3D6FB6", mid = "#F5F0E1", high = "#C44E52",
      midpoint = 0, limits = c(-rank_score_scale, rank_score_scale), guide = "none"
    ) +
    scale_x_continuous(
      breaks = c(1, 180, 362), labels = c("1", "180", "362"),
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      limits = c(-0.98, 1.02), breaks = c(0, 0.5, 1.0),
      labels = c("0", "0.5", "1.0"), expand = c(0.01, 0.01)
    ) +
    labs(x = "Rank in CPTAC-derived gene list", y = "Running ES\nRanked list", tag = "D") +
    theme_classic(base_family = "Arial", base_size = 9.3) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        family = "Arial", face = "bold", size = 9.0,
        color = "#222222", lineheight = 0.86
      ),
      axis.text = element_text(size = 9.0, color = "#222222"),
      axis.title = element_text(size = 9.8, color = "#222222"),
      legend.position = "none",
      panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.22),
      panel.grid.minor.x = element_blank(),
      panel.spacing = unit(0.08, "inches"),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(6, 8, 5, 6)
    )

  figure_6 <- (panel_a + panel_b + plot_layout(widths = c(1.05, 1.30))) /
    (panel_c + panel_d + plot_layout(widths = c(0.95, 1.35))) +
    plot_layout(heights = c(1.05, 1.25)) &
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998)
    )

  png_file <- file.path(output_dir, "Figure_6_geo_external_validation_large_text_400dpi.png")
  pdf_file <- file.path(output_dir, "Figure_6_geo_external_validation_large_text.pdf")
  svg_file <- file.path(output_dir, "Figure_6_geo_external_validation_large_text.svg")
  ggsave(png_file, figure_6, width = 9, height = 9.6, units = "in", dpi = 400, device = ragg::agg_png)
  ggsave(pdf_file, figure_6, width = 9, height = 9.6, units = "in", device = grDevices::cairo_pdf)
  ggsave(svg_file, figure_6, width = 9, height = 9.6, units = "in", device = svglite::svglite)
  message("Wrote Figure 6 to: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
}

main()
