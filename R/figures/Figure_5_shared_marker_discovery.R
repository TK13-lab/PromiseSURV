# Purpose: Reproduce manuscript Figure 5 (shared marker discovery and biology).
# Inputs: Task 6 shared marker discovery outputs and Task 6/6b biology tables.
# Usage: Rscript R/figures/Figure_5_shared_marker_discovery.R [--input-root=PATH] [--output-dir=PATH]

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  argument_value <- function(prefix, default) {
    hit <- args[startsWith(args, prefix)]
    if (length(hit) == 0) default else sub(prefix, "", hit[[1]], fixed = TRUE)
  }
  script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][[1]])
  package_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
  input_root <- normalizePath(argument_value("--input-root=", package_root), winslash = "/", mustWork = TRUE)

  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(grid)
    library(ComplexHeatmap)
  })

  # ComplexHeatmap captures Panel A on an off-screen PostScript device.
  # Register Arial-compatible metrics locally so its labels use the requested family.
  arial_type1_metrics <- grDevices::pdfFonts()[["ArialMT"]]
  arial_type1_metrics$family <- "Arial"
  grDevices::pdfFonts(Arial = arial_type1_metrics)

  figure_id <- "Figure_5_shared_marker_discovery_large_text"
  script_rel <- "R/figures/Figure_5_shared_marker_discovery.R"
  out_dir <- argument_value("--output-dir=", file.path(package_root, "outputs", "figures"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  shared_dir <- file.path(input_root, "outputs/task6/shared_markers")
  core_file <- file.path(shared_dir, "task6_shared_marker_core.csv")
  contrast_file <- file.path(shared_dir, "task6_shared_marker_model_contrast.csv")
  summary_file <- file.path(shared_dir, "task6_shared_marker_summary.csv")
  overlap_file <- file.path(shared_dir, "task6_shared_marker_overlap_membership.csv")
  selection_steps_file <- file.path(shared_dir, "task6_shared_marker_selection_steps.csv")
  enrichment_panel_file <- file.path(shared_dir, "task6_shared_marker_enrichment_panel_b.csv")
  enrichment_ranked_file <- file.path(shared_dir, "task6_shared_marker_enrichment_ranked_gene_scores.csv")
  enrichment_gsea_file <- file.path(shared_dir, "task6_shared_marker_enrichment_gsea.csv")
  enrichment_ora_file <- file.path(shared_dir, "task6_shared_marker_enrichment_ora.csv")
  enrichment_curve_file <- file.path(shared_dir, "task6_shared_marker_enrichment_curve_data.csv")
  program_summary_file <- file.path(shared_dir, "task6_shared_program_summary.csv")
  program_evidence_file <- file.path(shared_dir, "task6_shared_program_evidence.csv")
  application_file <- file.path(shared_dir, "task6_shared_marker_application_summary.csv")
  selection_config_file <- file.path(shared_dir, "task6_shared_marker_selection_config.json")

  input_files <- c(
    core_file,
    contrast_file,
    summary_file,
    overlap_file,
    selection_steps_file,
    enrichment_panel_file,
    enrichment_ranked_file,
    enrichment_gsea_file,
    enrichment_ora_file,
    enrichment_curve_file,
    program_summary_file,
    program_evidence_file,
    application_file,
    selection_config_file
  )
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0) {
    stop("Missing Figure 6 input files: ", paste(missing_inputs, collapse = ", "))
  }

  core <- read.csv(core_file, stringsAsFactors = FALSE)
  contrast <- read.csv(contrast_file, stringsAsFactors = FALSE)
  marker_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
  overlap <- read.csv(overlap_file, stringsAsFactors = FALSE)
  selection_steps <- read.csv(selection_steps_file, stringsAsFactors = FALSE)
  enrichment_panel <- read.csv(enrichment_panel_file, stringsAsFactors = FALSE)
  enrichment_curve <- read.csv(enrichment_curve_file, stringsAsFactors = FALSE)
  program_summary <- read.csv(program_summary_file, stringsAsFactors = FALSE)
  program_evidence <- read.csv(program_evidence_file, stringsAsFactors = FALSE)
  application <- read.csv(application_file, stringsAsFactors = FALSE)

  as_bool <- function(x) {
    if (is.logical(x)) {
      return(x)
    }
    tolower(as.character(x)) %in% c("true", "1", "yes")
  }
  for (bool_col in c(
    "HNSC_significant",
    "KIRC_significant",
    "LUAD_significant",
    "all_three_significant",
    "same_direction_all_three",
    "final_shared_core_marker",
    "promise_enriched_all_cancers"
  )) {
    overlap[[bool_col]] <- as_bool(overlap[[bool_col]])
  }

  cancer_levels <- c("HNSC", "KIRC", "LUAD")
  cancer_cols <- c(HNSC = "#C45A6A", KIRC = "#2F9C8B", LUAD = "#6A73B7")
  model_cols <- c("Clinical scaffold" = "#5C5C5C", "PROMISE student" = "#C44E52")
  layer_cols <- c(RNA = "#3F88C5", RPPA = "#D8902F")
  gain_cols <- c(`TRUE` = "#2A9D8F", `FALSE` = "#8E8E8E")
  selection_cols <- c("#4E79A7", "#59A14F", "#F28E2B", "#E15759", "#8A60A8")

  fmt2 <- function(x) sprintf("%.2f", x)
  fmt3 <- function(x) sprintf("%.3f", x)
  fmt_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  fmt_q <- function(x) {
    ifelse(x < 0.001, "<0.001", paste0("=", sprintf("%.3f", x)))
  }
  wrap_label <- function(x, width = 27) {
    vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", FUN.VALUE = character(1))
  }
  metric_value <- function(metric) {
    value <- application$value[application$metric == metric]
    if (length(value) == 0) {
      return(NA_real_)
    }
    as.numeric(value[[1]])
  }

  tagged_grid_panel <- function(expr, tag, width = 6, height = 4) {
    panel_grob <- grid.grabExpr(expr, width = width, height = height)
    tag_grob <- textGrob(
      tag,
      x = unit(0.01, "npc"),
      y = unit(0.99, "npc"),
      hjust = 0,
      vjust = 1,
      gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 24, col = "#111111")
    )
    wrap_elements(full = grobTree(panel_grob, tag_grob))
  }

  # Panel A: UpSet-style source of the final shared markers.
  significant_sets <- list(
    HNSC = overlap$gene_symbol[overlap$HNSC_significant],
    KIRC = overlap$gene_symbol[overlap$KIRC_significant],
    LUAD = overlap$gene_symbol[overlap$LUAD_significant]
  )
  panel_a_upset <- overlap[overlap$significant_count > 0, c(
    "gene_symbol", "HNSC_significant", "KIRC_significant", "LUAD_significant",
    "significant_count", "intersection_label", "all_three_significant",
    "same_direction_all_three", "final_shared_core_marker"
  )]
  comb_mat <- make_comb_mat(significant_sets, mode = "distinct")
  panel_a <- tagged_grid_panel(
    {
      ht <- UpSet(
        comb_mat,
        bg_col = "#F4F4F4",
        comb_col = "#2B2B2B",
        pt_size = unit(3.4, "mm"),
        lwd = 1.2,
        set_order = cancer_levels,
        top_annotation = upset_top_annotation(
          comb_mat,
          add_numbers = TRUE,
          gp = gpar(fill = "#B7B7B7", col = "white"),
          height = unit(34, "mm"),
          numbers_gp = gpar(fontfamily = "Arial", fontsize = 11.5, col = "#222222")
        ),
        right_annotation = upset_right_annotation(
          comb_mat,
          add_numbers = TRUE,
          gp = gpar(fill = unname(cancer_cols[cancer_levels]), col = "white"),
          width = unit(30, "mm"),
          numbers_gp = gpar(fontfamily = "Arial", fontsize = 11.5, col = "#222222")
        ),
        row_names_gp = gpar(fontfamily = "Arial", fontsize = 11.5)
      )
      draw(
        ht,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        padding = unit(c(5, 4, 5, 4), "mm")
      )
    },
    "A",
    width = 6.4,
    height = 4.7
  )

  # Panel B: MSigDB GSEA of the ranked CPTAC-derived RNA gene list.
  enrichment_panel <- enrichment_panel[order(enrichment_panel$panel_rank), , drop = FALSE]
  enrichment_panel$term_label_wrapped <- wrap_label(enrichment_panel$term_label, width = 24)
  enrichment_panel$term_label_wrapped <- factor(
    enrichment_panel$term_label_wrapped,
    levels = rev(enrichment_panel$term_label_wrapped)
  )
  enrichment_panel$fdr_label <- paste0("q", fmt_q(enrichment_panel$padj))
  enrichment_cols <- c(
    "Reactome" = "#4E79A7",
    "GO biological process" = "#59A14F",
    "Hallmark" = "#C45A6A"
  )
  gsea_plot_terms <- enrichment_panel[enrichment_panel$panel_rank <= 4, , drop = FALSE]
  curve_top <- enrichment_curve[enrichment_curve$panel_rank %in% gsea_plot_terms$panel_rank, , drop = FALSE]
  curve_top <- merge(
    curve_top,
    gsea_plot_terms[, c("pathway", "term_label_wrapped", "fdr_label", "neglog10_fdr")],
    by = "pathway",
    all.x = TRUE
  )
  curve_top$term_label_wrapped <- factor(curve_top$term_label_wrapped, levels = rev(levels(gsea_plot_terms$term_label_wrapped)))
  curve_hits <- curve_top[curve_top$hit, , drop = FALSE]
  # GseaVis convention: gseaNb(addGene=TRUE, markTopgene=TRUE, topGeneN=3)
  # marks the first ranked hit genes for positive-NES terms.
  top_gene_labels <- do.call(
    rbind,
    lapply(split(curve_hits, curve_hits$pathway), function(dat) {
      dat <- dat[order(dat$rank, -dat$rank_score, dat$gene_symbol), , drop = FALSE]
      dat <- dat[seq_len(min(3, nrow(dat))), , drop = FALSE]
      dat$top_gene_rank <- seq_len(nrow(dat))
      dat
    })
  )
  top_gene_labels$term_label_wrapped <- factor(top_gene_labels$term_label_wrapped, levels = levels(curve_top$term_label_wrapped))
  top_gene_labels$label_x <- ifelse(
    top_gene_labels$rank < 35,
    top_gene_labels$rank + 14 + (top_gene_labels$top_gene_rank - 1) * 28,
    top_gene_labels$rank + 5
  )
  top_gene_labels$label_x <- pmin(top_gene_labels$label_x, max(curve_top$rank, na.rm = TRUE) - 28)
  top_gene_labels$label_y <- c(0.78, 0.58, 0.38)[top_gene_labels$top_gene_rank]
  curve_labels <- unique(curve_top[, c("pathway", "term_label_wrapped", "collection_label", "NES", "fdr_label")])
  curve_y <- aggregate(running_es ~ pathway, data = curve_top, max)
  curve_labels <- merge(curve_labels, curve_y, by = "pathway", all.x = TRUE)
  curve_labels$x <- max(curve_top$rank, na.rm = TRUE) * 0.96
  curve_labels$y <- curve_labels$running_es * 0.88
  curve_labels$label <- paste0("NES ", fmt2(curve_labels$NES), "; ", curve_labels$fdr_label)
  rank_score_scale <- max(abs(curve_top$rank_score), na.rm = TRUE)
  rank_track <- unique(curve_top[, c("pathway", "term_label_wrapped", "rank", "rank_score")])
  rank_track$rank_base <- -0.68
  rank_track$rank_y <- rank_track$rank_base + 0.25 * rank_track$rank_score / rank_score_scale
  rank_track$term_label_wrapped <- factor(rank_track$term_label_wrapped, levels = levels(curve_top$term_label_wrapped))
  track_labels <- unique(curve_top[, c("pathway", "term_label_wrapped")])
  track_labels$x <- 6
  track_labels$hit_y <- -0.18
  track_labels$rank_y <- -0.68
  track_labels$term_label_wrapped <- factor(track_labels$term_label_wrapped, levels = levels(curve_top$term_label_wrapped))

  panel_b <- ggplot() +
    geom_hline(yintercept = 0, color = "#777777", linewidth = 0.25, linetype = "22") +
    geom_tile(
      data = rank_track,
      aes(x = rank, y = -0.34, fill = rank_score),
      width = 1.02,
      height = 0.10,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = curve_hits,
      aes(x = rank, xend = rank, y = -0.29, yend = -0.08),
      inherit.aes = FALSE,
      color = "#222222",
      linewidth = 0.16,
      alpha = 0.58
    ) +
    geom_segment(
      data = rank_track,
      aes(x = rank, xend = rank, y = rank_base, yend = rank_y),
      inherit.aes = FALSE,
      color = "#B8B8B8",
      linewidth = 0.22,
      alpha = 0.72
    ) +
    geom_hline(yintercept = -0.68, color = "#555555", linewidth = 0.22, linetype = "22") +
    geom_line(
      data = curve_top,
      aes(x = rank, y = running_es, color = collection_label, group = pathway),
      linewidth = 0.76,
      lineend = "round"
    ) +
    geom_segment(
      data = top_gene_labels,
      aes(x = rank, y = running_es, xend = label_x, yend = label_y),
      inherit.aes = FALSE,
      color = "#222222",
      linewidth = 0.22,
      arrow = arrow(length = unit(0.035, "inches"), type = "open")
    ) +
    geom_text(
      data = top_gene_labels,
      aes(x = label_x, y = label_y, label = gene_symbol),
      inherit.aes = FALSE,
      family = "Arial",
      color = "#111111",
      size = 3.30,
      fontface = "italic",
      hjust = 0
    ) +
    geom_text(
      data = curve_labels,
      aes(x = x, y = y, label = label, color = collection_label),
      inherit.aes = FALSE,
      family = "Arial",
      size = 3.30,
      hjust = 1,
      color = "#222222"
    ) +
    geom_text(
      data = track_labels,
      aes(x = x, y = hit_y, label = "hits"),
      inherit.aes = FALSE,
      family = "Arial",
      size = 3.30,
      hjust = 0,
      color = "#444444"
    ) +
    geom_text(
      data = track_labels,
      aes(x = x, y = rank_y, label = "ranked score"),
      inherit.aes = FALSE,
      family = "Arial",
      size = 3.30,
      hjust = 0,
      color = "#444444"
    ) +
    facet_wrap(~term_label_wrapped, ncol = 2, scales = "free_y") +
    scale_color_manual(values = enrichment_cols, name = NULL) +
    scale_fill_gradient2(
      low = "#3D6FB6",
      mid = "#F5F0E1",
      high = "#C44E52",
      midpoint = 0,
      limits = c(-rank_score_scale, rank_score_scale),
      guide = "none"
    ) +
    scale_x_continuous(breaks = c(1, 180, 362), labels = c("1", "180", "362"), expand = c(0.01, 0.01)) +
    scale_y_continuous(
      limits = c(-0.98, 1.02),
      breaks = c(0, 0.5, 1.0),
      labels = c("0", "0.5", "1.0"),
      expand = c(0.01, 0.01)
    ) +
    labs(x = "Rank in CPTAC-derived gene list", y = "Running ES\nRanked list", tag = "B") +
    theme_classic(base_family = "Arial", base_size = 10.5) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(family = "Arial", face = "bold", size = 9.5, color = "#222222", lineheight = 0.88),
      axis.text = element_text(size = 9.5, color = "#222222"),
      axis.title = element_text(size = 10.5, color = "#222222"),
      axis.line = element_line(color = "#222222", linewidth = 0.42),
      axis.ticks = element_line(color = "#222222", linewidth = 0.32),
      legend.position = "none",
      legend.key.size = unit(0.10, "inches"),
      panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.22),
      panel.grid.minor.x = element_blank(),
      panel.spacing = unit(0.08, "inches"),
      plot.title = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24, color = "#111111"),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(6, 8, 5, 6)
    )

  # Panel C: final same-direction shared marker core.
  marker_summary$gene_symbol <- factor(
    marker_summary$gene_symbol,
    levels = rev(marker_summary$gene_symbol)
  )
  core$gene_symbol <- factor(core$gene_symbol, levels = levels(marker_summary$gene_symbol))
  core$cancer_label <- factor(core$cancer_label, levels = cancer_levels)
  core$cancer_x <- match(core$cancer_label, cancer_levels)
  core$rho_label <- sprintf("%.2f", core$rna_spearman_rho_risk)

  panel_c <- ggplot(core, aes(x = cancer_x, y = gene_symbol, fill = rna_spearman_rho_risk)) +
    geom_tile(color = "white", linewidth = 0.9, width = 0.92, height = 0.84) +
    geom_text(aes(label = rho_label), family = "Arial", size = 3.80, color = "#111111") +
    scale_fill_gradient2(
      low = "#3D6FB6",
      mid = "#F5F0E1",
      high = "#C44E52",
      midpoint = 0,
      limits = c(-0.36, 0.36),
      name = "rho"
    ) +
    scale_x_continuous(
      breaks = seq_along(cancer_levels),
      labels = cancer_levels,
      expand = c(0.02, 0.02)
    ) +
    labs(x = "TCGA cancer", y = "Final shared RNA markers", tag = "C") +
    theme_classic(base_family = "Arial", base_size = 11.2) +
    theme(
      axis.text.x = element_text(face = "bold", size = 11.2),
      axis.text.y = element_text(face = "bold", size = 11.2),
      axis.title = element_text(size = 12.4),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11.2),
      legend.text = element_text(size = 11.2),
      legend.key.height = unit(0.48, "inches"),
      plot.title = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24, color = "#111111"),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(6, 7, 6, 6)
    )

  # Panel D: PROMISE marker gain versus the clinical scaffold.
  gain_summary <- marker_summary
  gain_summary$gene_symbol <- factor(gain_summary$gene_symbol, levels = levels(marker_summary$gene_symbol))
  gain_summary$gain_label <- paste0(
    ifelse(gain_summary$mean_abs_rho_gain_vs_clinical >= 0, "+", ""),
    fmt3(gain_summary$mean_abs_rho_gain_vs_clinical)
  )
  gain_summary$gain_all <- factor(gain_summary$promise_enriched_all_cancers, levels = c(TRUE, FALSE))

  gain_long <- rbind(
    data.frame(
      gene_symbol = gain_summary$gene_symbol,
      model_label = "Clinical scaffold",
      mean_abs_rho = gain_summary$mean_abs_clinical_rho,
      gain_all = gain_summary$gain_all,
      stringsAsFactors = FALSE
    ),
    data.frame(
      gene_symbol = gain_summary$gene_symbol,
      model_label = "PROMISE student",
      mean_abs_rho = gain_summary$mean_abs_student_rho,
      gain_all = gain_summary$gain_all,
      stringsAsFactors = FALSE
    )
  )
  gain_long$model_label <- factor(gain_long$model_label, levels = c("Clinical scaffold", "PROMISE student"))

  panel_d <- ggplot(gain_summary, aes(y = gene_symbol)) +
    geom_segment(
      aes(x = mean_abs_clinical_rho, xend = mean_abs_student_rho, yend = gene_symbol, color = gain_all),
      linewidth = 0.8,
      alpha = 0.78
    ) +
    geom_point(
      data = gain_long,
      aes(x = mean_abs_rho, y = gene_symbol, fill = model_label, shape = model_label),
      inherit.aes = FALSE,
      size = 2.9,
      stroke = 0.75,
      color = "#222222"
    ) +
    geom_text(
      aes(x = pmax(mean_abs_clinical_rho, mean_abs_student_rho) + 0.016, label = gain_label),
      family = "Arial",
      size = 3.40,
      hjust = 0,
      color = "#222222"
    ) +
    scale_color_manual(values = gain_cols, guide = "none") +
    scale_fill_manual(values = model_cols, name = NULL) +
    scale_shape_manual(values = c("Clinical scaffold" = 21, "PROMISE student" = 24), name = NULL) +
    guides(
      fill = guide_legend(override.aes = list(shape = c(21, 24), size = 2.6)),
      shape = "none"
    ) +
    scale_x_continuous(limits = c(0, 0.37), breaks = c(0, 0.10, 0.20, 0.30)) +
    labs(x = "Mean absolute rho across cancers", y = NULL, tag = "D") +
    theme_classic(base_family = "Arial", base_size = 10.5) +
    theme(
      axis.text = element_text(size = 9.8),
      axis.text.y = element_text(face = "bold", size = 9.8),
      axis.title = element_text(size = 10.5),
      axis.line = element_line(color = "#222222", linewidth = 0.42),
      axis.ticks = element_line(color = "#222222", linewidth = 0.32),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 9.5),
      legend.background = element_rect(fill = "white", color = NA),
      legend.margin = margin(1, 2, 1, 2),
      plot.title = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24, color = "#111111"),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(6, 7, 6, 6)
    )

  # Panel E: reusable RNA/RPPA program recurrence as a clean Sankey-style flow.
  program_plot <- program_summary[program_summary$total_supported_contexts > 0, , drop = FALSE]
  program_order <- unique(program_plot$program[order(-program_plot$total_supported_contexts, -program_plot$total_mean_abs_rho, program_plot$program)])
  program_plot$program <- factor(program_plot$program, levels = program_order)
  program_plot$layer <- factor(program_plot$layer, levels = c("RNA", "RPPA"))
  program_plot$support_label <- ifelse(
    program_plot$supported_cancers > 0,
    paste0(program_plot$supported_cancers, "/3"),
    ""
  )

  program_nodes <- unique(program_plot[, c("program", "total_supported_contexts", "total_mean_abs_rho")])
  program_nodes <- program_nodes[match(program_order, as.character(program_nodes$program)), ]
  program_nodes$program_y <- rev(seq(1, by = 1.22, length.out = nrow(program_nodes)))
  program_nodes$program_label <- paste0(program_nodes$program, " (", program_nodes$total_supported_contexts, ")")
  layer_nodes <- data.frame(
    layer = factor(c("RNA", "RPPA"), levels = c("RNA", "RPPA")),
    layer_y = c(max(program_nodes$program_y) * 0.72, max(program_nodes$program_y) * 0.28),
    stringsAsFactors = FALSE
  )
  program_flow <- merge(
    program_plot[program_plot$supported_cancers > 0, , drop = FALSE],
    program_nodes[, c("program", "program_y", "program_label")],
    by = "program",
    all.x = TRUE
  )
  program_flow <- merge(program_flow, layer_nodes, by = "layer", all.x = TRUE)
  panel_e <- ggplot() +
    annotate(
      "text",
      x = 0.10,
      y = max(program_nodes$program_y) + 0.82,
      label = "Layer",
      family = "Arial",
      fontface = "plain",
      size = 3.70,
      hjust = 0,
      color = "#222222"
    ) +
    annotate(
      "text",
      x = 0.82,
      y = max(program_nodes$program_y) + 0.82,
      label = "Program\n(supported contexts)",
      family = "Arial",
      fontface = "plain",
      size = 3.70,
      lineheight = 0.90,
      hjust = 0,
      color = "#222222"
    ) +
    geom_curve(
      data = program_flow,
      aes(
        x = 0.14,
        y = layer_y,
        xend = 0.82,
        yend = program_y,
        color = layer,
        linewidth = supported_cancers
      ),
      curvature = 0.34,
      alpha = 0.68,
      lineend = "round"
    ) +
    geom_point(
      data = layer_nodes,
      aes(x = 0.14, y = layer_y, fill = layer),
      shape = 21,
      size = 6.0,
      stroke = 0.55,
      color = "#222222"
    ) +
    geom_text(
      data = layer_nodes,
      aes(x = 0.055, y = layer_y, label = layer, color = layer),
      family = "Arial",
      fontface = "plain",
      size = 3.70,
      hjust = 0,
      show.legend = FALSE
    ) +
    geom_point(
      data = program_nodes,
      aes(x = 0.82, y = program_y),
      shape = 21,
      size = 3.2,
      stroke = 0.40,
      fill = "#F4F4F4",
      color = "#444444"
    ) +
    geom_text(
      data = program_nodes,
      aes(x = 0.86, y = program_y, label = program_label),
      family = "Arial",
      fontface = "plain",
      size = 3.70,
      hjust = 0,
      color = "#222222"
    ) +
    scale_color_manual(values = layer_cols, guide = "none") +
    scale_fill_manual(values = layer_cols, guide = "none") +
    scale_linewidth_continuous(range = c(1.1, 5.0), guide = "none") +
    coord_cartesian(
      xlim = c(0.03, 1.34),
      ylim = c(0.30, max(program_nodes$program_y) + 1.05),
      clip = "off"
    ) +
    labs(tag = "F") +
    theme_void(base_family = "Arial") +
    theme(
      plot.margin = margin(6, 24, 6, 7),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24, color = "#111111"),
      plot.tag.position = c(0.002, 0.998)
    )

  # Panel F: leading-edge genes supporting the top GSEA term.
  top_term <- enrichment_panel[order(enrichment_panel$panel_rank), ][1, , drop = FALSE]
  leading_genes <- strsplit(top_term$leading_edge_genes[[1]], ";")[[1]]
  leading_edge <- do.call(
    rbind,
    lapply(leading_genes, function(gene) {
      row <- overlap[overlap$gene_symbol == gene, , drop = FALSE]
      if (nrow(row) == 0) {
        return(NULL)
      }
      data.frame(
        gene_symbol = gene,
        cancer_label = cancer_levels,
        rho = as.numeric(unlist(row[1, paste0(cancer_levels, "_rho")], use.names = FALSE)),
        fdr = as.numeric(unlist(row[1, paste0(cancer_levels, "_fdr")], use.names = FALSE)),
        stringsAsFactors = FALSE
      )
    })
  )
  leading_edge$gene_symbol <- factor(leading_edge$gene_symbol, levels = rev(leading_genes))
  leading_edge$cancer_label <- factor(leading_edge$cancer_label, levels = cancer_levels)
  leading_edge$rho_label <- paste0(fmt2(leading_edge$rho), ifelse(leading_edge$fdr < 0.05, "*", ""))

  panel_f <- ggplot(leading_edge, aes(x = cancer_label, y = gene_symbol, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.78, width = 0.90, height = 0.82) +
    geom_text(aes(label = rho_label), family = "Arial", size = 3.70, color = "#111111") +
    scale_fill_gradient2(
      low = "#3D6FB6",
      mid = "#F5F0E1",
      high = "#C44E52",
      midpoint = 0,
      limits = c(-0.50, 0.50),
      name = "rho"
    ) +
    labs(
      x = "TCGA cancer",
      y = paste0(top_term$term_label[[1]], "\nleading edge"),
      tag = "E"
    ) +
    theme_classic(base_family = "Arial", base_size = 10.6) +
    theme(
      axis.text.x = element_text(face = "bold", size = 9.8),
      axis.text.y = element_text(face = "bold", size = 9.2),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10.8, margin = margin(t = 3)),
      axis.title.y = element_text(size = 10.8, margin = margin(r = 3), lineheight = 0.90),
      legend.position = "right",
      legend.key.height = unit(0.33, "inches"),
      legend.title = element_text(size = 10.0),
      legend.text = element_text(size = 9.0),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24, color = "#111111"),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(6, 7, 6, 5)
    )

  row1 <- wrap_plots(panel_a, panel_b, ncol = 2, widths = c(1.00, 1.45))
  row2 <- wrap_plots(panel_c, panel_d, ncol = 2, widths = c(0.95, 1.05))
  row3 <- wrap_plots(panel_f, panel_e, ncol = 2, widths = c(0.92, 1.08))
  main_figure <- wrap_plots(row1, row2, row3, ncol = 1, heights = c(1.28, 0.95, 1.12)) &
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.tag = element_text(family = "Arial", face = "bold", size = 24),
      plot.tag.position = c(0.002, 0.998)
    )

  bind_rows_fill <- function(tables) {
    all_names <- unique(unlist(lapply(tables, names), use.names = FALSE))
    normalized <- lapply(tables, function(dat) {
      missing <- setdiff(all_names, names(dat))
      for (name in missing) {
        dat[[name]] <- NA
      }
      dat[, all_names, drop = FALSE]
    })
    do.call(rbind, normalized)
  }

  png_file <- file.path(out_dir, paste0(figure_id, "_400dpi.png"))
  svg_file <- file.path(out_dir, paste0(figure_id, ".svg"))
  pdf_file <- file.path(out_dir, paste0(figure_id, ".pdf"))
  panel_data_file <- file.path(out_dir, paste0(figure_id, "_panel_data.csv"))
  panel_a_table_file <- file.path(out_dir, paste0(figure_id, "_panel_a_gene_overlap_upset.csv"))
  panel_b_table_file <- file.path(out_dir, paste0(figure_id, "_panel_b_gsea_enrichment.csv"))
  panel_b_gene_label_file <- file.path(out_dir, paste0(figure_id, "_panel_b_gsea_top3_genes.csv"))
  supp_selection_table_file <- file.path(out_dir, paste0(figure_id, "_supp_selection_funnel.csv"))
  panel_c_table_file <- file.path(out_dir, paste0(figure_id, "_panel_c_shared_marker_core.csv"))
  panel_d_table_file <- file.path(out_dir, paste0(figure_id, "_panel_d_marker_gain_vs_clinical.csv"))
  panel_e_table_file <- file.path(out_dir, paste0(figure_id, "_panel_e_gsea_leading_edge.csv"))
  panel_f_table_file <- file.path(out_dir, paste0(figure_id, "_panel_f_program_recurrence.csv"))
  manifest_file <- file.path(out_dir, paste0(figure_id, "_manifest.csv"))

  figure_width <- 9.0
  figure_height <- 11.6
  ggsave(png_file, main_figure,
    width = figure_width, height = figure_height, units = "in", dpi = 400,
    bg = "white", device = ragg::agg_png, limitsize = FALSE
  )
  ggsave(svg_file, main_figure,
    width = figure_width, height = figure_height, units = "in",
    bg = "white", device = svglite::svglite, limitsize = FALSE
  )
  ggsave(pdf_file, main_figure, width = figure_width, height = figure_height, bg = "white", device = grDevices::cairo_pdf)

  panel_data <- bind_rows_fill(list(
    transform(panel_a_upset, figure_panel = "A_gene_overlap_upset"),
    transform(enrichment_panel, figure_panel = "B_gsea_enrichment"),
    transform(top_gene_labels, figure_panel = "B_gsea_top_gene_labels"),
    transform(core, figure_panel = "C_shared_marker_core"),
    transform(gain_summary, figure_panel = "D_marker_gain_vs_clinical"),
    transform(leading_edge, figure_panel = "E_gsea_leading_edge"),
    transform(program_plot, figure_panel = "F_program_recurrence")
  ))
  write.csv(panel_data, panel_data_file, row.names = FALSE)
  write.csv(panel_a_upset[order(panel_a_upset$intersection_label, panel_a_upset$gene_symbol), ], panel_a_table_file, row.names = FALSE)
  write.csv(enrichment_panel[order(enrichment_panel$panel_rank), ], panel_b_table_file, row.names = FALSE)
  write.csv(top_gene_labels[order(top_gene_labels$panel_rank, top_gene_labels$top_gene_rank), ], panel_b_gene_label_file, row.names = FALSE)
  write.csv(selection_steps, supp_selection_table_file, row.names = FALSE)
  write.csv(core[order(core$gene_symbol, core$cancer_label), ], panel_c_table_file, row.names = FALSE)
  write.csv(gain_summary[order(gain_summary$gene_symbol), ], panel_d_table_file, row.names = FALSE)
  write.csv(leading_edge[order(leading_edge$gene_symbol, leading_edge$cancer_label), ], panel_e_table_file, row.names = FALSE)
  write.csv(program_plot[order(program_plot$program, program_plot$layer), ], panel_f_table_file, row.names = FALSE)
  manifest <- data.frame(
    figure_id = figure_id,
    panel_id = c("A", "B", "C", "D", "E", "F"),
    script = script_rel,
    input_table = c(
      overlap_file,
      paste(enrichment_panel_file, enrichment_ranked_file, enrichment_gsea_file, enrichment_ora_file, enrichment_curve_file, sep = "; "),
      core_file,
      paste(summary_file, contrast_file, sep = "; "),
      paste(enrichment_panel_file, enrichment_curve_file, overlap_file, sep = "; "),
      paste(program_summary_file, program_evidence_file, sep = "; ")
    ),
    output_file = png_file,
    width = figure_width,
    height = figure_height,
    dpi = 400,
    render_date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    notes = c(
      "UpSet-style overlap of TCGA RNA risk-associated genes from the CPTAC teacher high/low protein/phosphoprotein gene pool.",
      "MSigDB GSEA over the 362-gene CPTAC-derived source ranking; a 2x2 grid of compound mini-panels shows running ES, hit barcode/score strip, ranked-score distribution, and top three ranked hit genes per pathway.",
      "Final conserved RNA marker core: FDR < 0.05 and same Spearman direction in HNSC, KIRC, and LUAD.",
      "PROMISE-vs-clinical marker contrast: mean absolute RNA association with PROMISE student risk versus common clinical scaffold risk.",
      "Leading-edge gene support for the top GSEA term, shown as signed RNA-risk associations across HNSC/KIRC/LUAD.",
      "Reusable program recurrence as a clean Sankey-style flow: significant RNA cancer contexts use FDR < 0.05; RPPA contexts use FDR < 0.10."
    ),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, manifest_file, row.names = FALSE)

  message("Wrote Figure 5 to: ", png_file)
}

main()
