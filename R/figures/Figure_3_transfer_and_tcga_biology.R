# Purpose: Reproduce manuscript Figure 3 (risk transfer and TCGA biology).

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
  })

  FONT_FAMILY <- "Arial"
  figure_id <- "Figure_3_transfer_and_tcga_biology"
  out_dir <- argument_value("--output-dir=", file.path(package_root, "outputs", "figures"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  biology_audit_file <- file.path(input_root, "outputs/task6/task6_biology_audit.csv")
  cptac_scores_file <- file.path(input_root, "outputs/task6/cptac/task6_cptac_teacher_student_risk_scores.csv")
  rna_pathway_file <- file.path(input_root, "outputs/task6/tcga/task6_tcga_rna_pathway_risk_associations.csv")
  rna_pathway_scores_file <- file.path(input_root, "outputs/task6/tcga/task6_tcga_rna_pathway_scores.csv")
  rppa_module_file <- file.path(input_root, "outputs/task6b/rppa/task6b_rppa_module_risk_associations.csv")

  input_files <- c(
    biology_audit_file,
    cptac_scores_file,
    rna_pathway_file,
    rna_pathway_scores_file,
    rppa_module_file
  )
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0) {
    stop("Missing Figure 4 input files: ", paste(missing_inputs, collapse = ", "))
  }

  biology_audit <- read.csv(biology_audit_file, stringsAsFactors = FALSE)
  cptac_scores <- read.csv(cptac_scores_file, stringsAsFactors = FALSE)
  rna_pathways <- read.csv(rna_pathway_file, stringsAsFactors = FALSE)
  rna_pathway_scores <- read.csv(rna_pathway_scores_file, stringsAsFactors = FALSE)
  rppa_modules <- read.csv(rppa_module_file, stringsAsFactors = FALSE)

  cancer_levels <- c("HNSC", "KIRC", "LUAD")
  cancer_cols <- c(HNSC = "#B65C6A", KIRC = "#4E9F8E", LUAD = "#6C78A8")
  focus_programs <- c("PI3K/AKT", "Angiogenesis", "EMT/Stromal", "Proliferation", "Immune", "DNA damage")
  program_cols <- c(
    "PI3K/AKT" = "#E69F00",
    "Angiogenesis" = "#8C6D31",
    "EMT/Stromal" = "#CC79A7",
    "Proliferation" = "#D55E00",
    "Immune" = "#0072B2",
    "DNA damage" = "#009E73"
  )
  div_low <- "#4E79A7"
  div_mid <- "#F7F7F3"
  div_high <- "#B65C6A"
  neutral <- "#5C6470"
  facet_frame <- data.frame(
    cancer_label = factor(cancer_levels, levels = cancer_levels),
    border_col = unname(cancer_cols[cancer_levels]),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    stringsAsFactors = FALSE
  )

  fmt2 <- function(x) sprintf("%.2f", x)
  fmt3 <- function(x) sprintf("%.3f", x)
  fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  sig_mark <- function(fdr) ifelse(fdr < 0.05, "*", ifelse(fdr < 0.10, "+", ""))
  neglog_fdr <- function(fdr) pmin(-log10(pmax(fdr, 1e-50)), 12)
  z_by_group <- function(x, group) {
    ave(x, group, FUN = function(v) {
      s <- sd(v, na.rm = TRUE)
      if (!is.finite(s) || s < 1e-8) {
        rep(0, length(v))
      } else {
        (v - mean(v, na.rm = TRUE)) / s
      }
    })
  }

  teacher_transfer <- biology_audit[
    biology_audit$analysis == "teacher_vs_rna_clinical_student_risk" &
      biology_audit$cancer_label %in% cancer_levels,
  ]
  teacher_transfer$cancer_label <- factor(teacher_transfer$cancer_label, levels = cancer_levels)
  teacher_transfer$label <- paste0("rho=", fmt3(teacher_transfer$spearman_rho), "; n=", fmt_int(teacher_transfer$n))

  cptac_scores <- cptac_scores[cptac_scores$cancer_label %in% cancer_levels, ]
  cptac_scores$cancer_label <- factor(cptac_scores$cancer_label, levels = cancer_levels)
  cptac_scores$teacher_z <- z_by_group(cptac_scores$teacher_risk_score, cptac_scores$cancer_label)
  cptac_scores$student_z <- z_by_group(cptac_scores$rna_clinical_student_risk_score, cptac_scores$cancer_label)
  transfer_label <- merge(
    unique(cptac_scores[, c("cancer_label"), drop = FALSE]),
    teacher_transfer[, c("cancer_label", "label")],
    by = "cancer_label",
    all.x = TRUE
  )
  transfer_label$x <- -2.65
  transfer_label$y <- 2.35

  panel_a <- ggplot(cptac_scores, aes(x = teacher_z, y = student_z, color = cancer_label)) +
    geom_abline(slope = 1, intercept = 0, color = "#666666", linewidth = 0.35, linetype = "dashed") +
    geom_point(size = 1.45, alpha = 0.76) +
    geom_text(
      data = transfer_label,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      family = FONT_FAMILY,
      size = 3.90,
      hjust = 0,
      color = "#111111"
    ) +
    geom_rect(
      data = facet_frame,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      color = facet_frame$border_col,
      linewidth = 0.95
    ) +
    facet_wrap(~cancer_label, nrow = 1) +
    scale_color_manual(values = cancer_cols) +
    scale_x_continuous(breaks = c(-2, 0, 2), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(-2, 0, 2), expand = c(0, 0)) +
    coord_cartesian(xlim = c(-2.8, 2.8), ylim = c(-2.8, 2.8), expand = FALSE) +
    labs(x = "CPTAC multiomics teacher risk z-score", y = "RNA+clinical student risk z-score", tag = "A") +
    theme_classic(base_family = "Arial", base_size = 12.2) +
    theme(
      legend.position = "none",
      axis.line = element_line(color = "#222222", linewidth = 0.45),
      axis.ticks = element_line(color = "#222222", linewidth = 0.35),
      axis.text = element_text(family = "Arial", size = 12.2, color = "#222222"),
      axis.title = element_text(family = "Arial", size = 13.0, color = "#222222"),
      strip.text = element_text(family = "Arial", face = "bold", size = 12.2, color = "#222222"),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#B8B8B8", fill = NA, linewidth = 0.32),
      panel.spacing.x = unit(0.24, "inches"),
      plot.tag = element_text(family = "Arial", face = "bold", size = 22, color = "#111111"),
      plot.tag.position = c(0.060, 0.998),
      plot.margin = margin(5.5, 6.5, 5.5, 6.5)
    )

  rna_program_map <- c(
    rna_pi3k_akt_mtor_score = "PI3K/AKT",
    rna_angiogenesis_score = "Angiogenesis",
    rna_emt_score = "EMT/Stromal",
    rna_stromal_ecm_score = "EMT/Stromal",
    rna_dna_damage_score = "DNA damage",
    rna_proliferation_score = "Proliferation",
    rna_immune_cytotoxic_score = "Immune"
  )
  rna_pathways <- rna_pathways[rna_pathways$cancer_label %in% cancer_levels, ]
  rna_pathways$program <- unname(rna_program_map[rna_pathways$rna_signature])
  rna_pathways <- rna_pathways[!is.na(rna_pathways$program), ]
  rna_pathways$program <- factor(rna_pathways$program, levels = focus_programs)
  rna_pathways$cancer_label <- factor(rna_pathways$cancer_label, levels = rev(cancer_levels))
  rna_pathways$rho <- rna_pathways$spearman_rho_risk
  rna_pathways$label <- paste0(fmt2(rna_pathways$rho), sig_mark(rna_pathways$spearman_fdr))
  rna_pathways$neglog_fdr <- neglog_fdr(rna_pathways$spearman_fdr)

  gradient_rows <- list()
  rna_program_assoc_rows <- list()
  for (cancer in cancer_levels) {
    for (program in focus_programs) {
      signatures <- names(rna_program_map)[rna_program_map == program]
      signatures <- intersect(signatures, names(rna_pathway_scores))
      sub <- rna_pathway_scores[
        rna_pathway_scores$cancer_label == cancer,
        c("cancer_label", "patient_id", "risk_percentile", signatures),
        drop = FALSE
      ]
      if (length(signatures) == 0 || nrow(sub) == 0) {
        next
      }
      score <- rowMeans(sub[, signatures, drop = FALSE], na.rm = TRUE)
      score[is.nan(score)] <- NA_real_
      ok <- is.finite(score) & is.finite(sub$risk_percentile)
      assoc <- if (sum(ok) > 2) {
        suppressWarnings(cor.test(sub$risk_percentile[ok], score[ok], method = "spearman", exact = FALSE))
      } else {
        list(estimate = NA_real_, p.value = NA_real_)
      }
      rna_program_assoc_rows[[length(rna_program_assoc_rows) + 1]] <- data.frame(
        cancer_label = cancer,
        program = program,
        rna_program_n = sum(ok),
        rna_program_rho = unname(assoc$estimate),
        rna_program_p = assoc$p.value,
        stringsAsFactors = FALSE
      )
      sub$program_score <- score
      sub$risk_decile <- pmin(10, pmax(1, ceiling(sub$risk_percentile * 10)))
      agg <- aggregate(program_score ~ cancer_label + risk_decile, sub, mean, na.rm = TRUE)
      decile_n <- aggregate(program_score ~ cancer_label + risk_decile, sub, function(v) sum(is.finite(v)))
      names(decile_n)[names(decile_n) == "program_score"] <- "decile_n"
      agg <- merge(agg, decile_n, by = c("cancer_label", "risk_decile"), all.x = TRUE)
      agg$program <- program
      gradient_rows[[length(gradient_rows) + 1]] <- agg
    }
  }
  rna_program_assoc <- do.call(rbind, rna_program_assoc_rows)
  rna_program_assoc$rna_program_fdr <- p.adjust(rna_program_assoc$rna_program_p, method = "BH")
  risk_gradient <- do.call(rbind, gradient_rows)
  risk_gradient <- merge(risk_gradient, rna_program_assoc, by = c("cancer_label", "program"), all.x = TRUE)
  risk_gradient$cancer_label <- factor(risk_gradient$cancer_label, levels = cancer_levels)
  risk_gradient$program <- factor(risk_gradient$program, levels = focus_programs)

  panel_b <- ggplot(risk_gradient, aes(x = risk_decile, y = program_score, color = program, group = program)) +
    geom_hline(yintercept = 0, color = "#777777", linewidth = 0.30) +
    geom_line(linewidth = 0.72) +
    geom_point(size = 1.45) +
    geom_rect(
      data = facet_frame,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      color = facet_frame$border_col,
      linewidth = 0.95
    ) +
    facet_wrap(~cancer_label, nrow = 1) +
    scale_color_manual(values = program_cols, limits = focus_programs, drop = FALSE, name = "RNA programs") +
    scale_x_continuous(breaks = c(1, 5, 10), limits = c(1, 10), expand = c(0.02, 0)) +
    labs(x = "TCGA student-risk decile", y = "Mean RNA program score", tag = "B") +
    theme_classic(base_family = "Arial", base_size = 12.0) +
    theme(
      axis.line = element_line(color = "#222222", linewidth = 0.45),
      axis.ticks = element_line(color = "#222222", linewidth = 0.35),
      axis.text = element_text(family = "Arial", size = 12.0, color = "#222222"),
      axis.title = element_text(family = "Arial", size = 12.8, color = "#222222"),
      strip.text = element_text(family = "Arial", face = "bold", size = 12.0, color = "#222222"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(family = "Arial", size = 10.8, color = "#222222"),
      legend.text = element_text(family = "Arial", size = 10.3, color = "#222222"),
      legend.key.height = unit(0.10, "inches"),
      legend.key.width = unit(0.22, "inches"),
      legend.spacing.x = unit(0.05, "inches"),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#B8B8B8", fill = NA, linewidth = 0.32),
      panel.spacing.x = unit(0.24, "inches"),
      plot.tag = element_text(family = "Arial", face = "bold", size = 22, color = "#111111"),
      plot.tag.position = c(0.060, 0.998),
      plot.margin = margin(5.5, 6.5, 5.5, 6.5)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))

  rppa_program_map <- c(
    rppa_proliferation = "Proliferation",
    rppa_pi3k_akt_mtor = "PI3K/AKT",
    rppa_immune_stromal = "Immune",
    rppa_dna_damage = "DNA damage",
    rppa_emt_stromal = "EMT/Stromal"
  )
  rppa_modules <- rppa_modules[rppa_modules$cancer_label %in% cancer_levels, ]
  rppa_modules$program <- unname(rppa_program_map[rppa_modules$rppa_module])
  rppa_modules <- rppa_modules[!is.na(rppa_modules$program), ]
  rppa_modules$program <- factor(rppa_modules$program, levels = focus_programs)
  rppa_modules$cancer_label <- factor(rppa_modules$cancer_label, levels = rev(cancer_levels))
  rppa_modules$rho <- rppa_modules$spearman_rho_risk
  rppa_modules$label <- paste0(fmt2(rppa_modules$rho), sig_mark(rppa_modules$spearman_fdr))
  rppa_modules$neglog_fdr <- neglog_fdr(rppa_modules$spearman_fdr)

  rna_bubble <- data.frame(
    cancer_label = as.character(rna_pathways$cancer_label),
    program = as.character(rna_pathways$program),
    layer = "RNA",
    rho = rna_pathways$rho,
    fdr = rna_pathways$spearman_fdr,
    neglog_fdr = rna_pathways$neglog_fdr,
    stringsAsFactors = FALSE
  )
  rppa_bubble <- data.frame(
    cancer_label = as.character(rppa_modules$cancer_label),
    program = as.character(rppa_modules$program),
    layer = "RPPA",
    rho = rppa_modules$rho,
    fdr = rppa_modules$spearman_fdr,
    neglog_fdr = rppa_modules$neglog_fdr,
    stringsAsFactors = FALSE
  )
  layer_bubble <- rbind(rna_bubble, rppa_bubble)
  layer_bubble <- do.call(
    rbind,
    lapply(split(layer_bubble, paste(layer_bubble$cancer_label, layer_bubble$program, layer_bubble$layer, sep = "\r")), function(dat) {
      dat <- dat[order(dat$fdr, -abs(dat$rho)), ]
      dat[1, , drop = FALSE]
    })
  )
  layer_bubble <- layer_bubble[
    is.finite(layer_bubble$fdr) & (layer_bubble$fdr < 0.10 | abs(layer_bubble$rho) >= 0.15),
  ]
  layer_bubble$cancer_label <- factor(layer_bubble$cancer_label, levels = cancer_levels)
  layer_bubble$program <- factor(layer_bubble$program, levels = rev(focus_programs))
  layer_bubble$layer <- factor(layer_bubble$layer, levels = c("RNA", "RPPA"))
  layer_bubble$label <- ifelse(layer_bubble$fdr < 0.05, "*", ifelse(layer_bubble$fdr < 0.10, "+", ""))

  panel_c <- ggplot(layer_bubble, aes(x = cancer_label, y = program)) +
    geom_point(aes(size = neglog_fdr, fill = rho), shape = 21, color = "#444444", stroke = 0.25, alpha = 0.94) +
    geom_text(aes(label = label), family = "Arial", size = 3.50, color = "#111111", vjust = 0.35) +
    facet_grid(. ~ layer) +
    scale_fill_gradient2(
      low = div_low,
      mid = div_mid,
      high = div_high,
      midpoint = 0,
      limits = c(-0.32, 0.32),
      breaks = c(-0.30, 0, 0.30),
      name = "Spearman\nrho"
    ) +
    scale_size_continuous(range = c(1.8, 6.5), breaks = c(2, 5, 10), name = "-log10(FDR)") +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Cancer type", y = "Program", tag = "C") +
    theme_classic(base_family = "Arial", base_size = 11.6) +
    theme(
      axis.line = element_line(color = "#222222", linewidth = 0.45),
      axis.ticks = element_line(color = "#222222", linewidth = 0.35),
      axis.text.x = element_text(family = "Arial", face = "bold", size = 11.2, color = "#222222"),
      axis.text.y = element_text(family = "Arial", size = 11.2, color = "#222222"),
      axis.title = element_text(family = "Arial", size = 12.4, color = "#222222"),
      strip.text = element_text(family = "Arial", face = "bold", size = 11.2, color = "#222222"),
      legend.position = "bottom",
      legend.title = element_text(family = "Arial", size = 10.5, color = "#222222"),
      legend.text = element_text(family = "Arial", size = 10.0, color = "#222222"),
      legend.key.height = unit(0.12, "inches"),
      legend.key.width = unit(0.25, "inches"),
      legend.spacing.x = unit(0.05, "inches"),
      panel.grid.major = element_line(color = "#EFEFEF", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#B8B8B8", fill = NA, linewidth = 0.32),
      panel.spacing.x = unit(0.12, "inches"),
      plot.tag = element_text(family = "Arial", face = "bold", size = 22, color = "#111111"),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(5.5, 6.5, 5.5, 6.5)
    ) +
    guides(
      fill = guide_colorbar(order = 1, title.position = "top", barwidth = unit(0.72, "inches"), barheight = unit(0.10, "inches")),
      size = guide_legend(order = 2, title.position = "top", nrow = 1)
    )

  network_edges <- layer_bubble[order(layer_bubble$fdr, -abs(layer_bubble$rho)), ]
  network_edges <- network_edges[network_edges$fdr < 0.05 | (network_edges$cancer_label == "HNSC" & network_edges$program == "EMT/Stromal"), ]
  network_edges <- head(network_edges, 12)
  network_edges$cancer_label <- as.character(network_edges$cancer_label)
  network_edges$program <- as.character(network_edges$program)
  network_edges$layer <- as.character(network_edges$layer)

  cancer_nodes <- data.frame(
    node = cancer_levels,
    node_type = "Cancer",
    angle = seq(150, 210, length.out = length(cancer_levels)),
    stringsAsFactors = FALSE
  )
  program_nodes <- data.frame(
    node = focus_programs,
    node_type = "Program",
    angle = seq(-75, 105, length.out = length(focus_programs)),
    stringsAsFactors = FALSE
  )
  circle_nodes <- rbind(cancer_nodes, program_nodes)
  circle_nodes$theta <- circle_nodes$angle * pi / 180
  circle_nodes$x <- cos(circle_nodes$theta)
  circle_nodes$y <- sin(circle_nodes$theta)
  circle_nodes$label_x <- 1.17 * cos(circle_nodes$theta)
  circle_nodes$label_y <- 1.17 * sin(circle_nodes$theta)
  circle_nodes$hjust <- ifelse(cos(circle_nodes$theta) >= 0, 0, 1)
  circle_nodes$fontface <- ifelse(circle_nodes$node_type == "Cancer", "bold", "plain")
  circle_nodes$node_col <- ifelse(
    circle_nodes$node_type == "Cancer",
    unname(cancer_cols[circle_nodes$node]),
    unname(program_cols[circle_nodes$node])
  )

  arc_segments <- do.call(
    rbind,
    lapply(seq_len(nrow(circle_nodes)), function(i) {
      half_width <- if (circle_nodes$node_type[[i]] == "Cancer") 8 else 5.5
      theta <- seq(
        (circle_nodes$angle[[i]] - half_width) * pi / 180,
        (circle_nodes$angle[[i]] + half_width) * pi / 180,
        length.out = 36
      )
      data.frame(
        node = circle_nodes$node[[i]],
        node_type = circle_nodes$node_type[[i]],
        node_col = circle_nodes$node_col[[i]],
        x = cos(theta),
        y = sin(theta),
        stringsAsFactors = FALSE
      )
    })
  )

  cancer_xy <- circle_nodes[circle_nodes$node_type == "Cancer", c("node", "x", "y")]
  names(cancer_xy) <- c("cancer_label", "x_cancer", "y_cancer")
  program_xy <- circle_nodes[circle_nodes$node_type == "Program", c("node", "x", "y")]
  names(program_xy) <- c("program", "x_program", "y_program")
  network_edges <- merge(network_edges, cancer_xy, by = "cancer_label")
  network_edges <- merge(network_edges, program_xy, by = "program")
  network_edges$direction <- ifelse(network_edges$rho >= 0, "Higher with risk", "Lower with risk")
  network_edges$line_size <- pmax(0.35, pmin(1.60, network_edges$neglog_fdr / 5.2))
  network_edges$edge_id <- seq_len(nrow(network_edges))

  chord_paths <- do.call(
    rbind,
    lapply(seq_len(nrow(network_edges)), function(i) {
      edge <- network_edges[i, ]
      t <- seq(0, 1, length.out = 90)
      control_x <- 0.06 * sin(i * 1.7)
      control_y <- 0.06 * cos(i * 1.3)
      data.frame(
        edge_id = edge$edge_id,
        x = (1 - t)^2 * (0.92 * edge$x_cancer) + 2 * (1 - t) * t * control_x + t^2 * (0.92 * edge$x_program),
        y = (1 - t)^2 * (0.92 * edge$y_cancer) + 2 * (1 - t) * t * control_y + t^2 * (0.92 * edge$y_program),
        rho = edge$rho,
        line_size = edge$line_size,
        layer = edge$layer,
        stringsAsFactors = FALSE
      )
    })
  )

  panel_d <- ggplot() +
    geom_path(
      data = chord_paths,
      aes(x = x, y = y, group = edge_id, color = rho, linewidth = line_size, linetype = layer),
      alpha = 0.78,
      lineend = "round"
    ) +
    geom_point(
      data = circle_nodes[circle_nodes$node_type == "Cancer", ],
      aes(x = x, y = y),
      color = "#FFFFFF",
      fill = circle_nodes$node_col[circle_nodes$node_type == "Cancer"],
      shape = 21,
      stroke = 0.45,
      size = 5.0
    ) +
    geom_path(
      data = arc_segments[arc_segments$node_type == "Program", ],
      aes(x = x, y = y, group = node),
      color = arc_segments$node_col[arc_segments$node_type == "Program"],
      linewidth = 3.1,
      lineend = "butt"
    ) +
    geom_text(
      data = circle_nodes[circle_nodes$node_type == "Cancer", ],
      aes(x = label_x, y = label_y, label = node, hjust = hjust),
      family = "Arial",
      fontface = "bold",
      size = 3.80,
      color = "#4A4A4A",
      lineheight = 0.92
    ) +
    geom_text(
      data = circle_nodes[circle_nodes$node_type == "Program", ],
      aes(x = label_x, y = label_y, label = node, hjust = hjust, fontface = fontface),
      family = "Arial",
      size = 3.55,
      color = "#1D1D1D",
      lineheight = 0.92
    ) +
    scale_color_gradient2(
      low = div_low, mid = div_mid, high = div_high, midpoint = 0,
      limits = c(-0.32, 0.32), name = "rho"
    ) +
    scale_linewidth_identity() +
    scale_linetype_manual(values = c(RNA = "solid", RPPA = "22"), name = "Layer") +
    coord_fixed(xlim = c(-1.36, 1.48), ylim = c(-1.24, 1.22), expand = FALSE, clip = "off") +
    labs(tag = "D") +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      legend.position = "none",
      plot.tag = element_text(family = "Arial", face = "bold", size = 22, color = "#111111"),
      plot.tag.position = c(0.002, 0.998),
      plot.margin = margin(5.5, 11, 5.5, 11)
    )

  main_figure <- panel_a / panel_b / (panel_c + panel_d + plot_layout(widths = c(1.00, 1.05))) +
    plot_layout(heights = c(0.95, 1.08, 1.02)) &
    theme(plot.tag = element_text(family = "Arial", face = "bold", size = 22, color = "#111111"))

  png_file <- file.path(out_dir, "Figure_3_transfer_and_tcga_biology_large_text_400dpi.png")
  svg_file <- file.path(out_dir, "Figure_3_transfer_and_tcga_biology_large_text.svg")
  pdf_file <- file.path(out_dir, "Figure_3_transfer_and_tcga_biology_large_text.pdf")
  ggsave(png_file, main_figure,
    width = 9.0, height = 8.0, units = "in", dpi = 400,
    bg = "white", device = ragg::agg_png, limitsize = FALSE
  )
  ggsave(svg_file, main_figure,
    width = 9.0, height = 8.0, units = "in",
    bg = "white", device = svglite::svglite, limitsize = FALSE
  )
  ggsave(pdf_file, main_figure,
    width = 9.0, height = 8.0, units = "in",
    bg = "white", device = grDevices::cairo_pdf, limitsize = FALSE
  )

  write.csv(
    data.frame(
      figure = "Figure 3",
      panel = c("A", "B", "C", "D"),
      content = c(
        "CPTAC teacher-student transfer",
        "TCGA RNA-program gradients",
        "RNA/RPPA risk associations",
        "Cross-cancer program network"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Figure_3_panel_manifest.csv"),
    row.names = FALSE
  )
  message("Wrote ", png_file)
}

main()
