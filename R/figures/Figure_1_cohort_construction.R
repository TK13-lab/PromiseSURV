# Purpose: Reproduce manuscript Figure 1 (cohort construction and modality availability).
# Inputs: outputs/figures/figure2/figure2_cohort_construction_panel_data.csv
# Outputs: the exact 400-dpi Figure 1 PNG plus vector PDF/SVG.
# Usage: Rscript R/figures/Figure_1_cohort_construction.R [--input-root=PATH] [--output-dir=PATH]

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  argument_value <- function(prefix, default) {
    hit <- args[startsWith(args, prefix)]
    if (length(hit) == 0) default else sub(prefix, "", hit[[1]], fixed = TRUE)
  }
  script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][[1]])
  package_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
  input_root <- normalizePath(argument_value("--input-root=", package_root), winslash = "/", mustWork = TRUE)

  library(grid)

  FONT_FAMILY <- "Arial"
  input_file <- file.path(input_root, "outputs/figures/figure2/figure2_cohort_construction_panel_data.csv")
  out_dir <- argument_value("--output-dir=", file.path(package_root, "outputs", "figures"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_file <- file.path(out_dir, "Figure_1_manifest.csv")

  dat <- read.csv(input_file, stringsAsFactors = FALSE)
  required_cols <- c("figure_panel", "source", "cancer_label", "n")
  stopifnot(all(required_cols %in% colnames(dat)))

  cohort <- dat[dat$figure_panel == "B" & !is.na(dat$source), ]
  follow <- dat[dat$figure_panel == "C" & !is.na(dat$followup_years), ]
  sex <- dat[dat$figure_panel == "E" & !is.na(dat$sex), ]
  modalities <- dat[dat$figure_panel == "F" & !is.na(dat$modality), ]

  cancers <- c("LUAD", "KIRC", "HNSC")
  cancer_cols <- c(LUAD = "#79C6D3", KIRC = "#8DCB8E", HNSC = "#D88965")
  source_cols <- c(CPTAC = "#3C63A2", TCGA = "#9B7898")
  outcome_cols <- c(event = "#B24455", censored = "#C9D4E1")
  sex_cols <- c(Female = "#9B7898", Male = "#3C63A2", Missing = "#B9C0CA")
  modality_cols <- c(
    Clinical = "#E7EDF3",
    Proteomics = "#BBD7EA",
    Phosphoproteomics = "#F1C6C2",
    RNA = "#BFE0C8"
  )

  fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

  source_n <- function(source) {
    sum(cohort$n[cohort$source == source], na.rm = TRUE)
  }

  source_label <- function(source) {
    if (source == "CPTAC") {
      paste0("Discovery: CPTAC/PDC\n(n = ", fmt_int(source_n(source)), ")")
    } else {
      paste0("Deployment: TCGA\n(n = ", fmt_int(source_n(source)), ")")
    }
  }

  cancer_counts <- function(source) {
    out <- setNames(rep(0, length(cancers)), cancers)
    for (cancer in cancers) {
      value <- cohort$n[cohort$source == source & cohort$cancer_label == cancer]
      out[[cancer]] <- if (length(value) == 0) 0 else value[[1]]
    }
    out
  }

  sex_counts <- function(source) {
    labels <- c("Female", "Male", "Missing")
    out <- setNames(rep(0, length(labels)), labels)
    for (label in labels) {
      value <- sex$n[sex$source == source & sex$sex == label]
      out[[label]] <- sum(value, na.rm = TRUE)
    }
    out
  }

  sex_counts_by_cancer <- function(source, cancer) {
    labels <- c("Female", "Male", "Missing")
    out <- setNames(rep(0, length(labels)), labels)
    for (label in labels) {
      value <- sex$n[sex$source == source & sex$cancer_label == cancer & sex$sex == label]
      out[[label]] <- sum(value, na.rm = TRUE)
    }
    out
  }

  outcome_counts <- function(source) {
    rows <- cohort[cohort$source == source, ]
    c(
      event = sum(rows$events, na.rm = TRUE),
      censored = sum(rows$non_events, na.rm = TRUE)
    )
  }

  outcome_counts_by_cancer <- function(source, cancer) {
    rows <- cohort[cohort$source == source & cohort$cancer_label == cancer, ]
    c(
      event = sum(rows$events, na.rm = TRUE),
      censored = sum(rows$non_events, na.rm = TRUE)
    )
  }

  followup_values <- function(source) {
    follow$followup_years[follow$source == source & !is.na(follow$followup_years)]
  }

  followup_values_by_cancer <- function(source, cancer) {
    follow$followup_years[follow$source == source & follow$cancer_label == cancer & !is.na(follow$followup_years)]
  }

  modality_rows <- function(source, cancer) {
    if (source == "TCGA") {
      rows <- modalities[grepl("TCGA", modalities$source, ignore.case = TRUE) & modalities$cancer_label == cancer, ]
    } else {
      rows <- modalities[!grepl("TCGA", modalities$source, ignore.case = TRUE) & modalities$cancer_label == cancer, ]
    }
    rows
  }

  row_layout <- data.frame(
    label = c("Cancer\ncohort", "Data\nmodalities", "Sex", "Time to\nfollow-up", "Clinical\noutcomes"),
    y_top = c(0.88, 0.70, 0.52, 0.34, 0.16),
    y_bottom = c(0.70, 0.52, 0.34, 0.16, 0.02),
    stringsAsFactors = FALSE
  )

  label_x0 <- 0.015
  label_w <- 0.145
  cell_gap <- 0.014
  cptac_x0 <- label_x0 + label_w + cell_gap
  cell_w <- (0.985 - cptac_x0 - cell_gap) / 2
  tcga_x0 <- cptac_x0 + cell_w + cell_gap
  header_y0 <- 0.88
  header_h <- 0.10

  draw_cell_background <- function() {
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.rect(
      x = (cptac_x0 + cell_w / 2), y = header_y0 + header_h / 2,
      width = cell_w, height = header_h,
      gp = gpar(fill = "#D9D9D9", col = "#8C8C8C", lwd = 0.8)
    )
    grid.rect(
      x = (tcga_x0 + cell_w / 2), y = header_y0 + header_h / 2,
      width = cell_w, height = header_h,
      gp = gpar(fill = "#D9D9D9", col = "#8C8C8C", lwd = 0.8)
    )
    grid.text(source_label("CPTAC"),
      x = cptac_x0 + cell_w / 2, y = header_y0 + header_h / 2,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 13, fontface = "bold", col = "#222222")
    )
    grid.text(source_label("TCGA"),
      x = tcga_x0 + cell_w / 2, y = header_y0 + header_h / 2,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 13, fontface = "bold", col = "#222222")
    )

    line_x <- c(label_x0, cptac_x0 - cell_gap / 2, cptac_x0 + cell_w + cell_gap / 2, 0.985)
    for (x in line_x) {
      grid.lines(x = c(x, x), y = c(0.02, 0.98), gp = gpar(col = "#8D8D8D", lwd = 0.6))
    }
    for (y in c(0.98, row_layout$y_top, 0.02)) {
      grid.lines(x = c(label_x0, 0.985), y = c(y, y), gp = gpar(col = "#8D8D8D", lwd = 0.6))
    }
  }

  draw_row_label <- function(i) {
    y_mid <- (row_layout$y_top[i] + row_layout$y_bottom[i]) / 2
    grid.text(row_layout$label[i],
      x = label_x0 + 0.02, y = y_mid, just = "left",
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 12.5, col = "#111111")
    )
  }

  draw_icon <- function(kind, x0, y0, w, h) {
    cx <- x0 + w / 2
    cy <- y0 + h / 2
    if (kind == "cohort") {
      for (i in 1:3) {
        grid.circle(cx + (i - 2) * w * 0.18, cy + h * 0.08,
          r = unit(w * 0.055, "npc"),
          gp = gpar(fill = cancer_cols[[i]], col = NA)
        )
        grid.rect(cx + (i - 2) * w * 0.18, cy - h * 0.07,
          width = w * 0.09, height = h * 0.18,
          gp = gpar(fill = cancer_cols[[i]], col = NA)
        )
      }
    } else if (kind == "omics") {
      grid.rect(cx - w * 0.22, cy - h * 0.18,
        width = w * 0.16, height = h * 0.36,
        gp = gpar(fill = modality_cols[["Proteomics"]], col = "#777777")
      )
      grid.rect(cx, cy - h * 0.18,
        width = w * 0.16, height = h * 0.36,
        gp = gpar(fill = modality_cols[["Phosphoproteomics"]], col = "#777777")
      )
      grid.rect(cx + w * 0.22, cy - h * 0.18,
        width = w * 0.16, height = h * 0.36,
        gp = gpar(fill = modality_cols[["RNA"]], col = "#777777")
      )
      grid.lines(
        x = c(cx - w * 0.22, cx, cx + w * 0.22), y = c(cy + h * 0.24, cy + h * 0.34, cy + h * 0.24),
        gp = gpar(col = "#52616B", lwd = 2)
      )
    } else if (kind == "sex") {
      grid.text("\u2640", x = cx - w * 0.12, y = cy, gp = gpar(fontfamily = FONT_FAMILY, fontsize = 25, col = sex_cols[["Female"]]))
      grid.text("\u2642", x = cx + w * 0.14, y = cy, gp = gpar(fontfamily = FONT_FAMILY, fontsize = 25, col = sex_cols[["Male"]]))
    } else if (kind == "clock") {
      grid.circle(cx, cy, r = unit(min(w, h) * 0.22, "npc"), gp = gpar(fill = "#F7F7F7", col = "#8D8D8D", lwd = 1.4))
      grid.lines(x = c(cx, cx), y = c(cy, cy + h * 0.15), gp = gpar(col = "#777777", lwd = 1.4))
      grid.lines(x = c(cx, cx + w * 0.11), y = c(cy, cy - h * 0.05), gp = gpar(col = "#777777", lwd = 1.4))
    } else if (kind == "outcome") {
      xs <- seq(x0 + w * 0.18, x0 + w * 0.82, length.out = 5)
      ys <- c(y0 + h * 0.72, y0 + h * 0.66, y0 + h * 0.58, y0 + h * 0.50, y0 + h * 0.38)
      grid.lines(xs, ys, gp = gpar(col = outcome_cols[["event"]], lwd = 2))
      grid.text("OS", x = cx, y = y0 + h * 0.24, gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, col = "#333333"))
    }
  }

  draw_segment_bar <- function(x0, y0, w, h, values, colors, labels = names(values), cex = 10.5) {
    total <- sum(values, na.rm = TRUE)
    if (total <= 0) {
      return(invisible(NULL))
    }
    x <- x0
    for (i in seq_along(values)) {
      frac <- values[[i]] / total
      seg_w <- w * frac
      grid.rect(
        x = x + seg_w / 2, y = y0 + h / 2, width = seg_w, height = h,
        gp = gpar(fill = colors[[i]], col = "white", lwd = 0.8)
      )
      if (frac >= 0.09 && nzchar(labels[[i]])) {
        txt <- paste0(labels[[i]], "\n", fmt_int(values[[i]]), " (", round(100 * frac), "%)")
        grid.text(txt,
          x = x + seg_w / 2, y = y0 + h / 2,
          gp = gpar(fontfamily = FONT_FAMILY, fontsize = cex, col = "#111111", lineheight = 0.95)
        )
      }
      x <- x + seg_w
    }
  }

  draw_lung_icon <- function(cx, cy, s, col) {
    fill_col <- grDevices::adjustcolor(col, alpha.f = 0.78)
    edge_col <- "#4B5960"
    grid.rect(cx, cy + s * 0.18,
      width = s * 0.055, height = s * 0.23,
      gp = gpar(fill = edge_col, col = NA)
    )
    grid.lines(
      x = c(cx, cx - s * 0.14), y = c(cy + s * 0.06, cy - s * 0.08),
      gp = gpar(col = edge_col, lwd = 1.2)
    )
    grid.lines(
      x = c(cx, cx + s * 0.14), y = c(cy + s * 0.06, cy - s * 0.08),
      gp = gpar(col = edge_col, lwd = 1.2)
    )
    grid.circle(cx - s * 0.15, cy - s * 0.02,
      r = unit(s * 0.17, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.circle(cx + s * 0.15, cy - s * 0.02,
      r = unit(s * 0.17, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.circle(cx, cy - s * 0.02,
      r = unit(s * 0.08, "npc"),
      gp = gpar(fill = "white", col = NA)
    )
  }

  draw_kidney_icon <- function(cx, cy, s, col) {
    fill_col <- grDevices::adjustcolor(col, alpha.f = 0.82)
    edge_col <- "#4B5960"
    grid.circle(cx - s * 0.12, cy,
      r = unit(s * 0.16, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.circle(cx + s * 0.10, cy,
      r = unit(s * 0.16, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.circle(cx - s * 0.03, cy + s * 0.01,
      r = unit(s * 0.12, "npc"),
      gp = gpar(fill = "white", col = NA)
    )
    grid.lines(
      x = c(cx + s * 0.16, cx + s * 0.27), y = c(cy, cy + s * 0.05),
      gp = gpar(col = edge_col, lwd = 1.4)
    )
    grid.lines(
      x = c(cx + s * 0.16, cx + s * 0.27), y = c(cy - s * 0.02, cy - s * 0.08),
      gp = gpar(col = edge_col, lwd = 1.4)
    )
  }

  draw_uterus_icon <- function(cx, cy, s, col) {
    fill_col <- grDevices::adjustcolor(col, alpha.f = 0.82)
    edge_col <- "#4B5960"
    grid.circle(cx, cy - s * 0.03,
      r = unit(s * 0.15, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.rect(cx, cy - s * 0.20,
      width = s * 0.07, height = s * 0.18,
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.lines(
      x = c(cx - s * 0.10, cx - s * 0.26), y = c(cy + s * 0.06, cy + s * 0.15),
      gp = gpar(col = edge_col, lwd = 1.5)
    )
    grid.lines(
      x = c(cx + s * 0.10, cx + s * 0.26), y = c(cy + s * 0.06, cy + s * 0.15),
      gp = gpar(col = edge_col, lwd = 1.5)
    )
    grid.circle(cx - s * 0.30, cy + s * 0.16,
      r = unit(s * 0.045, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.6)
    )
    grid.circle(cx + s * 0.30, cy + s * 0.16,
      r = unit(s * 0.045, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.6)
    )
  }

  draw_head_neck_icon <- function(cx, cy, s, col) {
    fill_col <- grDevices::adjustcolor(col, alpha.f = 0.82)
    edge_col <- "#4B5960"
    grid.circle(cx, cy + s * 0.08,
      r = unit(s * 0.17, "npc"),
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.rect(cx - s * 0.02, cy - s * 0.18,
      width = s * 0.15, height = s * 0.24,
      gp = gpar(fill = fill_col, col = edge_col, lwd = 0.7)
    )
    grid.circle(cx + s * 0.12, cy + s * 0.10,
      r = unit(s * 0.035, "npc"),
      gp = gpar(fill = "white", col = edge_col, lwd = 0.5)
    )
    grid.lines(
      x = c(cx - s * 0.10, cx + s * 0.16), y = c(cy + s * 0.01, cy - s * 0.01),
      gp = gpar(col = edge_col, lwd = 1.1)
    )
    grid.lines(
      x = c(cx - s * 0.07, cx + s * 0.13), y = c(cy - s * 0.07, cy - s * 0.08),
      gp = gpar(col = edge_col, lwd = 1.0)
    )
  }

  draw_organ_icon <- function(cancer, cx, cy, s) {
    if (cancer == "LUAD") {
      draw_lung_icon(cx, cy, s, cancer_cols[[cancer]])
    } else if (cancer == "KIRC") {
      draw_kidney_icon(cx, cy, s, cancer_cols[[cancer]])
    } else if (cancer == "HNSC") {
      draw_head_neck_icon(cx, cy, s, cancer_cols[[cancer]])
    }
  }

  draw_cancer_composition <- function(x0, y0, w, h, source) {
    counts <- cancer_counts(source)
    total <- sum(counts, na.rm = TRUE)
    if (total <= 0) {
      return(invisible(NULL))
    }
    gap <- w * 0.025
    tile_w <- (w - gap * 2) / 3
    for (i in seq_along(cancers)) {
      cancer <- cancers[[i]]
      bx <- x0 + (i - 1) * (tile_w + gap)
      pct <- round(100 * counts[[cancer]] / total)
      draw_organ_icon(cancer, bx + tile_w / 2, y0 + h * 0.62, min(tile_w, h) * 0.62)
      grid.text(cancer,
        x = bx + tile_w / 2, y = y0 + h * 0.30,
        gp = gpar(fontfamily = FONT_FAMILY, fontsize = 11.2, fontface = "bold", col = "#111111")
      )
      grid.text(paste0(fmt_int(counts[[cancer]]), " (", pct, "%)"),
        x = bx + tile_w / 2, y = y0 + h * 0.13,
        gp = gpar(fontfamily = FONT_FAMILY, fontsize = 10, col = "#333333")
      )
    }
  }

  draw_modality_by_cancer <- function(x0, y0, w, h, source) {
    gap <- w * 0.035
    block_w <- (w - gap * 2) / 3
    modality_order <- if (source == "CPTAC") {
      c("clinical", "rna", "proteomics", "phosphoproteomics")
    } else {
      c("clinical", "rna")
    }
    modality_labels <- c(clinical = "Clin", proteomics = "Prot", phosphoproteomics = "Phos", rna = "RNA")
    modality_fill <- c(
      clinical = modality_cols[["Clinical"]],
      proteomics = modality_cols[["Proteomics"]],
      phosphoproteomics = modality_cols[["Phosphoproteomics"]],
      rna = modality_cols[["RNA"]]
    )
    for (i in seq_along(cancers)) {
      cancer <- cancers[[i]]
      bx <- x0 + (i - 1) * (block_w + gap)
      rows <- modality_rows(source, cancer)
      grid.text(cancer,
        x = bx + block_w / 2, y = y0 + h * 0.84,
        gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, fontface = "bold", col = cancer_cols[[cancer]])
      )
      tile_gap <- block_w * 0.035
      tile_w <- (block_w - tile_gap * (length(modality_order) - 1)) / length(modality_order)
      for (j in seq_along(modality_order)) {
        modality <- modality_order[[j]]
        row <- rows[rows$modality == modality, ]
        available <- if (nrow(row) > 0) row$available_count[[1]] else 0
        tx <- bx + (j - 1) * (tile_w + tile_gap)
        fill <- if (available > 0) modality_fill[[modality]] else "#F1F3F5"
        col <- if (available > 0) "#8A8A8A" else "#C4C7CC"
        grid.roundrect(
          x = tx + tile_w / 2, y = y0 + h * 0.45, width = tile_w, height = h * 0.38,
          r = unit(0.006, "npc"), gp = gpar(fill = fill, col = col, lwd = 0.5)
        )
        grid.text(paste0(modality_labels[[modality]], "\n", fmt_int(available)),
          x = tx + tile_w / 2, y = y0 + h * 0.45,
          gp = gpar(fontfamily = FONT_FAMILY, fontsize = 7.4, col = "#111111", lineheight = 0.92)
        )
      }
    }
  }

  draw_people <- function(x0, y0, w, h, counts) {
    total <- sum(counts, na.rm = TRUE)
    if (total <= 0) {
      return(invisible(NULL))
    }
    n_icons <- 20
    female_n <- round(n_icons * counts[["Female"]] / total)
    male_n <- n_icons - female_n
    cols <- c(rep(sex_cols[["Female"]], female_n), rep(sex_cols[["Male"]], male_n))
    cols <- cols[seq_len(n_icons)]
    cols[is.na(cols)] <- sex_cols[["Missing"]]
    xs <- seq(x0 + w * 0.06, x0 + w * 0.94, length.out = n_icons)
    for (i in seq_len(n_icons)) {
      grid.circle(xs[[i]], y0 + h * 0.65, r = unit(w * 0.012, "npc"), gp = gpar(fill = cols[[i]], col = NA))
      grid.rect(xs[[i]], y0 + h * 0.43, width = w * 0.018, height = h * 0.23, gp = gpar(fill = cols[[i]], col = NA))
      grid.lines(
        x = c(xs[[i]] - w * 0.018, xs[[i]] + w * 0.018), y = c(y0 + h * 0.51, y0 + h * 0.51),
        gp = gpar(col = cols[[i]], lwd = 1.5)
      )
      grid.lines(
        x = c(xs[[i]] - w * 0.01, xs[[i]] - w * 0.018), y = c(y0 + h * 0.31, y0 + h * 0.19),
        gp = gpar(col = cols[[i]], lwd = 1.5)
      )
      grid.lines(
        x = c(xs[[i]] + w * 0.01, xs[[i]] + w * 0.018), y = c(y0 + h * 0.31, y0 + h * 0.19),
        gp = gpar(col = cols[[i]], lwd = 1.5)
      )
    }
    grid.text(
      paste0("Female ", round(100 * counts[["Female"]] / total), "%  |  Male ", round(100 * counts[["Male"]] / total), "%"),
      x = x0 + w / 2, y = y0 + h * 0.08,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 10.5, col = "#222222")
    )
  }

  draw_sex_by_cancer <- function(x0, y0, w, h, source) {
    gap <- w * 0.035
    block_w <- (w - gap * 2) / 3
    for (i in seq_along(cancers)) {
      cancer <- cancers[[i]]
      counts <- sex_counts_by_cancer(source, cancer)
      total <- sum(counts, na.rm = TRUE)
      bx <- x0 + (i - 1) * (block_w + gap)
      grid.text(cancer,
        x = bx + block_w / 2, y = y0 + h * 0.86,
        gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, fontface = "bold", col = cancer_cols[[cancer]])
      )
      if (total > 0) {
        vals <- counts[c("Female", "Male")]
        cols <- sex_cols[c("Female", "Male")]
        draw_segment_bar(bx, y0 + h * 0.47, block_w, h * 0.22, vals, cols, rep("", length(vals)), cex = 8)
        grid.text(
          paste0("F ", round(100 * counts[["Female"]] / total), "% | M ", round(100 * counts[["Male"]] / total), "%"),
          x = bx + block_w / 2, y = y0 + h * 0.28,
          gp = gpar(fontfamily = FONT_FAMILY, fontsize = 8.6, col = "#222222")
        )
        grid.text(paste0("n=", fmt_int(total)),
          x = bx + block_w / 2, y = y0 + h * 0.12,
          gp = gpar(fontfamily = FONT_FAMILY, fontsize = 8.2, col = "#555555")
        )
      }
    }
  }

  draw_density_cell <- function(x0, y0, w, h, values, xlab, unit_label = "y") {
    values <- values[is.finite(values)]
    if (length(values) < 3) {
      return(invisible(NULL))
    }
    dens <- density(values, n = 100, from = 0, to = max(10, quantile(values, 0.98, na.rm = TRUE)))
    xs <- x0 + w * 0.08 + (dens$x - min(dens$x)) / (max(dens$x) - min(dens$x)) * w * 0.78
    ys <- y0 + h * 0.18 + dens$y / max(dens$y) * h * 0.62
    grid.lines(xs, ys, gp = gpar(col = "#575757", lwd = 2))
    grid.lines(x = c(x0 + w * 0.08, x0 + w * 0.88), y = c(y0 + h * 0.18, y0 + h * 0.18), gp = gpar(col = "#777777", lwd = 0.8))
    grid.lines(x = c(x0 + w * 0.08, x0 + w * 0.08), y = c(y0 + h * 0.18, y0 + h * 0.80), gp = gpar(col = "#777777", lwd = 0.8))
    med <- median(values, na.rm = TRUE)
    med_x <- x0 + w * 0.08 + (med - min(dens$x)) / (max(dens$x) - min(dens$x)) * w * 0.78
    grid.lines(x = c(med_x, med_x), y = c(y0 + h * 0.18, y0 + h * 0.76), gp = gpar(col = "#D33F49", lwd = 2))
    grid.text(paste0(round(med, 1), unit_label),
      x = med_x + w * 0.05, y = y0 + h * 0.72,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 13, col = "#111111")
    )
    grid.text(xlab, x = x0 + w * 0.50, y = y0 + h * 0.04, gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, col = "#333333"))
  }

  draw_mini_density <- function(x0, y0, w, h, values, cancer) {
    values <- values[is.finite(values)]
    if (length(values) < 3) {
      return(invisible(NULL))
    }
    upper <- max(10, quantile(values, 0.98, na.rm = TRUE))
    dens <- density(values, n = 80, from = 0, to = upper)
    xs <- x0 + w * 0.06 + (dens$x - min(dens$x)) / (max(dens$x) - min(dens$x)) * w * 0.86
    ys <- y0 + h * 0.23 + dens$y / max(dens$y) * h * 0.42
    grid.lines(xs, ys, gp = gpar(col = "#575757", lwd = 1.5))
    grid.lines(x = c(x0 + w * 0.06, x0 + w * 0.92), y = c(y0 + h * 0.23, y0 + h * 0.23), gp = gpar(col = "#777777", lwd = 0.6))
    med <- median(values, na.rm = TRUE)
    med_x <- x0 + w * 0.06 + (med - min(dens$x)) / (max(dens$x) - min(dens$x)) * w * 0.86
    grid.lines(x = c(med_x, med_x), y = c(y0 + h * 0.23, y0 + h * 0.66), gp = gpar(col = "#D33F49", lwd = 1.4))
    grid.text(cancer,
      x = x0 + w / 2, y = y0 + h * 0.84,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, fontface = "bold", col = cancer_cols[[cancer]])
    )
    grid.text(sprintf("%.1fy", med),
      x = x0 + w / 2, y = y0 + h * 0.08,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 8.8, col = "#222222")
    )
  }

  draw_followup_by_cancer <- function(x0, y0, w, h, source) {
    gap <- w * 0.035
    block_w <- (w - gap * 2) / 3
    for (i in seq_along(cancers)) {
      cancer <- cancers[[i]]
      bx <- x0 + (i - 1) * (block_w + gap)
      draw_mini_density(bx, y0, block_w, h, followup_values_by_cancer(source, cancer), cancer)
    }
  }

  draw_outcomes <- function(x0, y0, w, h, counts) {
    cols <- outcome_cols
    draw_segment_bar(x0 + w * 0.04, y0 + h * 0.36, w * 0.90, h * 0.34, counts, cols, rep("", length(counts)), cex = 11)
    total <- sum(counts)
    grid.text(
      paste0("OS events: ", fmt_int(counts[["event"]]), " / ", fmt_int(total)),
      x = x0 + w / 2, y = y0 + h * 0.16,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 11, col = "#222222")
    )
    grid.text("OS event",
      x = x0 + w * 0.18, y = y0 + h * 0.80,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, col = "#333333")
    )
    grid.text("Censored",
      x = x0 + w * 0.68, y = y0 + h * 0.80,
      gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, col = "#333333")
    )
  }

  draw_outcomes_by_cancer <- function(x0, y0, w, h, source) {
    gap <- w * 0.035
    block_w <- (w - gap * 2) / 3
    for (i in seq_along(cancers)) {
      cancer <- cancers[[i]]
      counts <- outcome_counts_by_cancer(source, cancer)
      total <- sum(counts, na.rm = TRUE)
      bx <- x0 + (i - 1) * (block_w + gap)
      grid.text(cancer,
        x = bx + block_w / 2, y = y0 + h * 0.82,
        gp = gpar(fontfamily = FONT_FAMILY, fontsize = 9.5, fontface = "bold", col = cancer_cols[[cancer]])
      )
      draw_segment_bar(bx, y0 + h * 0.48, block_w, h * 0.20, counts, outcome_cols, rep("", length(counts)), cex = 8)
      if (total > 0) {
        grid.text(
          paste0(fmt_int(counts[["event"]]), "/", fmt_int(total), " OS ev"),
          x = bx + block_w / 2, y = y0 + h * 0.26,
          gp = gpar(fontfamily = FONT_FAMILY, fontsize = 8.6, col = "#222222")
        )
        grid.text(
          paste0(round(100 * counts[["event"]] / total), "%"),
          x = bx + block_w / 2, y = y0 + h * 0.10,
          gp = gpar(fontfamily = FONT_FAMILY, fontsize = 8.2, col = "#555555")
        )
      }
    }
  }

  draw_landscape <- function() {
    grid.newpage()
    pushViewport(viewport(gp = gpar(fontfamily = FONT_FAMILY)))
    draw_cell_background()
    for (i in seq_len(nrow(row_layout))) {
      draw_row_label(i)
    }

    for (source in c("CPTAC", "TCGA")) {
      x0 <- if (source == "CPTAC") cptac_x0 else tcga_x0
      row1 <- row_layout[1, ]
      draw_cancer_composition(x0 + cell_w * 0.04, row1$y_bottom + 0.025, cell_w * 0.92, (row1$y_top - row1$y_bottom) - 0.05, source)

      row2 <- row_layout[2, ]
      draw_modality_by_cancer(x0 + cell_w * 0.04, row2$y_bottom + 0.020, cell_w * 0.92, (row2$y_top - row2$y_bottom) - 0.04, source)

      row3 <- row_layout[3, ]
      draw_sex_by_cancer(x0 + cell_w * 0.04, row3$y_bottom + 0.020, cell_w * 0.92, (row3$y_top - row3$y_bottom) - 0.04, source)

      row4 <- row_layout[4, ]
      draw_followup_by_cancer(x0 + cell_w * 0.04, row4$y_bottom + 0.020, cell_w * 0.92, (row4$y_top - row4$y_bottom) - 0.04, source)

      row5 <- row_layout[5, ]
      draw_outcomes_by_cancer(x0 + cell_w * 0.04, row5$y_bottom + 0.020, cell_w * 0.92, (row5$y_top - row5$y_bottom) - 0.04, source)
    }
    popViewport()
  }

  png_file <- file.path(out_dir, "Figure_1_cohort_construction_and_modality_availability_400dpi.png")
  svg_file <- file.path(out_dir, "Figure_1_cohort_construction_and_modality_availability.svg")
  pdf_file <- file.path(out_dir, "Figure_1_cohort_construction_and_modality_availability.pdf")

  png_width_in <- 11.5
  png_height_in <- 8.4
  png_dpi <- 400
  png_width_px <- png_width_in * png_dpi
  png_height_px <- png_height_in * png_dpi

  png(png_file, width = png_width_in, height = png_height_in, units = "in", res = png_dpi, type = "cairo", family = FONT_FAMILY)
  draw_landscape()
  dev.off()

  svg(svg_file, width = 11.5, height = 8.4, family = FONT_FAMILY)
  draw_landscape()
  dev.off()

  cairo_pdf(pdf_file, width = 11.5, height = 8.4, family = FONT_FAMILY)
  draw_landscape()
  dev.off()

  manifest <- data.frame(
    figure_id = "Figure 1",
    panel_id = c("A", "B", "C", "D", "E"),
    script = "R/figures/Figure_1_cohort_construction.R",
    input_table = "outputs/figures/figure2/figure2_cohort_construction_panel_data.csv",
    output_file = c(png_file, svg_file, pdf_file, png_file, svg_file),
    width = c(png_width_px, 11.5, 11.5, png_width_px, 11.5),
    height = c(png_height_px, 8.4, 8.4, png_height_px, 8.4),
    dpi = c(png_dpi, NA, NA, png_dpi, NA),
    date = as.character(Sys.Date()),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    font_family = FONT_FAMILY,
    notes = c(
      "Cancer cohort composition",
      "Data modalities",
      "Sex composition",
      "Follow-up density",
      "Clinical outcomes"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(manifest, manifest_file, row.names = FALSE)
  message("Wrote ", png_file)
  message("Wrote ", svg_file)
  message("Wrote ", pdf_file)
}

main()
