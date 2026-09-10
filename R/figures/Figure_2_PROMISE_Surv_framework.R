# Purpose: Reproduce manuscript Figure 2 (PROMISE-Surv workflow).


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  argument_value <- function(prefix, default) {
    hit <- args[startsWith(args, prefix)]
    if (length(hit) == 0) default else sub(prefix, "", hit[[1]], fixed = TRUE)
  }
  script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][[1]])
  package_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
  output_dir <- argument_value("--output-dir=", file.path(package_root, "outputs", "figures"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
  })

  nodes <- data.frame(
    id = c("disc", "teacher", "risk", "student", "deploy", "surv", "bio", "marker"),
    label = c(
      "A Discovery inputs\nClinical + RNA\nProteomics + phosphoproteomics\nCPTAC/PDC complete cases",
      "Protein-informed teacher\nLearns survival-risk structure\nfrom complete modalities",
      "B Risk-transfer bridge\nDistill teacher risk\nnot missing assays",
      "Deployable student\nClinical + RNA only",
      "C Partial-modality cohorts\nTCGA + public GEO\nproteomics absent",
      "External survival\nstress test\nOS + DFS/progression",
      "D Biological transfer\nRNA programs\nRPPA modules",
      "Shared marker discovery\nconserved RNA-risk\ncandidates"
    ),
    x = c(0.12, 0.35, 0.55, 0.73, 0.90, 0.62, 0.86, 0.74),
    y = c(0.67, 0.67, 0.67, 0.67, 0.67, 0.30, 0.30, 0.10),
    width = c(0.22, 0.22, 0.17, 0.17, 0.18, 0.22, 0.22, 0.30),
    height = c(0.23, 0.20, 0.18, 0.18, 0.18, 0.17, 0.17, 0.15),
    fill_group = c("discovery", "teacher", "bridge", "student", "deploy", "output", "output", "output"),
    stringsAsFactors = FALSE
  )
  nodes$xmin <- nodes$x - nodes$width / 2
  nodes$xmax <- nodes$x + nodes$width / 2
  nodes$ymin <- nodes$y - nodes$height / 2
  nodes$ymax <- nodes$y + nodes$height / 2

  edges <- data.frame(
    x = c(0.23, 0.46, 0.635, 0.815, 0.86, 0.92, 0.74),
    y = c(0.67, 0.67, 0.67, 0.67, 0.58, 0.58, 0.215),
    xend = c(0.24, 0.465, 0.645, 0.81, 0.66, 0.86, 0.74),
    yend = c(0.67, 0.67, 0.67, 0.67, 0.385, 0.385, 0.175),
    stringsAsFactors = FALSE
  )

  figure2 <- ggplot() +
    geom_segment(
      data = edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.55,
      color = "#4A4A4A",
      arrow = arrow(length = unit(0.11, "inches"), type = "closed")
    ) +
    geom_rect(
      data = nodes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_group),
      color = "#555555",
      linewidth = 0.40
    ) +
    geom_text(
      data = nodes,
      aes(x = x, y = y, label = label),
      family = "Arial",
      size = 2.85,
      lineheight = 0.90
    ) +
    annotate(
      "text",
      x = 0.12, y = 0.93,
      label = "PROMISE-Surv transfers protein-informed risk into deployable RNA+clinical cohorts",
      family = "Arial", fontface = "bold", size = 4.30, hjust = 0
    ) +
    annotate(
      "text",
      x = 0.12, y = 0.865,
      label = "Framework claim: transfer and interrogation of survival biology under missing-proteomics constraints",
      family = "Arial", size = 3.30, hjust = 0, color = "#333333"
    ) +
    scale_fill_manual(values = c(
      discovery = "#F6DEC7", teacher = "#F1CDD2", bridge = "#EFE7DD",
      student = "#D7ECE8", deploy = "#E5E9F1", output = "#F4E5BC"
    )) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off", expand = FALSE) +
    theme_void(base_family = "Arial") +
    theme(legend.position = "none", plot.margin = margin(12, 12, 12, 12))

  png_file <- file.path(output_dir, "Figure_2_PROMISE_SURV_framework_400dpi.png")
  svg_file <- file.path(output_dir, "Figure_2_PROMISE_SURV_framework.svg")
  pdf_file <- file.path(output_dir, "Figure_2_PROMISE_SURV_framework.pdf")
  ggsave(png_file, figure2,
    width = 11.2, height = 6.2, units = "in", dpi = 400,
    bg = "white", device = ragg::agg_png, limitsize = FALSE
  )
  ggsave(svg_file, figure2,
    width = 11.2, height = 6.2, units = "in",
    bg = "white", device = svglite::svglite, limitsize = FALSE
  )
  ggsave(pdf_file, figure2,
    width = 11.2, height = 6.2, units = "in",
    bg = "white", device = grDevices::cairo_pdf, limitsize = FALSE
  )
  message("Wrote ", png_file)
}

main()
