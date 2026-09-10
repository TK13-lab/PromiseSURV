main <- function() {
  # Purpose: Run reproducible enrichment/GSEA on the Figure 5/6 retained
  # CPTAC-derived TCGA RNA gene pool.
  # Inputs: outputs/task6/shared_markers/task6_shared_marker_overlap_membership.csv
  # Outputs: outputs/task6/shared_markers/task6_shared_marker_enrichment_*.csv
  # Usage: Rscript R/analysis/shared_marker_enrichment.R [--input-root=PATH]

  command_args <- commandArgs(trailingOnly = TRUE)
  argument_value <- function(prefix, default) {
    hit <- command_args[startsWith(command_args, prefix)]
    if (length(hit) == 0) default else sub(prefix, "", hit[[1]], fixed = TRUE)
  }
  args_file <- commandArgs(FALSE)
  file_arg <- args_file[grepl("^--file=", args_file)]
  script_path <- if (length(file_arg) > 0) {
    sub("^--file=", "", file_arg[[1]])
  } else {
    "R/analysis/shared_marker_enrichment.R"
  }
  package_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  project_dir <- normalizePath(
    argument_value("--input-root=", package_dir),
    winslash = "/", mustWork = TRUE
  )

  suppressPackageStartupMessages({
    library(msigdbr)
    library(fgsea)
    library(BiocParallel)
    library(clusterProfiler)
  })

  shared_dir <- file.path(project_dir, "outputs/task6/shared_markers")
  overlap_file <- file.path(shared_dir, "task6_shared_marker_overlap_membership.csv")
  if (!file.exists(overlap_file)) {
    stop("Missing overlap membership table: ", overlap_file)
  }

  overlap <- read.csv(overlap_file, stringsAsFactors = FALSE)
  as_bool <- function(x) {
    if (is.logical(x)) {
      return(x)
    }
    tolower(as.character(x)) %in% c("true", "1", "yes")
  }

  required_cols <- c(
    "gene_symbol",
    "tested_in_all_three",
    "significant_count",
    "all_three_significant",
    "same_direction_all_three",
    "promise_enriched_all_cancers",
    "HNSC_rho",
    "HNSC_fdr",
    "KIRC_rho",
    "KIRC_fdr",
    "LUAD_rho",
    "LUAD_fdr"
  )
  missing_cols <- setdiff(required_cols, names(overlap))
  if (length(missing_cols) > 0) {
    stop("Missing columns in overlap table: ", paste(missing_cols, collapse = ", "))
  }

  for (bool_col in c("tested_in_all_three", "all_three_significant", "same_direction_all_three", "promise_enriched_all_cancers")) {
    overlap[[bool_col]] <- as_bool(overlap[[bool_col]])
  }

  source <- overlap[overlap$tested_in_all_three, , drop = FALSE]
  universe <- unique(source$gene_symbol)
  target_sets <- list(
    risk_associated_any_cancer = unique(source$gene_symbol[source$significant_count > 0]),
    all_three_risk_associated = unique(source$gene_symbol[source$all_three_significant]),
    same_direction_core = unique(source$gene_symbol[source$same_direction_all_three]),
    promise_enriched_core = unique(source$gene_symbol[source$promise_enriched_all_cancers])
  )

  score_one <- function(rho, fdr) {
    fdr <- pmax(as.numeric(fdr), 1e-50)
    sign(as.numeric(rho)) * pmin(-log10(fdr), 20)
  }
  score_mat <- cbind(
    HNSC = score_one(source$HNSC_rho, source$HNSC_fdr),
    KIRC = score_one(source$KIRC_rho, source$KIRC_fdr),
    LUAD = score_one(source$LUAD_rho, source$LUAD_fdr)
  )
  source$signed_gsea_score <- rowSums(score_mat, na.rm = TRUE)
  source$mean_signed_gsea_score <- rowMeans(score_mat, na.rm = TRUE)
  source$absolute_gsea_evidence <- rowSums(abs(score_mat), na.rm = TRUE)
  source$direction_consistency <- abs(rowSums(sign(score_mat), na.rm = TRUE))

  ranked_scores <- source[, c(
    "gene_symbol",
    "signed_gsea_score",
    "mean_signed_gsea_score",
    "absolute_gsea_evidence",
    "direction_consistency",
    "significant_count",
    "all_three_significant",
    "same_direction_all_three",
    "promise_enriched_all_cancers"
  )]
  ranked_scores <- ranked_scores[order(-ranked_scores$signed_gsea_score, ranked_scores$gene_symbol), ]

  clean_term <- function(x) {
    out <- gsub("^HALLMARK_", "", x)
    out <- gsub("^REACTOME_", "", out)
    out <- gsub("^GOBP_", "", out)
    out <- gsub("_", " ", out)
    out <- tools::toTitleCase(tolower(out))
    out <- gsub("\\bRna\\b", "RNA", out)
    out <- gsub("\\bTrna\\b", "tRNA", out)
    out <- gsub("\\bRrna\\b", "rRNA", out)
    out <- gsub("\\bDna\\b", "DNA", out)
    out <- gsub("\\bE2f\\b", "E2F", out)
    out <- gsub("\\bMyc\\b", "MYC", out)
    out <- gsub("\\bMtorc1\\b", "mTORC1", out)
    out
  }

  load_msig <- function(collection, subcollection = NULL, label) {
    dat <- if (is.null(subcollection)) {
      msigdbr(species = "Homo sapiens", collection = collection)
    } else {
      msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    }
    dat <- unique(dat[, c("gs_name", "gene_symbol")])
    names(dat) <- c("term_id", "gene_symbol")
    dat$collection_label <- label
    dat
  }

  term2gene_all <- rbind(
    load_msig("H", label = "Hallmark"),
    load_msig("C2", "CP:REACTOME", "Reactome"),
    load_msig("C5", "GO:BP", "GO biological process")
  )
  term2gene_all <- unique(term2gene_all[term2gene_all$gene_symbol %in% universe, ])

  pathway_list <- split(term2gene_all$gene_symbol, term2gene_all$term_id)
  pathway_list <- lapply(pathway_list, unique)
  pathway_list <- pathway_list[lengths(pathway_list) >= 3]
  term_collection <- unique(term2gene_all[, c("term_id", "collection_label")])

  stats <- ranked_scores$signed_gsea_score
  names(stats) <- ranked_scores$gene_symbol
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  stats <- stats[!duplicated(names(stats))]

  set.seed(20260225)
  gsea <- fgsea(
    pathways = pathway_list,
    stats = stats,
    minSize = 3,
    maxSize = 200,
    eps = 0,
    BPPARAM = SerialParam()
  )
  gsea <- as.data.frame(gsea)
  gsea$leading_edge_genes <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  gsea$leadingEdge <- NULL
  gsea <- merge(gsea, term_collection, by.x = "pathway", by.y = "term_id", all.x = TRUE)
  gsea$term_label <- clean_term(gsea$pathway)
  gsea <- gsea[order(gsea$padj, -abs(gsea$NES), gsea$pathway), ]

  run_ora <- function(gene_set_name, genes, collection_label) {
    term2gene <- term2gene_all[term2gene_all$collection_label == collection_label, c("term_id", "gene_symbol")]
    if (length(intersect(genes, term2gene$gene_symbol)) < 2) {
      return(data.frame())
    }
    enr <- enricher(
      gene = genes,
      universe = universe,
      TERM2GENE = term2gene,
      pAdjustMethod = "BH",
      minGSSize = 2,
      maxGSSize = 200,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
    if (is.null(enr)) {
      return(data.frame())
    }
    out <- as.data.frame(enr)
    if (nrow(out) == 0) {
      return(out)
    }
    out$gene_set_name <- gene_set_name
    out$collection_label <- collection_label
    out$term_label <- clean_term(out$ID)
    out
  }

  ora <- do.call(
    rbind,
    unlist(
      lapply(names(target_sets), function(set_name) {
        lapply(unique(term2gene_all$collection_label), function(collection_label) {
          run_ora(set_name, target_sets[[set_name]], collection_label)
        })
      }),
      recursive = FALSE
    )
  )
  if (is.null(ora) || nrow(ora) == 0) {
    ora <- data.frame()
  } else {
    ora <- ora[order(ora$p.adjust, ora$pvalue, ora$ID), ]
  }

  preferred_terms <- c(
    "REACTOME_METABOLISM_OF_RNA",
    "GOBP_RIBOSOME_BIOGENESIS",
    "GOBP_RRNA_PROCESSING",
    "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
    "REACTOME_TRNA_PROCESSING",
    "REACTOME_TRNA_PROCESSING_IN_THE_NUCLEUS",
    "GOBP_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS",
    "HALLMARK_E2F_TARGETS"
  )
  panel_terms <- gsea[gsea$pathway %in% preferred_terms, , drop = FALSE]
  panel_terms$panel_rank <- match(panel_terms$pathway, preferred_terms)
  panel_terms <- panel_terms[order(panel_terms$panel_rank), ]
  panel_terms$neglog10_fdr <- -log10(pmax(panel_terms$padj, 1e-50))
  panel_terms$direction <- ifelse(panel_terms$NES >= 0, "Higher-risk edge", "Lower-risk edge")
  panel_terms$rank_metric <- "sum(sign(rho) * -log10(FDR)) across HNSC/KIRC/LUAD"

  make_curve <- function(term_row) {
    genes <- pathway_list[[term_row$pathway]]
    hit <- names(stats) %in% genes
    hit_weight <- abs(stats)
    hit_norm <- sum(hit_weight[hit])
    miss_norm <- sum(!hit)
    running <- cumsum(ifelse(hit, hit_weight / hit_norm, -1 / miss_norm))
    data.frame(
      pathway = term_row$pathway,
      term_label = term_row$term_label,
      collection_label = term_row$collection_label,
      panel_rank = term_row$panel_rank,
      NES = term_row$NES,
      padj = term_row$padj,
      rank = seq_along(stats),
      gene_symbol = names(stats),
      rank_score = as.numeric(stats),
      hit = hit,
      running_es = running,
      stringsAsFactors = FALSE
    )
  }
  curve_data <- do.call(rbind, lapply(seq_len(nrow(panel_terms)), function(i) make_curve(panel_terms[i, ])))

  ranked_scores_file <- file.path(shared_dir, "task6_shared_marker_enrichment_ranked_gene_scores.csv")
  gsea_file <- file.path(shared_dir, "task6_shared_marker_enrichment_gsea.csv")
  ora_file <- file.path(shared_dir, "task6_shared_marker_enrichment_ora.csv")
  panel_file <- file.path(shared_dir, "task6_shared_marker_enrichment_panel_b.csv")
  curve_file <- file.path(shared_dir, "task6_shared_marker_enrichment_curve_data.csv")

  write.csv(ranked_scores, ranked_scores_file, row.names = FALSE)
  write.csv(gsea, gsea_file, row.names = FALSE)
  write.csv(ora, ora_file, row.names = FALSE)
  write.csv(panel_terms, panel_file, row.names = FALSE)
  write.csv(curve_data, curve_file, row.names = FALSE)

  message("Wrote ranked scores to: ", ranked_scores_file)
  message("Wrote GSEA table to: ", gsea_file)
  message("Wrote ORA table to: ", ora_file)
  message("Wrote Figure 6 panel B enrichment table to: ", panel_file)
  message("Wrote Figure 6 GSEA curve data to: ", curve_file)
}

main()
