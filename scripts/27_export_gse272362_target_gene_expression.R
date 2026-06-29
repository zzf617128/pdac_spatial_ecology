#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required")
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"

target_genes <- unique(c(
  "SPP1", "TREM2", "APOE", "LGALS3", "COL1A1", "COL1A2", "FN1", "SPARC",
  "TGFB1", "TGFBI", "CTGF", "INHBA", "VIM", "ITGA5", "SERPINE1", "MMP14",
  "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "CXCL9", "CXCL10", "STAT1", "IRF1",
  "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1", "CXCL13", "BANK1",
  "CD3D", "CD3E", "CD8A", "TRAC", "PDCD1", "CTLA4", "LAG3", "TIGIT",
  "CD68", "CXCL11", "C1QA", "C1QB", "C1QC", "MARCO", "MRC1", "LST1", "FCGR3A", "ITGAM",
  "CD44", "ITGAV", "ITGB1", "ITGB5", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3",
  "COL3A1", "COL6A1", "COL6A2", "COL6A3", "POSTN",
  "IL6", "IL6R", "IL6ST", "LIF", "LIFR", "OSM", "OSMR", "JAK1", "STAT3", "SOCS3",
  "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "PDGFRA", "PDGFRB", "CXCL12", "CXCL14", "CFD", "DPT", "HAS1",
  "LAMP3", "CCR7", "CLEC10A", "FCER1A",
  "CD4", "IL7R", "NKG7", "GZMB",
  "EPCAM", "KRT8", "KRT18", "KRT19", "KRT17", "MSLN", "CEACAM6",
  "PECAM1", "VWF", "EMCN", "KDR",
  "SOX10", "S100B", "PLP1", "MPZ", "NRXN1",
  "PRSS1", "CPA1", "REG1A"
))

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
}

is_s4_with_slot <- function(object, slot_name) {
  slot_name %in% names(attributes(object))
}

get_slot_safe <- function(object, slot_name) {
  if (is_s4_with_slot(object, slot_name)) attributes(object)[[slot_name]] else NULL
}

extract_matrix_from_assay <- function(assay) {
  layers <- get_slot_safe(assay, "layers")
  if (!is.null(layers)) {
    for (nm in c("counts", "data")) {
      if (nm %in% names(layers)) return(list(matrix = layers[[nm]], layer = nm))
    }
  }
  for (slot_name in c("counts", "data")) {
    value <- get_slot_safe(assay, slot_name)
    if (!is.null(value) && length(dim(value)) == 2) return(list(matrix = value, layer = slot_name))
  }
  NULL
}

extract_expression_matrix <- function(object) {
  assays <- get_slot_safe(object, "assays")
  if (is.null(assays)) stop("Could not find @assays in the RDS object.")
  preferred <- intersect(c("Spatial", "RNA", "SCT"), names(assays))
  candidates <- c(preferred, setdiff(names(assays), preferred))
  for (assay_name in candidates) {
    got <- extract_matrix_from_assay(assays[[assay_name]])
    if (!is.null(got)) {
      matrix <- got$matrix
      if (!inherits(matrix, "sparseMatrix")) matrix <- Matrix::Matrix(matrix, sparse = TRUE)
      return(list(matrix = matrix, assay = assay_name, layer = got$layer))
    }
  }
  stop("Could not find a counts/data matrix in object assays.")
}

flatten_objects <- function(object) {
  if (is.list(object) && !is_s4_with_slot(object, "assays")) {
    hits <- list()
    object_names <- names(object)
    if (is.null(object_names)) object_names <- paste0("object_", seq_along(object))
    for (i in seq_along(object)) {
      nm <- object_names[[i]]
      item <- object[[i]]
      if (is_s4_with_slot(item, "assays")) hits[[nm]] <- item
    }
    if (length(hits)) return(hits)
  }
  list(PDAC_Updated = object)
}

extract_gene_frame <- function(object) {
  expr <- extract_expression_matrix(object)
  x <- expr$matrix
  genes <- rownames(x)
  barcodes <- colnames(x)
  if (is.null(genes) || is.null(barcodes)) stop("Expression matrix lacks rownames or colnames.")
  gene_upper <- toupper(genes)
  n_counts <- as.numeric(Matrix::colSums(x))
  scale <- rep(0, length(n_counts))
  scale[n_counts > 0] <- 10000 / n_counts[n_counts > 0]
  out <- data.frame(barcode = barcodes, stringsAsFactors = FALSE)
  for (gene in target_genes) {
    idx <- which(gene_upper == gene)
    if (length(idx) == 0) {
      out[[gene]] <- NA_real_
    } else {
      values <- as.numeric(Matrix::colSums(x[idx, , drop = FALSE]))
      out[[gene]] <- log1p(values * scale)
    }
  }
  out
}

main <- function() {
  rds_path <- file.path(project_root, "data", "external", "GSE272362_zenodo", "PDAC_Updated.rds")
  score_path <- file.path(project_root, "results", "tables", "gse272362_rds_spot_level_scores.csv")
  out_path <- file.path(project_root, "results", "tables", "gse272362_rds_target_gene_expression.csv")
  status_path <- file.path(project_root, "results", "logs", "stage_27_gse272362_target_gene_export_status.json")
  if (!file.exists(rds_path)) stop(paste("Missing RDS:", rds_path))
  if (!file.exists(score_path)) stop(paste("Missing score table:", score_path))

  message("Reading RDS: ", rds_path)
  object <- readRDS(rds_path)
  objects <- flatten_objects(object)
  frames <- lapply(objects, extract_gene_frame)
  gene_df <- do.call(rbind, frames)
  scores <- read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- c(
    "dataset_id", "sample_id", "specimen_type", "barcode", "x_pixel", "y_pixel",
    "score_caf_myeloid_barrier", "score_tumor_epithelial"
  )
  scores <- scores[, intersect(keep, names(scores)), drop = FALSE]
  merged <- merge(scores, gene_df, by = "barcode", all.x = FALSE, all.y = FALSE, sort = FALSE)
  # Preserve the Python table convention: barcode is not the first column.
  merged <- merged[, c(setdiff(names(merged), target_genes), target_genes), drop = FALSE]
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(merged, out_path, row.names = FALSE, fileEncoding = "UTF-8")
  dir.create(dirname(status_path), recursive = TRUE, showWarnings = FALSE)
  status <- sprintf(
    '{\n  "stage": "27_gse272362_target_gene_export",\n  "status": "success",\n  "n_spots": %d,\n  "n_genes_requested": %d,\n  "output": "%s"\n}\n',
    nrow(merged), length(target_genes), gsub("\\\\", "/", out_path)
  )
  writeLines(status, status_path, useBytes = TRUE)
  message("Wrote ", nrow(merged), " spots to ", out_path)
}

tryCatch(main(), error = function(e) {
  message("FAILED: ", conditionMessage(e))
  quit(status = 1)
})
