#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The R package 'Matrix' is required.")
  }
})

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
}
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"
if (!dir.exists(project_root)) {
  stop(paste("Project root does not exist:", project_root))
}

now_iso <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_status <- function(status, payload) {
  path <- file.path(project_root, "results", "logs", "stage_11_gse272362_rds_status.json")
  safe_dir_create(dirname(path))
  base <- list(
    stage = "11_gse272362_rds",
    status = status,
    timestamp_utc = now_iso()
  )
  payload <- c(base, payload)
  json <- paste0(
    "{\n",
    paste(sprintf('  "%s": %s', names(payload), vapply(payload, to_json_scalar, character(1))), collapse = ",\n"),
    "\n}\n"
  )
  writeLines(json, path, useBytes = TRUE)
}

to_json_scalar <- function(x) {
  if (length(x) == 0 || is.null(x)) return("null")
  if (length(x) > 1 || is.list(x)) {
    vals <- vapply(x, to_json_scalar, character(1))
    return(paste0("[", paste(vals, collapse = ", "), "]"))
  }
  if (is.numeric(x) || is.integer(x)) {
    if (is.na(x)) "null" else as.character(x)
  } else if (is.logical(x)) {
    if (is.na(x)) "null" else tolower(as.character(x))
  } else {
    paste0('"', gsub('"', '\\"', as.character(x), fixed = TRUE), '"')
  }
}

parse_signatures <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  signatures <- list()
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || !grepl(": \\[", line, fixed = FALSE)) next
    key <- sub(":.*$", "", line)
    genes <- sub("^.*\\[", "", line)
    genes <- sub("\\].*$", "", genes)
    genes <- toupper(trimws(strsplit(genes, ",", fixed = TRUE)[[1]]))
    genes <- genes[nzchar(genes)]
    signatures[[key]] <- genes
  }
  signatures
}

zscore <- function(values) {
  values <- as.numeric(values)
  s <- stats::sd(values, na.rm = TRUE)
  m <- mean(values, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(values)))
  (values - m) / s
}

bh_adjust <- function(p) {
  stats::p.adjust(p, method = "BH")
}

is_s4_with_slot <- function(object, slot_name) {
  slot_name %in% names(attributes(object))
}

get_slot_safe <- function(object, slot_name) {
  if (is_s4_with_slot(object, slot_name)) {
    return(attributes(object)[[slot_name]])
  }
  NULL
}

assay_candidates <- function(object) {
  assays <- get_slot_safe(object, "assays")
  names(assays) %||% character()
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

extract_metadata <- function(object, barcodes) {
  meta <- get_slot_safe(object, "meta.data")
  if (is.null(meta)) {
    meta <- data.frame(row.names = barcodes)
  }
  meta <- as.data.frame(meta)
  meta$barcode <- rownames(meta)
  missing <- setdiff(barcodes, meta$barcode)
  if (length(missing) > 0) {
    extra <- data.frame(barcode = missing, row.names = missing)
    meta <- rbind(meta, extra)
  }
  meta <- meta[match(barcodes, meta$barcode), , drop = FALSE]
  meta
}

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NA_character_
}

infer_sample_id <- function(meta, fallback) {
  col <- first_existing_col(meta, c("sample_id", "sample", "Sample", "orig.ident", "library_id", "slice", "image", "patient", "Patient"))
  if (is.na(col)) {
    rep(fallback, nrow(meta))
  } else {
    values <- as.character(meta[[col]])
    values[is.na(values) | !nzchar(values)] <- fallback
    values
  }
}

infer_specimen_type <- function(meta) {
  col <- first_existing_col(meta, c("specimen_type", "Origin", "origin", "tissue_type", "tissue", "Tissue", "site", "Site", "organ", "Organ", "condition", "Condition", "group", "Group"))
  if (is.na(col)) {
    rep("metadata_required", nrow(meta))
  } else {
    values <- tolower(as.character(meta[[col]]))
    values[grepl("normal", values)] <- "normal_pancreas"
    values[grepl("liver|hepatic", values)] <- "liver_metastasis"
    values[grepl("lymph|node|ln", values)] <- "lymph_node_metastasis"
    values[grepl("primary|tumou?r|pdac|^pancreas$", values)] <- "primary_tumor"
    values[is.na(values) | !nzchar(values)] <- "metadata_required"
    values
  }
}

extract_coordinates <- function(object, barcodes, meta) {
  coords <- data.frame(barcode = barcodes, x_pixel = NA_real_, y_pixel = NA_real_)
  images <- get_slot_safe(object, "images")
  if (!is.null(images) && length(images) > 0) {
    pieces <- list()
    for (image_name in names(images)) {
      coordinates <- get_slot_safe(images[[image_name]], "coordinates")
      if (is.null(coordinates)) next
      coordinates <- as.data.frame(coordinates)
      coordinates$barcode <- rownames(coordinates)
      x_col <- first_existing_col(coordinates, c("imagecol", "x", "col", "pxl_col_in_fullres", "x_pixel"))
      y_col <- first_existing_col(coordinates, c("imagerow", "y", "row", "pxl_row_in_fullres", "y_pixel"))
      if (is.na(x_col) || is.na(y_col)) next
      pieces[[image_name]] <- data.frame(
        barcode = coordinates$barcode,
        x_pixel = as.numeric(coordinates[[x_col]]),
        y_pixel = as.numeric(coordinates[[y_col]]),
        image_id = image_name
      )
    }
    if (length(pieces) > 0) {
      all_coords <- do.call(rbind, pieces)
      all_coords <- all_coords[!duplicated(all_coords$barcode), , drop = FALSE]
      coords <- merge(coords["barcode"], all_coords, by = "barcode", all.x = TRUE, sort = FALSE)
    }
  }
  if (all(is.na(coords$x_pixel)) || all(is.na(coords$y_pixel))) {
    x_col <- first_existing_col(meta, c("x_pixel", "pxl_col_in_fullres", "imagecol", "x", "col"))
    y_col <- first_existing_col(meta, c("y_pixel", "pxl_row_in_fullres", "imagerow", "y", "row"))
    if (!is.na(x_col) && !is.na(y_col)) {
      coords$x_pixel <- as.numeric(meta[[x_col]])
      coords$y_pixel <- as.numeric(meta[[y_col]])
    }
  }
  coords
}

score_one_object <- function(object, object_name, signatures) {
  expr <- extract_expression_matrix(object)
  x <- expr$matrix
  barcodes <- colnames(x)
  genes <- rownames(x)
  if (is.null(barcodes) || is.null(genes)) stop("Expression matrix needs rownames and colnames.")

  meta <- extract_metadata(object, barcodes)
  sample_id <- infer_sample_id(meta, object_name)
  specimen_type <- infer_specimen_type(meta)
  coords <- extract_coordinates(object, barcodes, meta)

  n_counts <- as.numeric(Matrix::colSums(x))
  n_genes <- as.numeric(Matrix::colSums(x > 0))
  scale <- rep(0, length(n_counts))
  scale[n_counts > 0] <- 10000 / n_counts[n_counts > 0]
  gene_upper <- toupper(genes)

  spot <- data.frame(
    dataset_id = "GSE272362",
    sample_id = sample_id,
    patient_id = if ("patient_id" %in% names(meta)) {
      as.character(meta$patient_id)
    } else if ("patient" %in% names(meta)) {
      as.character(meta$patient)
    } else {
      paste0("patient_unknown_", sample_id)
    },
    specimen_type = specimen_type,
    barcode = barcodes,
    x_pixel = coords$x_pixel,
    y_pixel = coords$y_pixel,
    n_counts = n_counts,
    n_genes = n_genes,
    stringsAsFactors = FALSE
  )

  coverage <- list()
  for (sig_name in names(signatures)) {
    sig_genes <- signatures[[sig_name]]
    idx <- which(gene_upper %in% sig_genes)
    present <- sort(unique(gene_upper[idx]))
    missing <- sort(setdiff(sig_genes, present))
    reliable <- length(present) >= 3
    if (reliable) {
      sub <- x[idx, , drop = FALSE]
      sub <- Matrix::t(Matrix::t(sub) * scale)
      values <- as.numeric(Matrix::colMeans(log1p(sub)))
    } else {
      values <- rep(NA_real_, ncol(x))
    }
    spot[[paste0("score_", sig_name)]] <- values
    spot[[paste0("z_", sig_name)]] <- ave(values, sample_id, FUN = zscore)
    coverage[[sig_name]] <- data.frame(
      dataset_id = "GSE272362",
      sample_id = object_name,
      signature = sig_name,
      n_genes_defined = length(sig_genes),
      n_genes_present = length(present),
      n_genes_missing = length(missing),
      present_genes = paste(present, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      reliable = reliable,
      stringsAsFactors = FALSE
    )
  }

  composites <- list(
    immune_hub_core = c("b_cell", "t_cell", "dc_apc", "tls_chemokine"),
    immune_hub_maturity = c("fdc_gc_like", "plasma_cell", "ifn_antigen_presentation"),
    caf_myeloid_barrier = c("pan_caf", "myeloid", "spp1_tam", "tgfb_pathway", "icaf"),
    tumor_aggressive = c("pdac_basal_like", "emt_invasion", "hypoxia", "proliferation")
  )
  for (comp_name in names(composites)) {
    cols <- paste0("z_", composites[[comp_name]])
    cols <- cols[cols %in% names(spot)]
    spot[[paste0("score_", comp_name)]] <- if (length(cols)) rowMeans(spot[, cols, drop = FALSE], na.rm = TRUE) else NA_real_
  }

  list(
    spot = spot,
    coverage = do.call(rbind, coverage),
    assay = expr$assay,
    layer = expr$layer,
    n_spots = ncol(x),
    n_genes = nrow(x),
    metadata_columns = names(meta)
  )
}

summarize_samples <- function(spot) {
  score_cols <- c(
    "score_immune_hub_core",
    "score_immune_hub_maturity",
    "score_caf_myeloid_barrier",
    "score_tumor_aggressive",
    "z_ifn_antigen_presentation"
  )
  rows <- list()
  groups <- split(spot, spot$sample_id, drop = TRUE)
  for (sample_name in names(groups)) {
    df <- groups[[sample_name]]
    row <- data.frame(
      dataset_id = "GSE272362",
      sample_id = sample_name,
      specimen_type = names(sort(table(df$specimen_type), decreasing = TRUE))[1],
      n_spots_qc = nrow(df),
      median_counts = median(df$n_counts, na.rm = TRUE),
      median_genes = median(df$n_genes, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (col in score_cols) {
      if (col %in% names(df)) {
        prefix <- sub("^score_", "", col)
        prefix <- sub("^z_", "", prefix)
        row[[paste0("mean_", prefix)]] <- mean(df[[col]], na.rm = TRUE)
        row[[paste0("p90_", prefix)]] <- as.numeric(stats::quantile(df[[col]], 0.9, na.rm = TRUE))
        row[[paste0("fraction_", prefix, "_z_gt1")]] <- mean(df[[col]] > 1, na.rm = TRUE)
      }
    }
    rows[[sample_name]] <- row
  }
  do.call(rbind, rows)
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

main <- function() {
  rds_path <- file.path(project_root, "data", "external", "GSE272362_zenodo", "PDAC_Updated.rds")
  signature_path <- file.path(project_root, "config", "signatures.yaml")
  out_dir <- file.path(project_root, "results", "tables")
  meta_dir <- file.path(project_root, "metadata")
  safe_dir_create(out_dir)
  safe_dir_create(meta_dir)

  if (!file.exists(rds_path)) stop(paste("Missing RDS:", rds_path))
  signatures <- parse_signatures(signature_path)

  message("Reading RDS: ", rds_path)
  object <- readRDS(rds_path)
  objects <- flatten_objects(object)
  message("Scoring ", length(objects), " object(s)")

  results <- list()
  for (nm in names(objects)) {
    message("Scoring object: ", nm)
    results[[nm]] <- score_one_object(objects[[nm]], nm, signatures)
  }

  spot <- do.call(rbind, lapply(results, `[[`, "spot"))
  coverage <- do.call(rbind, lapply(results, `[[`, "coverage"))
  sample <- summarize_samples(spot)
  metadata_columns <- unique(unlist(lapply(results, `[[`, "metadata_columns")))
  metadata_audit <- data.frame(column = metadata_columns, stringsAsFactors = FALSE)
  object_audit <- data.frame(
    object_name = names(results),
    assay = vapply(results, `[[`, character(1), "assay"),
    layer = vapply(results, `[[`, character(1), "layer"),
    n_genes = vapply(results, `[[`, numeric(1), "n_genes"),
    n_spots = vapply(results, `[[`, numeric(1), "n_spots"),
    stringsAsFactors = FALSE
  )

  write.csv(spot, file.path(out_dir, "gse272362_rds_spot_level_scores.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(sample, file.path(out_dir, "gse272362_rds_sample_level_scores.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(coverage, file.path(out_dir, "gse272362_rds_signature_gene_coverage.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(metadata_audit, file.path(meta_dir, "gse272362_rds_metadata_columns.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(object_audit, file.path(meta_dir, "gse272362_rds_object_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  write_status("success", list(
    n_objects = length(results),
    n_spots = nrow(spot),
    n_samples = length(unique(spot$sample_id)),
    output_spot_scores = file.path(out_dir, "gse272362_rds_spot_level_scores.csv"),
    output_sample_scores = file.path(out_dir, "gse272362_rds_sample_level_scores.csv")
  ))
  message("Done. Wrote GSE272362 RDS scores.")
}

tryCatch(
  main(),
  error = function(e) {
    write_status("failed", list(error = conditionMessage(e)))
    stop(e)
  }
)
