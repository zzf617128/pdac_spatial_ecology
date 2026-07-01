suppressPackageStartupMessages({
  library(data.table)
})

root <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search/cho_imc"), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "data")
table_dir <- file.path(root, "tables")
doc_dir <- file.path(root, "docs")
log_dir <- file.path(root, "logs")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

rds_path <- file.path(data_dir, "backup_output.rds")
log_path <- file.path(log_dir, "01_inspect_cho_imc_object.log")
sink(log_path, split = TRUE)
cat("Reading ", rds_path, "\n", sep = "")
obj <- readRDS(rds_path)
cat("Object class:", paste(class(obj), collapse = "/"), "\n")
cat("Object size:", format(object.size(obj), units = "auto"), "\n")

capture <- function(expr) paste(capture.output(print(expr)), collapse = "\n")

object_summary <- data.table(
  field = c("object_class", "object_size", "names", "slot_names"),
  value = c(
    paste(class(obj), collapse = "/"),
    format(object.size(obj), units = "auto"),
    if (!is.null(names(obj))) paste(head(names(obj), 100), collapse = ";") else NA_character_,
    if (isS4(obj)) paste(slotNames(obj), collapse = ";") else NA_character_
  )
)

metadata <- NULL
expr <- NULL
coords <- NULL
assay_notes <- character()

if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
  metadata <- as.data.table(as.data.frame(SummarizedExperiment::colData(obj)), keep.rownames = "cell_id")
  assay_names <- SummarizedExperiment::assayNames(obj)
  assay_notes <- c(assay_notes, paste0("assayNames=", paste(assay_names, collapse = ";")))
  if (length(assay_names) > 0) {
    expr <- SummarizedExperiment::assay(obj, assay_names[1])
  }
  if ("spatialCoords" %in% getNamespaceExports("SpatialExperiment") && inherits(obj, "SpatialExperiment")) {
    coords <- tryCatch(as.data.table(SpatialExperiment::spatialCoords(obj), keep.rownames = "cell_id"), error = function(e) NULL)
  }
} else if (inherits(obj, "Seurat")) {
  metadata <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  assay_notes <- c(assay_notes, paste0("assays=", paste(names(obj@assays), collapse = ";")))
  default_assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) names(obj@assays)[1])
  expr <- tryCatch(SeuratObject::GetAssayData(obj, assay = default_assay, slot = "data"), error = function(e) NULL)
  if (is.null(expr)) expr <- tryCatch(SeuratObject::GetAssayData(obj, assay = default_assay, slot = "counts"), error = function(e) NULL)
  if (length(obj@reductions) > 0) {
    red_names <- names(obj@reductions)
    assay_notes <- c(assay_notes, paste0("reductions=", paste(red_names, collapse = ";")))
    for (rn in red_names) {
      emb <- tryCatch(SeuratObject::Embeddings(obj, rn), error = function(e) NULL)
      if (!is.null(emb) && ncol(emb) >= 2) {
        coords <- as.data.table(emb[, 1:2, drop = FALSE], keep.rownames = "cell_id")
        setnames(coords, names(coords)[2:3], c("x", "y"))
        break
      }
    }
  }
} else if (is.list(obj)) {
  nm <- names(obj)
  cat("List names:\n")
  print(nm)
  possible_meta <- nm[grepl("meta|cell|phen|anno", nm, ignore.case = TRUE)]
  for (n in possible_meta) {
    x <- obj[[n]]
    if (is.data.frame(x) || is.data.table(x)) {
      metadata <- as.data.table(x, keep.rownames = "cell_id")
      assay_notes <- c(assay_notes, paste0("metadata_from=", n))
      break
    }
  }
  possible_expr <- nm[grepl("expr|assay|count|marker|data|mat", nm, ignore.case = TRUE)]
  for (n in possible_expr) {
    x <- obj[[n]]
    if (is.matrix(x) || inherits(x, "Matrix")) {
      expr <- x
      assay_notes <- c(assay_notes, paste0("matrix_from=", n))
      break
    }
  }
}

if (is.null(metadata)) {
  cat("No metadata table detected by generic rules.\n")
  metadata_summary <- data.table(field = "metadata_detected", value = "FALSE")
  marker_availability <- data.table(marker_query = character(), matched_marker = character(), available = logical())
  sample_roi_summary <- data.table()
} else {
  cat("Metadata dimensions:", nrow(metadata), ncol(metadata), "\n")
  cat("Metadata columns:\n")
  print(names(metadata))
  fwrite(metadata[, .SD[1:min(.N, 5)]], file.path(table_dir, "cho_imc_metadata_head.csv"))
  metadata_summary <- rbindlist(list(
    object_summary,
    data.table(field = "n_cells_metadata", value = as.character(nrow(metadata))),
    data.table(field = "n_metadata_columns", value = as.character(ncol(metadata))),
    data.table(field = "metadata_columns", value = paste(names(metadata), collapse = ";")),
    data.table(field = "assay_notes", value = paste(assay_notes, collapse = " | "))
  ), fill = TRUE)
}

if (!is.null(expr)) {
  cat("Expression class:", paste(class(expr), collapse = "/"), "\n")
  cat("Expression dim:", paste(dim(expr), collapse = " x "), "\n")
  marker_names <- rownames(expr)
  if (is.null(marker_names) && nrow(expr) < ncol(expr)) marker_names <- rownames(expr)
  if (is.null(marker_names)) marker_names <- character()
  metadata_summary <- rbind(metadata_summary, data.table(field = c("expression_class", "expression_dim", "n_markers", "marker_names"), value = c(paste(class(expr), collapse = "/"), paste(dim(expr), collapse = " x "), as.character(length(marker_names)), paste(marker_names, collapse = ";"))), fill = TRUE)
} else {
  marker_names <- character()
  metadata_summary <- rbind(metadata_summary, data.table(field = "expression_detected", value = "FALSE"), fill = TRUE)
}

queries <- list(
  SMA = c("SMA", "aSMA", "alphaSMA", "ACTA2", "SMA1"),
  Collagen = c("Collagen", "COL1A1", "COL1", "COL3", "COLLAGEN"),
  Vimentin = c("Vimentin", "VIM"),
  FAP = c("FAP"),
  PDPN = c("PDPN", "Podoplanin"),
  Fibronectin = c("Fibronectin", "FN1"),
  CD68 = c("CD68"),
  CD163 = c("CD163"),
  CD11b = c("CD11b", "ITGAM"),
  CD14 = c("CD14"),
  HLA_DR = c("HLA-DR", "HLADR", "HLA_DR", "HLA.DR", "HLA-DRA"),
  CD45 = c("CD45", "PTPRC"),
  CD206 = c("CD206", "MRC1"),
  panCK = c("panCK", "PanCK", "Pan-Cytokeratin", "Cytokeratin", "CK"),
  EPCAM = c("EPCAM", "EpCAM"),
  CD3 = c("CD3"),
  CD4 = c("CD4"),
  CD8 = c("CD8"),
  CD20 = c("CD20", "MS4A1"),
  FOXP3 = c("FOXP3", "FoxP3"),
  CD44 = c("CD44"),
  Integrin = c("Integrin", "ITGA", "ITGB"),
  SPP1 = c("SPP1", "Osteopontin")
)

match_marker <- function(qs, markers) {
  if (length(markers) == 0) return(character())
  hits <- unique(unlist(lapply(qs, function(q) markers[grepl(q, markers, ignore.case = TRUE, fixed = FALSE)])))
  hits
}

marker_availability <- rbindlist(lapply(names(queries), function(q) {
  hits <- match_marker(queries[[q]], marker_names)
  data.table(marker_query = q, matched_marker = if (length(hits)) paste(hits, collapse = ";") else NA_character_, available = length(hits) > 0)
}))

if (!is.null(metadata)) {
  col_lower <- tolower(names(metadata))
  candidate_fields <- list(
    patient_id = c("patient", "case", "donor"),
    sample_id = c("sample", "specimen"),
    roi_id = c("roi", "fov", "image", "region", "core"),
    site_status = c("primary", "metast", "tissue", "site", "organ"),
    x = c("^x$", "centroid.*x", "cell.*x", "xcoord", "x_coordinate", "cell_x"),
    y = c("^y$", "centroid.*y", "cell.*y", "ycoord", "y_coordinate", "cell_y"),
    cell_type = c("cell.*type", "celltype", "annotation", "phenotype", "cluster"),
    neighborhood = c("neighbor", "neighbour", "community", "niche")
  )
  field_matches <- rbindlist(lapply(names(candidate_fields), function(f) {
    pats <- candidate_fields[[f]]
    hits <- names(metadata)[Reduce(`|`, lapply(pats, function(p) grepl(p, col_lower, ignore.case = TRUE)))]
    data.table(field = f, matched_columns = paste(hits, collapse = ";"), n_matches = length(hits))
  }))
  fwrite(field_matches, file.path(table_dir, "cho_imc_metadata_field_matches.csv"))
  metadata_summary <- rbind(metadata_summary, data.table(field = paste0("field_match_", field_matches$field), value = field_matches$matched_columns), fill = TRUE)
  
  pick_first <- function(field) {
    z <- field_matches[field == !!field]
    NA_character_
  }
  patient_col <- names(metadata)[grepl("patient|case|donor", col_lower)][1]
  sample_col <- names(metadata)[grepl("sample|specimen", col_lower)][1]
  roi_col <- names(metadata)[grepl("roi|fov|image|region|core", col_lower)][1]
  type_col <- names(metadata)[grepl("cell.*type|celltype|annotation|phenotype", col_lower)][1]
  site_col <- names(metadata)[grepl("primary|metast|tissue|site|organ", col_lower)][1]
  x_col <- names(metadata)[grepl("^x$|centroid.*x|cell.*x|xcoord|x_coordinate|cell_x", col_lower)][1]
  y_col <- names(metadata)[grepl("^y$|centroid.*y|cell.*y|ycoord|y_coordinate|cell_y", col_lower)][1]
  
  sample_roi_summary <- data.table(summary_level = "all", n_cells = nrow(metadata))
  if (!is.na(patient_col)) sample_roi_summary[, n_patients := uniqueN(metadata[[patient_col]])]
  if (!is.na(sample_col)) sample_roi_summary[, n_samples := uniqueN(metadata[[sample_col]])]
  if (!is.na(roi_col)) sample_roi_summary[, n_rois := uniqueN(metadata[[roi_col]])]
  if (!is.na(type_col)) sample_roi_summary[, n_cell_types := uniqueN(metadata[[type_col]])]
  if (!is.na(site_col)) sample_roi_summary[, n_sites := uniqueN(metadata[[site_col]])]
  sample_roi_summary[, patient_col := patient_col]
  sample_roi_summary[, sample_col := sample_col]
  sample_roi_summary[, roi_col := roi_col]
  sample_roi_summary[, cell_type_col := type_col]
  sample_roi_summary[, site_col := site_col]
  sample_roi_summary[, x_col := x_col]
  sample_roi_summary[, y_col := y_col]
}

fwrite(metadata_summary, file.path(table_dir, "cho_imc_object_metadata_summary.csv"))
fwrite(marker_availability, file.path(table_dir, "cho_imc_marker_availability.csv"))
fwrite(sample_roi_summary, file.path(table_dir, "cho_imc_sample_roi_summary.csv"))

report <- c(
  "# Cho IMC Feasibility Report",
  "",
  paste0("Object class: ", paste(class(obj), collapse = "/")),
  paste0("Object size: ", format(object.size(obj), units = "auto")),
  "",
  "## Data structure",
  paste0("- Metadata detected: ", !is.null(metadata)),
  paste0("- Expression matrix detected: ", !is.null(expr)),
  paste0("- Assay notes: ", paste(assay_notes, collapse = " | ")),
  "",
  "## Marker availability",
  paste(capture.output(print(marker_availability)), collapse = "\n"),
  "",
  "## Sample/ROI summary",
  paste(capture.output(print(sample_roi_summary)), collapse = "\n"),
  "",
  "## Feasibility",
  "- Spatial-neighborhood analysis feasibility depends on detected x/y and ROI/sample columns above.",
  "- CAF/stromal, myeloid/macrophage, epithelial and immune definitions depend on marker and/or author annotation coverage above.",
  "- Full downstream analysis should proceed only if cell coordinates, cell/ROI identifiers and marker matrix align."
)
writeLines(report, file.path(doc_dir, "cho_imc_feasibility_report.md"), useBytes = TRUE)

sink()
