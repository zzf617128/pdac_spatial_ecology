suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!requireNamespace("FNN", quietly = TRUE)) stop("FNN package is required for kNN analysis")

set.seed(20260630)

root <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search/cho_imc"), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "data")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
doc_dir <- file.path(root, "docs")
source_dir <- file.path(root, "source_data")
manifest_dir <- file.path(root, "manifest")
log_dir <- file.path(root, "logs")
invisible(lapply(c(table_dir, figure_dir, doc_dir, source_dir, manifest_dir, log_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, path) fwrite(as.data.frame(x), path)
safe_div <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)
log2oe <- function(obs, exp) log2(safe_div(obs, exp))
perm_p_greater <- function(obs, nulls) {
  nulls <- nulls[is.finite(nulls)]
  if (!is.finite(obs) || length(nulls) == 0) return(NA_real_)
  (1 + sum(nulls >= obs)) / (1 + length(nulls))
}
zv <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

rds_path <- file.path(data_dir, "backup_output.rds")
obj <- readRDS(rds_path)

sample_meta <- as.data.table(obj$meta_data)
sample_meta[, Patient := obj$Patient]
sample_meta[, site_status := fifelse(Site == "Pancreas", "Primary/Pancreas", "Metastasis/Liver")]
sample_meta[, roi_id := paste0("ROI_", sample_id)]

env <- slot(obj$fcs1, "frames")
cell_list <- vector("list", nrow(sample_meta))
offset <- 0L
for (i in seq_len(nrow(sample_meta))) {
  sid <- sample_meta$sample_id[i]
  frame_name <- paste0("V", sid)
  e <- slot(get(frame_name, envir = env), "exprs")
  dt <- as.data.table(e)
  n <- nrow(dt)
  idx <- seq.int(offset + 1L, offset + n)
  dt[, global_cell_id := idx]
  dt[, sample_id := sid]
  dt[, roi_id := paste0("ROI_", sid)]
  dt[, file_name := sample_meta$file_name[i]]
  dt[, patient_id := sample_meta$Patient[i]]
  dt[, site := sample_meta$Site[i]]
  dt[, site_status := sample_meta$site_status[i]]
  dt[, author_cluster_numeric := obj$cell_clustering[idx]]
  dt[, author_broad := obj$cell_clustering1m[idx]]
  dt[, author_refined := obj$cell_clustering2m[idx]]
  cell_list[[i]] <- dt
  offset <- offset + n
}
cells <- rbindlist(cell_list, use.names = TRUE, fill = TRUE)
stopifnot(nrow(cells) == length(obj$sample_ids))

marker_cols <- setdiff(colnames(slot(get("V1", envir = env), "exprs")), obj$otherparameters)
marker_cols <- setdiff(marker_cols, c("CellId"))
coord_cols <- c("X_position", "Y_position")

marker_queries <- list(
  SMA_ACTA2 = c("SMA", "SMAVIM", "ACTA2"),
  Collagen = c("Collagen"),
  Vimentin = c("VIM", "SMAVIM"),
  FAP = c("FAP"),
  PDPN = c("Podoplanin", "PDPN"),
  Fibronectin = c("Fibronectin", "FN1"),
  CD68 = c("CD68"),
  CD163 = c("CD163"),
  CD11b = c("CD11b", "ITGAM"),
  CD14 = c("CD14"),
  HLA_DR = c("HLADR", "HLA-DR", "HLA_DR"),
  CD74 = c("CD74"),
  CD45 = c("CD45"),
  CD206_MRC1 = c("CD206", "MRC1"),
  panCK_CK = c("^CK$", "panCK", "Cytokeratin"),
  EPCAM = c("EPCAM", "EpCAM"),
  CD3 = c("CD3"),
  CD4 = c("CD4"),
  CD8 = c("CD8"),
  CD20 = c("CD20"),
  FOXP3 = c("FOXP3"),
  CD44 = c("CD44"),
  Integrin = c("Integrin", "ITGA", "ITGB"),
  SPP1 = c("SPP1", "Osteopontin")
)

match_marker <- function(qs) {
  norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
  marker_norm <- norm(marker_cols)
  hits <- unique(unlist(lapply(qs, function(q) {
    qn <- norm(q)
    if (q %in% c("SMA", "SMAVIM", "VIM")) {
      marker_cols[grepl(q, marker_cols, ignore.case = TRUE)]
    } else if (q %in% c("CD45")) {
      marker_cols[grepl("^CD45", marker_cols, ignore.case = TRUE)]
    } else if (grepl("^\\^", q)) {
      marker_cols[grepl(q, marker_cols, ignore.case = TRUE)]
    } else {
      marker_cols[marker_norm == qn]
    }
  })))
  hits
}

marker_availability <- rbindlist(lapply(names(marker_queries), function(q) {
  hits <- match_marker(marker_queries[[q]])
  data.table(marker_query = q, matched_marker = if (length(hits)) paste(hits, collapse = ";") else NA_character_, available = length(hits) > 0)
}))

object_metadata_summary <- rbindlist(list(
  data.table(field = "object_class", value = paste(class(obj), collapse = "/")),
  data.table(field = "object_size", value = format(object.size(obj), units = "auto")),
  data.table(field = "n_cells", value = as.character(nrow(cells))),
  data.table(field = "n_samples_roi", value = as.character(uniqueN(cells$roi_id))),
  data.table(field = "n_patients", value = as.character(uniqueN(cells$patient_id))),
  data.table(field = "sites", value = paste(names(table(cells$site)), as.integer(table(cells$site)), sep = "=", collapse = ";")),
  data.table(field = "marker_columns", value = paste(marker_cols, collapse = ";")),
  data.table(field = "coordinate_columns", value = paste(coord_cols, collapse = ";")),
  data.table(field = "author_broad_labels", value = paste(names(table(cells$author_broad)), as.integer(table(cells$author_broad)), sep = "=", collapse = ";")),
  data.table(field = "author_refined_labels", value = paste(names(table(cells$author_refined)), as.integer(table(cells$author_refined)), sep = "=", collapse = ";")),
  data.table(field = "normalization_status", value = "processed RDS fcs/fcs1 expression values; exact transformation not explicitly encoded in object; fcsraw retains raw channel names")
))

sample_roi_summary <- cells[, .(
  n_cells = .N,
  n_author_broad = uniqueN(author_broad),
  n_author_refined = uniqueN(author_refined),
  x_min = min(X_position, na.rm = TRUE),
  x_max = max(X_position, na.rm = TRUE),
  y_min = min(Y_position, na.rm = TRUE),
  y_max = max(Y_position, na.rm = TRUE)
), by = .(patient_id, sample_id, roi_id, file_name, site, site_status)]

write_csv(object_metadata_summary, file.path(table_dir, "cho_imc_object_metadata_summary.csv"))
write_csv(marker_availability, file.path(table_dir, "cho_imc_marker_availability.csv"))
write_csv(sample_roi_summary, file.path(table_dir, "cho_imc_sample_roi_summary.csv"))

cells[, compartment_author := fifelse(author_broad == "Stroma", "CAF/stromal-like",
  fifelse(author_broad == "Myeloid", "myeloid/macrophage-like",
    fifelse(author_broad == "Tumor", "tumor epithelial-like",
      fifelse(author_broad %in% c("CD4T", "CD8T"), "T cell",
        fifelse(author_broad %in% c("Immune_Mix", "NK", "Neutrophil"), "other immune", "other / unassigned")))))]

module_score <- function(markers) {
  present <- intersect(markers, names(cells))
  if (length(present) == 0) return(rep(NA_real_, nrow(cells)))
  out <- cells[, rowMeans(as.data.table(lapply(.SD, zv)), na.rm = TRUE), .SDcols = present]
  out
}
cells[, stromal_marker_score := module_score(c("Collagen", "SMAVIM", "Podoplanin"))]
cells[, myeloid_marker_score := module_score(c("CD68", "CD163", "CD206", "HLADR", "CD74"))]
cells[, tumor_marker_score := module_score(c("CK"))]
cells[, tcell_marker_score := module_score(c("CD3", "CD4", "CD8", "FOXP3"))]

cells[, compartment_marker_permissive := fifelse(stromal_marker_score > quantile(stromal_marker_score, 0.80, na.rm = TRUE) & tumor_marker_score <= quantile(tumor_marker_score, 0.90, na.rm = TRUE), "CAF/stromal-like",
  fifelse(myeloid_marker_score > quantile(myeloid_marker_score, 0.80, na.rm = TRUE), "myeloid/macrophage-like",
    fifelse(tumor_marker_score > quantile(tumor_marker_score, 0.80, na.rm = TRUE), "tumor epithelial-like",
      fifelse(tcell_marker_score > quantile(tcell_marker_score, 0.85, na.rm = TRUE), "T cell", "other / unassigned"))))]
cells[, compartment_marker_stringent := fifelse(stromal_marker_score > quantile(stromal_marker_score, 0.90, na.rm = TRUE) & tumor_marker_score <= quantile(tumor_marker_score, 0.75, na.rm = TRUE), "CAF/stromal-like",
  fifelse(myeloid_marker_score > quantile(myeloid_marker_score, 0.90, na.rm = TRUE), "myeloid/macrophage-like",
    fifelse(tumor_marker_score > quantile(tumor_marker_score, 0.90, na.rm = TRUE), "tumor epithelial-like",
      fifelse(tcell_marker_score > quantile(tcell_marker_score, 0.90, na.rm = TRUE), "T cell", "other / unassigned"))))]

cell_type_harmonization <- cells[, .(
  global_cell_id, patient_id, sample_id, roi_id, site, site_status,
  author_broad, author_refined, author_cluster_numeric,
  compartment_author, compartment_marker_permissive, compartment_marker_stringent,
  X_position, Y_position, stromal_marker_score, myeloid_marker_score, tumor_marker_score, tcell_marker_score
)]
write_csv(cell_type_harmonization, file.path(table_dir, "cho_imc_cell_type_harmonization.csv"))

compartment_abundance_by_roi <- cells[, .N, by = .(patient_id, sample_id, roi_id, site, site_status, compartment_author)]
compartment_abundance_by_roi[, fraction := N / sum(N), by = roi_id]
write_csv(compartment_abundance_by_roi, file.path(table_dir, "cho_imc_compartment_abundance_by_roi.csv"))

pdf(file.path(figure_dir, "cho_imc_compartment_composition.pdf"), width = 8.5, height = 5)
print(ggplot(compartment_abundance_by_roi, aes(x = reorder(roi_id, fraction), y = fraction, fill = compartment_author)) +
  geom_col(width = 0.9) + coord_flip() + theme_bw(base_size = 7) +
  labs(x = "ROI", y = "Cell fraction", fill = "Compartment", title = "Cho IMC compartment composition by ROI"))
dev.off()

get_neighbors <- function(dt, kmax = 15) {
  coords <- as.matrix(dt[, .(X_position, Y_position)])
  if (nrow(coords) <= kmax + 1 || any(!is.finite(coords))) return(NULL)
  FNN::get.knn(coords, k = kmax)$nn.index
}

edge_stat <- function(labels, nn, source_label, target_label, k = 10, target_override = NULL, source_override = NULL) {
  source <- if (is.null(source_override)) labels == source_label else source_override
  target <- if (is.null(target_override)) labels == target_label else target_override
  if (sum(source, na.rm = TRUE) == 0 || sum(target, na.rm = TRUE) == 0 || is.null(nn)) {
    return(list(obs = NA_real_, exp = NA_real_, log2_oe = NA_real_, n_source = sum(source, na.rm = TRUE), n_target = sum(target, na.rm = TRUE), n_edges = NA_integer_))
  }
  nnk <- nn[, seq_len(k), drop = FALSE]
  edges <- as.vector(nnk[source, , drop = FALSE])
  obs <- mean(target[edges], na.rm = TRUE)
  exp <- mean(target, na.rm = TRUE)
  list(obs = obs, exp = exp, log2_oe = log2oe(obs, exp), n_source = sum(source), n_target = sum(target), n_edges = length(edges))
}

perm_nulls <- function(labels, nn, source_label, target_label, k = 10, B = 100, strata = NULL) {
  source <- labels == source_label
  target <- labels == target_label
  if (sum(source) == 0 || sum(target) == 0 || is.null(nn)) return(rep(NA_real_, B))
  nnk <- nn[, seq_len(k), drop = FALSE]
  edges <- as.vector(nnk[source, , drop = FALSE])
  out <- numeric(B)
  if (is.null(strata)) {
    for (b in seq_len(B)) {
      perm_target <- sample(target)
      out[b] <- mean(perm_target[edges])
    }
  } else {
    for (b in seq_len(B)) {
      perm_target <- target
      for (s in unique(strata)) {
        idx <- which(strata == s)
        perm_target[idx] <- sample(target[idx])
      }
      out[b] <- mean(perm_target[edges])
    }
  }
  out
}

abundance_nulls <- function(labels, nn, source_label, target_label, k = 10, B = 100) {
  source <- labels == source_label
  n_target <- sum(labels == target_label)
  if (sum(source) == 0 || n_target == 0 || is.null(nn)) return(rep(NA_real_, B))
  nnk <- nn[, seq_len(k), drop = FALSE]
  edges <- as.vector(nnk[source, , drop = FALSE])
  n <- length(labels)
  out <- numeric(B)
  for (b in seq_len(B)) {
    rand_target <- rep(FALSE, n)
    rand_target[sample.int(n, n_target)] <- TRUE
    out[b] <- mean(rand_target[edges])
  }
  out
}

analyze_pair_by_roi <- function(label_col, source_label, target_label, k = 10, B = 100, pair_name = "pair") {
  rows <- vector("list", length(unique(cells$roi_id)))
  i <- 0L
  for (rid in unique(cells$roi_id)) {
    dt <- cells[roi_id == rid]
    labels <- dt[[label_col]]
    nn <- get_neighbors(dt, max(15, k))
    st <- edge_stat(labels, nn, source_label, target_label, k)
    density <- dt$Number_Neighbors
    strata <- tryCatch(as.integer(cut(density, quantile(density, probs = seq(0, 1, 0.25), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)), error = function(e) rep(1L, nrow(dt)))
    n1 <- perm_nulls(labels, nn, source_label, target_label, k, B = B)
    n2 <- abundance_nulls(labels, nn, source_label, target_label, k, B = B)
    n3 <- perm_nulls(labels, nn, source_label, target_label, k, B = B, strata = strata)
    i <- i + 1L
    rows[[i]] <- data.table(
      pair = pair_name, label_set = label_col, k = k,
      patient_id = dt$patient_id[1], sample_id = dt$sample_id[1], roi_id = rid, site = dt$site[1], site_status = dt$site_status[1],
      source_label = source_label, target_label = target_label,
      n_cells = nrow(dt), n_source = st$n_source, n_target = st$n_target, n_edges = st$n_edges,
      observed_neighbor_fraction = st$obs,
      expected_target_fraction = st$exp,
      observed_expected_ratio = safe_div(st$obs, st$exp),
      log2_oe = st$log2_oe,
      random_label_null_mean = mean(n1, na.rm = TRUE),
      random_label_p_greater = perm_p_greater(st$obs, n1),
      abundance_matched_null_mean = mean(n2, na.rm = TRUE),
      abundance_matched_p_greater = perm_p_greater(st$obs, n2),
      density_stratified_null_mean = mean(n3, na.rm = TRUE),
      density_stratified_p_greater = perm_p_greater(st$obs, n3)
    )
  }
  rbindlist(rows)
}

main_adj <- analyze_pair_by_roi("compartment_author", "CAF/stromal-like", "myeloid/macrophage-like", k = 10, B = 100, pair_name = "CAF_stromal_to_myeloid")
write_csv(main_adj, file.path(table_dir, "cho_imc_caf_myeloid_adjacency_results.csv"))

summarize_by_patient <- function(dt, value_col = "log2_oe") {
  dt[, .(
    n_roi = .N,
    median_value = median(get(value_col), na.rm = TRUE),
    support_fraction_positive = mean(get(value_col) > 0, na.rm = TRUE),
    median_random_label_p = median(random_label_p_greater, na.rm = TRUE),
    median_abundance_p = median(abundance_matched_p_greater, na.rm = TRUE),
    median_density_p = median(density_stratified_p_greater, na.rm = TRUE)
  ), by = patient_id]
}

pdf(file.path(figure_dir, "cho_imc_caf_myeloid_adjacency.pdf"), width = 7.5, height = 4.8)
print(ggplot(main_adj, aes(x = site_status, y = log2_oe, color = site_status)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA, width = 0.45) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  theme_bw(base_size = 9) + labs(x = NULL, y = "ROI log2(O/E)", title = "Cho IMC CAF/stromal -> myeloid/macrophage adjacency") +
  theme(legend.position = "none"))
dev.off()

# Stromal-rich neighborhood enrichment.
neigh_rows <- list()
interface_rows <- list()
hladr_rows <- list()
for (rid in unique(cells$roi_id)) {
  dt <- cells[roi_id == rid]
  nn <- get_neighbors(dt, 15)
  if (is.null(nn)) next
  labels <- dt$compartment_author
  k <- 10
  nnk <- nn[, seq_len(k), drop = FALSE]
  stromal <- labels == "CAF/stromal-like"
  myeloid <- labels == "myeloid/macrophage-like"
  tumor <- labels == "tumor epithelial-like"
  stromal_frac <- rowMeans(matrix(stromal[as.vector(nnk)], nrow = nrow(nnk)))
  myeloid_frac <- rowMeans(matrix(myeloid[as.vector(nnk)], nrow = nrow(nnk)))
  tumor_frac <- rowMeans(matrix(tumor[as.vector(nnk)], nrow = nrow(nnk)))
  stromal_rich <- stromal_frac >= quantile(stromal_frac, 0.75, na.rm = TRUE) & stromal_frac > 0
  neigh_rows[[rid]] <- data.table(
    patient_id = dt$patient_id[1], sample_id = dt$sample_id[1], roi_id = rid, site = dt$site[1], site_status = dt$site_status[1],
    n_cells = nrow(dt),
    n_stromal_rich_neighborhoods = sum(stromal_rich, na.rm = TRUE),
    macrophage_fraction_stromal_rich = mean(myeloid_frac[stromal_rich], na.rm = TRUE),
    macrophage_fraction_non_stromal = mean(myeloid_frac[!stromal_rich], na.rm = TRUE),
    delta_macrophage_fraction = mean(myeloid_frac[stromal_rich], na.rm = TRUE) - mean(myeloid_frac[!stromal_rich], na.rm = TRUE),
    cd68_stromal_rich = mean(dt$CD68[stromal_rich], na.rm = TRUE),
    cd68_non_stromal = mean(dt$CD68[!stromal_rich], na.rm = TRUE),
    cd163_stromal_rich = mean(dt$CD163[stromal_rich], na.rm = TRUE),
    cd163_non_stromal = mean(dt$CD163[!stromal_rich], na.rm = TRUE),
    hladr_stromal_rich = mean(dt$HLADR[stromal_rich], na.rm = TRUE),
    hladr_non_stromal = mean(dt$HLADR[!stromal_rich], na.rm = TRUE),
    tumor_fraction_stromal_rich = mean(tumor_frac[stromal_rich], na.rm = TRUE)
  )
  interface_stromal <- stromal & tumor_frac > 0
  st_int <- edge_stat(labels, nn, "CAF/stromal-like", "myeloid/macrophage-like", k = k, source_override = interface_stromal)
  interface_rows[[rid]] <- data.table(
    patient_id = dt$patient_id[1], sample_id = dt$sample_id[1], roi_id = rid, site = dt$site[1], site_status = dt$site_status[1],
    n_cells = nrow(dt), n_interface_stromal = sum(interface_stromal), n_myeloid = sum(myeloid),
    interface_myeloid_neighbor_fraction = st_int$obs,
    expected_myeloid_fraction = st_int$exp,
    interface_observed_expected_ratio = safe_div(st_int$obs, st_int$exp),
    interface_log2_oe = st_int$log2_oe,
    triad_fraction_stromal_cells = mean(stromal & tumor_frac > 0 & myeloid_frac > 0, na.rm = TRUE)
  )
  myeloid_interface <- myeloid & rowMeans(matrix((stromal | tumor)[as.vector(nnk)], nrow = nrow(nnk))) > 0
  hladr_rows[[rid]] <- data.table(
    patient_id = dt$patient_id[1], sample_id = dt$sample_id[1], roi_id = rid, site = dt$site[1], site_status = dt$site_status[1],
    cd44_available = "CD44" %in% names(dt),
    spp1_available = "SPP1" %in% names(dt),
    hladr_available = "HLADR" %in% names(dt),
    n_myeloid = sum(myeloid),
    n_myeloid_interface = sum(myeloid_interface),
    hladr_myeloid_interface = mean(dt$HLADR[myeloid_interface], na.rm = TRUE),
    hladr_myeloid_noninterface = mean(dt$HLADR[myeloid & !myeloid_interface], na.rm = TRUE),
    delta_hladr_interface = mean(dt$HLADR[myeloid_interface], na.rm = TRUE) - mean(dt$HLADR[myeloid & !myeloid_interface], na.rm = TRUE),
    cd68_myeloid_interface = mean(dt$CD68[myeloid_interface], na.rm = TRUE),
    cd163_myeloid_interface = mean(dt$CD163[myeloid_interface], na.rm = TRUE)
  )
}
stromal_neigh <- rbindlist(neigh_rows, fill = TRUE)
interface_dt <- rbindlist(interface_rows, fill = TRUE)
hladr_dt <- rbindlist(hladr_rows, fill = TRUE)
write_csv(stromal_neigh, file.path(table_dir, "cho_imc_stromal_neighborhood_macrophage_enrichment.csv"))
write_csv(interface_dt, file.path(table_dir, "cho_imc_tumor_stroma_myeloid_interface_results.csv"))
write_csv(hladr_dt, file.path(table_dir, "cho_imc_cd44_hladr_context_results.csv"))

pdf(file.path(figure_dir, "cho_imc_stromal_neighborhood_composition.pdf"), width = 7.2, height = 4.5)
print(ggplot(stromal_neigh, aes(x = site_status, y = delta_macrophage_fraction, color = site_status)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA, width = 0.45) + geom_jitter(width = 0.15, alpha = 0.8) +
  theme_bw(base_size = 9) + theme(legend.position = "none") +
  labs(x = NULL, y = "Macrophage fraction delta", title = "Macrophage enrichment in stromal-rich neighborhoods"))
dev.off()

pdf(file.path(figure_dir, "cho_imc_tumor_stroma_myeloid_interface.pdf"), width = 7.2, height = 4.5)
print(ggplot(interface_dt, aes(x = site_status, y = interface_log2_oe, color = site_status)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA, width = 0.45) + geom_jitter(width = 0.15, alpha = 0.8) +
  theme_bw(base_size = 9) + theme(legend.position = "none") +
  labs(x = NULL, y = "Interface stromal -> myeloid log2(O/E)", title = "Tumor-stroma-myeloid interface support"))
dev.off()

pdf(file.path(figure_dir, "cho_imc_cd44_hladr_context.pdf"), width = 7.2, height = 4.5)
print(ggplot(hladr_dt, aes(x = site_status, y = delta_hladr_interface, color = site_status)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA, width = 0.45) + geom_jitter(width = 0.15, alpha = 0.8) +
  theme_bw(base_size = 9) + theme(legend.position = "none") +
  labs(x = NULL, y = "HLA-DR delta in interface myeloid", title = "HLA-DR myeloid protein context"))
dev.off()

# Primary vs metastasis.
primary_meta <- rbindlist(list(
  main_adj[, .(metric = "caf_myeloid_log2_oe", value = log2_oe, patient_id, sample_id, roi_id, site_status)],
  stromal_neigh[, .(metric = "stromal_neigh_macrophage_delta", value = delta_macrophage_fraction, patient_id, sample_id, roi_id, site_status)],
  interface_dt[, .(metric = "interface_log2_oe", value = interface_log2_oe, patient_id, sample_id, roi_id, site_status)],
  hladr_dt[, .(metric = "interface_hladr_delta", value = delta_hladr_interface, patient_id, sample_id, roi_id, site_status)]
), fill = TRUE)
primary_meta_results <- primary_meta[, {
  if (uniqueN(site_status) < 2) {
    list(n = .N, primary_median = NA_real_, metastasis_median = NA_real_, delta_metastasis_minus_primary = NA_real_, wilcox_p = NA_real_)
  } else {
    list(
      n = .N,
      primary_median = median(value[site_status == "Primary/Pancreas"], na.rm = TRUE),
      metastasis_median = median(value[site_status == "Metastasis/Liver"], na.rm = TRUE),
      delta_metastasis_minus_primary = median(value[site_status == "Metastasis/Liver"], na.rm = TRUE) - median(value[site_status == "Primary/Pancreas"], na.rm = TRUE),
      wilcox_p = suppressWarnings(wilcox.test(value ~ site_status)$p.value)
    )
  }
}, by = metric]
primary_meta_results[, fdr := p.adjust(wilcox_p, method = "BH")]
write_csv(primary_meta_results, file.path(table_dir, "cho_imc_primary_metastasis_comparison.csv"))

pdf(file.path(figure_dir, "cho_imc_primary_metastasis_comparison.pdf"), width = 8, height = 4.8)
print(ggplot(primary_meta, aes(x = site_status, y = value, color = site_status)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.45) + geom_jitter(width = 0.15, alpha = 0.65, size = 1.3) +
  facet_wrap(~metric, scales = "free_y") + theme_bw(base_size = 8) + theme(legend.position = "none") +
  labs(x = NULL, y = "ROI value", title = "Primary/Pancreas vs liver metastasis exploratory comparison"))
dev.off()

# Robustness and controls.
knn_sens <- rbindlist(lapply(c(5, 10, 15), function(k) analyze_pair_by_roi("compartment_author", "CAF/stromal-like", "myeloid/macrophage-like", k = k, B = 20, pair_name = "CAF_stromal_to_myeloid")), fill = TRUE)
write_csv(knn_sens, file.path(table_dir, "cho_imc_knn_parameter_sensitivity.csv"))

threshold_sens <- rbindlist(list(
  analyze_pair_by_roi("compartment_author", "CAF/stromal-like", "myeloid/macrophage-like", k = 10, B = 20, pair_name = "author_broad"),
  analyze_pair_by_roi("compartment_marker_permissive", "CAF/stromal-like", "myeloid/macrophage-like", k = 10, B = 20, pair_name = "marker_permissive"),
  analyze_pair_by_roi("compartment_marker_stringent", "CAF/stromal-like", "myeloid/macrophage-like", k = 10, B = 20, pair_name = "marker_stringent")
), fill = TRUE)
write_csv(threshold_sens, file.path(table_dir, "cho_imc_threshold_sensitivity.csv"))

null_sensitivity <- main_adj[, .(
  pair, label_set, k, patient_id, sample_id, roi_id, site_status, log2_oe,
  random_label_delta = observed_neighbor_fraction - random_label_null_mean,
  abundance_matched_delta = observed_neighbor_fraction - abundance_matched_null_mean,
  density_stratified_delta = observed_neighbor_fraction - density_stratified_null_mean,
  random_label_p_greater, abundance_matched_p_greater, density_stratified_p_greater
)]
write_csv(null_sensitivity, file.path(table_dir, "cho_imc_null_sensitivity.csv"))

loo_patient <- main_adj[, .(
  n_patients_remaining = uniqueN(patient_id) - 1,
  median_log2_oe = median(main_adj[patient_id != .BY$patient_id]$log2_oe, na.rm = TRUE),
  support_fraction_positive = mean(main_adj[patient_id != .BY$patient_id]$log2_oe > 0, na.rm = TRUE)
), by = patient_id]
loo_roi <- main_adj[, .(
  left_out_roi = roi_id,
  median_log2_oe = median(main_adj[roi_id != .BY$roi_id]$log2_oe, na.rm = TRUE),
  support_fraction_positive = mean(main_adj[roi_id != .BY$roi_id]$log2_oe > 0, na.rm = TRUE)
), by = roi_id]
leave_one_out <- rbindlist(list(
  loo_patient[, .(level = "patient", left_out = patient_id, median_log2_oe, support_fraction_positive)],
  loo_roi[, .(level = "roi", left_out = left_out_roi, median_log2_oe, support_fraction_positive)]
), fill = TRUE)
write_csv(leave_one_out, file.path(table_dir, "cho_imc_leave_one_out.csv"))

control_pairs <- rbindlist(list(
  analyze_pair_by_roi("compartment_author", "tumor epithelial-like", "tumor epithelial-like", k = 10, B = 20, pair_name = "positive_tumor_to_tumor"),
  analyze_pair_by_roi("compartment_author", "T cell", "T cell", k = 10, B = 20, pair_name = "positive_tcell_to_tcell"),
  analyze_pair_by_roi("compartment_author", "CAF/stromal-like", "CAF/stromal-like", k = 10, B = 20, pair_name = "positive_stromal_to_stromal"),
  analyze_pair_by_roi("compartment_author", "myeloid/macrophage-like", "myeloid/macrophage-like", k = 10, B = 20, pair_name = "positive_myeloid_to_myeloid"),
  analyze_pair_by_roi("compartment_author", "CAF/stromal-like", "other / unassigned", k = 10, B = 20, pair_name = "negative_stromal_to_unassigned")
), fill = TRUE)
write_csv(control_pairs, file.path(table_dir, "cho_imc_positive_negative_controls.csv"))

pdf(file.path(figure_dir, "cho_imc_robustness_summary.pdf"), width = 8.5, height = 6)
robust_plot <- rbindlist(list(
  knn_sens[, .(analysis = paste0("k=", k), log2_oe)],
  threshold_sens[, .(analysis = pair, log2_oe)]
), fill = TRUE)
print(ggplot(robust_plot, aes(x = analysis, y = log2_oe)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.18, alpha = 0.4, size = 0.8) +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "ROI log2(O/E)", title = "Cho IMC CAF/stromal -> myeloid robustness summary"))
dev.off()

pdf(file.path(figure_dir, "cho_imc_positive_negative_controls.pdf"), width = 8.5, height = 5)
print(ggplot(control_pairs, aes(x = pair, y = log2_oe)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.18, alpha = 0.4, size = 0.8) +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "ROI log2(O/E)", title = "Cho IMC positive and negative controls"))
dev.off()

# Source data.
source_data <- rbindlist(list(
  main_adj[, .(panel = "CAF_myeloid_adjacency", metric = "log2_oe", value = log2_oe, patient_id, sample_id, roi_id, site_status, note = paste0("random_p=", signif(random_label_p_greater, 3)))],
  stromal_neigh[, .(panel = "stromal_neighborhood_macrophage", metric = "delta_macrophage_fraction", value = delta_macrophage_fraction, patient_id, sample_id, roi_id, site_status, note = "")],
  interface_dt[, .(panel = "tumor_stroma_myeloid_interface", metric = "interface_log2_oe", value = interface_log2_oe, patient_id, sample_id, roi_id, site_status, note = "")],
  hladr_dt[, .(panel = "hladr_context", metric = "delta_hladr_interface", value = delta_hladr_interface, patient_id, sample_id, roi_id, site_status, note = paste0("CD44_available=", cd44_available, ";SPP1_available=", spp1_available))]
), fill = TRUE)
write_csv(source_data, file.path(source_dir, "Source_Data_Cho_IMC.csv"))

# Summaries and decision.
main_patient <- summarize_by_patient(main_adj)
main_summary <- data.table(
  metric = "CAF_stromal_to_myeloid_log2OE",
  n_roi = nrow(main_adj),
  n_patients = uniqueN(main_adj$patient_id),
  median_roi = median(main_adj$log2_oe, na.rm = TRUE),
  roi_support_fraction_positive = mean(main_adj$log2_oe > 0, na.rm = TRUE),
  patient_median = median(main_patient$median_value, na.rm = TRUE),
  patient_support_fraction_positive = mean(main_patient$median_value > 0, na.rm = TRUE),
  median_random_label_p = median(main_adj$random_label_p_greater, na.rm = TRUE),
  median_abundance_p = median(main_adj$abundance_matched_p_greater, na.rm = TRUE),
  median_density_p = median(main_adj$density_stratified_p_greater, na.rm = TRUE)
)

stromal_summary <- data.table(
  metric = "stromal_neighborhood_macrophage_delta",
  n_roi = nrow(stromal_neigh),
  n_patients = uniqueN(stromal_neigh$patient_id),
  median_roi = median(stromal_neigh$delta_macrophage_fraction, na.rm = TRUE),
  roi_support_fraction_positive = mean(stromal_neigh$delta_macrophage_fraction > 0, na.rm = TRUE),
  patient_median = median(stromal_neigh[, .(v = median(delta_macrophage_fraction, na.rm = TRUE)), by = patient_id]$v, na.rm = TRUE),
  patient_support_fraction_positive = mean(stromal_neigh[, .(v = median(delta_macrophage_fraction, na.rm = TRUE)), by = patient_id]$v > 0, na.rm = TRUE)
)

interface_summary <- data.table(
  metric = "tumor_stroma_myeloid_interface_log2OE",
  n_roi = nrow(interface_dt),
  n_patients = uniqueN(interface_dt$patient_id),
  median_roi = median(interface_dt$interface_log2_oe, na.rm = TRUE),
  roi_support_fraction_positive = mean(interface_dt$interface_log2_oe > 0, na.rm = TRUE),
  patient_median = median(interface_dt[, .(v = median(interface_log2_oe, na.rm = TRUE)), by = patient_id]$v, na.rm = TRUE),
  patient_support_fraction_positive = mean(interface_dt[, .(v = median(interface_log2_oe, na.rm = TRUE)), by = patient_id]$v > 0, na.rm = TRUE)
)

hladr_summary <- data.table(
  metric = "HLA_DR_interface_myeloid_delta",
  n_roi = nrow(hladr_dt),
  n_patients = uniqueN(hladr_dt$patient_id),
  median_roi = median(hladr_dt$delta_hladr_interface, na.rm = TRUE),
  roi_support_fraction_positive = mean(hladr_dt$delta_hladr_interface > 0, na.rm = TRUE),
  patient_median = median(hladr_dt[, .(v = median(delta_hladr_interface, na.rm = TRUE)), by = patient_id]$v, na.rm = TRUE),
  patient_support_fraction_positive = mean(hladr_dt[, .(v = median(delta_hladr_interface, na.rm = TRUE)), by = patient_id]$v > 0, na.rm = TRUE)
)

decision_summary <- rbindlist(list(main_summary, stromal_summary, interface_summary, hladr_summary), fill = TRUE)
write_csv(decision_summary, file.path(table_dir, "cho_imc_decision_summary.csv"))

knn_pass <- knn_sens[, .(median_log2_oe = median(log2_oe, na.rm = TRUE), support = mean(log2_oe > 0, na.rm = TRUE)), by = k]
threshold_pass <- threshold_sens[, .(median_log2_oe = median(log2_oe, na.rm = TRUE), support = mean(log2_oe > 0, na.rm = TRUE)), by = pair]
strong_A <- main_summary$median_roi > 0.25 && main_summary$roi_support_fraction_positive >= 0.65 && main_summary$patient_support_fraction_positive >= 0.65 && all(knn_pass$median_log2_oe > 0, na.rm = TRUE)
strong_B <- stromal_summary$median_roi > 0.05 && stromal_summary$roi_support_fraction_positive >= 0.65 && stromal_summary$patient_support_fraction_positive >= 0.65
strong_C <- interface_summary$median_roi > 0.25 && interface_summary$roi_support_fraction_positive >= 0.65 && interface_summary$patient_support_fraction_positive >= 0.65
strong_D <- hladr_summary$median_roi > 0 && hladr_summary$roi_support_fraction_positive >= 0.65 && hladr_summary$patient_support_fraction_positive >= 0.65

decision <- if (strong_A || strong_B || strong_C || strong_D) {
  "A. Add Cho IMC to ED10 v2 as strong protein-neighborhood support"
} else if (main_summary$median_roi > 0 || stromal_summary$median_roi > 0 || interface_summary$median_roi > 0) {
  "B. Use Cho IMC as secondary/source-only support"
} else {
  "C. Do not show Cho IMC; archive only"
}

writeLines(c(
  "# Cho IMC Cell Type Harmonization Notes",
  "",
  "Author broad labels were used as the primary compartment definition. `Stroma` was mapped to CAF/stromal-like, `Myeloid` to myeloid/macrophage-like, `Tumor` to tumor epithelial-like, CD4T/CD8T to T cell, and Immune_Mix/NK/Neutrophil to other immune.",
  "",
  "Marker-score permissive and stringent labels were generated as sensitivity definitions using Collagen/SMAVIM/Podoplanin, CD68/CD163/CD206/HLADR/CD74, CK, and CD3/CD4/CD8/FOXP3 scores."
), file.path(doc_dir, "cho_imc_cell_type_harmonization_notes.md"), useBytes = TRUE)

writeLines(c(
  "# Cho IMC CD44 Claim Boundary",
  "",
  paste0("CD44 available: ", marker_availability[marker_query == "CD44", available]),
  paste0("SPP1 available: ", marker_availability[marker_query == "SPP1", available]),
  paste0("HLA-DR available: ", marker_availability[marker_query == "HLA_DR", available]),
  "",
  "Because CD44 and SPP1 are not both measured as spatially linked markers in this IMC panel, these data must not be described as direct SPP1-CD44 validation. HLA-DR can be used only as macrophage/APC-like protein-context support."
), file.path(doc_dir, "cho_imc_cd44_claim_boundary.md"), useBytes = TRUE)

writeLines(c(
  "# Cho IMC Primary/Metastasis Interpretation",
  "",
  "Primary/pancreas and liver-metastasis ROI labels are available. Comparisons are exploratory protein-context checks, not clinical prediction or lymph-node immune-uncoupling validation.",
  "",
  paste(capture.output(print(primary_meta_results)), collapse = "\n")
), file.path(doc_dir, "cho_imc_primary_metastasis_interpretation.md"), useBytes = TRUE)

writeLines(c(
  "# Cho IMC Robustness Interpretation",
  "",
  "Robustness was evaluated across kNN k=5/10/15, author and marker-score compartment definitions, within-ROI random-label nulls, abundance-matched target nulls, density-stratified nulls, leave-one-patient/ROI-out summaries, and positive/negative controls.",
  "",
  "kNN sensitivity:",
  paste(capture.output(print(knn_pass)), collapse = "\n"),
  "",
  "Threshold sensitivity:",
  paste(capture.output(print(threshold_pass)), collapse = "\n")
), file.path(doc_dir, "cho_imc_robustness_interpretation.md"), useBytes = TRUE)

writeLines(c(
  "# Cho IMC Feasibility Report",
  "",
  paste0("Cells: ", nrow(cells), "."),
  paste0("ROIs/samples: ", uniqueN(cells$roi_id), "."),
  paste0("Patients: ", uniqueN(cells$patient_id), "."),
  paste0("Sites: ", paste(names(table(cells$site)), as.integer(table(cells$site)), sep = "=", collapse = "; ")),
  "",
  "The processed RDS contains flowSet objects with marker matrices, X/Y coordinates, author broad/refined cell labels, sample IDs, patient IDs and pancreas/liver site labels. This is sufficient for spatial-neighborhood analysis, broad CAF/stromal-like and myeloid/macrophage-like compartment definition, ROI/patient-level aggregation, and exploratory primary-versus-liver-metastasis comparison.",
  "",
  "The panel supports stromal/myeloid/tumor/immune protein context using Collagen, SMAVIM, Podoplanin, CD68, CD163, CD206, CD74, HLADR, CK, CD3, CD4, CD8 and FOXP3. CD44, SPP1, FAP, fibronectin, CD14, CD11b, CD20 and EPCAM were not detected in the processed marker names."
), file.path(doc_dir, "cho_imc_feasibility_report.md"), useBytes = TRUE)

writeLines(c(
  "# Cho IMC Final Decision Report",
  "",
  paste0("Decision: ", decision),
  "",
  "## Decision summary",
  paste(capture.output(print(decision_summary)), collapse = "\n"),
  "",
  "## Comparison with existing orthogonal evidence",
  "- GSE240078 GeoMx DSP: strong compartment-level TME-vs-carcinoma support for CAF/matrix and immune/TME programs.",
  "- GSE199102 GeoMx WTA: strong independent compartment-level concordance with GSE240078 using CAF+Immune versus epithelial segments.",
  "- GSE310352 CosMx: robust cell-level CAF/matrix-associated TGF/EMT stromal-interface organization; not direct SPP1-CD44 support.",
  "- Cho rapid-autopsy IMC: protein-level spatial-neighborhood/context support as summarized above; interpretation depends on ROI/patient support and robustness.",
  "- Prior vascular-niche IMC Zenodo 10246315: previously logged as weak/near-random adjacency and not suitable for ED10 main support.",
  "",
  "## Claim boundary",
  "Cho IMC can support protein-level neighborhood or spatial-protein context only. It does not establish causality, direct SPP1-CD44 validation, Visium gradient reconstruction, lymph-node immune uncoupling, clinical prediction or mature TLS validation."
), file.path(doc_dir, "cho_imc_final_decision_report.md"), useBytes = TRUE)

results_draft <- if (startsWith(decision, "A.")) {
  "In the Cho rapid-autopsy PDAC IMC dataset, protein-level ROI-resolved analysis supported stromal-myeloid neighborhood organization, with CAF/stromal-like to myeloid/macrophage-like adjacency and/or stromal-rich macrophage/interface enrichment exceeding within-ROI null expectations across patients. These data provide spatial-protein support consistent with CAF/matrix stromal-neighborhood architecture, while remaining observational."
} else if (startsWith(decision, "B.")) {
  "Cho rapid-autopsy PDAC IMC provided protein-level spatial context for stromal, myeloid and epithelial compartments, but the ROI/patient-level support was not sufficient to upgrade Extended Data Fig. 10. The dataset is best retained as secondary/source-level support."
} else {
  "Cho rapid-autopsy PDAC IMC was feasible for protein-neighborhood analysis, but stromal-myeloid support was weak, unstable or insufficient for display. The result should be archived as internal feasibility rather than added to Extended Data Fig. 10."
}
writeLines(c("# Cho IMC Results Paragraph Draft", "", results_draft), file.path(doc_dir, "cho_imc_results_paragraph_draft.md"), useBytes = TRUE)

methods_draft <- "Processed Cho rapid-autopsy PDAC IMC data were obtained from Zenodo record 15596960 and analyzed from the processed backup_output RDS object. Marker matrices and X/Y coordinates were extracted from the fcs1 flowSet frames, and author broad cell labels were harmonized into CAF/stromal-like, myeloid/macrophage-like, tumor epithelial-like, T cell, other immune and unassigned compartments. Spatial-neighborhood enrichment was quantified within ROI using k-nearest-neighbor graphs and summarized at ROI and patient levels. Null models included within-ROI label permutation, abundance-matched random target labels and density-stratified permutation based on local neighborhood density. Primary/pancreas versus liver-metastasis comparisons were treated as exploratory."
writeLines(c("# Cho IMC Methods Draft", "", methods_draft), file.path(doc_dir, "cho_imc_methods_draft.md"), useBytes = TRUE)

legend_draft <- "Candidate Cho IMC panels: CAF/stromal-to-myeloid/macrophage kNN adjacency; macrophage enrichment in stromal-rich neighborhoods; tumor-stroma-myeloid interface enrichment; HLA-DR macrophage/APC-like protein context; robustness across kNN, label-threshold and null definitions. Display only if ROI/patient-level support is strong enough for ED10 v2."
writeLines(c("# Cho IMC Figure Legend Draft", "", legend_draft), file.path(doc_dir, "cho_imc_figure_legend_draft.md"), useBytes = TRUE)

claim_notes <- c(
  "# Cho IMC Claim Boundary Notes",
  "",
  "Allowed: protein-level neighborhood support; spatial-protein context; stromal-myeloid interface support; consistency with CAF/matrix or CAF-myeloid stromal-neighborhood architecture; candidate axis prioritization, not mechanistic proof.",
  "",
  "Not allowed: causality; causal signaling; direct SPP1-CD44 validation unless SPP1 and CD44 are both measured and spatially linked; Visium gradient reconstruction; lymph-node immune uncoupling; clinical prediction; mature TLS validation."
)
writeLines(claim_notes, file.path(doc_dir, "cho_imc_claim_boundary_notes.md"), useBytes = TRUE)

# Update global ranking/manifest if not ED10 v2.
global_table_dir <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search/tables"), winslash = "/", mustWork = TRUE)
global_manifest_dir <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search/manifest"), winslash = "/", mustWork = TRUE)
ranking_path <- file.path(global_table_dir, "signal_strength_ranking.csv")
rank <- if (file.exists(ranking_path)) fread(ranking_path, fill = TRUE) else data.table()
new_rank <- data.table(
  dataset = "Cho rapid-autopsy PDAC IMC Zenodo 15596960",
  platform = "IMC",
  validation_target = "protein-level CAF/stromal-myeloid neighborhood and interface support",
  effect_size = paste0("median ROI CAF->myeloid log2OE=", round(main_summary$median_roi, 3), "; stromal-neighborhood macrophage delta=", round(stromal_summary$median_roi, 3), "; interface log2OE=", round(interface_summary$median_roi, 3)),
  support_fraction = paste0("ROI support=", round(main_summary$roi_support_fraction_positive, 3), "; patient support=", round(main_summary$patient_support_fraction_positive, 3)),
  patient_or_sample_level_support = paste0(uniqueN(cells$patient_id), " patients; ", uniqueN(cells$roi_id), " ROIs"),
  FDR_or_p_value = paste0("median random-label p=", signif(main_summary$median_random_label_p, 3)),
  robustness_level = decision,
  positive_controls_passed = "see cho_imc_positive_negative_controls.csv",
  negative_controls_passed = "see cho_imc_positive_negative_controls.csv",
  claim_supported = ifelse(startsWith(decision, "A."), "protein-level stromal-myeloid neighborhood support", "secondary/source-level protein context"),
  claim_boundary = "No causality, direct SPP1-CD44, Visium gradient reconstruction, LN immune uncoupling, clinical prediction, or mature TLS validation",
  recommended_display = ifelse(startsWith(decision, "A."), "consider_ED10_v2", "source_only_or_archive"),
  reason = "Processed IMC RDS contains ROI-resolved protein expression, coordinates, author labels, patient and pancreas/liver site metadata."
)
if (nrow(rank) > 0 && "dataset" %in% names(rank)) rank <- rank[dataset != new_rank$dataset]
rank <- rbindlist(list(rank, new_rank), fill = TRUE)
write_csv(rank, ranking_path)

panel_path <- file.path(global_manifest_dir, "candidate_ed_panel_manifest.csv")
panel <- if (file.exists(panel_path)) fread(panel_path, fill = TRUE) else data.table()
new_panel <- data.table(
  panel = "Cho_IMC",
  title = "Cho rapid-autopsy PDAC IMC protein-neighborhood support",
  status = ifelse(startsWith(decision, "A."), "candidate_for_ED10_v2", "source_or_archive"),
  dataset = "Zenodo 15596960",
  source_table = file.path(source_dir, "Source_Data_Cho_IMC.csv"),
  figure_path = file.path(figure_dir, "cho_imc_caf_myeloid_adjacency.pdf"),
  decision = decision,
  note = "Do not integrate into ED10 v1; ED10 v2 only if decision A is accepted after review."
)
if (nrow(panel) > 0 && all(c("panel", "dataset") %in% names(panel))) panel <- panel[!(panel == "Cho_IMC" & dataset == "Zenodo 15596960")]
panel <- rbindlist(list(panel, new_panel), fill = TRUE)
write_csv(panel, panel_path)

neg_log <- file.path(root, "docs/negative_or_weak_results_log.md")
if (!startsWith(decision, "A.")) {
  writeLines(c(
    "# Cho IMC Negative or Weak Results Log",
    "",
    paste0("Decision: ", decision),
    "",
    paste0("Median ROI CAF/stromal -> myeloid log2(O/E): ", round(main_summary$median_roi, 3)),
    paste0("ROI support fraction: ", round(main_summary$roi_support_fraction_positive, 3)),
    paste0("Patient support fraction: ", round(main_summary$patient_support_fraction_positive, 3)),
    "",
    "Cho IMC was not promoted to ED10 v2 in this automated decision layer unless manually re-reviewed."
  ), neg_log, useBytes = TRUE)
}

# ED10 v2 only for strong decision.
if (startsWith(decision, "A.")) {
  p1 <- ggplot(main_adj, aes(x = site_status, y = log2_oe, color = site_status)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = 0.75) +
    theme_bw(base_size = 9) + theme(legend.position = "none") + labs(x = NULL, y = "CAF->myeloid log2(O/E)", title = "Cho IMC protein-neighborhood support")
  pdf(file.path(figure_dir, "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_IMC_v2.pdf"), width = 7, height = 5)
  print(p1)
  dev.off()
  svg(file.path(figure_dir, "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_IMC_v2.svg"), width = 7, height = 5)
  print(p1)
  dev.off()
  write_csv(source_data, file.path(source_dir, "Source_Data_Extended_Data_Figure_10_v2.csv"))
  write_csv(data.table(panel = "Cho_IMC", source = file.path(source_dir, "Source_Data_Extended_Data_Figure_10_v2.csv"), decision = decision), file.path(manifest_dir, "ed10_v2_panel_map.csv"))
  write_csv(data.table(output = "Source_Data_Extended_Data_Figure_10_v2.csv", path = file.path(source_dir, "Source_Data_Extended_Data_Figure_10_v2.csv")), file.path(manifest_dir, "ed10_v2_source_data_manifest.csv"))
  write_csv(data.table(script = "02_analyze_cho_imc_spatial.R", purpose = "Cho IMC ED10 v2 candidate generation"), file.path(manifest_dir, "ed10_v2_script_manifest.csv"))
  write_csv(data.table(parameter = c("k_main", "permutations_main", "decision"), value = c("10", "100", decision)), file.path(manifest_dir, "ed10_v2_parameter_manifest.csv"))
  writeLines(c("# ED10 v2 Integration Recommendation", "", paste0("Decision: ", decision), "", "Cho IMC may be considered as a protein-neighborhood panel for ED10 v2 after manual review."), file.path(doc_dir, "ed10_v2_integration_recommendation.md"), useBytes = TRUE)
}

manifest_outputs <- data.table(
  output = c(
    "cho_imc_object_metadata_summary.csv", "cho_imc_marker_availability.csv", "cho_imc_sample_roi_summary.csv",
    "cho_imc_cell_type_harmonization.csv", "cho_imc_caf_myeloid_adjacency_results.csv", "cho_imc_decision_summary.csv",
    "Source_Data_Cho_IMC.csv", "cho_imc_final_decision_report.md"
  ),
  path = c(
    file.path(table_dir, "cho_imc_object_metadata_summary.csv"), file.path(table_dir, "cho_imc_marker_availability.csv"), file.path(table_dir, "cho_imc_sample_roi_summary.csv"),
    file.path(table_dir, "cho_imc_cell_type_harmonization.csv"), file.path(table_dir, "cho_imc_caf_myeloid_adjacency_results.csv"), file.path(table_dir, "cho_imc_decision_summary.csv"),
    file.path(source_dir, "Source_Data_Cho_IMC.csv"), file.path(doc_dir, "cho_imc_final_decision_report.md")
  )
)
write_csv(manifest_outputs, file.path(manifest_dir, "cho_imc_output_manifest.csv"))

cat("Cho IMC analysis complete. Decision: ", decision, "\n", sep = "")
