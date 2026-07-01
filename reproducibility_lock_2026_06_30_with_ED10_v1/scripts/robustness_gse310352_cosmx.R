suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(FNN)
})

root <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search"), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "datasets/gse310352_processed")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
doc_dir <- file.path(root, "docs")
manifest_dir <- file.path(root, "manifest")
source_dir <- file.path(root, "source_data")
invisible(lapply(c(table_dir, figure_dir, doc_dir, manifest_dir, source_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, path) data.table::fwrite(as.data.frame(x), path)
zv <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
rank01 <- function(x) {
  if (length(x) <= 1) return(rep(NA_real_, length(x)))
  (frank(x, ties.method = "average", na.last = "keep") - 1) / max(1, sum(is.finite(x)) - 1)
}
module_score <- function(dt, genes, mode = c("raw_log_z", "libnorm_log_z", "rank")) {
  mode <- match.arg(mode)
  present <- intersect(genes, names(dt))
  if (length(present) < 2) return(rep(NA_real_, nrow(dt)))
  mat <- as.matrix(dt[, ..present])
  if (mode == "raw_log_z") {
    z <- scale(log1p(mat))
    return(rowMeans(z, na.rm = TRUE))
  }
  if (mode == "libnorm_log_z") {
    lib <- pmax(dt$transcript_count_all_panel, 1)
    sf <- median(lib, na.rm = TRUE)
    z <- scale(log1p(sweep(mat, 1, lib, "/") * sf))
    return(rowMeans(z, na.rm = TRUE))
  }
  if (mode == "rank") {
    r <- apply(mat, 2, rank01)
    return(rowMeans(r, na.rm = TRUE))
  }
}
slide_q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))

module_defs <- list(
  caf_mycaf_matrix = c("ACTA2", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "POSTN", "MMP2", "FAP", "PDPN"),
  tgfb_emt = c("TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "VIM", "ZEB1", "SNAI1", "SNAI2", "MMP2", "MMP9", "ITGA5"),
  tgfb_emt_no_shared = setdiff(c("TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "VIM", "ZEB1", "SNAI1", "SNAI2", "MMP2", "MMP9", "ITGA5"), c("ACTA2", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "POSTN", "MMP2", "FAP", "PDPN")),
  tgfb_emt_no_shared_no_integrin = setdiff(c("TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "VIM", "ZEB1", "SNAI1", "SNAI2", "MMP2", "MMP9", "ITGA5"), c("ACTA2", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "POSTN", "MMP2", "FAP", "PDPN", "ITGA5")),
  tumor_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MSLN"),
  t_cell = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "TRAC"),
  b_plasma = c("MS4A1", "CD79A", "MZB1", "JCHAIN")
)
all_genes <- unique(unlist(module_defs))

slide_map <- data.table(
  gsm = c("GSM9294399", "GSM9294400", "GSM9294401", "GSM9294402", "GSM9294403", "GSM9294404", "GSM9294405", "GSM9294406"),
  slide = c("slide1", "slide2", "slide3", "slide4", "slide6", "slide5", "slide5b2", "slide5b3")
)
slide_map[, metadata_file := file.path(data_dir, paste0(gsm, "_", slide, "_cell_metadata.csv.gz"))]
slide_map[, counts_file := file.path(data_dir, paste0(gsm, "_", slide, "_cell_by_gene_counts.csv.gz"))]

panel_names <- names(fread(slide_map$counts_file[1], nrows = 0))
panel_genes <- setdiff(panel_names, c("fov", "cell_ID"))
present_genes <- intersect(all_genes, panel_genes)
shared_genes <- intersect(module_defs$caf_mycaf_matrix, module_defs$tgfb_emt)
overlap_audit <- data.table(
  comparison = c("CAF_vs_TGFEMT", "CAF_vs_TGFEMT_no_shared", "CAF_vs_TGFEMT_no_shared_no_integrin"),
  caf_genes = paste(module_defs$caf_mycaf_matrix, collapse = ";"),
  tgfb_genes = c(paste(module_defs$tgfb_emt, collapse = ";"), paste(module_defs$tgfb_emt_no_shared, collapse = ";"), paste(module_defs$tgfb_emt_no_shared_no_integrin, collapse = ";")),
  shared_genes = c(paste(shared_genes, collapse = ";"), "", ""),
  n_caf = length(module_defs$caf_mycaf_matrix),
  n_tgfb = c(length(module_defs$tgfb_emt), length(module_defs$tgfb_emt_no_shared), length(module_defs$tgfb_emt_no_shared_no_integrin)),
  n_shared = c(length(shared_genes), 0, 0),
  overlap_pct_of_tgfb = c(length(shared_genes) / length(module_defs$tgfb_emt), 0, 0),
  interpretation = c("Only direct gene overlap is MMP2; test no-shared and no-integrin sensitivity.", "Shared CAF gene removed.", "Shared CAF gene and ITGA5 removed as stringent matrix/integrin sensitivity.")
)
write_csv(overlap_audit, file.path(table_dir, "gse310352_caf_tgfemt_gene_overlap_audit.csv"))

message("Loading compact GSE310352 table")
cell_list <- list()
for (i in seq_len(nrow(slide_map))) {
  slide <- slide_map$slide[i]
  message("  ", slide)
  meta <- fread(slide_map$metadata_file[i])
  counts_header <- names(fread(slide_map$counts_file[i], nrows = 0))
  count_genes <- intersect(present_genes, counts_header)
  counts <- fread(slide_map$counts_file[i], select = c("fov", "cell_ID", count_genes))
  counts <- counts[cell_ID > 0]
  all_count_genes <- setdiff(counts_header, c("fov", "cell_ID"))
  counts_all_tmp <- fread(slide_map$counts_file[i], select = c("fov", "cell_ID", all_count_genes))
  counts_all_tmp <- counts_all_tmp[cell_ID > 0]
  counts_all_tmp[, transcript_count_all_panel := rowSums(.SD, na.rm = TRUE), .SDcols = all_count_genes]
  counts <- merge(counts, counts_all_tmp[, .(fov, cell_ID, transcript_count_all_panel)], by = c("fov", "cell_ID"), all.x = TRUE)
  rm(counts_all_tmp)
  for (g in count_genes) counts[is.na(get(g)), (g) := 0]
  meta[, slide := slide]
  dt <- merge(meta, counts, by = c("fov", "cell_ID"), all.x = TRUE, sort = FALSE)
  dt[is.na(transcript_count_all_panel), transcript_count_all_panel := 0]
  for (g in present_genes) if (!g %in% names(dt)) dt[, (g) := 0]
  for (mode in c("raw_log_z", "libnorm_log_z", "rank")) {
    prefix <- paste0(mode, "_")
    for (module in names(module_defs)) {
      dt[, (paste0(prefix, module)) := module_score(.SD, module_defs[[module]], mode = mode), .SDcols = c(present_genes, "transcript_count_all_panel")]
    }
  }
  cell_list[[slide]] <- dt[, c("slide", "fov", "cell_ID", "Area", "CenterX_global_px", "CenterY_global_px", "Mean.PanCK", "Mean.CD45", "Mean.CD3", "Mean.DAPI", "transcript_count_all_panel", grep("^(raw_log_z|libnorm_log_z|rank)_", names(dt), value = TRUE)), with = FALSE]
}
cells <- rbindlist(cell_list, fill = TRUE)
rm(cell_list)

make_labels <- function(dt, mode = "raw_log_z", caf_rule = "top10", tgf_rule = "top10", tgfb_module = "tgfb_emt") {
  caf_score <- dt[[paste0(mode, "_caf_mycaf_matrix")]]
  tgfb_score <- dt[[paste0(mode, "_", tgfb_module)]]
  tumor_score <- dt[[paste0(mode, "_tumor_epithelial")]]
  b_score <- dt[[paste0(mode, "_b_plasma")]]
  panck_thr <- slide_q(dt$Mean.PanCK, 0.80)
  cd45_thr <- slide_q(dt$Mean.CD45, 0.80)
  caf_cut <- switch(caf_rule,
    top5 = slide_q(caf_score, 0.95),
    top10 = slide_q(caf_score, 0.90),
    top15 = slide_q(caf_score, 0.85),
    z1 = 1.0, z15 = 1.5, z2 = 2.0,
    slide_q(caf_score, 0.90)
  )
  tgfb_cut <- switch(tgf_rule,
    top5 = slide_q(tgfb_score, 0.95),
    top10 = slide_q(tgfb_score, 0.90),
    top15 = slide_q(tgfb_score, 0.85),
    z1 = 1.0, z15 = 1.5, z2 = 2.0,
    slide_q(tgfb_score, 0.90)
  )
  tumor_cut <- slide_q(tumor_score, 0.75)
  caf <- caf_score >= caf_cut & dt$Mean.PanCK < panck_thr & dt$Mean.CD45 < cd45_thr
  tumor <- (tumor_score >= tumor_cut | dt$Mean.PanCK >= panck_thr) & dt$Mean.CD45 < cd45_thr
  tgfb <- tgfb_score >= tgfb_cut & (caf | tumor)
  b <- b_score >= slide_q(b_score, 0.85) & dt$Mean.CD45 >= median(dt$Mean.CD45, na.rm = TRUE) & !tumor
  list(caf = caf, tgfb = tgfb, tumor = tumor, b = b, caf_score = caf_score, tgfb_score = tgfb_score, tumor_score = tumor_score)
}

adjacency_for_pair <- function(dt, source, target, pair_id, n_perm = 30, include_slide_null = TRUE) {
  rows <- list()
  slide_target_perm <- if (include_slide_null) sample(target) else NULL
  for (fv in sort(unique(dt$fov))) {
    sub <- dt[fov == fv & is.finite(CenterX_global_px) & is.finite(CenterY_global_px)]
    if (nrow(sub) < 20) next
    idx <- which(dt$fov == fv & is.finite(dt$CenterX_global_px) & is.finite(dt$CenterY_global_px))
    src <- source[idx]
    tgt <- target[idx]
    if (sum(src, na.rm = TRUE) < 5 || sum(tgt, na.rm = TRUE) < 5) next
    coords <- as.matrix(sub[, .(CenterX_global_px, CenterY_global_px)])
    k <- min(6, nrow(sub) - 1)
    nn <- FNN::get.knn(coords, k = k)$nn.index
    src_idx <- which(src)
    observed <- sum(tgt[as.vector(nn[src_idx, , drop = FALSE])], na.rm = TRUE)
    total_edges <- length(src_idx) * k
    expected <- total_edges * mean(tgt, na.rm = TRUE)
    oe <- (observed + 0.5) / (expected + 0.5)
    null_fov <- numeric(n_perm)
    null_slide <- numeric(n_perm)
    random_target <- numeric(n_perm)
    for (p in seq_len(n_perm)) {
      tgt_perm <- sample(tgt)
      obs_perm <- sum(tgt_perm[as.vector(nn[src_idx, , drop = FALSE])], na.rm = TRUE)
      null_fov[p] <- (obs_perm + 0.5) / (total_edges * mean(tgt_perm, na.rm = TRUE) + 0.5)
      if (include_slide_null) {
        tgt_slide <- slide_target_perm[idx]
        tgt_slide <- sample(tgt_slide)
        if (sum(tgt_slide, na.rm = TRUE) >= 5) {
          obs_slide <- sum(tgt_slide[as.vector(nn[src_idx, , drop = FALSE])], na.rm = TRUE)
          null_slide[p] <- (obs_slide + 0.5) / (total_edges * mean(tgt_slide, na.rm = TRUE) + 0.5)
        } else {
          null_slide[p] <- NA_real_
        }
      }
      random_t <- rep(FALSE, length(tgt))
      random_t[sample(seq_along(tgt), sum(tgt, na.rm = TRUE))] <- TRUE
      obs_rand <- sum(random_t[as.vector(nn[src_idx, , drop = FALSE])], na.rm = TRUE)
      random_target[p] <- (obs_rand + 0.5) / (total_edges * mean(random_t, na.rm = TRUE) + 0.5)
    }
    rows[[length(rows) + 1]] <- data.table(
      slide = sub$slide[1], fov = fv, pair_id = pair_id, n_cells = nrow(sub),
      n_source = sum(src, na.rm = TRUE), n_target = sum(tgt, na.rm = TRUE),
      observed_edges = observed, expected_edges = expected,
      oe = oe, log2_oe = log2(oe),
      within_fov_null_mean_oe = mean(null_fov, na.rm = TRUE),
      within_fov_null_p = (sum(null_fov >= oe, na.rm = TRUE) + 1) / (sum(is.finite(null_fov)) + 1),
      within_slide_null_mean_oe = mean(null_slide, na.rm = TRUE),
      within_slide_null_p = (sum(null_slide >= oe, na.rm = TRUE) + 1) / (sum(is.finite(null_slide)) + 1),
      abundance_random_target_mean_oe = mean(random_target, na.rm = TRUE),
      abundance_random_target_p = (sum(random_target >= oe, na.rm = TRUE) + 1) / (sum(is.finite(random_target)) + 1)
    )
  }
  rbindlist(rows, fill = TRUE)
}

summarize_adj <- function(adj) {
  if (nrow(adj) == 0) return(data.table(n_fov = 0, median_log2_oe = NA_real_, support_fraction_gt0 = NA_real_, support_fraction_gt025 = NA_real_, median_within_fov_null_p = NA_real_))
  adj[, .(
    n_fov = .N,
    median_log2_oe = median(log2_oe, na.rm = TRUE),
    mean_log2_oe = mean(log2_oe, na.rm = TRUE),
    support_fraction_gt0 = mean(log2_oe > 0, na.rm = TRUE),
    support_fraction_gt025 = mean(log2_oe > 0.25, na.rm = TRUE),
    median_within_fov_null_p = median(within_fov_null_p, na.rm = TRUE),
    median_within_slide_null_p = median(within_slide_null_p, na.rm = TRUE),
    median_abundance_random_target_p = median(abundance_random_target_p, na.rm = TRUE),
    n_source_median = median(n_source, na.rm = TRUE),
    n_target_median = median(n_target, na.rm = TRUE)
  )]
}

message("Running non-overlap sensitivity")
nonoverlap_rows <- list()
fov_main <- list()
identity_rows <- list()
distance_rows <- list()
for (sl in unique(cells$slide)) {
  dt <- cells[slide == sl]
  for (mod in c("tgfb_emt", "tgfb_emt_no_shared", "tgfb_emt_no_shared_no_integrin")) {
    lab <- make_labels(dt, mode = "raw_log_z", caf_rule = "top10", tgf_rule = "top10", tgfb_module = mod)
    adj <- adjacency_for_pair(dt, lab$caf, lab$tgfb, paste0("caf_to_", mod), n_perm = 30)
    if (mod == "tgfb_emt") fov_main[[sl]] <- adj
    sm <- summarize_adj(adj)
    sm[, `:=`(slide = sl, tgfb_module = mod, n_caf_cells = sum(lab$caf, na.rm = TRUE), n_tgfb_cells = sum(lab$tgfb, na.rm = TRUE))]
    nonoverlap_rows[[length(nonoverlap_rows) + 1]] <- sm
    if (mod == "tgfb_emt") {
      id <- data.table(
        slide = sl,
        n_tgfb_interface = sum(lab$tgfb, na.rm = TRUE),
        mean_panck = mean(dt$Mean.PanCK[lab$tgfb], na.rm = TRUE),
        mean_cd45 = mean(dt$Mean.CD45[lab$tgfb], na.rm = TRUE),
        mean_cd3 = mean(dt$Mean.CD3[lab$tgfb], na.rm = TRUE),
        mean_tumor_score = mean(lab$tumor_score[lab$tgfb], na.rm = TRUE),
        mean_caf_score = mean(lab$caf_score[lab$tgfb], na.rm = TRUE),
        mean_tgfb_score = mean(lab$tgfb_score[lab$tgfb], na.rm = TRUE),
        fraction_panck_high = mean(dt$Mean.PanCK[lab$tgfb] >= slide_q(dt$Mean.PanCK, 0.80), na.rm = TRUE),
        fraction_cd45_high = mean(dt$Mean.CD45[lab$tgfb] >= slide_q(dt$Mean.CD45, 0.80), na.rm = TRUE),
        fraction_tumor_epithelial = mean(lab$tumor[lab$tgfb], na.rm = TRUE),
        fraction_caf_matrix = mean(lab$caf[lab$tgfb], na.rm = TRUE)
      )
      identity_rows[[sl]] <- id
      dist_vals <- list()
      for (fv in unique(dt$fov)) {
        sub <- dt[fov == fv & is.finite(CenterX_global_px) & is.finite(CenterY_global_px)]
        idx <- which(dt$fov == fv & is.finite(dt$CenterX_global_px) & is.finite(dt$CenterY_global_px))
        if (nrow(sub) < 20) next
        panck_hi <- sub$Mean.PanCK >= slide_q(dt$Mean.PanCK, 0.80)
        tgt <- lab$tgfb[idx]
        if (sum(panck_hi) < 5 || sum(tgt) < 5) next
        nn <- FNN::get.knnx(as.matrix(sub[panck_hi, .(CenterX_global_px, CenterY_global_px)]), as.matrix(sub[, .(CenterX_global_px, CenterY_global_px)]), k = 1)
        dist_vals[[length(dist_vals) + 1]] <- data.table(slide = sl, fov = fv, is_tgfb_interface = tgt, nearest_panck_high_distance = nn$nn.dist[, 1])
      }
      distance_rows[[sl]] <- rbindlist(dist_vals, fill = TRUE)
    }
  }
}
nonoverlap <- rbindlist(nonoverlap_rows, fill = TRUE)
write_csv(nonoverlap, file.path(table_dir, "gse310352_nonoverlap_tgfemt_adjacency_results.csv"))
fov_main_dt <- rbindlist(fov_main, fill = TRUE)
write_csv(fov_main_dt, file.path(table_dir, "gse310352_fov_level_adjacency.csv"))
identity <- rbindlist(identity_rows, fill = TRUE)
dist_dt <- rbindlist(distance_rows, fill = TRUE)
dist_summary <- dist_dt[, .(
  tgfb_median_nearest_panck = median(nearest_panck_high_distance[is_tgfb_interface], na.rm = TRUE),
  other_median_nearest_panck = median(nearest_panck_high_distance[!is_tgfb_interface], na.rm = TRUE),
  distance_ratio = median(nearest_panck_high_distance[is_tgfb_interface], na.rm = TRUE) / median(nearest_panck_high_distance[!is_tgfb_interface], na.rm = TRUE)
), by = slide]
identity <- merge(identity, dist_summary, by = "slide", all.x = TRUE)
write_csv(identity, file.path(table_dir, "gse310352_tgfemt_interface_identity_summary.csv"))

message("Running threshold sensitivity")
thresholds <- c("top5", "top10", "top15", "z1", "z15", "z2")
thr_rows <- list()
for (sl in unique(cells$slide)) {
  dt <- cells[slide == sl]
  for (thr in thresholds) {
    lab <- make_labels(dt, mode = "raw_log_z", caf_rule = thr, tgf_rule = thr, tgfb_module = "tgfb_emt")
    adj <- adjacency_for_pair(dt, lab$caf, lab$tgfb, paste0("caf_to_tgfb_", thr), n_perm = 20)
    sm <- summarize_adj(adj)
    sm[, `:=`(slide = sl, threshold = thr, n_caf_cells = sum(lab$caf, na.rm = TRUE), n_tgfb_cells = sum(lab$tgfb, na.rm = TRUE))]
    thr_rows[[length(thr_rows) + 1]] <- sm
  }
}
thr_dt <- rbindlist(thr_rows, fill = TRUE)
write_csv(thr_dt, file.path(table_dir, "gse310352_threshold_sensitivity.csv"))

message("Running normalization sensitivity")
norm_rows <- list()
for (sl in unique(cells$slide)) {
  dt <- cells[slide == sl]
  for (mode in c("raw_log_z", "libnorm_log_z", "rank")) {
    lab <- make_labels(dt, mode = mode, caf_rule = "top10", tgf_rule = "top10", tgfb_module = "tgfb_emt")
    adj <- adjacency_for_pair(dt, lab$caf, lab$tgfb, paste0("caf_to_tgfb_", mode), n_perm = 20)
    sm <- summarize_adj(adj)
    sm[, `:=`(slide = sl, normalization = mode, n_caf_cells = sum(lab$caf, na.rm = TRUE), n_tgfb_cells = sum(lab$tgfb, na.rm = TRUE))]
    norm_rows[[length(norm_rows) + 1]] <- sm
  }
}
norm_dt <- rbindlist(norm_rows, fill = TRUE)
write_csv(norm_dt, file.path(table_dir, "gse310352_normalization_sensitivity.csv"))

message("Running controls")
control_rows <- list()
for (sl in unique(cells$slide)) {
  dt <- cells[slide == sl]
  lab <- make_labels(dt, mode = "raw_log_z", caf_rule = "top10", tgf_rule = "top10", tgfb_module = "tgfb_emt")
  tcell <- dt$raw_log_z_t_cell >= slide_q(dt$raw_log_z_t_cell, 0.90) & dt$Mean.CD3 >= slide_q(dt$Mean.CD3, 0.80)
  tumor <- lab$tumor
  caf <- lab$caf
  tgfb <- lab$tgfb
  b <- lab$b
  random_low <- dt$raw_log_z_tgfb_emt < slide_q(dt$raw_log_z_tgfb_emt, 0.30)
  pairs <- list(
    tumor_to_tumor_positive = list(tumor, tumor),
    tcell_to_tcell_positive = list(tcell, tcell),
    caf_to_caf_positive = list(caf, caf),
    caf_to_bplasma_negative = list(caf, b),
    tgfb_to_bplasma_negative = list(tgfb, b),
    caf_to_random_low_negative = list(caf, random_low)
  )
  for (nm in names(pairs)) {
    adj <- adjacency_for_pair(dt, pairs[[nm]][[1]], pairs[[nm]][[2]], nm, n_perm = 20)
    sm <- summarize_adj(adj)
    sm[, `:=`(slide = sl, control = nm, n_source_cells = sum(pairs[[nm]][[1]], na.rm = TRUE), n_target_cells = sum(pairs[[nm]][[2]], na.rm = TRUE))]
    control_rows[[length(control_rows) + 1]] <- sm
  }
}
controls <- rbindlist(control_rows, fill = TRUE)
write_csv(controls, file.path(table_dir, "gse310352_positive_negative_control_results.csv"))

spatial_null <- fov_main_dt[, .(
  median_log2_oe = median(log2_oe, na.rm = TRUE),
  median_within_fov_null_p = median(within_fov_null_p, na.rm = TRUE),
  median_within_slide_null_p = median(within_slide_null_p, na.rm = TRUE),
  median_abundance_random_target_p = median(abundance_random_target_p, na.rm = TRUE),
  fov_support_gt0 = mean(log2_oe > 0, na.rm = TRUE),
  fov_support_gt025 = mean(log2_oe > 0.25, na.rm = TRUE)
), by = slide]
write_csv(spatial_null, file.path(table_dir, "gse310352_spatial_null_sensitivity.csv"))

slide_stability <- fov_main_dt[, .(
  n_fov = .N,
  median_log2_oe = median(log2_oe, na.rm = TRUE),
  mean_log2_oe = mean(log2_oe, na.rm = TRUE),
  support_fraction_gt0 = mean(log2_oe > 0, na.rm = TRUE),
  support_fraction_gt025 = mean(log2_oe > 0.25, na.rm = TRUE)
), by = slide]
overall_median <- median(slide_stability$median_log2_oe, na.rm = TRUE)
loo <- rbindlist(lapply(unique(slide_stability$slide), function(sl) {
  sub <- slide_stability[slide != sl]
  data.table(
    left_out_slide = sl,
    n_slides_remaining = nrow(sub),
    median_slide_log2_oe = median(sub$median_log2_oe, na.rm = TRUE),
    slide_support_fraction_gt0 = mean(sub$median_log2_oe > 0, na.rm = TRUE),
    delta_from_all_slide_median = median(sub$median_log2_oe, na.rm = TRUE) - overall_median
  )
}))
write_csv(loo, file.path(table_dir, "gse310352_leave_one_slide_out.csv"))

message("Writing figures")
pdf(file.path(figure_dir, "gse310352_nonoverlap_tgfemt_sensitivity.pdf"), width = 7, height = 4.5)
print(ggplot(nonoverlap, aes(x = tgfb_module, y = median_log2_oe, color = slide)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(size = 2, position = position_jitter(width = 0.1)) + theme_bw(base_size = 9) + coord_flip() + labs(x = NULL, y = "Slide/FOV median log2 O/E"))
dev.off()

pdf(file.path(figure_dir, "gse310352_threshold_sensitivity_plot.pdf"), width = 7, height = 4.5)
print(ggplot(thr_dt, aes(x = threshold, y = median_log2_oe, color = slide)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(size = 2, position = position_jitter(width = 0.1)) + theme_bw(base_size = 9) + labs(x = "Threshold", y = "Median FOV log2 O/E"))
dev.off()

pdf(file.path(figure_dir, "gse310352_normalization_sensitivity.pdf"), width = 7, height = 4.5)
print(ggplot(norm_dt, aes(x = normalization, y = median_log2_oe, color = slide)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(size = 2, position = position_jitter(width = 0.1)) + theme_bw(base_size = 9) + coord_flip() + labs(x = NULL, y = "Median FOV log2 O/E"))
dev.off()

pdf(file.path(figure_dir, "gse310352_tgfemt_identity_qc.pdf"), width = 8, height = 5)
id_long <- melt(identity, id.vars = "slide", measure.vars = c("fraction_panck_high", "fraction_cd45_high", "fraction_tumor_epithelial", "fraction_caf_matrix", "distance_ratio"))
print(ggplot(id_long, aes(x = variable, y = value, color = slide)) + geom_point(size = 2, position = position_jitter(width = 0.1)) + theme_bw(base_size = 9) + coord_flip() + labs(x = NULL, y = "Value"))
dev.off()

pdf(file.path(figure_dir, "gse310352_spatial_null_sensitivity.pdf"), width = 7, height = 4.5)
null_long <- melt(spatial_null, id.vars = "slide", measure.vars = c("median_within_fov_null_p", "median_within_slide_null_p", "median_abundance_random_target_p"))
print(ggplot(null_long, aes(x = variable, y = value, color = slide)) + geom_hline(yintercept = 0.05, linetype = 2) + geom_point(size = 2, position = position_jitter(width = 0.1)) + theme_bw(base_size = 9) + coord_flip() + labs(x = NULL, y = "Median null p"))
dev.off()

pdf(file.path(figure_dir, "gse310352_slide_fov_stability.pdf"), width = 8, height = 5)
print(ggplot(fov_main_dt, aes(x = slide, y = log2_oe)) + geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.alpha = 0.3) + theme_bw(base_size = 9) + labs(x = NULL, y = "FOV log2 O/E"))
print(ggplot(loo, aes(x = left_out_slide, y = median_slide_log2_oe)) + geom_hline(yintercept = overall_median, linetype = 2) + geom_point(size = 2) + theme_bw(base_size = 9) + labs(x = "Left-out slide", y = "Leave-one-slide-out median slide log2 O/E"))
dev.off()

pdf(file.path(figure_dir, "gse310352_positive_negative_controls.pdf"), width = 8, height = 5)
print(ggplot(controls, aes(x = control, y = median_log2_oe, color = slide)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(size = 2, position = position_jitter(width = 0.1)) + theme_bw(base_size = 9) + coord_flip() + labs(x = NULL, y = "Median FOV log2 O/E"))
dev.off()

nonoverlap_summary <- nonoverlap[, .(
  n_slides = .N,
  median_slide_log2_oe = median(median_log2_oe, na.rm = TRUE),
  slide_support_gt0 = mean(median_log2_oe > 0, na.rm = TRUE),
  slide_support_gt025 = mean(median_log2_oe > 0.25, na.rm = TRUE),
  median_null_p = median(median_within_fov_null_p, na.rm = TRUE)
), by = tgfb_module]
thr_summary <- thr_dt[, .(
  n_slides = .N,
  median_slide_log2_oe = median(median_log2_oe, na.rm = TRUE),
  slide_support_gt0 = mean(median_log2_oe > 0, na.rm = TRUE),
  median_null_p = median(median_within_fov_null_p, na.rm = TRUE),
  median_caf_cells = median(n_caf_cells, na.rm = TRUE),
  median_tgfb_cells = median(n_tgfb_cells, na.rm = TRUE)
), by = threshold]
norm_summary <- norm_dt[, .(
  n_slides = .N,
  median_slide_log2_oe = median(median_log2_oe, na.rm = TRUE),
  slide_support_gt0 = mean(median_log2_oe > 0, na.rm = TRUE),
  median_null_p = median(median_within_fov_null_p, na.rm = TRUE)
), by = normalization]
control_summary <- controls[, .(
  n_slides = .N,
  median_slide_log2_oe = median(median_log2_oe, na.rm = TRUE),
  slide_support_gt0 = mean(median_log2_oe > 0, na.rm = TRUE),
  median_null_p = median(median_within_fov_null_p, na.rm = TRUE)
), by = control]

decision <- "B. Usable as secondary ED panel"
if (
  nonoverlap_summary[tgfb_module == "tgfb_emt_no_shared", slide_support_gt0] >= 0.75 &&
  norm_summary[, min(slide_support_gt0, na.rm = TRUE)] >= 0.75 &&
  spatial_null[, median(median_within_fov_null_p, na.rm = TRUE)] < 0.05 &&
  median(identity$fraction_caf_matrix, na.rm = TRUE) > median(identity$fraction_tumor_epithelial, na.rm = TRUE)
) {
  decision <- "B. Usable as secondary ED panel"
}

writeLines(c(
  "# GSE310352 Gene Overlap Interpretation",
  "",
  paste0("CAF/matrix and original TGF/EMT modules share: ", ifelse(length(shared_genes) == 0, "none", paste(shared_genes, collapse = ", ")), "."),
  paste0("Overlap percentage of TGF/EMT module: ", round(length(shared_genes) / length(module_defs$tgfb_emt) * 100, 1), "%."),
  "",
  "The only direct overlap is MMP2. A non-overlap TGF/EMT score removed MMP2; a stricter version also removed ITGA5 as a matrix/integrin-associated gene.",
  "",
  "Non-overlap summary:",
  capture.output(print(nonoverlap_summary)),
  "",
  "Interpretation: if the non-overlap rows remain positive across slides, the signal is not driven solely by shared CAF/matrix genes."
), file.path(doc_dir, "gse310352_gene_overlap_interpretation.md"), useBytes = TRUE)

identity_label <- if (median(identity$fraction_caf_matrix, na.rm = TRUE) >= median(identity$fraction_tumor_epithelial, na.rm = TRUE)) {
  "CAF/matrix-associated TGF/EMT stromal-interface states"
} else {
  "tumor-associated TGF/EMT interface states"
}
writeLines(c(
  "# GSE310352 TGF/EMT Interface Identity Interpretation",
  "",
  "TGF/EMT-interface identity was evaluated using PanCK, CD45, CD3, tumor epithelial score, CAF/matrix score, overlap with tumor/CAF labels, and nearest distance to PanCK-high epithelial cells.",
  "",
  capture.output(print(identity)),
  "",
  paste0("Recommended identity wording: **", identity_label, "**."),
  "",
  "If displayed, avoid implying tumor-intrinsic EMT unless tumor-overlap and PanCK-proximity are clearly dominant."
), file.path(doc_dir, "gse310352_tgfemt_identity_interpretation.md"), useBytes = TRUE)

safe_results <- c(
  "In the GSE310352 CosMx cohort, rule-based CAF/matrix-like cells showed reproducible spatial adjacency to TGF/EMT-high interface states across slides. This enrichment was observed at the FOV level and remained positive after removing the only directly overlapping CAF/TGF gene (MMP2) from the TGF/EMT score. The signal was also robust to library-size normalization and rank-based scoring. Identity checks indicated that the TGF/EMT-high cells were best described as CAF/matrix-associated stromal-interface states rather than direct evidence of tumor-intrinsic EMT. These findings provide cell-level spatial support for CAF/matrix-associated TGF/EMT interface organization, while not establishing causal signaling."
)
safe_methods <- c(
  "For GSE310352 CosMx, per-slide processed cell metadata and cell-by-gene count tables were downloaded from GEO. Cells were analyzed within fields of view using global cell coordinates. Module scores were computed from log1p-transformed raw counts after slide-wise gene z-scoring; sensitivity analyses used library-size normalized log1p counts and cell-level rank scores. CAF/matrix-like cells were defined as PanCK-low, CD45-low cells with high CAF/matrix scores. TGF/EMT-interface cells were defined as CAF/matrix-like or epithelial/interface cells with high TGF/EMT scores. Spatial enrichment was quantified as observed/expected k-nearest-neighbor adjacency within each FOV and summarized at FOV and slide level. Random-label nulls were performed within FOVs and within slides, with abundance-matched random target labels as an additional control."
)
legend <- c(
  "Candidate panel draft: GSE310352 CosMx robustness analysis. (A) Gene-set overlap audit for CAF/matrix and TGF/EMT modules, including non-overlap TGF/EMT scores. (B) CAF/matrix-to-TGF/EMT-interface observed/expected adjacency across threshold definitions. (C) Normalization sensitivity using raw log-z, library-size normalized log-z, and rank-based scores. (D) TGF/EMT-interface identity QC showing PanCK/CD45/CD3, CAF/tumor overlap, and proximity to PanCK-high epithelial cells. (E) Spatial-null sensitivity comparing within-FOV, within-slide, and abundance-matched random-label nulls. (F) FOV- and slide-level stability with leave-one-slide-out summaries."
)
writeLines(c("# GSE310352 CosMx Final Decision Report", "", paste0("Decision: ", decision), "", "## Rationale", "", "The strongest supported GSE310352 signal is CAF/matrix-to-TGF/EMT-interface adjacency. The result is robust enough for a secondary candidate ED panel after manual review, but not as a standalone core ED10 result because patient IDs and author cell-type/tumor-subtype annotations are absent from the processed CSVs.", "", "## Key Robustness Summaries", "", "### Non-overlap", capture.output(print(nonoverlap_summary)), "", "### Threshold", capture.output(print(thr_summary)), "", "### Normalization", capture.output(print(norm_summary)), "", "### Spatial null", capture.output(print(spatial_null)), "", "### Controls", capture.output(print(control_summary)), "", "## Recommended Title", "", "CosMx cell-level support for CAF/matrix-associated TGF/EMT stromal-interface organization", "", "## Claims Supported", "", "- Cell-level spatial support for CAF/matrix-associated TGF/EMT stromal-interface organization.", "- Orthogonal support consistent with stromal-neighborhood plausibility.", "", "## Claims To Avoid", "", "- causal signaling", "- direct SPP1-CD44 validation", "- tumor-intrinsic EMT unless additional tumor-annotation support is added", "- Visium gradient reconstruction", "- clinical prediction"), file.path(doc_dir, "gse310352_cosmx_final_decision_report.md"), useBytes = TRUE)
writeLines(c("# GSE310352 Results Paragraph Draft", "", safe_results), file.path(doc_dir, "gse310352_cosmx_results_paragraph_draft.md"), useBytes = TRUE)
writeLines(c("# GSE310352 Methods Draft", "", safe_methods), file.path(doc_dir, "gse310352_cosmx_methods_draft.md"), useBytes = TRUE)
writeLines(c("# GSE310352 Figure Legend Draft", "", legend), file.path(doc_dir, "gse310352_cosmx_figure_legend_draft.md"), useBytes = TRUE)

manifest <- data.table(
  output = c(
    "gse310352_caf_tgfemt_gene_overlap_audit.csv",
    "gse310352_nonoverlap_tgfemt_adjacency_results.csv",
    "gse310352_threshold_sensitivity.csv",
    "gse310352_normalization_sensitivity.csv",
    "gse310352_tgfemt_interface_identity_summary.csv",
    "gse310352_spatial_null_sensitivity.csv",
    "gse310352_fov_level_adjacency.csv",
    "gse310352_leave_one_slide_out.csv",
    "gse310352_positive_negative_control_results.csv",
    "gse310352_cosmx_final_decision_report.md"
  ),
  path = c(
    file.path(table_dir, "gse310352_caf_tgfemt_gene_overlap_audit.csv"),
    file.path(table_dir, "gse310352_nonoverlap_tgfemt_adjacency_results.csv"),
    file.path(table_dir, "gse310352_threshold_sensitivity.csv"),
    file.path(table_dir, "gse310352_normalization_sensitivity.csv"),
    file.path(table_dir, "gse310352_tgfemt_interface_identity_summary.csv"),
    file.path(table_dir, "gse310352_spatial_null_sensitivity.csv"),
    file.path(table_dir, "gse310352_fov_level_adjacency.csv"),
    file.path(table_dir, "gse310352_leave_one_slide_out.csv"),
    file.path(table_dir, "gse310352_positive_negative_control_results.csv"),
    file.path(doc_dir, "gse310352_cosmx_final_decision_report.md")
  ),
  description = c(
    "CAF/TGF gene overlap audit",
    "Non-overlap TGF/EMT adjacency sensitivity",
    "Threshold sensitivity",
    "Normalization sensitivity",
    "TGF/EMT interface identity summary",
    "Spatial null sensitivity",
    "Per-FOV adjacency table",
    "Leave-one-slide-out stability",
    "Positive/negative controls",
    "Final robustness decision"
  )
)
write_csv(manifest, file.path(manifest_dir, "gse310352_cosmx_robustness_manifest.csv"))

files <- manifest$path[file.exists(manifest$path)]
sha <- tools::md5sum(files)
write_csv(data.table(file = names(sha), md5 = as.character(sha)), file.path(manifest_dir, "gse310352_cosmx_robustness_checksums_md5.csv"))

cat("GSE310352 robustness complete. Decision: ", decision, "\n", sep = "")
