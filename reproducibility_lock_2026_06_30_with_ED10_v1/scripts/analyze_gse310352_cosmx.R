suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(FNN)
})

root <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search"), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "datasets/gse310352_processed")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
source_dir <- file.path(root, "source_data")
doc_dir <- file.path(root, "docs")
manifest_dir <- file.path(root, "manifest")
invisible(lapply(c(table_dir, figure_dir, source_dir, doc_dir, manifest_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, path) fwrite(as.data.frame(x), path)

module_defs <- list(
  caf_mycaf_matrix = c("ACTA2", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "POSTN", "MMP2", "FAP", "PDPN"),
  myeloid_macrophage_tam = c("LYZ", "CD68", "CD163", "CD14", "FCGR3A", "C1QA", "C1QB", "C1QC", "APOE", "MRC1", "MARCO", "MSR1"),
  spp1_cd44_integrin = c("SPP1", "CD44", "ITGA5", "ITGAV", "ITGB1", "ITGB5", "FN1", "COL1A1"),
  tgfb_emt = c("TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "VIM", "ZEB1", "SNAI1", "SNAI2", "MMP2", "MMP9", "ITGA5"),
  ifn_apc = c("HLA-DRA", "HLA-DRB1", "CD74", "B2M", "IFIT1", "IFIT3", "CXCL9", "CXCL10", "STAT1", "IRF1"),
  t_cell = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "TRAC"),
  b_plasma = c("MS4A1", "CD79A", "MZB1", "JCHAIN"),
  tumor_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MSLN"),
  tumor_aggressive = c("KRT17", "KRT5", "KRT6A", "S100A2", "LAMC2", "ITGA3")
)

slide_map <- data.table(
  gsm = c("GSM9294399", "GSM9294400", "GSM9294401", "GSM9294402", "GSM9294403", "GSM9294404", "GSM9294405", "GSM9294406"),
  slide = c("slide1", "slide2", "slide3", "slide4", "slide6", "slide5", "slide5b2", "slide5b3")
)
slide_map[, metadata_file := file.path(data_dir, paste0(gsm, "_", slide, "_cell_metadata.csv.gz"))]
slide_map[, counts_file := file.path(data_dir, paste0(gsm, "_", slide, "_cell_by_gene_counts.csv.gz"))]
slide_map[, metadata_exists := file.exists(metadata_file)]
slide_map[, counts_exists := file.exists(counts_file)]

if (!all(slide_map$metadata_exists & slide_map$counts_exists)) {
  missing <- slide_map[!(metadata_exists & counts_exists)]
  stop("Missing GSE310352 processed files: ", paste(missing$slide, collapse = ", "))
}

first_counts <- fread(slide_map$counts_file[1], nrows = 0)
panel_genes <- setdiff(names(first_counts), c("fov", "cell_ID"))
all_requested_genes <- unique(unlist(module_defs, use.names = FALSE))
present_requested <- intersect(all_requested_genes, panel_genes)
selected_cols <- c("fov", "cell_ID", present_requested)

coverage <- rbindlist(lapply(names(module_defs), function(module) {
  genes <- module_defs[[module]]
  present <- intersect(genes, panel_genes)
  data.table(
    dataset = "GSE310352",
    platform = "CosMx",
    module = module,
    requested_genes = paste(genes, collapse = ";"),
    available_genes = paste(present, collapse = ";"),
    n_requested = length(genes),
    n_available = length(present),
    coverage = length(present) / length(genes),
    usable = length(present) >= 2
  )
}))
write_csv(coverage, file.path(table_dir, "gse310352_cosmx_gene_coverage.csv"))

score_module <- function(dt, genes) {
  present <- intersect(genes, names(dt))
  if (length(present) < 2) return(rep(NA_real_, nrow(dt)))
  mat <- as.matrix(dt[, ..present])
  z <- scale(log1p(mat))
  rowMeans(z, na.rm = TRUE)
}

q_by_slide <- function(x, p) {
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

slide_summaries <- list()
cell_scores_list <- list()
controls_list <- list()

for (i in seq_len(nrow(slide_map))) {
  gsm <- slide_map$gsm[i]
  slide <- slide_map$slide[i]
  message("Loading ", slide)
  meta <- fread(slide_map$metadata_file[i])
  meta[, slide := slide]
  meta[, gsm := gsm]
  meta[, cell_key := paste(fov, cell_ID, sep = "_")]
  count_header <- names(fread(slide_map$counts_file[i], nrows = 0))
  gene_cols <- setdiff(count_header, c("fov", "cell_ID"))
  counts_all <- fread(slide_map$counts_file[i])
  counts_all <- counts_all[cell_ID > 0]
  counts_all[, transcript_count_all_panel := rowSums(.SD, na.rm = TRUE), .SDcols = gene_cols]
  keep_cols <- intersect(c("fov", "cell_ID", "transcript_count_all_panel", present_requested), names(counts_all))
  counts <- counts_all[, ..keep_cols]
  rm(counts_all)
  counts[, cell_key := paste(fov, cell_ID, sep = "_")]
  dt <- merge(meta, counts, by = c("cell_key", "fov", "cell_ID"), all.x = TRUE, sort = FALSE)
  dt[is.na(transcript_count_all_panel), transcript_count_all_panel := 0]
  for (g in present_requested) {
    if (!g %in% names(dt)) dt[, (g) := 0]
    dt[is.na(get(g)), (g) := 0]
  }

  for (module in names(module_defs)) {
    dt[, (paste0("score_", module)) := score_module(.SD, module_defs[[module]]), .SDcols = present_requested]
  }

  panck_thr <- q_by_slide(dt$Mean.PanCK, 0.80)
  cd45_thr <- q_by_slide(dt$Mean.CD45, 0.80)
  cd3_thr <- q_by_slide(dt$Mean.CD3, 0.80)
  caf_thr <- q_by_slide(dt$score_caf_mycaf_matrix, 0.80)
  myeloid_thr <- q_by_slide(dt$score_myeloid_macrophage_tam, 0.80)
  tumor_thr <- q_by_slide(dt$score_tumor_epithelial, 0.75)
  tcell_thr <- q_by_slide(dt$score_t_cell, 0.80)
  bcell_thr <- q_by_slide(dt$score_b_plasma, 0.85)
  aggressive_thr <- q_by_slide(dt$score_tumor_aggressive, 0.75)
  tgfb_thr <- q_by_slide(dt$score_tgfb_emt, 0.80)
  cd44_thr <- if ("CD44" %in% names(dt)) q_by_slide(log1p(dt$CD44), 0.80) else Inf
  spp1_thr <- if ("SPP1" %in% names(dt)) q_by_slide(log1p(dt$SPP1), 0.80) else Inf

  dt[, panck_low := Mean.PanCK < panck_thr]
  dt[, cd45_low := Mean.CD45 < cd45_thr]
  dt[, cd45_detectable := Mean.CD45 >= median(Mean.CD45, na.rm = TRUE)]
  dt[, state_tumor_epithelial := (score_tumor_epithelial >= tumor_thr | Mean.PanCK >= panck_thr) & cd45_low]
  dt[, state_t_cell := (score_t_cell >= tcell_thr | Mean.CD3 >= cd3_thr) & cd45_detectable & !state_tumor_epithelial]
  dt[, state_b_plasma := score_b_plasma >= bcell_thr & cd45_detectable & !state_tumor_epithelial & !state_t_cell]
  dt[, state_myeloid_macrophage := score_myeloid_macrophage_tam >= myeloid_thr & panck_low & cd45_detectable & !state_tumor_epithelial & !state_t_cell & !state_b_plasma]
  dt[, state_caf_matrix := score_caf_mycaf_matrix >= caf_thr & panck_low & cd45_low & !state_tumor_epithelial & !state_t_cell & !state_b_plasma]
  dt[, state_spp1_macrophage := state_myeloid_macrophage & "SPP1" %in% names(dt) & log1p(SPP1) >= spp1_thr]
  dt[, state_cd44_tumor_stromal := "CD44" %in% names(dt) & log1p(CD44) >= cd44_thr & (state_tumor_epithelial | state_caf_matrix)]
  dt[, state_tumor_aggressive_epi := state_tumor_epithelial & score_tumor_aggressive >= aggressive_thr]
  dt[, state_tgfb_emt_interface := (state_tumor_epithelial | state_caf_matrix) & score_tgfb_emt >= tgfb_thr]
  dt[, state_other := !(state_tumor_epithelial | state_t_cell | state_b_plasma | state_myeloid_macrophage | state_caf_matrix)]

  state_cols <- grep("^state_", names(dt), value = TRUE)
  score_cols <- grep("^score_", names(dt), value = TRUE)
  if_cols <- grep("^(Mean|Max)\\.", names(dt), value = TRUE)
  cell_scores <- dt[, c("slide", "gsm", "fov", "cell_ID", "cell_key", "Area", "CenterX_global_px", "CenterY_global_px", "transcript_count_all_panel", if_cols, score_cols, state_cols), with = FALSE]
  cell_scores_list[[slide]] <- cell_scores

  slide_summaries[[slide]] <- data.table(
    slide = slide,
    gsm = gsm,
    n_cells_metadata = nrow(meta),
    n_cells_joined = nrow(dt),
    n_fov = uniqueN(dt$fov),
    patient_id_available = FALSE,
    patient_id_field = NA_character_,
    cell_type_annotation_exists = FALSE,
    tumor_subtype_annotation_exists = FALSE,
    coordinate_fields = paste(intersect(c("CenterX_local_px", "CenterY_local_px", "CenterX_global_px", "CenterY_global_px"), names(dt)), collapse = ";"),
    IF_marker_fields = paste(if_cols, collapse = ";"),
    metadata_fields = paste(names(meta), collapse = ";"),
    count_gene_columns = length(gene_cols),
    raw_counts_loaded = TRUE,
    normalized_counts_loaded = FALSE,
    seurat_object_loaded = FALSE,
    median_transcript_count_all_panel = median(dt$transcript_count_all_panel, na.rm = TRUE),
    caf_matrix_fraction = mean(dt$state_caf_matrix, na.rm = TRUE),
    myeloid_macrophage_fraction = mean(dt$state_myeloid_macrophage, na.rm = TRUE),
    spp1_macrophage_fraction = mean(dt$state_spp1_macrophage, na.rm = TRUE),
    tumor_epithelial_fraction = mean(dt$state_tumor_epithelial, na.rm = TRUE),
    tumor_aggressive_epi_fraction = mean(dt$state_tumor_aggressive_epi, na.rm = TRUE),
    t_cell_fraction = mean(dt$state_t_cell, na.rm = TRUE),
    b_plasma_fraction = mean(dt$state_b_plasma, na.rm = TRUE)
  )

  controls_list[[slide]] <- rbindlist(list(
    data.table(slide = slide, control = "tumor_epithelial_vs_PanCK", state = "state_tumor_epithelial", marker_or_score = "Mean.PanCK", state_mean = mean(dt[state_tumor_epithelial == TRUE]$Mean.PanCK, na.rm = TRUE), other_mean = mean(dt[state_tumor_epithelial != TRUE]$Mean.PanCK, na.rm = TRUE)),
    data.table(slide = slide, control = "t_cell_vs_CD3_IF", state = "state_t_cell", marker_or_score = "Mean.CD3", state_mean = mean(dt[state_t_cell == TRUE]$Mean.CD3, na.rm = TRUE), other_mean = mean(dt[state_t_cell != TRUE]$Mean.CD3, na.rm = TRUE)),
    data.table(slide = slide, control = "immune_vs_CD45_IF", state = "state_myeloid_or_t", marker_or_score = "Mean.CD45", state_mean = mean(dt[state_myeloid_macrophage == TRUE | state_t_cell == TRUE]$Mean.CD45, na.rm = TRUE), other_mean = mean(dt[!(state_myeloid_macrophage == TRUE | state_t_cell == TRUE)]$Mean.CD45, na.rm = TRUE)),
    data.table(slide = slide, control = "CAF_state_vs_CAF_score", state = "state_caf_matrix", marker_or_score = "score_caf_mycaf_matrix", state_mean = mean(dt[state_caf_matrix == TRUE]$score_caf_mycaf_matrix, na.rm = TRUE), other_mean = mean(dt[state_caf_matrix != TRUE]$score_caf_mycaf_matrix, na.rm = TRUE)),
    data.table(slide = slide, control = "myeloid_state_vs_myeloid_score", state = "state_myeloid_macrophage", marker_or_score = "score_myeloid_macrophage_tam", state_mean = mean(dt[state_myeloid_macrophage == TRUE]$score_myeloid_macrophage_tam, na.rm = TRUE), other_mean = mean(dt[state_myeloid_macrophage != TRUE]$score_myeloid_macrophage_tam, na.rm = TRUE))
  ))
}

metadata_summary <- rbindlist(slide_summaries, fill = TRUE)
write_csv(metadata_summary, file.path(table_dir, "gse310352_cosmx_metadata_summary.csv"))

cell_scores_all <- rbindlist(cell_scores_list, fill = TRUE)
write_csv(cell_scores_all, file.path(table_dir, "gse310352_cosmx_cell_state_scores.csv"))

controls <- rbindlist(controls_list, fill = TRUE)
controls[, fold_change := state_mean / pmax(other_mean, 1e-9)]
controls[, log2_fold_change := log2(fold_change)]
controls[, positive_control_passed := log2_fold_change > 0.5]
write_csv(controls, file.path(table_dir, "gse310352_cosmx_positive_negative_controls.csv"))

pairs <- data.table(
  pair_id = c("caf_to_myeloid", "caf_to_spp1_macrophage", "spp1_macrophage_to_cd44_tumor_stromal", "caf_to_tumor_aggressive_epi", "caf_to_tgfb_emt_interface", "tumor_epi_to_tumor_epi_positive_control", "tcell_to_tcell_positive_control", "caf_to_b_plasma_negative_control"),
  source_state = c("state_caf_matrix", "state_caf_matrix", "state_spp1_macrophage", "state_caf_matrix", "state_caf_matrix", "state_tumor_epithelial", "state_t_cell", "state_caf_matrix"),
  target_state = c("state_myeloid_macrophage", "state_spp1_macrophage", "state_cd44_tumor_stromal", "state_tumor_aggressive_epi", "state_tgfb_emt_interface", "state_tumor_epithelial", "state_t_cell", "state_b_plasma"),
  expected_direction = c("positive", "positive", "positive", "positive", "positive", "positive_control", "positive_control", "negative_control")
)

adj_rows <- list()
set.seed(20260630)
for (slide in unique(cell_scores_all$slide)) {
  slide_name <- slide
  dt_slide <- cell_scores_all[slide == slide_name]
  for (fv in sort(unique(dt_slide$fov))) {
    dt <- dt_slide[fov == fv & is.finite(CenterX_global_px) & is.finite(CenterY_global_px)]
    if (nrow(dt) < 20) next
    coords <- as.matrix(dt[, .(CenterX_global_px, CenterY_global_px)])
    k <- min(6, nrow(dt) - 1)
    if (k < 1) next
    nn <- FNN::get.knn(coords, k = k)$nn.index
    for (j in seq_len(nrow(pairs))) {
      src <- dt[[pairs$source_state[j]]] == TRUE
      tgt <- dt[[pairs$target_state[j]]] == TRUE
      n_src <- sum(src, na.rm = TRUE)
      n_tgt <- sum(tgt, na.rm = TRUE)
      if (n_src < 5 || n_tgt < 5) {
        adj_rows[[length(adj_rows) + 1]] <- data.table(slide = slide, fov = fv, pair_id = pairs$pair_id[j], source_state = pairs$source_state[j], target_state = pairs$target_state[j], n_cells = nrow(dt), n_source = n_src, n_target = n_tgt, observed_edges = NA_real_, expected_edges = NA_real_, oe = NA_real_, log2_oe = NA_real_, null_mean_oe = NA_real_, null_p_greater = NA_real_, analyzable = FALSE)
        next
      }
      src_idx <- which(src)
      observed <- sum(tgt[as.vector(nn[src_idx, , drop = FALSE])], na.rm = TRUE)
      total_edges <- length(src_idx) * k
      expected <- total_edges * mean(tgt, na.rm = TRUE)
      oe <- (observed + 0.5) / (expected + 0.5)
      null_oe <- rep(NA_real_, 50)
      for (perm in seq_along(null_oe)) {
        tgt_perm <- sample(tgt)
        obs_perm <- sum(tgt_perm[as.vector(nn[src_idx, , drop = FALSE])], na.rm = TRUE)
        exp_perm <- total_edges * mean(tgt_perm, na.rm = TRUE)
        null_oe[perm] <- (obs_perm + 0.5) / (exp_perm + 0.5)
      }
      adj_rows[[length(adj_rows) + 1]] <- data.table(
        slide = slide, fov = fv, pair_id = pairs$pair_id[j],
        source_state = pairs$source_state[j], target_state = pairs$target_state[j],
        n_cells = nrow(dt), n_source = n_src, n_target = n_tgt,
        observed_edges = observed, expected_edges = expected, oe = oe, log2_oe = log2(oe),
        null_mean_oe = mean(null_oe, na.rm = TRUE),
        null_p_greater = (sum(null_oe >= oe, na.rm = TRUE) + 1) / (sum(is.finite(null_oe)) + 1),
        analyzable = TRUE
      )
    }
  }
}

adjacency <- rbindlist(adj_rows, fill = TRUE)
adjacency <- merge(adjacency, pairs[, .(pair_id, expected_direction)], by = "pair_id", all.x = TRUE)
write_csv(adjacency, file.path(table_dir, "gse310352_cosmx_spatial_adjacency_results.csv"))

slide_support <- adjacency[analyzable == TRUE, .(
  n_fov_analyzable = .N,
  median_log2_oe = median(log2_oe, na.rm = TRUE),
  mean_log2_oe = mean(log2_oe, na.rm = TRUE),
  support_fraction_log2oe_gt0 = mean(log2_oe > 0, na.rm = TRUE),
  support_fraction_log2oe_gt025 = mean(log2_oe > 0.25, na.rm = TRUE),
  median_null_p = median(null_p_greater, na.rm = TRUE)
), by = .(slide, pair_id, expected_direction)]

overall_support <- slide_support[, .(
  n_slides = .N,
  median_slide_log2_oe = median(median_log2_oe, na.rm = TRUE),
  mean_slide_log2_oe = mean(median_log2_oe, na.rm = TRUE),
  slide_support_fraction_gt0 = mean(median_log2_oe > 0, na.rm = TRUE),
  slide_support_fraction_gt025 = mean(median_log2_oe > 0.25, na.rm = TRUE),
  median_fov_support_fraction_gt0 = median(support_fraction_log2oe_gt0, na.rm = TRUE),
  median_null_p = median(median_null_p, na.rm = TRUE)
), by = .(pair_id, expected_direction)]
overall_support[, patient_level_available := FALSE]
overall_support[, primary_support_unit := "slide/FOV; patient IDs unavailable in processed CSV"]

write_csv(slide_support, file.path(table_dir, "gse310352_cosmx_slide_level_support.csv"))
write_csv(overall_support, file.path(table_dir, "gse310352_cosmx_patient_level_support.csv"))

source_data <- merge(
  slide_support,
  overall_support[, .(pair_id, overall_median_slide_log2_oe = median_slide_log2_oe, overall_slide_support_fraction_gt0 = slide_support_fraction_gt0, overall_slide_support_fraction_gt025 = slide_support_fraction_gt025)],
  by = "pair_id",
  all.x = TRUE
)
write_csv(source_data, file.path(source_dir, "Source_Data_GSE310352_CosMx.csv"))

candidate_pairs <- overall_support[pair_id %in% c("caf_to_myeloid", "caf_to_spp1_macrophage", "spp1_macrophage_to_cd44_tumor_stromal", "caf_to_tumor_aggressive_epi", "caf_to_tgfb_emt_interface")]
success <- candidate_pairs[
  slide_support_fraction_gt0 >= 0.60 &
    (median_slide_log2_oe > 0.25 | slide_support_fraction_gt025 >= 0.50)
]
recommendation <- if (nrow(success) > 0) "extended_data_candidate_after_review" else "source_only_or_do_not_show"

panel_pdf <- file.path(figure_dir, "gse310352_cosmx_candidate_panels.pdf")
pdf(panel_pdf, width = 8.5, height = 6.5, onefile = TRUE)
print(
  ggplot(coverage, aes(x = reorder(module, coverage), y = coverage, fill = usable)) +
    geom_col(width = 0.72) +
    coord_flip() +
    theme_bw(base_size = 9) +
    scale_fill_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#3182BD")) +
    labs(x = NULL, y = "Gene coverage", fill = "Usable")
)
print(
  ggplot(overall_support, aes(x = reorder(pair_id, median_slide_log2_oe), y = median_slide_log2_oe, fill = expected_direction)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_col(width = 0.72) +
    coord_flip() +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "Median slide log2 observed/expected", fill = "Pair type")
)
print(
  ggplot(slide_support[pair_id %in% candidate_pairs$pair_id], aes(x = pair_id, y = median_log2_oe, color = slide)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_point(size = 2.2, alpha = 0.85, position = position_jitter(width = 0.12, height = 0)) +
    coord_flip() +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "Slide-level median log2 observed/expected", color = "Slide")
)
print(
  ggplot(controls, aes(x = control, y = log2_fold_change, color = slide)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_point(size = 2, alpha = 0.85, position = position_jitter(width = 0.16, height = 0)) +
    coord_flip() +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "State marker/score log2 fold-change", color = "Slide")
)
dev.off()

feasibility_lines <- c(
  "# GSE310352 CosMx Feasibility Report",
  "",
  "## Loaded Data",
  paste0("- Loaded processed metadata/counts for ", nrow(metadata_summary), " slides: ", paste(metadata_summary$slide, collapse = ", "), "."),
  paste0("- Total joined cells: ", format(sum(metadata_summary$n_cells_joined), big.mark = ","), "."),
  paste0("- Total FOVs: ", sum(metadata_summary$n_fov), "."),
  "",
  "## Metadata",
  paste0("- Coordinate fields: ", metadata_summary$coordinate_fields[1]),
  paste0("- IF marker fields include: ", metadata_summary$IF_marker_fields[1]),
  "- No patient ID, author cell-type annotation, or tumor-subtype field was present in the processed CSV metadata inspected here.",
  "- Patient-level support is therefore not available from these processed CSVs; the analysis uses slide and FOV as the primary independent units.",
  "",
  "## Count Data",
  paste0("- Count matrix gene columns: ", metadata_summary$count_gene_columns[1], "."),
  "- Raw counts were loaded from `cell_by_gene_counts.csv.gz`; no normalized counts or Seurat objects were loaded for this deep pass.",
  "",
  "## Coverage",
  paste0("- SPP1 available: ", "SPP1" %in% panel_genes, "; CD44 available: ", "CD44" %in% panel_genes, "."),
  paste0("- See `tables/gse310352_cosmx_gene_coverage.csv`.")
)
writeLines(feasibility_lines, file.path(doc_dir, "gse310352_cosmx_feasibility_report.md"), useBytes = TRUE)

interpretation_lines <- c(
  "# GSE310352 CosMx Interpretation Report",
  "",
  "## Analysis Scope",
  "Rule-based cell states were defined from CosMx RNA module scores and IF PanCK/CD45/CD3 fields. Spatial adjacency was computed within FOVs using k-nearest neighbors and summarized by FOV and slide.",
  "",
  "## Main Signal Summary",
  capture.output(print(overall_support[order(-median_slide_log2_oe)])),
  "",
  "## Candidate ED Decision",
  paste0("- Recommendation: ", recommendation),
  if (nrow(success) > 0) paste0("- Passing pairs: ", paste(success$pair_id, collapse = ", ")) else "- No candidate pair met the pre-defined strong-display criteria in this automated pass.",
  "",
  "## Claim Boundary",
  "- This analysis can support cell-level spatial consistency only if displayed after review.",
  "- It does not prove causal signaling.",
  "- Direct SPP1-CD44 validation should be claimed only for the SPP1 macrophage to CD44 tumor/stromal pair if that pair remains robust after manual review and patient mapping is resolved."
)
writeLines(interpretation_lines, file.path(doc_dir, "gse310352_cosmx_interpretation_report.md"), useBytes = TRUE)

ranking_path <- file.path(table_dir, "signal_strength_ranking.csv")
if (file.exists(ranking_path)) {
  old_rank <- fread(ranking_path, fill = TRUE)
} else {
  old_rank <- data.table()
}
gse_row <- data.table(
  dataset = "GSE310352 in situ multi-modal PDAC CosMx",
  platform = "CosMx",
  validation_target = "cell-level CAF/myeloid/SPP1-CD44/matrix-interface spatial adjacency",
  effect_size = paste0("best candidate median slide log2OE=", round(max(candidate_pairs$median_slide_log2_oe, na.rm = TRUE), 3)),
  support_fraction = paste0("best candidate slide support >0=", round(max(candidate_pairs$slide_support_fraction_gt0, na.rm = TRUE), 3)),
  patient_or_sample_level_support = "slide/FOV-level only; patient IDs unavailable in processed CSV",
  FDR_or_p_value = "random-label null summarized per FOV; see adjacency table",
  robustness_level = ifelse(recommendation == "extended_data_candidate_after_review", "moderate automated", "insufficient for ED display"),
  positive_controls_passed = paste0(sum(controls$positive_control_passed, na.rm = TRUE), "/", nrow(controls), " slide-control rows"),
  negative_controls_passed = "see caf_to_b_plasma_negative_control in adjacency table",
  claim_supported = ifelse(recommendation == "extended_data_candidate_after_review", "potential cell-level spatial support", "source-level exploratory only"),
  claim_boundary = "No causality; no Visium gradient reconstruction; direct SPP1-CD44 only if spatial pair remains robust",
  recommended_display = ifelse(recommendation == "extended_data_candidate_after_review", "extended_data_after_review", "source_only_or_do_not_show"),
  reason = ifelse(recommendation == "extended_data_candidate_after_review", "At least one pre-defined candidate pair met slide-level support criteria.", "Automated slide/FOV support did not meet strong-display criteria or lacks patient mapping.")
)
updated_rank <- if (nrow(old_rank) > 0 && "dataset" %in% names(old_rank)) {
  old_rank <- old_rank[dataset != gse_row$dataset]
  rbindlist(list(gse_row, old_rank), fill = TRUE)
} else {
  gse_row
}
write_csv(updated_rank, ranking_path)
write_csv(updated_rank, file.path(manifest_dir, "signal_strength_ranking.csv"))

panel_manifest <- data.table(
  panel = c("A", "B", "C", "D"),
  title = c("GSE310352 gene/module coverage", "GSE310352 adjacency O/E summary", "GSE310352 slide-level support", "GSE310352 positive controls"),
  status = ifelse(recommendation == "extended_data_candidate_after_review", "candidate_after_review", "source_only"),
  figure_path = panel_pdf,
  source_table = c(
    file.path(table_dir, "gse310352_cosmx_gene_coverage.csv"),
    file.path(table_dir, "gse310352_cosmx_spatial_adjacency_results.csv"),
    file.path(table_dir, "gse310352_cosmx_patient_level_support.csv"),
    file.path(table_dir, "gse310352_cosmx_positive_negative_controls.csv")
  ),
  note = "Automated Phase 3 GSE310352 deep analysis; final ED10 generation deferred."
)
write_csv(panel_manifest, file.path(manifest_dir, "gse310352_cosmx_panel_manifest.csv"))

files_for_checksum <- c(
  file.path(table_dir, "gse310352_cosmx_metadata_summary.csv"),
  file.path(table_dir, "gse310352_cosmx_gene_coverage.csv"),
  file.path(table_dir, "gse310352_cosmx_cell_state_scores.csv"),
  file.path(table_dir, "gse310352_cosmx_spatial_adjacency_results.csv"),
  file.path(table_dir, "gse310352_cosmx_patient_level_support.csv"),
  file.path(table_dir, "gse310352_cosmx_positive_negative_controls.csv"),
  file.path(source_dir, "Source_Data_GSE310352_CosMx.csv"),
  panel_pdf,
  file.path(doc_dir, "gse310352_cosmx_feasibility_report.md"),
  file.path(doc_dir, "gse310352_cosmx_interpretation_report.md")
)
sha <- tools::md5sum(files_for_checksum[file.exists(files_for_checksum)])
write_csv(data.table(file = names(sha), md5 = as.character(sha)), file.path(manifest_dir, "gse310352_cosmx_checksums_md5.csv"))

cat("GSE310352 CosMx analysis complete. Recommendation: ", recommendation, "\n", sep = "")
