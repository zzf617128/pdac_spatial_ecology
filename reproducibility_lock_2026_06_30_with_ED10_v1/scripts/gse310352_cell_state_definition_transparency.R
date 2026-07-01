suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(FNN)
  library(gridExtra)
  library(patchwork)
})

set.seed(20260630)

root <- normalizePath(
  file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search"),
  winslash = "/",
  mustWork = TRUE
)
data_dir <- file.path(root, "datasets/gse310352_processed")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
doc_dir <- file.path(root, "docs")
script_dir <- file.path(root, "scripts")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) fwrite(as.data.frame(x), path)

state_file <- file.path(table_dir, "gse310352_cosmx_cell_state_scores.csv")
if (!file.exists(state_file)) stop("Missing cell-state score table: ", state_file)

marker_sets <- list(
  caf_matrix = c("ACTA2", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "POSTN", "MMP2", "PDPN"),
  tgfb_emt = c("TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "VIM", "ZEB1", "SNAI1", "SNAI2", "MMP2", "MMP9", "ITGA5"),
  tumor_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MSLN", "KRT17", "KRT5", "KRT6A", "S100A2", "LAMC2"),
  immune_myeloid = c("PTPRC", "LYZ", "CD68", "CD163", "CD14", "C1QA", "C1QB", "APOE", "HLA-DRA", "CD74"),
  t_cell = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "TRAC")
)

module_features <- data.table(
  marker_family = c("caf_matrix", "tgfb_emt", "tumor_epithelial", "tumor_epithelial", "immune_myeloid", "immune_myeloid", "t_cell"),
  marker = c("score_caf_mycaf_matrix", "score_tgfb_emt", "score_tumor_epithelial", "score_tumor_aggressive", "score_myeloid_macrophage_tam", "score_ifn_apc", "score_t_cell"),
  display_marker = c("CAF/matrix module", "TGF/EMT module", "tumor epithelial module", "tumor-aggressive module", "myeloid/macrophage module", "IFN/APC module", "T cell module")
)

all_requested_genes <- unique(unlist(marker_sets, use.names = FALSE))

message("Loading GSE310352 cell-state scores")
state_cols <- c(
  "slide", "gsm", "fov", "cell_ID", "cell_key",
  "CenterX_global_px", "CenterY_global_px",
  "Mean.MembraneStain", "Max.MembraneStain",
  "Mean.PanCK", "Max.PanCK", "Mean.CD45", "Max.CD45", "Mean.CD3", "Max.CD3",
  "Mean.DAPI", "Max.DAPI",
  "score_caf_mycaf_matrix", "score_myeloid_macrophage_tam", "score_tgfb_emt",
  "score_ifn_apc", "score_t_cell", "score_b_plasma", "score_tumor_epithelial", "score_tumor_aggressive",
  "state_tumor_epithelial", "state_t_cell", "state_b_plasma", "state_myeloid_macrophage",
  "state_caf_matrix", "state_tgfb_emt_interface", "state_other"
)
state_dt <- fread(state_file, select = state_cols)

logical_cols <- grep("^state_", names(state_dt), value = TRUE)
for (cc in logical_cols) {
  if (!is.logical(state_dt[[cc]])) state_dt[, (cc) := as.logical(get(cc))]
  state_dt[is.na(get(cc)), (cc) := FALSE]
}
state_dt[, state_immune_like := state_t_cell | state_b_plasma | state_myeloid_macrophage]
state_dt[, state_other_background := !(state_caf_matrix | state_tgfb_emt_interface | state_tumor_epithelial | state_immune_like)]

cell_state_defs <- data.table(
  cell_state = c("CAF/matrix-like", "TGF/EMT stromal-interface", "tumor epithelial-like", "immune-like", "other/background"),
  flag_col = c("state_caf_matrix", "state_tgfb_emt_interface", "state_tumor_epithelial", "state_immune_like", "state_other_background")
)

slide_map <- unique(state_dt[, .(slide, gsm)])
slide_map[, metadata_file := file.path(data_dir, paste0(gsm, "_", slide, "_cell_metadata.csv.gz"))]
slide_map[, counts_file := file.path(data_dir, paste0(gsm, "_", slide, "_cell_by_gene_counts.csv.gz"))]
if (!all(file.exists(slide_map$counts_file))) {
  stop("Missing count files for slides: ", paste(slide_map[!file.exists(counts_file)]$slide, collapse = ", "))
}

first_header <- names(fread(slide_map$counts_file[1], nrows = 0))
available_genes <- intersect(all_requested_genes, first_header)
gene_availability <- rbindlist(lapply(names(marker_sets), function(ms) {
  data.table(
    marker_family = ms,
    marker = marker_sets[[ms]],
    available = marker_sets[[ms]] %in% available_genes
  )
}))

message("Loading selected marker counts")
expr_list <- vector("list", nrow(slide_map))
for (i in seq_len(nrow(slide_map))) {
  sl <- slide_map$slide[i]
  message("  ", sl)
  keep <- c("fov", "cell_ID", available_genes)
  counts <- fread(slide_map$counts_file[i], select = keep)
  counts[, slide := sl]
  expr_list[[i]] <- counts
}
expr_dt <- rbindlist(expr_list, fill = TRUE)
rm(expr_list)
dt <- merge(state_dt, expr_dt, by = c("slide", "fov", "cell_ID"), all.x = TRUE, sort = FALSE)
rm(state_dt, expr_dt)
for (g in available_genes) dt[is.na(get(g)), (g) := 0]

slide_levels <- c(sort(unique(dt$slide)), "all")

summarize_feature <- function(data, flag, values, is_gene = TRUE) {
  inside <- data[[flag]] == TRUE
  outside <- !inside
  n_in <- sum(inside, na.rm = TRUE)
  n_out <- sum(outside, na.rm = TRUE)
  if (n_in == 0 || n_out == 0 || all(is.na(values))) {
    return(list(
      n_cells = n_in, n_other = n_out, mean_value = NA_real_, median_value = NA_real_,
      pct_detected = NA_real_, other_mean_value = NA_real_, delta_mean_value = NA_real_,
      log2_fc_raw_vs_other = NA_real_
    ))
  }
  if (is_gene) {
    raw <- values
    val <- log1p(raw)
    raw_in <- raw[inside]
    raw_out <- raw[outside]
    list(
      n_cells = n_in,
      n_other = n_out,
      mean_value = mean(val[inside], na.rm = TRUE),
      median_value = median(val[inside], na.rm = TRUE),
      pct_detected = mean(raw_in > 0, na.rm = TRUE),
      other_mean_value = mean(val[outside], na.rm = TRUE),
      delta_mean_value = mean(val[inside], na.rm = TRUE) - mean(val[outside], na.rm = TRUE),
      log2_fc_raw_vs_other = log2((mean(raw_in, na.rm = TRUE) + 0.05) / (mean(raw_out, na.rm = TRUE) + 0.05))
    )
  } else {
    list(
      n_cells = n_in,
      n_other = n_out,
      mean_value = mean(values[inside], na.rm = TRUE),
      median_value = median(values[inside], na.rm = TRUE),
      pct_detected = NA_real_,
      other_mean_value = mean(values[outside], na.rm = TRUE),
      delta_mean_value = mean(values[inside], na.rm = TRUE) - mean(values[outside], na.rm = TRUE),
      log2_fc_raw_vs_other = NA_real_
    )
  }
}

message("Summarizing marker/module enrichment")
enrich_rows <- list()
for (sl in slide_levels) {
  dsub <- if (sl == "all") dt else dt[slide == sl]
  summary_level <- if (sl == "all") "overall" else "slide"
  for (state_i in seq_len(nrow(cell_state_defs))) {
    cell_state <- cell_state_defs$cell_state[state_i]
    flag <- cell_state_defs$flag_col[state_i]
    for (ms in names(marker_sets)) {
      for (gene in marker_sets[[ms]]) {
        avail <- gene %in% available_genes
        vals <- if (avail) dsub[[gene]] else rep(NA_real_, nrow(dsub))
        ss <- summarize_feature(dsub, flag, vals, is_gene = TRUE)
        enrich_rows[[length(enrich_rows) + 1]] <- data.table(
          summary_level = summary_level, slide = sl, cell_state = cell_state,
          feature_type = "gene", marker_family = ms, marker = gene, display_marker = gene,
          available = avail,
          n_cells = ss$n_cells, n_other = ss$n_other,
          mean_value = ss$mean_value, median_value = ss$median_value,
          pct_detected = ss$pct_detected, other_mean_value = ss$other_mean_value,
          delta_mean_value = ss$delta_mean_value, log2_fc_raw_vs_other = ss$log2_fc_raw_vs_other
        )
      }
    }
    for (j in seq_len(nrow(module_features))) {
      feat <- module_features$marker[j]
      ss <- summarize_feature(dsub, flag, dsub[[feat]], is_gene = FALSE)
      enrich_rows[[length(enrich_rows) + 1]] <- data.table(
        summary_level = summary_level, slide = sl, cell_state = cell_state,
        feature_type = "module_score", marker_family = module_features$marker_family[j],
        marker = feat, display_marker = module_features$display_marker[j],
        available = TRUE,
        n_cells = ss$n_cells, n_other = ss$n_other,
        mean_value = ss$mean_value, median_value = ss$median_value,
        pct_detected = ss$pct_detected, other_mean_value = ss$other_mean_value,
        delta_mean_value = ss$delta_mean_value, log2_fc_raw_vs_other = ss$log2_fc_raw_vs_other
      )
    }
  }
}
marker_enrichment <- rbindlist(enrich_rows, fill = TRUE)
write_csv(marker_enrichment, file.path(table_dir, "gse310352_cell_state_definition_marker_enrichment.csv"))

heat_dt <- marker_enrichment[
  summary_level == "overall" & available == TRUE &
    (feature_type == "module_score" | marker %in% all_requested_genes)
]
heat_dt[, feature_label := ifelse(feature_type == "module_score", display_marker, marker)]
family_order <- c("caf_matrix", "tgfb_emt", "tumor_epithelial", "immune_myeloid", "t_cell")
heat_dt[, marker_family := factor(marker_family, levels = family_order)]
state_order <- c("CAF/matrix-like", "TGF/EMT stromal-interface", "tumor epithelial-like", "immune-like", "other/background")
heat_dt[, cell_state := factor(cell_state, levels = state_order)]
heat_dt[, feature_label := factor(feature_label, levels = rev(unique(feature_label[order(marker_family, feature_type, feature_label)])))]

p_heat <- ggplot(heat_dt, aes(x = cell_state, y = feature_label, fill = delta_mean_value)) +
  geom_tile(color = "white", linewidth = 0.18) +
  facet_grid(marker_family ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#B5453C", midpoint = 0, na.value = "grey90") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text.y = element_text(angle = 0, size = 7),
    panel.spacing.y = unit(0.08, "lines")
  ) +
  labs(x = NULL, y = NULL, fill = "Mean log1p/score\ndelta vs rest")
ggsave(file.path(figure_dir, "gse310352_cell_state_marker_enrichment_heatmap.pdf"), p_heat, width = 8, height = 11, useDingbats = FALSE)

message("Summarizing IF marker QC")
if_fields <- c("Mean.PanCK", "Mean.CD45", "Mean.CD3", "Mean.DAPI", "Mean.MembraneStain")
if_long_thresholds <- rbindlist(lapply(if_fields, function(fld) {
  dt[, .(q80 = as.numeric(quantile(get(fld), 0.80, na.rm = TRUE))), by = slide][, if_marker := fld]
}))

if_rows <- list()
for (sl in slide_levels) {
  dsub <- if (sl == "all") dt else dt[slide == sl]
  summary_level <- if (sl == "all") "overall" else "slide"
  for (state_i in seq_len(nrow(cell_state_defs))) {
    cell_state <- cell_state_defs$cell_state[state_i]
    flag <- cell_state_defs$flag_col[state_i]
    for (fld in if_fields) {
      vals <- dsub[[fld]]
      if (sl == "all") {
        q_dt <- if_long_thresholds[if_marker == fld, .(slide, q80)]
        local_q <- q_dt$q80[match(dsub$slide, q_dt$slide)]
      } else {
        local_q <- rep(if_long_thresholds[slide == sl & if_marker == fld]$q80, nrow(dsub))
      }
      inside <- dsub[[flag]] == TRUE
      if_rows[[length(if_rows) + 1]] <- data.table(
        summary_level = summary_level, slide = sl, cell_state = cell_state, if_marker = fld,
        n_cells = sum(inside, na.rm = TRUE),
        mean_value = mean(vals[inside], na.rm = TRUE),
        median_value = median(vals[inside], na.rm = TRUE),
        q25 = as.numeric(quantile(vals[inside], 0.25, na.rm = TRUE)),
        q75 = as.numeric(quantile(vals[inside], 0.75, na.rm = TRUE)),
        pct_above_slide_q80 = mean(vals[inside] >= local_q[inside], na.rm = TRUE),
        other_mean_value = mean(vals[!inside], na.rm = TRUE),
        delta_log1p_vs_rest = mean(log1p(vals[inside]), na.rm = TRUE) - mean(log1p(vals[!inside]), na.rm = TRUE)
      )
    }
  }
}
if_qc <- rbindlist(if_rows, fill = TRUE)
write_csv(if_qc, file.path(table_dir, "gse310352_cell_state_if_marker_qc.csv"))

sample_if <- rbindlist(lapply(seq_len(nrow(cell_state_defs)), function(i) {
  flag <- cell_state_defs$flag_col[i]
  rows <- which(dt[[flag]] == TRUE)
  if (length(rows) > 8000) rows <- sample(rows, 8000)
  out <- dt[rows, c("slide", if_fields), with = FALSE]
  out[, cell_state := cell_state_defs$cell_state[i]]
  out
}), fill = TRUE)
if_plot_dt <- melt(sample_if, id.vars = c("slide", "cell_state"), measure.vars = if_fields, variable.name = "if_marker", value.name = "value")
if_plot_dt[, cell_state := factor(cell_state, levels = state_order)]
p_if <- ggplot(if_plot_dt, aes(x = cell_state, y = log1p(value), fill = cell_state)) +
  geom_boxplot(outlier.shape = NA, width = 0.72, linewidth = 0.25) +
  facet_wrap(~ if_marker, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("#5B8C5A", "#B38B00", "#A04545", "#4D70A8", "#BDBDBD")) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "log1p IF intensity")
ggsave(file.path(figure_dir, "gse310352_cell_state_if_marker_qc.pdf"), p_if, width = 10.5, height = 3.4, useDingbats = FALSE)

message("Summarizing cell-state overlaps")
composition_for_tgf <- function(dsub) {
  d <- dsub[state_tgfb_emt_interface == TRUE]
  data.table(
    tgf_composition = c("CAF/matrix-like only", "tumor epithelial-like only", "CAF and tumor overlap", "immune-like", "neither CAF/tumor/immune"),
    n = c(
      sum(d$state_caf_matrix & !d$state_tumor_epithelial, na.rm = TRUE),
      sum(d$state_tumor_epithelial & !d$state_caf_matrix, na.rm = TRUE),
      sum(d$state_caf_matrix & d$state_tumor_epithelial, na.rm = TRUE),
      sum(d$state_immune_like & !d$state_caf_matrix & !d$state_tumor_epithelial, na.rm = TRUE),
      sum(!d$state_caf_matrix & !d$state_tumor_epithelial & !d$state_immune_like, na.rm = TRUE)
    ),
    denominator = nrow(d)
  )[, fraction := n / pmax(denominator, 1)]
}

overlap_rows <- list()
for (sl in slide_levels) {
  dsub <- if (sl == "all") dt else dt[slide == sl]
  summary_level <- if (sl == "all") "overall" else "slide"
  n_total <- nrow(dsub)
  n_caf <- sum(dsub$state_caf_matrix, na.rm = TRUE)
  n_tgf <- sum(dsub$state_tgfb_emt_interface, na.rm = TRUE)
  n_tumor <- sum(dsub$state_tumor_epithelial, na.rm = TRUE)
  n_immune <- sum(dsub$state_immune_like, na.rm = TRUE)
  metrics <- data.table(
    metric = c(
      "CAF and TGF/EMT overlap as fraction of CAF",
      "CAF and TGF/EMT overlap as fraction of TGF/EMT",
      "TGF/EMT and tumor epithelial overlap as fraction of TGF/EMT",
      "TGF/EMT and immune overlap as fraction of TGF/EMT",
      "TGF/EMT cells with CAF/matrix-high state",
      "TGF/EMT cells with tumor epithelial-like state"
    ),
    numerator = c(
      sum(dsub$state_caf_matrix & dsub$state_tgfb_emt_interface, na.rm = TRUE),
      sum(dsub$state_caf_matrix & dsub$state_tgfb_emt_interface, na.rm = TRUE),
      sum(dsub$state_tgfb_emt_interface & dsub$state_tumor_epithelial, na.rm = TRUE),
      sum(dsub$state_tgfb_emt_interface & dsub$state_immune_like, na.rm = TRUE),
      sum(dsub$state_tgfb_emt_interface & dsub$state_caf_matrix, na.rm = TRUE),
      sum(dsub$state_tgfb_emt_interface & dsub$state_tumor_epithelial, na.rm = TRUE)
    ),
    denominator = c(n_caf, n_tgf, n_tgf, n_tgf, n_tgf, n_tgf)
  )
  metrics[, fraction := numerator / pmax(denominator, 1)]
  metrics[, `:=`(summary_level = summary_level, slide = sl, n_total = n_total, n_caf = n_caf, n_tgf = n_tgf, n_tumor = n_tumor, n_immune = n_immune)]
  comp <- composition_for_tgf(dsub)
  comp[, `:=`(
    metric = paste0("TGF/EMT composition: ", tgf_composition),
    numerator = n,
    summary_level = summary_level, slide = sl, n_total = n_total, n_caf = n_caf, n_tgf = n_tgf, n_tumor = n_tumor, n_immune = n_immune
  )]
  overlap_rows[[length(overlap_rows) + 1]] <- rbindlist(list(metrics, comp[, names(metrics), with = FALSE]), fill = TRUE)
}
overlap_audit <- rbindlist(overlap_rows, fill = TRUE)
overlap_audit[, interpretation := fifelse(
  grepl("immune", metric) & fraction < 0.05, "minimal immune overlap",
  fifelse(grepl("tumor epithelial", metric) & fraction > 0.25, "TGF/EMT interface includes epithelial-proximal component; avoid tumor-intrinsic EMT claim",
  fifelse(grepl("CAF", metric) & fraction > 0.25, "TGF/EMT interface includes CAF/matrix-high component", "descriptive overlap metric"))
)]
write_csv(overlap_audit, file.path(table_dir, "gse310352_cell_state_overlap_audit.csv"))

comp_plot <- overlap_audit[grepl("^TGF/EMT composition:", metric)]
comp_plot[, tgf_composition := sub("^TGF/EMT composition: ", "", metric)]
comp_plot <- comp_plot[summary_level == "slide"]
p_overlap <- ggplot(comp_plot, aes(x = slide, y = fraction, fill = tgf_composition)) +
  geom_col(width = 0.76, color = "white", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "CAF/matrix-like only" = "#5B8C5A",
    "tumor epithelial-like only" = "#A04545",
    "CAF and tumor overlap" = "#7B5FA3",
    "immune-like" = "#4D70A8",
    "neither CAF/tumor/immune" = "#BDBDBD"
  )) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right") +
  labs(x = NULL, y = "Fraction of TGF/EMT interface cells", fill = "TGF/EMT composition")
ggsave(file.path(figure_dir, "gse310352_cell_state_overlap_upset_or_barplot.pdf"), p_overlap, width = 7.8, height = 4.2, useDingbats = FALSE)

message("Computing spatial identity context")
dt[, panck_q80_slide := as.numeric(quantile(Mean.PanCK, 0.80, na.rm = TRUE)), by = slide]
dt[, state_panck_high_tumor := state_tumor_epithelial & Mean.PanCK >= panck_q80_slide]

nearest_distance <- function(query_coords, target_coords) {
  if (nrow(query_coords) == 0 || nrow(target_coords) == 0) return(rep(NA_real_, nrow(query_coords)))
  FNN::get.knnx(target_coords, query_coords, k = 1)$nn.dist[, 1]
}

spatial_rows <- list()
target_defs <- data.table(
  target_type = c("CAF/matrix-like", "PanCK-high tumor epithelial-like"),
  target_col = c("state_caf_matrix", "state_panck_high_tumor")
)

for (sl in sort(unique(dt$slide))) {
  for (fv in sort(unique(dt[slide == sl]$fov))) {
    d <- dt[slide == sl & fov == fv & is.finite(CenterX_global_px) & is.finite(CenterY_global_px)]
    if (nrow(d) < 50) next
    query_idx <- which(d$state_tgfb_emt_interface == TRUE)
    if (length(query_idx) < 5) next
    if (length(query_idx) > 1000) query_idx <- sample(query_idx, 1000)
    query_coords <- as.matrix(d[query_idx, .(CenterX_global_px, CenterY_global_px)])
    random_query_pool <- which(d$state_tgfb_emt_interface != TRUE)
    if (length(random_query_pool) >= length(query_idx)) {
      random_query_idx <- sample(random_query_pool, length(query_idx))
    } else {
      random_query_idx <- random_query_pool
    }
    random_query_coords <- as.matrix(d[random_query_idx, .(CenterX_global_px, CenterY_global_px)])
    for (j in seq_len(nrow(target_defs))) {
      target_col <- target_defs$target_col[j]
      target_idx <- which(d[[target_col]] == TRUE)
      if (length(target_idx) < 5) next
      target_coords <- as.matrix(d[target_idx, .(CenterX_global_px, CenterY_global_px)])
      random_target_pool <- which(d[[target_col]] != TRUE)
      if (length(random_target_pool) >= 5) {
        random_target_n <- min(length(target_idx), length(random_target_pool), 3000)
        random_target_idx <- sample(random_target_pool, random_target_n)
        random_target_coords <- as.matrix(d[random_target_idx, .(CenterX_global_px, CenterY_global_px)])
      } else {
        random_target_coords <- matrix(numeric(), ncol = 2)
      }
      actual_dist <- nearest_distance(query_coords, target_coords)
      random_query_dist <- nearest_distance(random_query_coords, target_coords)
      random_target_dist <- nearest_distance(query_coords, random_target_coords)
      spatial_rows[[length(spatial_rows) + 1]] <- data.table(
        summary_level = "fov", slide = sl, fov = fv,
        target_type = target_defs$target_type[j],
        n_cells = nrow(d), n_tgf_query = length(query_idx), n_target = length(target_idx),
        median_actual_distance = median(actual_dist, na.rm = TRUE),
        median_random_query_distance = median(random_query_dist, na.rm = TRUE),
        median_random_target_distance = median(random_target_dist, na.rm = TRUE),
        ratio_actual_to_random_query = median(actual_dist, na.rm = TRUE) / median(random_query_dist, na.rm = TRUE),
        ratio_actual_to_random_target = median(actual_dist, na.rm = TRUE) / median(random_target_dist, na.rm = TRUE)
      )
    }
  }
}
spatial_fov <- rbindlist(spatial_rows, fill = TRUE)
spatial_slide <- spatial_fov[, .(
  n_fov = .N,
  median_actual_distance = median(median_actual_distance, na.rm = TRUE),
  median_random_query_distance = median(median_random_query_distance, na.rm = TRUE),
  median_random_target_distance = median(median_random_target_distance, na.rm = TRUE),
  ratio_actual_to_random_query = median(ratio_actual_to_random_query, na.rm = TRUE),
  ratio_actual_to_random_target = median(ratio_actual_to_random_target, na.rm = TRUE)
), by = .(slide, target_type)]
spatial_slide[, `:=`(summary_level = "slide", fov = NA_integer_, n_cells = NA_integer_, n_tgf_query = NA_integer_, n_target = NA_integer_)]
spatial_overall <- spatial_slide[, .(
  n_fov = sum(n_fov, na.rm = TRUE),
  median_actual_distance = median(median_actual_distance, na.rm = TRUE),
  median_random_query_distance = median(median_random_query_distance, na.rm = TRUE),
  median_random_target_distance = median(median_random_target_distance, na.rm = TRUE),
  ratio_actual_to_random_query = median(ratio_actual_to_random_query, na.rm = TRUE),
  ratio_actual_to_random_target = median(ratio_actual_to_random_target, na.rm = TRUE)
), by = target_type]
spatial_overall[, `:=`(summary_level = "overall", slide = "all", fov = NA_integer_, n_cells = NA_integer_, n_tgf_query = NA_integer_, n_target = NA_integer_)]
spatial_context <- rbindlist(list(spatial_fov, spatial_slide, spatial_overall), fill = TRUE)
write_csv(spatial_context, file.path(table_dir, "gse310352_tgfemt_spatial_identity_context.csv"))

spatial_plot <- melt(
  spatial_slide,
  id.vars = c("slide", "target_type"),
  measure.vars = c("median_actual_distance", "median_random_query_distance", "median_random_target_distance"),
  variable.name = "comparison",
  value.name = "median_distance"
)
spatial_plot[, comparison := factor(comparison, levels = c("median_actual_distance", "median_random_query_distance", "median_random_target_distance"), labels = c("TGF/EMT cells", "random query cells", "random target cells"))]
p_spatial <- ggplot(spatial_plot, aes(x = comparison, y = median_distance, color = slide, group = slide)) +
  geom_line(alpha = 0.55, linewidth = 0.35) +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_wrap(~ target_type, scales = "free_y") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right") +
  labs(x = NULL, y = "Median nearest distance (global px)", color = "Slide")
ggsave(file.path(figure_dir, "gse310352_tgfemt_spatial_identity_context.pdf"), p_spatial, width = 8.2, height = 4.2, useDingbats = FALSE)

message("Building compact QC panel")
module_heat <- marker_enrichment[
  summary_level == "overall" & feature_type == "module_score",
  .(cell_state, display_marker, delta_mean_value)
]
module_heat[, cell_state := factor(cell_state, levels = state_order)]
module_heat[, display_marker := factor(display_marker, levels = rev(module_features$display_marker))]
p_module <- ggplot(module_heat, aes(x = cell_state, y = display_marker, fill = delta_mean_value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#B5453C", midpoint = 0) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right") +
  labs(x = NULL, y = NULL, fill = "Delta")

if_comp <- if_qc[summary_level == "overall" & if_marker %in% c("Mean.PanCK", "Mean.CD45", "Mean.CD3")]
if_comp[, cell_state := factor(cell_state, levels = state_order)]
p_if_comp <- ggplot(if_comp, aes(x = cell_state, y = delta_log1p_vs_rest, fill = if_marker)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right") +
  labs(x = NULL, y = "IF log1p delta vs rest", fill = "IF marker")

p_qc <- (p_module | p_if_comp) / (p_overlap + theme(legend.position = "bottom") | p_spatial + theme(legend.position = "bottom")) +
  plot_annotation(title = "GSE310352 CosMx cell-state definition transparency QC")
ggsave(file.path(figure_dir, "gse310352_cell_state_definition_transparency_qc.pdf"), p_qc, width = 12, height = 9, useDingbats = FALSE)
svg(file.path(figure_dir, "gse310352_cell_state_definition_transparency_qc.svg"), width = 12, height = 9)
print(p_qc)
dev.off()

message("Writing transparency report")
overall_marker <- marker_enrichment[summary_level == "overall"]
get_delta <- function(state, target_marker) {
  val <- overall_marker[cell_state == state & display_marker == target_marker & feature_type == "module_score"]$delta_mean_value
  if (length(val) == 0) NA_real_ else val[1]
}
if_overall <- if_qc[summary_level == "overall"]
get_if <- function(state, marker, field = "pct_above_slide_q80") {
  val <- if_overall[cell_state == state & if_marker == marker][[field]]
  if (length(val) == 0) NA_real_ else val[1]
}
over_all <- overlap_audit[summary_level == "overall"]
frac_metric <- function(pattern) {
  val <- over_all[grepl(pattern, metric)]$fraction
  if (length(val) == 0) NA_real_ else val[1]
}
nonoverlap_path <- file.path(table_dir, "gse310352_nonoverlap_tgfemt_adjacency_results.csv")
nonoverlap_note <- "Non-overlap robustness file was not found."
if (file.exists(nonoverlap_path)) {
  non <- fread(nonoverlap_path)
  non_sum <- non[, .(
    n_slides = uniqueN(slide),
    positive_slides = sum(median_log2_oe > 0, na.rm = TRUE),
    median_slide_log2_oe = median(median_log2_oe, na.rm = TRUE)
  ), by = tgfb_module]
  nonoverlap_note <- paste(capture.output(print(non_sum)), collapse = "\n")
}

sp_over <- spatial_context[summary_level == "overall"]
spatial_note <- paste(capture.output(print(sp_over[, .(target_type, median_actual_distance, median_random_query_distance, median_random_target_distance, ratio_actual_to_random_query, ratio_actual_to_random_target)])), collapse = "\n")

report <- c(
  "# GSE310352 Cell-State Definition Transparency Report",
  "",
  "## Scope",
  "This QC analysis evaluates whether the existing rule-based GSE310352 CosMx cell-state labels are marker- and context-consistent. It does not redefine labels, change adjacency results, modify ED10, or make patient-level claims.",
  "",
  "## Data Boundary",
  "- Public metadata did not allow reliable recovery of patient, specimen or tissue-block identifiers for the processed CosMx slides.",
  "- Processed CSV files did not include author cell-type annotations.",
  "- All support remains slide-level/FOV-level and rule-based.",
  "",
  "## CAF/matrix-like Label",
  paste0("- CAF/matrix-like cells show positive overall CAF/matrix module enrichment versus the rest of cells (delta=", round(get_delta("CAF/matrix-like", "CAF/matrix module"), 3), ")."),
  paste0("- The fraction above each slide's 80th percentile is low for PanCK (", round(get_if("CAF/matrix-like", "Mean.PanCK"), 3), "), CD45 (", round(get_if("CAF/matrix-like", "Mean.CD45"), 3), ") and CD3 (", round(get_if("CAF/matrix-like", "Mean.CD3"), 3), "), consistent with a non-epithelial, non-immune stromal-like assignment."),
  "",
  "## TGF/EMT Stromal-Interface Label",
  paste0("- TGF/EMT stromal-interface cells show positive TGF/EMT module enrichment versus the rest of cells (delta=", round(get_delta("TGF/EMT stromal-interface", "TGF/EMT module"), 3), ")."),
  paste0("- Among TGF/EMT interface cells, the fraction overlapping CAF/matrix-like cells is ", round(frac_metric("CAF/matrix-high"), 3), " and the fraction overlapping tumor epithelial-like cells is ", round(frac_metric("tumor epithelial-like state"), 3), "."),
  paste0("- The fraction overlapping immune-like cells is ", round(frac_metric("immune overlap"), 3), ", supporting the interpretation that this label is not primarily an immune state."),
  "- Because the TGF/EMT interface state contains both CAF/matrix-like and epithelial-proximal components, the safest label remains stromal-interface rather than tumor-intrinsic EMT.",
  "",
  "## Tumor And Immune IF Checks",
  paste0("- Tumor epithelial-like cells show high PanCK enrichment: fraction above slide q80=", round(get_if("tumor epithelial-like", "Mean.PanCK"), 3), "."),
  paste0("- Immune-like cells show elevated CD45/CD3 context: CD45 fraction above slide q80=", round(get_if("immune-like", "Mean.CD45"), 3), "; CD3 fraction above slide q80=", round(get_if("immune-like", "Mean.CD3"), 3), "."),
  "",
  "## Spatial Identity Context",
  "Nearest-distance QC was computed within FOVs using sampled TGF/EMT interface query cells where needed for tractability. Overall summary:",
  "```",
  spatial_note,
  "```",
  "These distances are QC context only and are not a new biological claim or a replacement for the original adjacency analysis.",
  "",
  "## CAF/TGF Gene-Overlap Boundary",
  "The direct CAF versus TGF/EMT gene overlap was previously limited to MMP2. Prior robustness tables retained positive CAF/matrix to TGF/EMT adjacency after removing MMP2 and after removing MMP2 plus ITGA5:",
  "```",
  nonoverlap_note,
  "```",
  "This supports the view that the displayed adjacency is not solely a single shared-gene artifact, while still remaining observational and rule-based.",
  "",
  "## Safest Current Claim",
  "GSE310352 CosMx provides slide-level/FOV-level transparency support that rule-based CAF/matrix-like and TGF/EMT stromal-interface labels are marker-consistent enough for use as orthogonal stromal-interface support.",
  "",
  "## Claims Not Supported",
  "- Do not claim patient-level or specimen-level validation.",
  "- Do not claim author-annotated cell types.",
  "- Do not claim tumor-intrinsic EMT.",
  "- Do not claim causal EMT induction or CAF-to-tumor signaling.",
  "- Do not claim direct Visium distance-gradient reconstruction.",
  "",
  "## Outputs",
  "- `tables/gse310352_cell_state_definition_marker_enrichment.csv`",
  "- `tables/gse310352_cell_state_if_marker_qc.csv`",
  "- `tables/gse310352_cell_state_overlap_audit.csv`",
  "- `tables/gse310352_tgfemt_spatial_identity_context.csv`",
  "- `figures/gse310352_cell_state_marker_enrichment_heatmap.pdf`",
  "- `figures/gse310352_cell_state_if_marker_qc.pdf`",
  "- `figures/gse310352_cell_state_overlap_upset_or_barplot.pdf`",
  "- `figures/gse310352_tgfemt_spatial_identity_context.pdf`",
  "- `figures/gse310352_cell_state_definition_transparency_qc.pdf`",
  "- `figures/gse310352_cell_state_definition_transparency_qc.svg`"
)
writeLines(report, file.path(doc_dir, "gse310352_cell_state_definition_transparency_report.md"), useBytes = TRUE)

cat("GSE310352 cell-state definition transparency QC complete.\n")
