#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(FNN)
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"

program_cols <- c(
  mycaf = "z_mycaf",
  icaf = "z_icaf",
  apcaf = "z_apcaf",
  pan_caf = "z_pan_caf",
  myeloid = "z_myeloid",
  spp1_tam = "z_spp1_tam",
  tgfb = "z_tgfb_pathway",
  emt_invasion = "z_emt_invasion",
  hypoxia = "z_hypoxia",
  basal_like = "z_pdac_basal_like",
  classical_like = "z_pdac_classical_like",
  ifn_mhc = "z_ifn_antigen_presentation",
  immune_core = "score_immune_hub_core",
  tumor_aggressive = "score_tumor_aggressive",
  t_cell = "z_t_cell",
  b_cell = "z_b_cell",
  dc_apc = "z_dc_apc",
  plasma_cell = "z_plasma_cell",
  neural_schwann = "z_neural_schwann"
)

program_labels <- c(
  mycaf = "myCAF",
  icaf = "iCAF",
  apcaf = "apCAF",
  pan_caf = "pan-CAF",
  myeloid = "Myeloid",
  spp1_tam = "SPP1/TREM2 TAM",
  tgfb = "TGF-beta",
  emt_invasion = "EMT/invasion",
  hypoxia = "Hypoxia",
  basal_like = "Basal-like tumor",
  classical_like = "Classical-like tumor",
  ifn_mhc = "IFN/MHC",
  immune_core = "Immune core",
  tumor_aggressive = "Tumor-aggressive",
  t_cell = "T cell",
  b_cell = "B cell",
  dc_apc = "DC/APC",
  plasma_cell = "Plasma cell",
  neural_schwann = "Neural/Schwann"
)

site_labels <- c(
  primary_tumor = "Primary tumor",
  liver_metastasis = "Liver metastasis",
  lymph_node_metastasis = "Lymph node metastasis",
  normal_pancreas = "Normal pancreas"
)

site_colors <- c(
  "Primary tumor" = "#1B9E77",
  "Liver metastasis" = "#D95F02",
  "Lymph node metastasis" = "#7570B3",
  "Normal pancreas" = "#666666"
)

normalize_specimen_type <- function(x) {
  y <- tolower(as.character(x))
  out <- rep("metadata_required", length(y))
  out[grepl("primary|tumou?r|pdac|pancreas$", y)] <- "primary_tumor"
  out[grepl("liver|hepatic", y)] <- "liver_metastasis"
  out[grepl("lymph|node|ln", y)] <- "lymph_node_metastasis"
  out[grepl("normal", y)] <- "normal_pancreas"
  out
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 30) return(NA_real_)
  if (length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

median_neighbor_distance <- function(coords) {
  if (nrow(coords) < 3) return(NA_real_)
  sample_idx <- seq_len(nrow(coords))
  if (nrow(coords) > 1500) sample_idx <- sample(sample_idx, 1500)
  sampled <- coords[sample_idx, , drop = FALSE]
  median(FNN::get.knnx(data = sampled, query = sampled, k = 2)$nn.dist[, 2], na.rm = TRUE)
}

nearest_distance <- function(coords, core_coords) {
  FNN::get.knnx(data = core_coords, query = coords, k = 1)$nn.dist[, 1]
}

compute_sample_summary <- function(dt) {
  rows <- list()
  for (program in names(program_cols)) {
    values <- dt[[program_cols[[program]]]]
    rows[[length(rows) + 1]] <- data.table(
      dataset_id = dt$dataset_id[[1]],
      patient_id = dt$patient_id[[1]],
      sample_id = dt$sample_id[[1]],
      specimen_type = dt$specimen_type[[1]],
      program = program,
      program_label = program_labels[[program]],
      n_spots = nrow(dt),
      mean_score = mean(values, na.rm = TRUE),
      median_score = median(values, na.rm = TRUE),
      fraction_score_gt0 = mean(values > 0, na.rm = TRUE),
      fraction_score_gt1 = mean(values > 1, na.rm = TRUE)
    )
  }
  rbindlist(rows)
}

compute_core_gradients <- function(dt) {
  dt <- dt[is.finite(x_pixel) & is.finite(y_pixel)]
  if (nrow(dt) < 100) return(NULL)
  coords <- as.matrix(dt[, .(x_pixel, y_pixel)])
  nn <- median_neighbor_distance(coords)
  if (!is.finite(nn) || nn <= 0) return(NULL)
  caf <- dt$score_caf_myeloid_barrier
  threshold <- quantile(caf, 0.9, na.rm = TRUE, names = FALSE)
  core_idx <- which(caf >= threshold)
  if (length(core_idx) < 10) return(NULL)
  dist <- nearest_distance(coords, coords[core_idx, , drop = FALSE]) / nn
  rbindlist(lapply(names(program_cols), function(program) {
    data.table(
      dataset_id = dt$dataset_id[[1]],
      patient_id = dt$patient_id[[1]],
      sample_id = dt$sample_id[[1]],
      specimen_type = dt$specimen_type[[1]],
      program = program,
      program_label = program_labels[[program]],
      n_spots = nrow(dt),
      n_core_spots = length(core_idx),
      rho_distance_to_caf_core = safe_spearman(dist, dt[[program_cols[[program]]]])
    )
  }))
}

make_patient_deltas <- function(sample_summary, gradients) {
  gradients_wide <- gradients[
    ,
    .(patient_id, sample_id, specimen_type, program, program_label, rho_distance_to_caf_core)
  ]
  metric_dt <- rbindlist(list(
    sample_summary[, .(patient_id, sample_id, specimen_type, program, program_label, metric = "mean_score", value = mean_score)],
    gradients_wide[, .(patient_id, sample_id, specimen_type, program, program_label, metric = "rho_distance_to_caf_core", value = rho_distance_to_caf_core)]
  ), fill = TRUE)

  primary <- metric_dt[specimen_type == "primary_tumor", .(
    patient_id,
    program,
    program_label,
    metric,
    primary_sample_id = sample_id,
    primary_value = value
  )]
  comparison <- metric_dt[specimen_type != "primary_tumor", .(
    patient_id,
    program,
    program_label,
    metric,
    comparison_specimen_type = specimen_type,
    comparison_sample_id = sample_id,
    comparison_value = value
  )]
  deltas <- merge(comparison, primary, by = c("patient_id", "program", "program_label", "metric"), allow.cartesian = TRUE)
  deltas[, delta_vs_primary := comparison_value - primary_value]
  deltas[]
}

summarize_deltas <- function(deltas) {
  deltas[
    ,
    .(
      n_pairs = sum(is.finite(delta_vs_primary)),
      median_delta_vs_primary = median(delta_vs_primary, na.rm = TRUE),
      iqr_delta_vs_primary = IQR(delta_vs_primary, na.rm = TRUE),
      n_increased_vs_primary = sum(delta_vs_primary > 0, na.rm = TRUE),
      n_decreased_vs_primary = sum(delta_vs_primary < 0, na.rm = TRUE)
    ),
    by = .(comparison_specimen_type, program, program_label, metric)
  ][order(metric, comparison_specimen_type, program)]
}

summarize_gradients <- function(gradients) {
  gradients[
    ,
    .(
      n_samples = sum(is.finite(rho_distance_to_caf_core)),
      median_rho = median(rho_distance_to_caf_core, na.rm = TRUE),
      iqr_rho = IQR(rho_distance_to_caf_core, na.rm = TRUE),
      n_negative = sum(rho_distance_to_caf_core < 0, na.rm = TRUE),
      n_positive = sum(rho_distance_to_caf_core > 0, na.rm = TRUE)
    ),
    by = .(specimen_type, program, program_label)
  ][order(specimen_type, median_rho)]
}

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

plot_figure3 <- function(delta_summary, gradient_summary, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  selected_programs <- c(
    "mycaf", "icaf", "apcaf", "myeloid", "spp1_tam", "tgfb",
    "emt_invasion", "basal_like", "ifn_mhc", "immune_core", "tumor_aggressive",
    "t_cell", "b_cell", "dc_apc"
  )
  heat_dt <- gradient_summary[program %in% selected_programs & specimen_type != "normal_pancreas"]
  heat_dt[, specimen_label := factor(site_labels[specimen_type], levels = site_labels[c("primary_tumor", "liver_metastasis", "lymph_node_metastasis")])]
  heat_dt[, program_label := factor(program_label, levels = rev(program_labels[selected_programs]))]
  p_heat <- ggplot(heat_dt, aes(x = specimen_label, y = program_label, fill = median_rho)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = paste0(n_negative, "/", n_samples)), size = 2.3, color = "black") +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-0.55, 0.55), oob = scales::squish) +
    labs(
      title = "CAF-core subprogram decomposition",
      subtitle = "Tile color: median rho between distance to CAF core and program score; labels: negative samples/total",
      x = NULL,
      y = NULL,
      fill = "Median rho"
    ) +
    theme_paper(8.5) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))

  delta_programs <- c("ifn_mhc", "immune_core", "tumor_aggressive", "spp1_tam", "tgfb", "emt_invasion", "basal_like")
  delta_dt <- delta_summary[
    metric == "rho_distance_to_caf_core" &
      comparison_specimen_type %in% c("liver_metastasis", "lymph_node_metastasis") &
      program %in% delta_programs
  ]
  delta_dt[, specimen_label := factor(site_labels[comparison_specimen_type], levels = site_labels[c("liver_metastasis", "lymph_node_metastasis")])]
  delta_dt[, program_label := factor(program_label, levels = rev(program_labels[delta_programs]))]
  p_delta <- ggplot(delta_dt, aes(x = median_delta_vs_primary, y = program_label, fill = specimen_label)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.35) +
    geom_col(position = position_dodge2(width = 0.75), width = 0.65, color = "white", linewidth = 0.15) +
    geom_text(aes(label = paste0(n_decreased_vs_primary, "/", n_pairs, " more core-like")),
              position = position_dodge2(width = 0.75), hjust = ifelse(delta_dt$median_delta_vs_primary < 0, 1.05, -0.05),
              size = 2.1, color = "grey20") +
    scale_fill_manual(values = site_colors, drop = FALSE) +
    labs(
      title = "Patient-matched metastatic remodeling",
      subtitle = "Delta is metastasis minus matched primary for distance-to-CAF-core rho",
      x = "Median delta vs matched primary",
      y = NULL,
      fill = NULL
    ) +
    theme_paper(8.5)

  paired_dt <- fread(file.path(project_root, "results", "tables", "gse272362_patient_matched_site_deltas.csv"))
  paired_dt <- paired_dt[
    metric == "rho_distance_to_caf_core" &
      comparison_specimen_type %in% c("liver_metastasis", "lymph_node_metastasis") &
      program %in% c("ifn_mhc", "immune_core", "tumor_aggressive")
  ]
  paired_dt[, specimen_label := factor(site_labels[comparison_specimen_type], levels = site_labels[c("liver_metastasis", "lymph_node_metastasis")])]
  paired_dt[, program_label := factor(program_label, levels = program_labels[c("ifn_mhc", "immune_core", "tumor_aggressive")])]
  p_pairs <- ggplot(paired_dt, aes(x = primary_value, y = comparison_value, color = specimen_label)) +
    geom_abline(slope = 1, intercept = 0, color = "grey55", linewidth = 0.35) +
    geom_point(size = 2.0, alpha = 0.9) +
    facet_wrap(~program_label, scales = "free") +
    scale_color_manual(values = site_colors, drop = FALSE) +
    labs(
      title = "Matched-patient core-gradient shifts",
      subtitle = "Each point compares one metastasis with the matched primary from the same patient",
      x = "Primary rho",
      y = "Metastasis rho",
      color = NULL
    ) +
    theme_paper(8.5)

  figure <- arrangeGrob(
    p_heat,
    arrangeGrob(p_delta, p_pairs, ncol = 2, widths = c(1.05, 1)),
    nrow = 2,
    heights = c(1.2, 1)
  )

  ggsave(file.path(out_dir, "figure3_gse272362_matched_decomposition.pdf"), figure, width = 12, height = 10, bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, "figure3_gse272362_matched_decomposition.png"), figure, width = 12, height = 10, dpi = 260, bg = "white", limitsize = FALSE)
}

write_status <- function(status, payload) {
  path <- file.path(project_root, "results", "logs", "stage_20_gse272362_patient_matched_decomposition_status.json")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  payload <- c(list(stage = "20_gse272362_patient_matched_decomposition", status = status), payload)
  json <- paste0(
    "{\n",
    paste(sprintf('  "%s": "%s"', names(payload), gsub('"', '\\"', as.character(payload), fixed = TRUE)), collapse = ",\n"),
    "\n}\n"
  )
  writeLines(json, path, useBytes = TRUE)
}

main <- function() {
  set.seed(20260624)
  spot_path <- file.path(project_root, "results", "tables", "gse272362_rds_spot_level_scores.csv")
  if (!file.exists(spot_path)) stop("Missing ", spot_path)
  spot <- fread(spot_path)
  spot[, specimen_type := normalize_specimen_type(specimen_type)]
  missing <- setdiff(c("patient_id", "sample_id", "specimen_type", "x_pixel", "y_pixel", "score_caf_myeloid_barrier", program_cols), names(spot))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

  by_sample <- split(spot, by = c("patient_id", "sample_id"), keep.by = TRUE)
  message("Computing sample-level program summaries")
  sample_summary <- rbindlist(lapply(by_sample, compute_sample_summary), fill = TRUE)
  message("Computing CAF-core subprogram gradients")
  gradients <- rbindlist(lapply(by_sample, compute_core_gradients), fill = TRUE)
  gradient_summary <- summarize_gradients(gradients)
  deltas <- make_patient_deltas(sample_summary, gradients)
  delta_summary <- summarize_deltas(deltas)

  table_dir <- file.path(project_root, "results", "tables")
  fwrite(sample_summary, file.path(table_dir, "gse272362_sample_site_program_summary.csv"))
  fwrite(gradients, file.path(table_dir, "gse272362_caf_core_subprogram_gradients.csv"))
  fwrite(gradient_summary, file.path(table_dir, "gse272362_caf_core_subprogram_gradient_summary.csv"))
  fwrite(deltas, file.path(table_dir, "gse272362_patient_matched_site_deltas.csv"))
  fwrite(delta_summary, file.path(table_dir, "gse272362_patient_matched_site_delta_summary.csv"))

  plot_figure3(delta_summary, gradient_summary, file.path(project_root, "results", "figures", "main"))

  write_status("success", list(
    n_samples = length(by_sample),
    n_gradient_rows = nrow(gradients),
    n_delta_rows = nrow(deltas)
  ))
  message("Done.")
}

tryCatch(
  main(),
  error = function(e) {
    write_status("failed", list(error = conditionMessage(e)))
    stop(e)
  }
)
