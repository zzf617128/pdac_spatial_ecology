#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(FNN)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"

target_scores <- c(
  immune_core = "score_immune_hub_core",
  ifn_mhc = "z_ifn_antigen_presentation",
  tumor_aggressive = "score_tumor_aggressive",
  immune_maturity = "score_immune_hub_maturity"
)
thresholds <- data.table(core_label = c("top_15", "top_10", "top_5"), quantile = c(0.85, 0.90, 0.95))

now_iso <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}

write_status <- function(status, payload) {
  path <- file.path(project_root, "results", "logs", "stage_18_caf_core_threshold_sensitivity_status.json")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  payload <- c(list(stage = "18_caf_core_threshold_sensitivity", status = status, timestamp_utc = now_iso()), payload)
  to_json <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) return("null")
    if (length(x) > 1) return(paste0("[", paste(vapply(x, to_json, character(1)), collapse = ", "), "]"))
    if (is.numeric(x) || is.integer(x)) return(as.character(x))
    if (is.logical(x)) return(tolower(as.character(x)))
    paste0('"', gsub('"', '\\"', as.character(x), fixed = TRUE), '"')
  }
  lines <- c("{", paste(sprintf('  "%s": %s', names(payload), vapply(payload, to_json, character(1))), collapse = ",\n"), "}")
  writeLines(lines, path, useBytes = TRUE)
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

normalize_specimen_type <- function(x) {
  y <- tolower(as.character(x))
  out <- rep("metadata_required", length(y))
  out[grepl("primary|tumou?r|pdac|pancreas$", y)] <- "primary_tumor"
  out[grepl("liver|hepatic", y)] <- "liver_metastasis"
  out[grepl("lymph|node|ln", y)] <- "lymph_node_metastasis"
  out[grepl("normal", y)] <- "normal_pancreas"
  out
}

analyze_sample <- function(dt, cohort_label, thresholds) {
  if ("edge_or_background_risk" %in% names(dt)) {
    dt <- dt[!(as.logical(edge_or_background_risk) %in% TRUE)]
  }
  dt <- dt[is.finite(x_pixel) & is.finite(y_pixel)]
  if (nrow(dt) < 100) return(NULL)
  coords <- as.matrix(dt[, .(x_pixel, y_pixel)])
  nn <- median_neighbor_distance(coords)
  if (!is.finite(nn) || nn <= 0) return(NULL)
  caf <- dt$score_caf_myeloid_barrier

  rows <- list()
  for (i in seq_len(nrow(thresholds))) {
    threshold <- quantile(caf, thresholds$quantile[[i]], na.rm = TRUE, names = FALSE)
    core_idx <- which(caf >= threshold)
    if (length(core_idx) < 10) next
    dist <- nearest_distance(coords, coords[core_idx, , drop = FALSE]) / nn
    for (target in names(target_scores)) {
      rows[[length(rows) + 1]] <- data.table(
        cohort_label = cohort_label,
        dataset_id = dt$dataset_id[[1]],
        sample_id = dt$sample_id[[1]],
        specimen_type = if ("specimen_type" %in% names(dt)) dt$specimen_type[[1]] else "metadata_required",
        core_label = thresholds$core_label[[i]],
        core_quantile = thresholds$quantile[[i]],
        target = target,
        n_spots = nrow(dt),
        n_core_spots = length(core_idx),
        rho_distance_to_core = safe_spearman(dist, dt[[target_scores[[target]]]])
      )
    }
  }
  rbindlist(rows)
}

summarize_sensitivity <- function(sample_stats) {
  sample_stats[
    ,
    .(
      n_samples = sum(is.finite(rho_distance_to_core)),
      median_rho = median(rho_distance_to_core, na.rm = TRUE),
      iqr_rho = IQR(rho_distance_to_core, na.rm = TRUE),
      n_negative = sum(rho_distance_to_core < 0, na.rm = TRUE),
      n_positive = sum(rho_distance_to_core > 0, na.rm = TRUE)
    ),
    by = .(cohort_label, dataset_id, specimen_type, core_label, core_quantile, target)
  ][order(cohort_label, dataset_id, specimen_type, target, core_quantile)]
}

plot_sensitivity <- function(summary_stats) {
  out_dir <- file.path(project_root, "results", "figures", "mvp", "sensitivity")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dt <- copy(summary_stats)
  plot_dt[, support_fraction := n_negative / n_samples]
  plot_dt[, group := fifelse(cohort_label == "gse272362", specimen_type, dataset_id)]
  plot_dt <- plot_dt[group %in% c("GSE282302", "GSE274103", "primary_tumor", "liver_metastasis", "lymph_node_metastasis")]
  plot_dt[, core_label := factor(core_label, levels = c("top_15", "top_10", "top_5"))]
  p <- ggplot(plot_dt, aes(x = core_label, y = median_rho, group = group, color = group)) +
    geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2) +
    facet_wrap(~target, scales = "free_y") +
    labs(
      title = "CAF-core threshold sensitivity",
      subtitle = "Distance-to-core gradients remain directionally stable across top 15%, 10%, and 5% CAF-myeloid cores",
      x = "CAF-myeloid core definition",
      y = "Median Spearman rho: distance to CAF core vs target score",
      color = "cohort/site"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(out_dir, "caf_core_threshold_sensitivity.png"), p, width = 11, height = 6, dpi = 180)
  ggsave(file.path(out_dir, "caf_core_threshold_sensitivity.pdf"), p, width = 11, height = 6)
}

main <- function() {
  set.seed(20260624)
  mvp_path <- file.path(project_root, "results", "tables", "mvp_spot_level_scores_with_edge_qc.csv")
  gse_path <- file.path(project_root, "results", "tables", "gse272362_rds_spot_level_scores.csv")
  if (!file.exists(mvp_path)) stop("Missing ", mvp_path)
  if (!file.exists(gse_path)) stop("Missing ", gse_path)

  mvp <- fread(mvp_path)
  gse <- fread(gse_path)
  gse[, specimen_type := normalize_specimen_type(specimen_type)]

  message("Analyzing MVP cohorts")
  mvp_results <- rbindlist(lapply(split(mvp, by = c("dataset_id", "sample_id"), keep.by = TRUE), analyze_sample, cohort_label = "mvp", thresholds = thresholds), fill = TRUE)
  message("Analyzing GSE272362")
  gse_results <- rbindlist(lapply(split(gse, by = c("dataset_id", "sample_id"), keep.by = TRUE), analyze_sample, cohort_label = "gse272362", thresholds = thresholds), fill = TRUE)
  sample_stats <- rbindlist(list(mvp_results, gse_results), fill = TRUE)
  summary_stats <- summarize_sensitivity(sample_stats)

  out_dir <- file.path(project_root, "results", "tables")
  fwrite(sample_stats, file.path(out_dir, "caf_core_threshold_sensitivity_sample_stats.csv"))
  fwrite(summary_stats, file.path(out_dir, "caf_core_threshold_sensitivity_summary.csv"))
  plot_sensitivity(summary_stats)

  write_status("success", list(
    n_sample_rows = nrow(sample_stats),
    n_summary_rows = nrow(summary_stats),
    output_summary = file.path(out_dir, "caf_core_threshold_sensitivity_summary.csv")
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
