#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(FNN)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"
n_perm <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L
if (is.na(n_perm) || n_perm < 100L) stop("n_perm must be at least 100.")

target_scores <- c(
  immune_core = "score_immune_hub_core",
  ifn_mhc = "z_ifn_antigen_presentation",
  tumor_aggressive = "score_tumor_aggressive",
  immune_maturity = "score_immune_hub_maturity"
)

now_iso <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}

write_status <- function(status, payload) {
  path <- file.path(project_root, "results", "logs", "stage_15_mvp_random_core_permutation_status.json")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  payload <- c(list(stage = "15_mvp_random_core_permutation", status = status, timestamp_utc = now_iso()), payload)
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

nearest_distance <- function(coords, core_coords) {
  FNN::get.knnx(data = core_coords, query = coords, k = 1)$nn.dist[, 1]
}

median_neighbor_distance <- function(coords) {
  if (nrow(coords) < 3) return(NA_real_)
  sample_idx <- seq_len(nrow(coords))
  if (nrow(coords) > 1500) sample_idx <- sample(sample_idx, 1500)
  sampled <- coords[sample_idx, , drop = FALSE]
  median(FNN::get.knnx(data = sampled, query = sampled, k = 2)$nn.dist[, 2], na.rm = TRUE)
}

analyze_one_sample <- function(dt, n_perm) {
  if ("edge_or_background_risk" %in% names(dt)) {
    dt <- dt[!(as.logical(edge_or_background_risk) %in% TRUE)]
  }
  dt <- dt[is.finite(x_pixel) & is.finite(y_pixel)]
  if (nrow(dt) < 100) return(NULL)
  coords <- as.matrix(dt[, .(x_pixel, y_pixel)])
  nn <- median_neighbor_distance(coords)
  if (!is.finite(nn) || nn <= 0) return(NULL)
  caf <- dt$score_caf_myeloid_barrier
  threshold <- quantile(caf, 0.9, na.rm = TRUE, names = FALSE)
  core_idx <- which(caf >= threshold)
  if (length(core_idx) < 10) return(NULL)

  real_dist <- nearest_distance(coords, coords[core_idx, , drop = FALSE]) / nn
  target_values <- lapply(target_scores, function(col) dt[[col]])
  real_rho <- vapply(target_values, function(values) safe_spearman(real_dist, values), numeric(1))

  null_rho <- matrix(NA_real_, nrow = n_perm, ncol = length(target_scores))
  colnames(null_rho) <- names(target_scores)
  all_idx <- seq_len(nrow(dt))
  n_core <- length(core_idx)
  for (perm_id in seq_len(n_perm)) {
    random_idx <- sample(all_idx, n_core, replace = FALSE)
    random_dist <- nearest_distance(coords, coords[random_idx, , drop = FALSE]) / nn
    null_rho[perm_id, ] <- vapply(target_values, function(values) safe_spearman(random_dist, values), numeric(1))
  }

  dataset_id <- dt$dataset_id[[1]]
  sample_id <- dt$sample_id[[1]]
  specimen_type <- if ("specimen_type" %in% names(dt)) dt$specimen_type[[1]] else "metadata_required"
  sample_rows <- rbindlist(lapply(names(target_scores), function(target) {
    null_values <- null_rho[, target]
    observed <- real_rho[[target]]
    n_finite <- sum(is.finite(null_values))
    data.table(
      dataset_id = dataset_id,
      sample_id = sample_id,
      specimen_type = specimen_type,
      target = target,
      n_safe_spots = nrow(dt),
      n_core_spots = n_core,
      n_perm = n_perm,
      observed_rho = observed,
      null_median_rho = median(null_values, na.rm = TRUE),
      null_p05_rho = quantile(null_values, 0.05, na.rm = TRUE, names = FALSE),
      null_p95_rho = quantile(null_values, 0.95, na.rm = TRUE, names = FALSE),
      delta_vs_null_median = observed - median(null_values, na.rm = TRUE),
      empirical_p_more_negative = (1 + sum(null_values <= observed, na.rm = TRUE)) / (1 + n_finite)
    )
  }))
  null_rows <- rbindlist(lapply(names(target_scores), function(target) {
    data.table(
      dataset_id = dataset_id,
      sample_id = sample_id,
      specimen_type = specimen_type,
      target = target,
      perm_id = seq_len(n_perm),
      null_rho = null_rho[, target]
    )
  }))
  list(sample = sample_rows, null = null_rows)
}

summarize_results <- function(sample_results) {
  sample_results[, q_value_empirical := p.adjust(empirical_p_more_negative, method = "BH")]
  sample_results[
    ,
    .(
      n_samples = .N,
      median_observed_rho = median(observed_rho, na.rm = TRUE),
      median_null_rho = median(null_median_rho, na.rm = TRUE),
      median_delta_vs_null = median(delta_vs_null_median, na.rm = TRUE),
      n_observed_more_negative_than_null = sum(delta_vs_null_median < 0, na.rm = TRUE),
      n_empirical_p_lt_0_05 = sum(empirical_p_more_negative < 0.05, na.rm = TRUE),
      median_empirical_p_more_negative = median(empirical_p_more_negative, na.rm = TRUE),
      min_q_value_empirical = min(q_value_empirical, na.rm = TRUE)
    ),
    by = .(dataset_id, target)
  ][order(dataset_id, target)]
}

plot_results <- function(sample_results, summary_results) {
  out_dir <- file.path(project_root, "results", "figures", "mvp", "random_core_permutation")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  p1 <- ggplot(sample_results, aes(x = dataset_id, y = delta_vs_null_median, fill = dataset_id)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey30") +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.1, alpha = 0.45) +
    facet_wrap(~target, scales = "free_y") +
    labs(
      title = "MVP cohorts random-core control",
      subtitle = "Edge-QC-safe spots; negative delta means true CAF-myeloid core beats random same-size cores",
      x = NULL,
      y = "Observed rho minus random-core median rho"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  ggsave(file.path(out_dir, "mvp_random_core_permutation_delta.png"), p1, width = 10.5, height = 6, dpi = 180)
  ggsave(file.path(out_dir, "mvp_random_core_permutation_delta.pdf"), p1, width = 10.5, height = 6)

  p2 <- ggplot(summary_results, aes(x = dataset_id, y = n_observed_more_negative_than_null, fill = target)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = paste0(n_observed_more_negative_than_null, "/", n_samples)), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    labs(title = "Samples where true CAF-myeloid core beats random cores", x = NULL, y = "Number of samples") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(out_dir, "mvp_random_core_permutation_counts.png"), p2, width = 10.5, height = 5.5, dpi = 180)
  ggsave(file.path(out_dir, "mvp_random_core_permutation_counts.pdf"), p2, width = 10.5, height = 5.5)
}

main <- function() {
  set.seed(20260624)
  spot_path <- file.path(project_root, "results", "tables", "mvp_spot_level_scores_with_edge_qc.csv")
  if (!file.exists(spot_path)) stop("Missing edge-QC spot-level score table: ", spot_path)
  message("Reading ", spot_path)
  spot <- fread(spot_path)
  required <- c("dataset_id", "sample_id", "x_pixel", "y_pixel", "score_caf_myeloid_barrier", target_scores)
  missing <- setdiff(required, names(spot))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

  sample_keys <- unique(spot[, .(dataset_id, sample_id)])
  message("Running random-core permutation for ", nrow(sample_keys), " samples; n_perm=", n_perm)
  results <- vector("list", nrow(sample_keys))
  for (i in seq_len(nrow(sample_keys))) {
    dataset <- sample_keys$dataset_id[[i]]
    sid <- sample_keys$sample_id[[i]]
    message(sprintf("[%d/%d] %s %s", i, nrow(sample_keys), dataset, sid))
    results[[i]] <- analyze_one_sample(spot[dataset_id == dataset & sample_id == sid], n_perm = n_perm)
  }
  results <- results[!vapply(results, is.null, logical(1))]
  sample_results <- rbindlist(lapply(results, `[[`, "sample"))
  null_results <- rbindlist(lapply(results, `[[`, "null"))
  summary_results <- summarize_results(sample_results)

  out_dir <- file.path(project_root, "results", "tables")
  fwrite(sample_results, file.path(out_dir, "mvp_random_core_permutation_sample_stats.csv"))
  fwrite(summary_results, file.path(out_dir, "mvp_random_core_permutation_summary.csv"))
  fwrite(null_results, file.path(out_dir, "mvp_random_core_permutation_null_rhos.csv"))
  plot_results(sample_results, summary_results)

  write_status("success", list(
    n_samples = length(unique(sample_results$sample_id)),
    n_perm = n_perm,
    output_sample_stats = file.path(out_dir, "mvp_random_core_permutation_sample_stats.csv"),
    output_summary = file.path(out_dir, "mvp_random_core_permutation_summary.csv")
  ))
  message("Done.")
}

tryCatch(
  main(),
  error = function(e) {
    write_status("failed", list(error = conditionMessage(e), n_perm = n_perm))
    stop(e)
  }
)

