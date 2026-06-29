#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"

specimen_order <- c("normal_pancreas", "primary_tumor", "liver_metastasis", "lymph_node_metastasis")
specimen_labels <- c(
  normal_pancreas = "Normal\npancreas",
  primary_tumor = "Primary\ntumor",
  liver_metastasis = "Liver\nmetastasis",
  lymph_node_metastasis = "Lymph node\nmetastasis"
)
specimen_colors <- c(
  normal_pancreas = "#4C78A8",
  primary_tumor = "#F58518",
  liver_metastasis = "#54A24B",
  lymph_node_metastasis = "#B279A2"
)

read_table <- function(name) {
  read.csv(file.path(project_root, "results", "tables", name), stringsAsFactors = FALSE)
}

sample_summary <- read_table("gse272362_rds_sample_specimen_summary.csv")
gradients <- read_table("gse272362_rds_caf_myeloid_gradient_stats.csv")
bins <- read_table("gse272362_rds_caf_myeloid_distance_bins.csv")

plot_story <- function() {
  par(mfrow = c(1, 3), mar = c(7, 4.5, 3, 1), oma = c(0, 0, 2, 0))

  present <- specimen_order[specimen_order %in% sample_summary$specimen_type]
  data <- lapply(present, function(x) sample_summary$p90_caf_myeloid_barrier[sample_summary$specimen_type == x])
  boxplot(
    data,
    names = specimen_labels[present],
    las = 1,
    col = specimen_colors[present],
    outline = FALSE,
    ylab = "Sample p90 score",
    main = "CAF-myeloid high-score tail"
  )
  stripchart(data, vertical = TRUE, method = "jitter", pch = 16, col = "#22222299", add = TRUE)

  grad <- gradients[gradients$target %in% c("immune_core", "ifn_mhc", "tumor_aggressive"), ]
  data <- lapply(present, function(x) grad$rho_distance_to_caf_core[grad$specimen_type == x])
  boxplot(
    data,
    names = specimen_labels[present],
    las = 1,
    col = specimen_colors[present],
    outline = FALSE,
    ylab = "Spearman rho",
    main = "Distance to CAF-myeloid core"
  )
  abline(h = 0, lwd = 1, col = "#333333")
  mtext("Negative = co-localized with CAF-myeloid core", side = 3, line = 0.1, cex = 0.75)

  bin_order <- c("core_0_1.5", "near_1.5_3", "mid_3_6", "far_gt6")
  plot(
    seq_along(bin_order),
    rep(NA_real_, length(bin_order)),
    ylim = range(bins$mean_score[bins$target == "ifn_mhc"], na.rm = TRUE),
    xaxt = "n",
    xlab = "",
    ylab = "Median sample-level mean score",
    main = "IFN/MHC around CAF-myeloid core"
  )
  axis(1, at = seq_along(bin_order), labels = c("Core", "Near", "Mid", "Far"), las = 1)
  for (specimen in present) {
    sub <- bins[bins$specimen_type == specimen & bins$target == "ifn_mhc", ]
    y <- sapply(bin_order, function(bin) {
      values <- sub$mean_score[sub$distance_bin == bin]
      if (length(values) == 0) NA_real_ else median(values, na.rm = TRUE)
    })
    lines(seq_along(bin_order), y, type = "b", pch = 16, lwd = 2, col = specimen_colors[specimen])
  }
  legend("topright", legend = gsub("_", " ", present), col = specimen_colors[present], lwd = 2, pch = 16, bty = "n", cex = 0.75)

  mtext("GSE272362 PDAC_Updated.rds spatial ecology", outer = TRUE, cex = 1.2, font = 2)
}

out_dir <- file.path(project_root, "results", "figures", "mvp", "gse272362_rds")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(out_dir, "gse272362_specimen_ecology.png"), width = 2400, height = 800, res = 180)
plot_story()
dev.off()

pdf(file.path(out_dir, "gse272362_specimen_ecology.pdf"), width = 13.5, height = 4.5)
plot_story()
dev.off()

cat("Wrote GSE272362 story figures\n")

