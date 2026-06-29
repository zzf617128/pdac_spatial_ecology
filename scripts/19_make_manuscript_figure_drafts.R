#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(png)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"

target_labels <- c(
  ifn_mhc = "IFN/MHC",
  immune_core = "Immune core",
  immune_maturity = "Immune maturity-like",
  tumor_aggressive = "Tumor-aggressive"
)
site_labels <- c(
  primary_tumor = "Primary tumor",
  liver_metastasis = "Liver metastasis",
  lymph_node_metastasis = "Lymph node metastasis",
  normal_pancreas = "Normal pancreas"
)
site_palette <- c(
  "Primary tumor" = "#1B9E77",
  "Liver metastasis" = "#D95F02",
  "Lymph node metastasis" = "#7570B3",
  "Normal pancreas" = "#666666"
)
cohort_palette <- c("GSE282302" = "#1F78B4", "GSE274103" = "#33A02C")

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
      axis.title = element_text(size = base_size),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

label_targets <- function(x) {
  factor(target_labels[x], levels = target_labels[c("ifn_mhc", "immune_core", "tumor_aggressive", "immune_maturity")])
}

panel_label <- function(label, grob) {
  arrangeGrob(
    textGrob(label, x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
             gp = gpar(fontface = "bold", fontsize = 14)),
    grob,
    heights = unit.c(unit(0.28, "in"), unit(1, "null"))
  )
}

raster_panel <- function(path, title, subtitle = NULL) {
  if (!file.exists(path)) {
    return(textGrob(paste("Missing image:", basename(path)), gp = gpar(col = "red3")))
  }
  img <- readPNG(path)
  image_grob <- rasterGrob(img, interpolate = TRUE)
  title_grob <- textGrob(title, x = 0, just = "left", gp = gpar(fontface = "bold", fontsize = 10))
  if (is.null(subtitle)) {
    arrangeGrob(title_grob, image_grob, heights = unit.c(unit(0.24, "in"), unit(1, "null")))
  } else {
    arrangeGrob(
      title_grob,
      textGrob(subtitle, x = 0, just = "left", gp = gpar(fontsize = 8.5, col = "grey25")),
      image_grob,
      heights = unit.c(unit(0.22, "in"), unit(0.20, "in"), unit(1, "null"))
    )
  }
}

cohort_overview_panel <- function() {
  dt <- data.table(
    cohort = c("GSE282302", "GSE274103", "GSE272362 RDS"),
    context = c("Post-neoadjuvant PDAC ST-H&E", "Treatment-naive PDAC ST-H&E", "PDAC primary, liver, lymph node and normal pancreas ST-H&E"),
    samples = c("108 sections", "5 sections", "30 samples"),
    role = c("Discovery / main evidence", "Treatment-naive support", "Independent site validation")
  )
  dt[, y := rev(seq_len(.N))]
  ggplot(dt, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = y - 0.42, ymax = y + 0.42),
              fill = "#F4F4F1", color = "#BDBDBD", linewidth = 0.35) +
    geom_text(aes(x = 0.04, label = cohort), hjust = 0, fontface = "bold", size = 3.2) +
    geom_text(aes(x = 0.27, label = samples), hjust = 0, size = 2.8, color = "#333333") +
    geom_text(aes(x = 0.45, label = role), hjust = 0, size = 2.8, color = "#333333") +
    geom_text(aes(x = 0.04, y = y - 0.20, label = context), hjust = 0, size = 2.5, color = "#555555") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 3.55), expand = FALSE) +
    labs(title = "Study design", subtitle = "CAF-myeloid core gradients tested across cohorts and tissue sites") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey25"))
}

random_delta_panel <- function(dt, title, subtitle, color_field = "dataset_id") {
  plot_dt <- copy(dt)
  plot_dt <- plot_dt[target %in% names(target_labels)]
  plot_dt[, target_label := label_targets(target)]
  plot_dt[, support_label := paste0(n_observed_more_negative_than_null, "/", n_samples)]
  plot_dt[, group_label := if ("specimen_type" %in% names(plot_dt)) site_labels[specimen_type] else dataset_id]
  plot_dt[is.na(group_label), group_label := dataset_id]
  ggplot(plot_dt, aes(x = median_delta_vs_null, y = target_label, fill = group_label)) +
    geom_vline(xintercept = 0, color = "grey35", linewidth = 0.35) +
    geom_col(position = position_dodge2(width = 0.75), width = 0.62, color = "white", linewidth = 0.15) +
    geom_text(
      data = plot_dt[median_delta_vs_null < 0],
      aes(x = 0, label = support_label),
      position = position_dodge2(width = 0.75), hjust = 1.12, color = "white", size = 2.3
    ) +
    geom_text(
      data = plot_dt[median_delta_vs_null >= 0],
      aes(label = support_label),
      position = position_dodge2(width = 0.75), hjust = -0.12, color = "grey25", size = 2.3
    ) +
    scale_fill_manual(values = c(cohort_palette, site_palette), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    labs(title = title, subtitle = subtitle, x = "Median delta vs random core", y = NULL, fill = NULL) +
    theme_paper()
}

sensitivity_panel <- function(dt, keep_groups, title, subtitle) {
  plot_dt <- copy(dt)
  plot_dt <- plot_dt[target %in% names(target_labels)]
  plot_dt[, target_label := label_targets(target)]
  plot_dt[, group_label := ifelse(cohort_label == "gse272362", site_labels[specimen_type], dataset_id)]
  plot_dt <- plot_dt[group_label %in% keep_groups]
  plot_dt[, core_label := factor(core_label, levels = c("top_15", "top_10", "top_5"),
                                 labels = c("Top 15%", "Top 10%", "Top 5%"))]
  ggplot(plot_dt, aes(x = core_label, y = median_rho, group = group_label, color = group_label)) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.35) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.7) +
    facet_wrap(~target_label, scales = "free_y", nrow = 2) +
    scale_color_manual(values = c(cohort_palette, site_palette), drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = "CAF-myeloid core definition", y = "Median rho: distance to core vs score", color = NULL) +
    theme_paper(8.5)
}

site_counts_panel <- function(gse_scores) {
  dt <- unique(gse_scores[, .(sample_id, specimen_type)])
  spot_dt <- gse_scores[, .(n_spots = .N), by = .(specimen_type)]
  sample_dt <- dt[, .(n_samples = .N), by = specimen_type]
  plot_dt <- merge(sample_dt, spot_dt, by = "specimen_type")
  plot_dt[, site := site_labels[specimen_type]]
  plot_dt[, site := factor(site, levels = site_labels[c("primary_tumor", "liver_metastasis", "lymph_node_metastasis", "normal_pancreas")])]
  plot_dt[, label := paste0(n_samples, " samples\n", format(n_spots, big.mark = ","), " spots")]
  ggplot(plot_dt, aes(x = site, y = n_spots, fill = site)) +
    geom_col(width = 0.62, color = "white", linewidth = 0.2) +
    geom_text(aes(label = label), vjust = -0.15, size = 2.6) +
    scale_fill_manual(values = site_palette, drop = FALSE) +
    labs(title = "Independent RDS cohort", subtitle = "GSE272362 contains matched spatial profiles across PDAC tissue sites", x = NULL, y = "Spots") +
    theme_paper(9) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 18, hjust = 1))
}

write_sources <- function(mvp_random, gse_random, sensitivity, gse_scores, out_dir) {
  figure1_source <- rbindlist(list(
    data.table(figure = "Figure 1", panel = "A", source_file = "manual cohort summary",
               metric = "cohort/sample counts", note = "GSE282302: 108 sections; GSE274103: 5 sections; GSE272362 RDS: 30 samples"),
    mvp_random[, .(figure = "Figure 1", panel = "B", source_file = "mvp_random_core_permutation_summary.csv",
                   metric = paste0(target, " delta_vs_null"), group = dataset_id,
                   n_samples, value = median_delta_vs_null, support = paste0(n_observed_more_negative_than_null, "/", n_samples))],
    sensitivity[cohort_label == "mvp", .(figure = "Figure 1", panel = "C", source_file = "caf_core_threshold_sensitivity_summary.csv",
                   metric = paste0(target, " median_rho"), group = dataset_id, core_label, n_samples, value = median_rho)]
  ), fill = TRUE)

  site_counts <- merge(
    unique(gse_scores[, .(sample_id, specimen_type)])[, .(n_samples = .N), by = specimen_type],
    gse_scores[, .(n_spots = .N), by = specimen_type],
    by = "specimen_type"
  )
  figure2_source <- rbindlist(list(
    site_counts[, .(figure = "Figure 2", panel = "A", source_file = "gse272362_rds_spot_level_scores.csv",
                    metric = "site sample and spot counts", group = specimen_type, n_samples, value = n_spots)],
    gse_random[specimen_type != "normal_pancreas", .(figure = "Figure 2", panel = "B", source_file = "gse272362_rds_random_core_permutation_summary.csv",
                    metric = paste0(target, " delta_vs_null"), group = specimen_type,
                    n_samples, value = median_delta_vs_null, support = paste0(n_observed_more_negative_than_null, "/", n_samples))],
    sensitivity[cohort_label == "gse272362" & specimen_type != "normal_pancreas",
                .(figure = "Figure 2", panel = "C", source_file = "caf_core_threshold_sensitivity_summary.csv",
                  metric = paste0(target, " median_rho"), group = specimen_type, core_label, n_samples, value = median_rho)]
  ), fill = TRUE)
  fwrite(figure1_source, file.path(out_dir, "figure1_source.csv"))
  fwrite(figure2_source, file.path(out_dir, "figure2_source.csv"))
}

main <- function() {
  fig_dir <- file.path(project_root, "results", "figures", "main")
  table_dir <- file.path(project_root, "results", "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  mvp_random <- fread(file.path(table_dir, "mvp_random_core_permutation_summary.csv"))
  gse_random <- fread(file.path(table_dir, "gse272362_rds_random_core_permutation_summary.csv"))
  sensitivity <- fread(file.path(table_dir, "caf_core_threshold_sensitivity_summary.csv"))
  gse_scores <- fread(file.path(table_dir, "gse272362_rds_spot_level_scores.csv"))
  gse_scores[, specimen_type := fifelse(grepl("normal", tolower(specimen_type)), "normal_pancreas",
                                fifelse(grepl("liver", tolower(specimen_type)), "liver_metastasis",
                                fifelse(grepl("lymph|node|ln", tolower(specimen_type)), "lymph_node_metastasis", "primary_tumor")))]

  write_sources(mvp_random, gse_random, sensitivity, gse_scores, table_dir)

  fig1 <- arrangeGrob(
    panel_label("A", ggplotGrob(cohort_overview_panel())),
    panel_label("B", ggplotGrob(random_delta_panel(
      mvp_random,
      "Main cohorts exceed random cores",
      "Labels show samples where the observed CAF core is more negative than matched random cores"
    ))),
    panel_label("C", ggplotGrob(sensitivity_panel(
      sensitivity,
      keep_groups = c("GSE282302", "GSE274103"),
      "CAF-core threshold sensitivity",
      "Distance-to-core gradients are directionally stable across top 15%, 10% and 5% definitions"
    ))),
    panel_label("D", raster_panel(
      file.path(project_root, "results", "figures", "mvp", "he_overlays", "GSE282302", "GSM8641105_C3_D8_ROI3_he_overlay.png"),
      "Representative post-neoadjuvant section",
      "Spatial overlays summarize CAF core and target programs on H&E"
    )),
    layout_matrix = rbind(c(1, 2), c(3, 3), c(4, 4)),
    widths = c(1.05, 1.15),
    heights = c(0.70, 1.10, 0.70)
  )

  ggsave(file.path(fig_dir, "figure1_draft.pdf"), fig1, width = 12, height = 11, limitsize = FALSE, bg = "white")
  ggsave(file.path(fig_dir, "figure1_draft.png"), fig1, width = 12, height = 11, dpi = 220, limitsize = FALSE, bg = "white")
  ggsave(file.path(fig_dir, "figure1_main.pdf"), fig1, width = 12, height = 11, limitsize = FALSE, bg = "white")
  ggsave(file.path(fig_dir, "figure1_main.png"), fig1, width = 12, height = 11, dpi = 300, limitsize = FALSE, bg = "white")

  gse_plot_dt <- gse_random[specimen_type != "normal_pancreas"]
  fig2 <- arrangeGrob(
    panel_label("A", ggplotGrob(site_counts_panel(gse_scores))),
    panel_label("B", ggplotGrob(random_delta_panel(
      gse_plot_dt,
      "RDS validation by tissue site",
      "Primary and liver reproduce immune/IFN gradients; lymph node retains tumor-aggressive CAF association"
    ))),
    panel_label("C", raster_panel(
      file.path(project_root, "results", "figures", "mvp", "gse272362_rds", "he_overlays", "IU_PDA_HM10_liver_strong_ifn_caf_core_he_overlay.png"),
      "Liver metastasis example",
      "Strong IFN/MHC gradient around CAF-myeloid core"
    )),
    panel_label("D", raster_panel(
      file.path(project_root, "results", "figures", "mvp", "gse272362_rds", "he_overlays", "IU_PDA_T1_primary_strong_ifn_caf_core_he_overlay.png"),
      "Primary tumor example",
      "CAF-core-centered IFN/MHC organization"
    )),
    panel_label("E", raster_panel(
      file.path(project_root, "results", "figures", "mvp", "gse272362_rds", "he_overlays", "IU_PDA_LNM7_lnm_immune_divergent_he_overlay.png"),
      "Lymph node divergent example",
      "Immune-core signal is not consistently CAF-core centered"
    )),
    layout_matrix = rbind(c(1, 2), c(3, 3), c(4, 4), c(5, 5)),
    widths = c(1, 1.65),
    heights = c(0.74, 0.64, 0.64, 0.64)
  )

  ggsave(file.path(fig_dir, "figure2_draft.pdf"), fig2, width = 14, height = 15, limitsize = FALSE, bg = "white")
  ggsave(file.path(fig_dir, "figure2_draft.png"), fig2, width = 14, height = 15, dpi = 220, limitsize = FALSE, bg = "white")
  ggsave(file.path(fig_dir, "figure2_main.pdf"), fig2, width = 14, height = 15, limitsize = FALSE, bg = "white")
  ggsave(file.path(fig_dir, "figure2_main.png"), fig2, width = 14, height = 15, dpi = 300, limitsize = FALSE, bg = "white")

  manifest <- data.table(
    figure = c("Figure 1", "Figure 1", "Figure 1", "Figure 1", "Figure 2", "Figure 2", "Figure 2", "Figure 2"),
    artifact = c(
      "figure1_draft.pdf", "figure1_draft.png", "figure1_main.pdf", "figure1_main.png",
      "figure2_draft.pdf", "figure2_draft.png", "figure2_main.pdf", "figure2_main.png"
    ),
    path = file.path(fig_dir, c(
      "figure1_draft.pdf", "figure1_draft.png", "figure1_main.pdf", "figure1_main.png",
      "figure2_draft.pdf", "figure2_draft.png", "figure2_main.pdf", "figure2_main.png"
    ))
  )
  fwrite(manifest, file.path(table_dir, "manuscript_figure_draft_manifest.csv"))
  message("Figure drafts written to ", fig_dir)
}

main()
