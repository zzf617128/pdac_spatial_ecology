suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

root <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search"), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
doc_dir <- file.path(root, "docs")
source_dir <- file.path(root, "source_data")
manifest_dir <- file.path(root, "manifest")
script_path <- file.path(root, "scripts/make_extended_data_fig10_v1.R")
invisible(lapply(c(table_dir, figure_dir, doc_dir, source_dir, manifest_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, path) fwrite(as.data.frame(x), path)
fmt <- function(x, digits = 2) ifelse(is.na(x), "NA", formatC(as.numeric(x), digits = digits, format = "f"))
cap01 <- function(x) pmax(0, pmin(1, x))

theme_ed <- function(base_size = 7.4) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.3, colour = "grey15"),
      plot.title = element_text(size = base_size + 1.0, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "grey30", hjust = 0),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.size = unit(0.35, "lines"),
      plot.margin = margin(3, 5, 3, 5)
    )
}

panel_tag <- function(p, tag) {
  p + plot_annotation(tag_levels = list(tag)) &
    theme(plot.tag = element_text(face = "bold", size = 11), plot.tag.position = c(0, 1))
}

module_order <- c(
  "caf_mycaf", "pancaf_matrix", "myeloid_macrophage", "immune_core",
  "ifn_apc", "tgfb_emt", "tumor_epithelial", "tumor_aggressive"
)
module_labels <- c(
  caf_mycaf = "CAF/myCAF",
  pancaf_matrix = "panCAF/matrix",
  myeloid_macrophage = "Myeloid/macrophage",
  immune_core = "Immune-core",
  ifn_apc = "IFN/APC",
  tgfb_emt = "TGF/EMT",
  tumor_epithelial = "Tumor epithelial",
  tumor_aggressive = "Tumor-aggressive",
  t_cell = "T cell",
  spp1_tam = "SPP1/TAM"
)

g240_path <- file.path(getwd(), "pdac_spatial_ecology/outputs/orthogonal_validation_2026_06_30/tables/geomx_gse240078_tme_vs_carcinoma_module_tests.csv")
g199_controls_path <- file.path(table_dir, "gse199102_geomx_positive_negative_control_results.csv")
g199_paired_path <- file.path(table_dir, "gse199102_geomx_paired_segment_results.csv")
concord_path <- file.path(table_dir, "gse199102_gse240078_direction_concordance.csv")
cosmx_nonoverlap_path <- file.path(table_dir, "gse310352_nonoverlap_tgfemt_adjacency_results.csv")
cosmx_threshold_path <- file.path(table_dir, "gse310352_threshold_sensitivity.csv")
cosmx_norm_path <- file.path(table_dir, "gse310352_normalization_sensitivity.csv")
cosmx_null_path <- file.path(table_dir, "gse310352_spatial_null_sensitivity.csv")
cosmx_loso_path <- file.path(table_dir, "gse310352_leave_one_slide_out.csv")

required <- c(g240_path, g199_controls_path, g199_paired_path, concord_path, cosmx_nonoverlap_path, cosmx_threshold_path, cosmx_norm_path, cosmx_null_path, cosmx_loso_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) stop("Missing required input files: ", paste(missing, collapse = "; "))

g240 <- fread(g240_path)
g199 <- fread(g199_controls_path)
paired <- fread(g199_paired_path)
concord <- fread(concord_path)
cosmx_nonoverlap <- fread(cosmx_nonoverlap_path)
cosmx_threshold <- fread(cosmx_threshold_path)
cosmx_norm <- fread(cosmx_norm_path)
cosmx_null <- fread(cosmx_null_path)
cosmx_loso <- fread(cosmx_loso_path)

g240_plot <- g240[module %in% module_order]
g240_plot[, module_label := factor(module_labels[module], levels = rev(module_labels[module_order]))]
g240_plot[, expected := fifelse(module %in% c("tumor_epithelial", "tumor_aggressive"), "Tumor-high control", "TME-high program")]

g199_plot <- g199[module %in% module_order]
g199_plot[, module_label := factor(module_labels[module], levels = rev(module_labels[module_order]))]
g199_plot[, expected := fifelse(module %in% c("tumor_epithelial", "tumor_aggressive"), "Tumor-high control", "TME-high program")]

g240_spp1 <- g240[module == "spp1_tam", .(module, delta = delta_stroma_minus_tumor, fdr)]
g199_spp1 <- g199[module == "spp1_tam", .(module, delta = median_tme_minus_tumor, support = patient_support_fraction_tme_gt_tumor)]

concord_keep <- concord[module %in% c(module_order, "t_cell", "spp1_tam")]
concord_keep[, module_label := factor(module_labels[module], levels = rev(module_labels[c(module_order, "t_cell", "spp1_tam")]))]
concord_long <- melt(
  concord_keep,
  id.vars = c("module", "module_label", "direction_concordant"),
  measure.vars = c("gse240078_delta_stroma_minus_tumor", "gse199102_median_tme_minus_tumor"),
  variable.name = "dataset",
  value.name = "delta"
)
concord_long[, dataset := fifelse(dataset == "gse240078_delta_stroma_minus_tumor", "GSE240078\nDSP", "GSE199102\nWTA")]
concord_long[, direction := fifelse(delta > 0, "TME/stroma high", "Tumor high")]
concord_long[module == "spp1_tam", direction := "Tumor high\n(boundary)"]

paired_plot <- paired[comparison == "paired_CAF_minus_Epithelial" & module %in% c("caf_mycaf", "pancaf_matrix", "tgfb_emt", "immune_core", "ifn_apc")]
paired_plot[, module_label := factor(module_labels[module], levels = rev(module_labels[c("caf_mycaf", "pancaf_matrix", "tgfb_emt", "immune_core", "ifn_apc")]))]

cosmx_f <- copy(cosmx_nonoverlap)
cosmx_f[, tgfb_label := fifelse(tgfb_module == "tgfb_emt", "Original TGF/EMT",
  fifelse(tgfb_module == "tgfb_emt_no_shared", "Remove MMP2", "Remove MMP2 + ITGA5"))]
cosmx_f[, tgfb_label := factor(tgfb_label, levels = c("Original TGF/EMT", "Remove MMP2", "Remove MMP2 + ITGA5"))]
cosmx_f[, slide := factor(slide, levels = sort(unique(as.character(slide))))]
cosmx_f_summary <- cosmx_f[, .(
  median_slide_log2_oe = median(median_log2_oe, na.rm = TRUE),
  slide_support_gt0 = sum(median_log2_oe > 0, na.rm = TRUE),
  n_slides = .N
), by = tgfb_label]

threshold_sum <- cosmx_threshold[, .(
  support = mean(median_log2_oe > 0, na.rm = TRUE),
  median_log2_oe = median(median_log2_oe, na.rm = TRUE),
  n_tests = .N
), by = .(setting = threshold)]
threshold_sum[, group := "Threshold"]
norm_sum <- cosmx_norm[, .(
  support = mean(median_log2_oe > 0, na.rm = TRUE),
  median_log2_oe = median(median_log2_oe, na.rm = TRUE),
  n_tests = .N
), by = .(setting = normalization)]
norm_sum[, group := "Normalization"]
null_sum <- data.table(
  group = "Spatial null",
  setting = c("Within-FOV", "Within-slide", "Abundance-matched"),
  support = c(
    mean(cosmx_null$median_within_fov_null_p <= 0.05, na.rm = TRUE),
    mean(cosmx_null$median_within_slide_null_p <= 0.05, na.rm = TRUE),
    mean(cosmx_null$median_abundance_random_target_p <= 0.05, na.rm = TRUE)
  ),
  median_log2_oe = median(cosmx_null$median_log2_oe, na.rm = TRUE),
  n_tests = nrow(cosmx_null)
)
loso_sum <- data.table(
  group = "Leave-one-slide-out",
  setting = "All omissions",
  support = mean(cosmx_loso$slide_support_fraction_gt0 == 1, na.rm = TRUE),
  median_log2_oe = median(cosmx_loso$median_slide_log2_oe, na.rm = TRUE),
  n_tests = nrow(cosmx_loso)
)
robust_sum <- rbindlist(list(threshold_sum, norm_sum, null_sum, loso_sum), fill = TRUE)
robust_sum[, label := paste0(round(100 * support), "%")]
robust_sum[, group := factor(group, levels = c("Threshold", "Normalization", "Spatial null", "Leave-one-slide-out"))]

source_rows <- rbindlist(list(
  data.table(panel = "A", dataset = "Conceptual schematic", module = "orthogonal_validation_layers", module_label = "Discovery to GeoMx/CosMx", metric = "schematic", value = NA_real_, support = NA_real_, note = "Visium/Xenium discovery model; GeoMx compartment-level validation; CosMx stromal-interface validation"),
  g240_plot[, .(panel = "B", dataset = "GSE240078 GeoMx DSP", module, module_label = as.character(module_label), metric = "delta_TME_minus_carcinoma", value = delta_stroma_minus_tumor, support = NA_real_, note = expected)],
  g199_plot[, .(panel = "C", dataset = "GSE199102 GeoMx WTA", module, module_label = as.character(module_label), metric = "delta_TME_like_minus_epithelial", value = median_tme_minus_tumor, support = patient_support_fraction_tme_gt_tumor, note = expected)],
  concord_long[, .(panel = "D", dataset, module, module_label = as.character(module_label), metric = "delta_for_direction_concordance", value = delta, support = NA_real_, note = direction)],
  paired_plot[, .(panel = "E", dataset = "GSE199102 GeoMx WTA", module, module_label = as.character(module_label), metric = "paired_CAF_minus_epithelial_positive_fraction", value = support_fraction_positive, support = support_fraction_positive, note = paste0("n_roi=", n_roi, "; median_delta=", fmt(median_delta, 3)))],
  cosmx_f[, .(panel = "F", dataset = "GSE310352 CosMx", module = as.character(tgfb_module), module_label = as.character(tgfb_label), metric = "slide_level_log2_observed_expected", value = median_log2_oe, support = support_fraction_gt0, note = as.character(slide))],
  robust_sum[, .(panel = "G", dataset = "GSE310352 CosMx", module = as.character(setting), module_label = as.character(group), metric = "robustness_support_fraction", value = support, support = support, note = paste0("median_log2OE=", fmt(median_log2_oe, 3), "; n=", n_tests))],
  data.table(panel = "H", dataset = "Inference scope", module = "claim_boundaries", module_label = "Supported and not claimed", metric = "scope_statement", value = NA_real_, support = NA_real_, note = "Compartment-level GeoMx and CosMx stromal-interface support only; no causality, direct SPP1-CD44, tumor-intrinsic EMT, Visium gradient reconstruction, or LN immune uncoupling"),
  data.table(panel = "D", dataset = "Boundary note", module = "spp1_tam", module_label = "SPP1/TAM", metric = "not_stromal_support", value = g199_spp1$delta, support = g199_spp1$support, note = paste0("GSE240078 delta=", fmt(g240_spp1$delta, 3), "; GSE199102 delta=", fmt(g199_spp1$delta, 3), "; report separately"))
), fill = TRUE)
write_csv(source_rows, file.path(source_dir, "Source_Data_Extended_Data_Figure_10_v1.csv"))

pal <- c("TME-high program" = "#26706B", "Tumor-high control" = "#B64B3C")
dir_pal <- c("TME/stroma high" = "#26706B", "Tumor high" = "#B64B3C", "Tumor high\n(boundary)" = "#A06A34")

pA <- ggplot() +
  annotate("rect", xmin = 0.04, xmax = 0.30, ymin = 0.42, ymax = 0.80, fill = "#E7EEF8", colour = "#4B647F", linewidth = 0.35) +
  annotate("rect", xmin = 0.38, xmax = 0.64, ymin = 0.42, ymax = 0.80, fill = "#E8F4EF", colour = "#2F6B58", linewidth = 0.35) +
  annotate("rect", xmin = 0.72, xmax = 0.98, ymin = 0.42, ymax = 0.80, fill = "#F8EFE6", colour = "#925E33", linewidth = 0.35) +
  annotate("segment", x = 0.31, xend = 0.37, y = 0.61, yend = 0.61, arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35, colour = "grey35") +
  annotate("segment", x = 0.65, xend = 0.71, y = 0.61, yend = 0.61, arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35, colour = "grey35") +
  annotate("text", x = 0.17, y = 0.67, label = "Discovery\nVisium/Xenium", size = 2.5, fontface = "bold") +
  annotate("text", x = 0.17, y = 0.51, label = "CAF-core spatial\nprograms", size = 2.1) +
  annotate("text", x = 0.51, y = 0.67, label = "GeoMx\nDSP + WTA", size = 2.5, fontface = "bold") +
  annotate("text", x = 0.51, y = 0.51, label = "Compartment-level\nconcordance", size = 2.1) +
  annotate("text", x = 0.85, y = 0.67, label = "CosMx\ncell-level", size = 2.5, fontface = "bold") +
  annotate("text", x = 0.85, y = 0.51, label = "CAF/matrix to\nTGF/EMT interface", size = 2.1) +
  annotate("text", x = 0.51, y = 0.23, label = "Orthogonal support only: no causal or direct SPP1-CD44 claim", size = 2.15, colour = "grey25") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  labs(title = "Orthogonal validation logic")

pB <- ggplot(g240_plot, aes(x = delta_stroma_minus_tumor, y = module_label, colour = expected)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_segment(aes(x = 0, xend = delta_stroma_minus_tumor, yend = module_label), linewidth = 0.55) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = pal) +
  theme_ed() +
  labs(title = "GSE240078 GeoMx DSP", subtitle = "TME AOI - carcinoma AOI", x = "Module score delta", y = NULL, colour = NULL) +
  theme(legend.position = "none")

pC <- ggplot(g199_plot, aes(x = median_tme_minus_tumor, y = module_label, colour = expected)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_segment(aes(x = 0, xend = median_tme_minus_tumor, yend = module_label), linewidth = 0.55) +
  geom_point(aes(size = patient_support_fraction_tme_gt_tumor), alpha = 0.95) +
  scale_colour_manual(values = pal) +
  scale_size_continuous(range = c(1.2, 2.8), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_ed() +
  labs(title = "GSE199102 GeoMx WTA", subtitle = "CAF+Immune segments - epithelial segments", x = "Patient-level median delta", y = NULL, colour = NULL, size = "Patient\nsupport") +
  theme(legend.position = "none")

pD <- ggplot(concord_long, aes(x = dataset, y = module_label, fill = direction)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = fmt(delta, 2)), size = 2.0, colour = "white", fontface = "bold") +
  scale_fill_manual(values = dir_pal, na.value = "grey80") +
  theme_ed() +
  theme(axis.text.x = element_text(size = 6.2), legend.position = "none") +
  labs(title = "Cross-GeoMx direction", subtitle = "Green=TME/stroma high; SPP1/TAM boundary", x = NULL, y = NULL, fill = NULL)

pE <- ggplot(paired_plot, aes(x = support_fraction_positive, y = module_label)) +
  geom_col(fill = "#3C7D72", width = 0.62) +
  geom_text(aes(label = paste0(round(100 * support_fraction_positive), "%")), hjust = -0.08, size = 2.15) +
  scale_x_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1), labels = c("0", "50", "100")) +
  theme_ed() +
  labs(title = "GSE199102 paired ROI", subtitle = "CAF segment > epithelial segment", x = "Positive paired ROIs (%)", y = NULL)

pF <- ggplot(cosmx_f, aes(x = tgfb_label, y = median_log2_oe)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_point(aes(colour = slide), size = 1.45, alpha = 0.9, position = position_jitter(width = 0.08, height = 0)) +
  geom_crossbar(data = cosmx_f_summary, aes(x = tgfb_label, y = median_slide_log2_oe, ymin = median_slide_log2_oe, ymax = median_slide_log2_oe), inherit.aes = FALSE, width = 0.55, linewidth = 0.35, colour = "black") +
  geom_text(data = cosmx_f_summary, aes(x = tgfb_label, y = median_slide_log2_oe + 0.16, label = paste0(slide_support_gt0, "/", n_slides)), inherit.aes = FALSE, size = 2.05) +
  scale_colour_viridis_d(option = "D", end = 0.82) +
  theme_ed() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none") +
  labs(title = "GSE310352 CosMx", subtitle = "CAF/matrix -> TGF/EMT interface adjacency", x = NULL, y = "Slide log2(O/E)")

pG <- ggplot(robust_sum, aes(x = setting, y = group, fill = support)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 2.0, colour = "white", fontface = "bold") +
  scale_fill_gradient(low = "#D6E4E1", high = "#26706B", limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
  theme_ed() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.9), legend.position = "none") +
  labs(title = "CosMx robustness summary", subtitle = "Support across sensitivity checks", x = NULL, y = NULL, fill = "Support")

scope_text <- paste(
  "Inference scope",
  "GeoMx: compartment-level CAF/matrix and immune/TME programs.",
  "CosMx: CAF/matrix-associated TGF/EMT stromal-interface organization.",
  "No causal signaling.",
  "No direct SPP1-CD44 validation.",
  "No tumor-intrinsic EMT claim.",
  "No Visium distance-gradient reconstruction.",
  "No LN immune-uncoupling validation.",
  sep = "\n"
)
pH <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.05, ymax = 0.95, fill = "#F6F6F4", colour = "#B8B8B2", linewidth = 0.35) +
  annotate("text", x = 0.07, y = 0.88, label = "Inference scope", hjust = 0, size = 2.65, fontface = "bold") +
  annotate("text", x = 0.07, y = 0.77, label = "Supported", hjust = 0, size = 2.25, fontface = "bold", colour = "#26706B") +
  annotate("text", x = 0.07, y = 0.67, label = "GeoMx compartment-level CAF/matrix and immune/TME programs\nCosMx CAF/matrix-associated TGF/EMT stromal-interface organization", hjust = 0, size = 1.95, lineheight = 0.95) +
  annotate("text", x = 0.07, y = 0.45, label = "Not claimed", hjust = 0, size = 2.25, fontface = "bold", colour = "#B64B3C") +
  annotate("text", x = 0.07, y = 0.27, label = "Causality; direct SPP1-CD44 validation; tumor-intrinsic EMT;\nVisium distance-gradient reconstruction; LN immune uncoupling.", hjust = 0, size = 1.95, lineheight = 0.95) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void()

fig <- (
  (pA | pB | pC) /
    (pD | pE | pF) /
    (pG | pH)
) +
  plot_layout(heights = c(1.0, 1.05, 0.92), widths = c(1.05, 1.05, 1.08)) +
  plot_annotation(
    title = "Extended Data Fig. 10 | Orthogonal support for CAF/matrix compartment and TGF/EMT stromal-interface organization",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
      plot.tag = element_text(size = 11, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

pdf_path <- file.path(figure_dir, "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.pdf")
svg_path <- file.path(figure_dir, "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.svg")
ggsave(pdf_path, fig, width = 13.5, height = 9.2, units = "in", device = cairo_pdf)
svg(svg_path, width = 13.5, height = 9.2, onefile = TRUE)
print(fig)
dev.off()

legend_lines <- c(
  "# Extended Data Fig. 10 Legend Draft",
  "",
  "Extended Data Fig. 10 | Orthogonal support for CAF/matrix compartment and TGF/EMT stromal-interface organization.",
  "",
  "(A) Schematic of the orthogonal validation strategy, linking the discovery Visium/Xenium CAF-core model to independent GeoMx compartment-level and CosMx cell-level validation layers.",
  "(B) GSE240078 GeoMx DSP TME AOI versus carcinoma AOI module enrichment. Positive values indicate TME-enriched modules; tumor epithelial and tumor-aggressive modules serve as tumor-high controls. SPP1/TAM is not included as stromal support.",
  "(C) GSE199102 GeoMx WTA TME-like segment enrichment, comparing CAF plus Immune segments against epithelial segments at the patient level. Dot size indicates patient support fraction.",
  "(D) Direction concordance between GSE240078 and GSE199102 for shared modules. SPP1/TAM is shown separately as tumor-high/boundary evidence and is not interpreted as stromal SPP1 support.",
  "(E) Paired ROI validation in GSE199102, showing the fraction of paired ROIs in which CAF segments exceed epithelial segments for selected modules.",
  "(F) GSE310352 CosMx cell-level adjacency of CAF/matrix-like cells to TGF/EMT stromal-interface states across slides, including the original TGF/EMT module and non-overlap sensitivity versions removing MMP2 or MMP2 plus ITGA5. Labels indicate positive slides out of eight.",
  "(G) Summary of CosMx sensitivity analyses across threshold, normalization, spatial-null and leave-one-slide-out checks.",
  "(H) Inference scope. The figure supports orthogonal compartment-level and stromal-interface plausibility, not causal signaling, direct SPP1-CD44 validation, tumor-intrinsic EMT, Visium distance-gradient reconstruction or lymph-node immune-uncoupling validation."
)
writeLines(legend_lines, file.path(doc_dir, "ed10_v1_figure_legend_draft.md"), useBytes = TRUE)

results_lines <- c(
  "# ED10 v1 Results Paragraph Draft",
  "",
  "Independent GeoMx datasets provided compartment-level support for the stromal programs identified in the discovery analysis. In GSE240078 DSP, TME AOIs were enriched for CAF/myCAF, panCAF/matrix, myeloid/macrophage, immune-core, IFN/APC and TGF/EMT modules, whereas tumor epithelial and tumor-aggressive controls were carcinoma-enriched. In GSE199102 WTA, CAF plus Immune segments showed concordant TME-like enrichment for the same stromal and immune programs across patients, with paired ROI support for CAF/myCAF, panCAF/matrix and TGF/EMT modules. Separately, GSE310352 CosMx supported a cell-level organization in which CAF/matrix-like cells neighbor TGF/EMT stromal-interface states, and this pattern remained stable after gene-overlap, threshold, normalization, spatial-null and leave-one-slide-out sensitivity checks. These analyses provide orthogonal spatial-protein and spatial-transcriptomic support for CAF/matrix compartment and stromal-interface organization, while remaining compartmental and observational."
)
writeLines(results_lines, file.path(doc_dir, "ed10_v1_results_paragraph_draft.md"), useBytes = TRUE)

methods_lines <- c(
  "# ED10 v1 Methods Draft",
  "",
  "Extended Data Fig. 10 was generated from completed orthogonal-validation output tables without modifying manuscript files or prior figure panels. For GSE240078 GeoMx DSP, module-level TME AOI versus carcinoma AOI deltas were read from the existing normalized module-test table. For GSE199102 GeoMx WTA, patient-level median TME-like deltas were computed by comparing CAF plus Immune segments with epithelial segments; paired ROI support used CAF-minus-epithelial paired segment results. Module scores had been computed by z-scoring available genes across segments and averaging genes within each module. Cross-GeoMx concordance used the sign of the GSE240078 TME-minus-carcinoma delta and the GSE199102 TME-like-minus-epithelial delta. For GSE310352 CosMx, slide-level CAF/matrix to TGF/EMT adjacency log2 observed/expected values and robustness summaries were read from the finalized robustness tables. Source data for all plotted quantities were exported with panel, dataset, module, metric and note fields."
)
writeLines(methods_lines, file.path(doc_dir, "ed10_v1_methods_draft.md"), useBytes = TRUE)

claim_lines <- c(
  "# ED10 v1 Claim Boundary Notes",
  "",
  "Allowed wording:",
  "- independent GeoMx compartment-level support;",
  "- cross-dataset compartment-level concordance;",
  "- stromal/TME enrichment of CAF/matrix and related immune/interface programs;",
  "- CosMx support for CAF/matrix-associated TGF/EMT stromal-interface organization;",
  "- consistent with CAF/matrix stromal-neighborhood architecture.",
  "",
  "Do not claim:",
  "- causal signaling;",
  "- direct SPP1-CD44 validation;",
  "- tumor-intrinsic EMT;",
  "- reconstruction of the Visium CAF-core distance gradient;",
  "- direct CAF-myeloid adjacency validation from GeoMx;",
  "- lymph-node immune-uncoupling validation;",
  "- clinical prediction.",
  "",
  "SPP1/TAM handling: report as a boundary/control observation only. In these GeoMx summaries it is tumor-high or non-stromal and should not be framed as stromal SPP1-CD44 support."
)
writeLines(claim_lines, file.path(doc_dir, "ed10_v1_claim_boundary_notes.md"), useBytes = TRUE)

panel_map <- data.table(
  panel = LETTERS[1:8],
  title = c(
    "Orthogonal validation schematic",
    "GSE240078 GeoMx DSP TME vs carcinoma enrichment",
    "GSE199102 GeoMx WTA TME-like vs epithelial enrichment",
    "Cross-GeoMx direction concordance",
    "GSE199102 paired ROI CAF-vs-epithelial validation",
    "GSE310352 CosMx CAF/matrix to TGF/EMT interface adjacency",
    "GSE310352 CosMx robustness summary",
    "Inference scope"
  ),
  source = c(
    "conceptual schematic",
    g240_path,
    g199_controls_path,
    concord_path,
    g199_paired_path,
    cosmx_nonoverlap_path,
    paste(c(cosmx_threshold_path, cosmx_norm_path, cosmx_null_path, cosmx_loso_path), collapse = "; "),
    "claim boundary notes"
  ),
  primary_interpretation = c(
    "validation layers only",
    "independent compartment-level TME support",
    "independent compartment-level TME-like support",
    "directional concordance across GeoMx datasets",
    "within-ROI CAF compartment support",
    "cell-level stromal-interface organization",
    "robustness of CosMx interface signal",
    "explicit claim limits"
  )
)
write_csv(panel_map, file.path(manifest_dir, "ed10_v1_panel_map.csv"))

source_manifest <- data.table(
  output = c("figure_pdf", "figure_svg", "source_data", "legend", "results_paragraph", "methods", "claim_boundaries"),
  path = c(pdf_path, svg_path, file.path(source_dir, "Source_Data_Extended_Data_Figure_10_v1.csv"), file.path(doc_dir, "ed10_v1_figure_legend_draft.md"), file.path(doc_dir, "ed10_v1_results_paragraph_draft.md"), file.path(doc_dir, "ed10_v1_methods_draft.md"), file.path(doc_dir, "ed10_v1_claim_boundary_notes.md")),
  status = "generated"
)
write_csv(source_manifest, file.path(manifest_dir, "ed10_v1_source_data_manifest.csv"))

script_manifest <- data.table(
  script = script_path,
  purpose = "Generate Extended Data Figure 10 v1 and associated source/provenance files.",
  inputs = paste(required, collapse = "; "),
  outputs = paste(c(pdf_path, svg_path, file.path(source_dir, "Source_Data_Extended_Data_Figure_10_v1.csv")), collapse = "; "),
  notes = "Does not modify manuscript, existing figures, or reproducibility lock."
)
write_csv(script_manifest, file.path(manifest_dir, "ed10_v1_script_manifest.csv"))

parameter_manifest <- data.table(
  parameter = c(
    "figure_width_in", "figure_height_in", "geomx_main_modules", "gse199102_tme_like_definition",
    "gse310352_cosmx_modules", "cosmx_robustness_groups", "spp1_tam_policy", "claim_policy"
  ),
  value = c(
    "13.5", "9.2", paste(module_order, collapse = ";"),
    "CAF plus Immune segments minus Epithelial segments at patient level",
    "tgfb_emt; tgfb_emt_no_shared; tgfb_emt_no_shared_no_integrin",
    "threshold; normalization; spatial null; leave-one-slide-out",
    "shown only as boundary/tumor-high observation; not stromal support",
    "orthogonal compartment/interface support only; no causality or direct SPP1-CD44"
  )
)
write_csv(parameter_manifest, file.path(manifest_dir, "ed10_v1_parameter_manifest.csv"))

cat("Extended Data Fig. 10 v1 generated:\n", pdf_path, "\n", svg_path, "\n", sep = "")
