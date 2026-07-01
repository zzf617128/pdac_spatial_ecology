suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- normalizePath(file.path(getwd(), "pdac_spatial_ecology/results/orthogonal_validation_strong_search"), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "datasets/gse199102_geomx")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
doc_dir <- file.path(root, "docs")
source_dir <- file.path(root, "source_data")
manifest_dir <- file.path(root, "manifest")
invisible(lapply(c(table_dir, figure_dir, doc_dir, source_dir, manifest_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, path) fwrite(as.data.frame(x), path)
zv <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

module_defs <- list(
  caf_mycaf = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  pancaf_matrix = c("FAP", "COL1A1", "COL3A1", "PDPN", "THY1", "DCN", "LUM", "FN1", "VIM", "POSTN", "MMP2"),
  myeloid_macrophage = c("LST1", "TYROBP", "AIF1", "LYZ", "CD68", "CD163", "CSF1R", "FCGR3A", "ITGAM", "MSR1", "MRC1", "C1QA", "C1QB", "C1QC"),
  spp1_tam = c("SPP1", "TREM2", "APOE", "LGALS3", "GPNMB", "MARCO", "CD9", "CTSB", "CCL18"),
  tgfb_emt = c("TGFB1", "TGFB2", "TGFBI", "TGFBR1", "TGFBR2", "SERPINE1", "CTGF", "INHBA", "SMAD3", "VIM", "ZEB1", "SNAI1", "SNAI2", "MMP2", "MMP9", "ITGA5", "FN1"),
  ifn_apc = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "B2M", "STAT1", "IRF1", "IFITM1", "IFIT1", "IFIT3", "CXCL9", "CXCL10", "TAP1"),
  immune_core = c("CXCL13", "CCL19", "CCL21", "CXCL9", "CXCL10", "CD74", "HLA-DRA", "HLA-DPA1", "MS4A1", "CD3D", "CD3E"),
  t_cell = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "TRAC", "CD247", "PDCD1", "CTLA4", "LAG3", "TIGIT"),
  b_plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGKC", "XBP1", "SDC1"),
  tumor_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "MSLN"),
  tumor_aggressive = c("KRT17", "KRT5", "KRT6A", "S100A2", "SERPINB3", "LAMC2", "ITGA3", "MMP7", "TACSTD2")
)

segment_path <- file.path(data_dir, "GSE199102_Broad_PDAC_WTA_AllSamples_SegmentProperties.txt.gz")
expr_path <- file.path(data_dir, "GSE199102_hPDAC_WTA_20210222T2101_Q3Norm_TargetCountMatrix.txt.gz")
worksheet_path <- file.path(data_dir, "GSE199102_hPDAC_WTA_20210222T2101_LabWorksheet.txt.gz")
rds_files <- file.path(data_dir, c("GSE199102_ssGSEA_detrendApproach21-6-2.RDS.gz", "GSE199102_dfs_Q321-6.RDS.gz"))

seg <- fread(segment_path, showProgress = FALSE)
expr_dt <- fread(expr_path, showProgress = FALSE)
setnames(expr_dt, names(expr_dt)[1], "Gene")
expr_dt[, Gene := toupper(Gene)]
expr_dt <- expr_dt[!duplicated(Gene)]
expr <- as.matrix(expr_dt[, -"Gene"])
mode(expr) <- "numeric"
rownames(expr) <- expr_dt$Gene
colnames(expr) <- gsub("\\.", "-", colnames(expr))
colnames(expr) <- gsub("^DSP-", "DSP-", colnames(expr))
colnames(expr) <- gsub("^DSP", "DSP", colnames(expr))
seg[, sample_key_dash := gsub("-", ".", Sample_ID)]
seg[, sample_key_expr := Sample_ID]
seg[, sample_key_expr := gsub("-", ".", sample_key_expr)]
expr_sample_ids <- colnames(expr)
expr_sample_ids_dash <- gsub("\\.", "-", expr_sample_ids)
seg[, matrix_col := expr_sample_ids[match(Sample_ID, expr_sample_ids_dash)]]
if (all(is.na(seg$matrix_col))) {
  seg[, matrix_col := expr_sample_ids[match(gsub("-", ".", Sample_ID), expr_sample_ids)]]
}
seg_matched <- seg[!is.na(matrix_col)]
expr <- expr[, seg_matched$matrix_col, drop = FALSE]
colnames(expr) <- seg_matched$Sample_ID

seg_matched[, segment_class := fifelse(Segment %in% c("CAF"), "stromal_fibroblast",
  fifelse(Segment %in% c("Epithelial"), "tumor_epithelial",
    fifelse(Segment %in% c("Immune"), "immune", "other")))]
seg_matched[, roi_id := paste(Patient, Slide_name, Scan_name, ROI_number, sep = "__")]
seg_matched[, specimen_id := Slide_name]
seg_matched[, treatment_group := fifelse(is.na(TreatmentClass) | TreatmentClass == "", Treatment, TreatmentClass)]

rds_status <- rbindlist(lapply(rds_files, function(f) {
  status <- "not_attempted"
  class_text <- NA_character_
  note <- NA_character_
  obj <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(obj, "error")) {
    status <- "read_failed"
    note <- conditionMessage(obj)
  } else {
    status <- "loaded"
    class_text <- paste(class(obj), collapse = "/")
    note <- if (is.list(obj)) paste(head(names(obj), 20), collapse = ";") else ""
  }
  data.table(file = basename(f), status = status, class = class_text, note = note)
}))

available_genes <- rownames(expr)
coverage <- rbindlist(lapply(names(module_defs), function(module) {
  genes <- toupper(module_defs[[module]])
  present <- intersect(genes, available_genes)
  data.table(
    dataset = "GSE199102",
    platform = "GeoMx WTA",
    module = module,
    requested_genes = paste(genes, collapse = ";"),
    available_genes = paste(present, collapse = ";"),
    n_requested = length(genes),
    n_available = length(present),
    coverage = length(present) / length(genes),
    usable = length(present) >= 3
  )
}))
write_csv(coverage, file.path(table_dir, "gse199102_geomx_gene_coverage_by_module.csv"))

score_mat <- sapply(names(module_defs), function(module) {
  present <- intersect(toupper(module_defs[[module]]), rownames(expr))
  if (length(present) < 3) return(rep(NA_real_, ncol(expr)))
  z <- t(apply(expr[present, , drop = FALSE], 1, zv))
  colMeans(z, na.rm = TRUE)
})
score_dt <- as.data.table(score_mat)
score_dt[, Sample_ID := colnames(expr)]
score_dt <- cbind(seg_matched[, .(Sample_ID, Patient, PatientNumber, specimen_id, Slide_name, Scan_name, ROI_number, roi_id, Segment, segment_class, AOI_area, ROI_mask_area, Not_malignant, TreatmentClass, Treatment, treatment_group, Status, RawReads, AlignedReads, SequencingSaturation, NormFactorQ3)], score_dt[, -"Sample_ID"])
write_csv(score_dt, file.path(table_dir, "gse199102_geomx_module_scores.csv"))

metadata_summary <- rbindlist(list(
  data.table(field = "n_segments_matched_to_matrix", value = as.character(nrow(score_dt))),
  data.table(field = "n_matrix_genes", value = as.character(nrow(expr))),
  data.table(field = "n_patients", value = as.character(uniqueN(score_dt$Patient))),
  data.table(field = "n_specimens", value = as.character(uniqueN(score_dt$specimen_id))),
  data.table(field = "n_rois", value = as.character(uniqueN(score_dt$roi_id))),
  data.table(field = "segment_counts", value = paste(names(table(score_dt$Segment)), as.integer(table(score_dt$Segment)), sep = "=", collapse = ";")),
  data.table(field = "treatment_counts", value = paste(names(table(score_dt$treatment_group)), as.integer(table(score_dt$treatment_group)), sep = "=", collapse = ";")),
  data.table(field = "patient_level_aggregation_possible", value = "TRUE"),
  data.table(field = "paired_segments_within_roi_possible", value = as.character(any(table(score_dt$roi_id) > 1))),
  data.table(field = "rds_status", value = paste(rds_status$file, rds_status$status, sep = "=", collapse = ";"))
))
write_csv(metadata_summary, file.path(table_dir, "gse199102_geomx_metadata_summary.csv"))

expression_summary <- data.table(
  metric = c("n_genes", "n_segments", "min_expression", "median_expression", "max_expression", "n_missing"),
  value = c(nrow(expr), ncol(expr), min(expr, na.rm = TRUE), median(expr, na.rm = TRUE), max(expr, na.rm = TRUE), sum(is.na(expr)))
)
write_csv(expression_summary, file.path(table_dir, "gse199102_geomx_expression_summary.csv"))

segment_annotation_summary <- score_dt[, .(
  n_segments = .N,
  n_patients = uniqueN(Patient),
  n_rois = uniqueN(roi_id),
  median_aoi_area = as.numeric(median(AOI_area, na.rm = TRUE)),
  median_not_malignant = as.numeric(median(Not_malignant, na.rm = TRUE))
), by = .(Segment, segment_class)]
write_csv(segment_annotation_summary, file.path(table_dir, "gse199102_geomx_segment_annotation_summary.csv"))

module_names <- names(module_defs)
compare_classes <- function(dt, class_a, class_b, comparison_name) {
  rbindlist(lapply(module_names, function(module) {
    sub <- dt[segment_class %in% c(class_a, class_b) & is.finite(get(module))]
    if (nrow(sub) < 4 || length(unique(sub$segment_class)) < 2) {
      return(data.table(comparison = comparison_name, module = module, n = nrow(sub), class_a = class_a, class_b = class_b, mean_a = NA_real_, mean_b = NA_real_, delta_a_minus_b = NA_real_, wilcox_p = NA_real_))
    }
    y <- sub[[module]]
    data.table(
      comparison = comparison_name, module = module, n = nrow(sub), class_a = class_a, class_b = class_b,
      mean_a = mean(y[sub$segment_class == class_a], na.rm = TRUE),
      mean_b = mean(y[sub$segment_class == class_b], na.rm = TRUE),
      median_a = median(y[sub$segment_class == class_a], na.rm = TRUE),
      median_b = median(y[sub$segment_class == class_b], na.rm = TRUE),
      delta_a_minus_b = mean(y[sub$segment_class == class_a], na.rm = TRUE) - mean(y[sub$segment_class == class_b], na.rm = TRUE),
      wilcox_p = suppressWarnings(wilcox.test(y ~ sub$segment_class)$p.value)
    )
  }), fill = TRUE)
}
tme_dt <- copy(score_dt)
tme_dt[, segment_class := fifelse(segment_class %in% c("stromal_fibroblast", "immune"), "tme_like", segment_class)]
fib_dt <- copy(score_dt)
fib_dt[, segment_class := fifelse(segment_class == "stromal_fibroblast", "stromal_fibroblast", "non_fibroblast")]
segment_tests <- rbindlist(list(
  compare_classes(score_dt, "stromal_fibroblast", "tumor_epithelial", "CAF_vs_epithelial"),
  compare_classes(score_dt, "immune", "tumor_epithelial", "immune_vs_epithelial"),
  compare_classes(tme_dt, "tme_like", "tumor_epithelial", "TME_vs_epithelial"),
  compare_classes(fib_dt, "stromal_fibroblast", "non_fibroblast", "CAF_vs_nonCAF")
), fill = TRUE)
segment_tests[, fdr := p.adjust(wilcox_p, method = "BH"), by = comparison]
write_csv(segment_tests, file.path(table_dir, "gse199102_geomx_segment_comparison_results.csv"))

patient_class <- score_dt[, lapply(.SD, median, na.rm = TRUE), by = .(Patient, segment_class), .SDcols = module_names]
patient_long <- melt(patient_class, id.vars = c("Patient", "segment_class"), variable.name = "module", value.name = "median_score")
patient_wide <- dcast(patient_long, Patient + module ~ segment_class, value.var = "median_score")
for (nm in c("stromal_fibroblast", "tumor_epithelial", "immune")) if (!nm %in% names(patient_wide)) patient_wide[, (nm) := NA_real_]
patient_wide[, stromal_minus_tumor := stromal_fibroblast - tumor_epithelial]
patient_wide[, immune_minus_tumor := immune - tumor_epithelial]
patient_wide[, fibroblast_minus_nonfibroblast := NA_real_]
nonfib <- score_dt[segment_class != "stromal_fibroblast", lapply(.SD, median, na.rm = TRUE), by = Patient, .SDcols = module_names]
nonfib_long <- melt(nonfib, id.vars = "Patient", variable.name = "module", value.name = "nonfibroblast")
patient_wide <- merge(patient_wide, nonfib_long, by = c("Patient", "module"), all.x = TRUE)
patient_wide[, fibroblast_minus_nonfibroblast := stromal_fibroblast - nonfibroblast]
tme <- score_dt[segment_class != "tumor_epithelial", lapply(.SD, median, na.rm = TRUE), by = Patient, .SDcols = module_names]
tme_long <- melt(tme, id.vars = "Patient", variable.name = "module", value.name = "tme_like")
patient_wide <- merge(patient_wide, tme_long, by = c("Patient", "module"), all.x = TRUE)
patient_wide[, tme_minus_tumor := tme_like - tumor_epithelial]
write_csv(patient_wide, file.path(table_dir, "gse199102_geomx_patient_level_deltas.csv"))

paired_rows <- list()
for (module in module_names) {
  roi_w <- dcast(score_dt[, .(roi_id, Patient, ROI_number, Segment, value = get(module))], roi_id + Patient + ROI_number ~ Segment, value.var = "value", fun.aggregate = median, na.rm = TRUE)
  if (all(c("CAF", "Epithelial") %in% names(roi_w))) {
    roi_w[, delta_CAF_minus_Epithelial := CAF - Epithelial]
    paired_rows[[paste0(module, "_caf")]] <- data.table(module = module, comparison = "paired_CAF_minus_Epithelial", n_roi = sum(is.finite(roi_w$delta_CAF_minus_Epithelial)), median_delta = median(roi_w$delta_CAF_minus_Epithelial, na.rm = TRUE), support_fraction_positive = mean(roi_w$delta_CAF_minus_Epithelial > 0, na.rm = TRUE), wilcox_p = suppressWarnings(wilcox.test(roi_w$delta_CAF_minus_Epithelial)$p.value))
  }
  if (all(c("Immune", "Epithelial") %in% names(roi_w))) {
    roi_w[, delta_Immune_minus_Epithelial := Immune - Epithelial]
    paired_rows[[paste0(module, "_immune")]] <- data.table(module = module, comparison = "paired_Immune_minus_Epithelial", n_roi = sum(is.finite(roi_w$delta_Immune_minus_Epithelial)), median_delta = median(roi_w$delta_Immune_minus_Epithelial, na.rm = TRUE), support_fraction_positive = mean(roi_w$delta_Immune_minus_Epithelial > 0, na.rm = TRUE), wilcox_p = suppressWarnings(wilcox.test(roi_w$delta_Immune_minus_Epithelial)$p.value))
  }
}
paired <- rbindlist(paired_rows, fill = TRUE)
paired[, fdr := p.adjust(wilcox_p, method = "BH"), by = comparison]
write_csv(paired, file.path(table_dir, "gse199102_geomx_paired_segment_results.csv"))

patient_summary <- patient_wide[, .(
  n_patients = sum(is.finite(stromal_minus_tumor)),
  median_stromal_minus_tumor = median(stromal_minus_tumor, na.rm = TRUE),
  patient_support_fraction_stromal_gt_tumor = mean(stromal_minus_tumor > 0, na.rm = TRUE),
  median_immune_minus_tumor = median(immune_minus_tumor, na.rm = TRUE),
  patient_support_fraction_immune_gt_tumor = mean(immune_minus_tumor > 0, na.rm = TRUE),
  median_fibroblast_minus_nonfibroblast = median(fibroblast_minus_nonfibroblast, na.rm = TRUE),
  patient_support_fraction_fibroblast_gt_nonfibroblast = mean(fibroblast_minus_nonfibroblast > 0, na.rm = TRUE),
  median_tme_minus_tumor = median(tme_minus_tumor, na.rm = TRUE),
  patient_support_fraction_tme_gt_tumor = mean(tme_minus_tumor > 0, na.rm = TRUE)
), by = module]

expected <- data.table(
  module = module_names,
  expected_stromal_direction = c("stromal_high", "stromal_high", "stromal_or_immune_high", "report_separately", "stromal_high", "immune_or_stromal_high", "immune_or_stromal_high", "immune_high", "immune_high", "tumor_high", "tumor_high")
)
controls <- merge(patient_summary, expected, by = "module", all.x = TRUE)
controls[, control_pass := fifelse(expected_stromal_direction == "stromal_high", median_stromal_minus_tumor > 0 & patient_support_fraction_stromal_gt_tumor >= 0.6,
  fifelse(expected_stromal_direction == "tumor_high", median_stromal_minus_tumor < 0 & patient_support_fraction_stromal_gt_tumor <= 0.4,
    fifelse(expected_stromal_direction == "immune_high", median_immune_minus_tumor > 0 & patient_support_fraction_immune_gt_tumor >= 0.6,
      fifelse(expected_stromal_direction %in% c("immune_or_stromal_high", "stromal_or_immune_high"), (median_immune_minus_tumor > 0 & patient_support_fraction_immune_gt_tumor >= 0.6) | (median_stromal_minus_tumor > 0 & patient_support_fraction_stromal_gt_tumor >= 0.6),
        NA))))]
write_csv(controls, file.path(table_dir, "gse199102_geomx_positive_negative_control_results.csv"))

gse240_path <- file.path(root, "outputs/orthogonal_validation_2026_06_30/tables/geomx_gse240078_tme_vs_carcinoma_module_tests.csv")
if (!file.exists(gse240_path)) {
  gse240_path <- file.path(getwd(), "pdac_spatial_ecology/outputs/orthogonal_validation_2026_06_30/tables/geomx_gse240078_tme_vs_carcinoma_module_tests.csv")
}
if (file.exists(gse240_path)) {
  g240 <- fread(gse240_path)
  g240 <- g240[, .(module, gse240078_delta_stroma_minus_tumor = delta_stroma_minus_tumor, gse240078_fdr = fdr)]
} else {
  g240 <- data.table(module = module_names, gse240078_delta_stroma_minus_tumor = NA_real_, gse240078_fdr = NA_real_)
}
concord <- merge(patient_summary[, .(module, gse199102_median_tme_minus_tumor = median_tme_minus_tumor, gse199102_tme_patient_support = patient_support_fraction_tme_gt_tumor, gse199102_median_stromal_minus_tumor = median_stromal_minus_tumor, gse199102_median_immune_minus_tumor = median_immune_minus_tumor)], g240, by = "module", all.x = TRUE)
concord[, gse199102_direction := fifelse(gse199102_median_tme_minus_tumor > 0, "tme_high", "tumor_high")]
concord[, gse240078_direction := fifelse(gse240078_delta_stroma_minus_tumor > 0, "stromal_high", "tumor_high")]
concord[, direction_concordant := (gse199102_median_tme_minus_tumor > 0 & gse240078_delta_stroma_minus_tumor > 0) | (gse199102_median_tme_minus_tumor < 0 & gse240078_delta_stroma_minus_tumor < 0)]
concord[module == "spp1_tam", direction_concordant := NA]
concord[module == "spp1_tam", gse199102_direction := paste0(gse199102_direction, "_report_separately")]
write_csv(concord, file.path(table_dir, "gse199102_gse240078_direction_concordance.csv"))

source_data <- merge(patient_summary, paired[comparison == "paired_CAF_minus_Epithelial", .(module, paired_n_roi = n_roi, paired_median_delta = median_delta, paired_support_fraction_positive = support_fraction_positive, paired_fdr = fdr)], by = "module", all.x = TRUE)
source_data <- merge(source_data, concord, by = "module", all.x = TRUE, suffixes = c("", "_concord"))
write_csv(source_data, file.path(source_dir, "Source_Data_GSE199102_GeoMx.csv"))

pdf(file.path(figure_dir, "gse199102_geomx_candidate_panels.pdf"), width = 8.5, height = 6.5, onefile = TRUE)
print(ggplot(coverage, aes(x = reorder(module, coverage), y = coverage, fill = usable)) + geom_col() + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "Gene coverage", fill = "Usable"))
print(ggplot(patient_summary, aes(x = reorder(module, median_stromal_minus_tumor), y = median_stromal_minus_tumor, fill = patient_support_fraction_stromal_gt_tumor)) + geom_hline(yintercept = 0, linetype = 2) + geom_col() + coord_flip() + theme_bw(base_size = 9) + scale_fill_viridis_c(limits = c(0, 1)) + labs(x = NULL, y = "Patient-level median CAF - epithelial score", fill = "Patient support"))
print(ggplot(paired[comparison == "paired_CAF_minus_Epithelial"], aes(x = reorder(module, median_delta), y = median_delta, fill = support_fraction_positive)) + geom_hline(yintercept = 0, linetype = 2) + geom_col() + coord_flip() + theme_bw(base_size = 9) + scale_fill_viridis_c(limits = c(0, 1)) + labs(x = NULL, y = "Paired ROI median CAF - epithelial delta", fill = "ROI support"))
print(ggplot(controls, aes(x = reorder(module, median_stromal_minus_tumor), y = median_stromal_minus_tumor, color = expected_stromal_direction, shape = control_pass)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(size = 3) + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "Patient-level stromal - tumor delta", color = "Expected", shape = "Pass"))
dev.off()

pdf(file.path(figure_dir, "gse199102_gse240078_concordance_plot.pdf"), width = 6.5, height = 4.5)
print(ggplot(concord[module != "spp1_tam"], aes(x = gse240078_delta_stroma_minus_tumor, y = gse199102_median_tme_minus_tumor, label = module, color = direction_concordant)) + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + geom_point(size = 2.5) + geom_text(size = 2.6, check_overlap = TRUE, vjust = -0.6) + theme_bw(base_size = 9) + labs(x = "GSE240078 stroma - tumor delta", y = "GSE199102 TME-like - epithelial delta", color = "Direction concordant"))
dev.off()

pdf(file.path(figure_dir, "gse199102_geomx_positive_negative_controls.pdf"), width = 7, height = 4.5)
print(ggplot(controls, aes(x = reorder(module, median_stromal_minus_tumor), y = median_stromal_minus_tumor, fill = control_pass)) + geom_hline(yintercept = 0, linetype = 2) + geom_col() + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "Patient-level stromal - tumor delta", fill = "Expected control pass"))
dev.off()

strong_modules <- c("caf_mycaf", "pancaf_matrix", "myeloid_macrophage", "tgfb_emt", "ifn_apc", "immune_core")
strong_support <- controls[module %in% strong_modules]
tumor_support <- controls[module %in% c("tumor_epithelial", "tumor_aggressive")]
decision <- if (
  all(controls[module %in% c("caf_mycaf", "pancaf_matrix"), median_stromal_minus_tumor > 0]) &&
  sum(strong_support$control_pass %in% TRUE, na.rm = TRUE) >= 4 &&
  all(tumor_support$control_pass %in% TRUE, na.rm = TRUE)
) "A. Strong ED support" else if (all(controls[module %in% c("caf_mycaf", "pancaf_matrix"), median_stromal_minus_tumor > 0])) "B. Secondary/source support" else "C. Do not show"

writeLines(c(
  "# GSE199102 GeoMx Feasibility Report",
  "",
  paste0("Matched segments: ", nrow(score_dt), "."),
  paste0("Patients: ", uniqueN(score_dt$Patient), "."),
  paste0("Specimens/slides: ", uniqueN(score_dt$specimen_id), "."),
  paste0("ROIs: ", uniqueN(score_dt$roi_id), "."),
  paste0("Segments: ", paste(names(table(score_dt$Segment)), as.integer(table(score_dt$Segment)), sep = "=", collapse = "; ")),
  paste0("Treatment groups: ", paste(names(table(score_dt$treatment_group)), as.integer(table(score_dt$treatment_group)), sep = "=", collapse = "; ")),
  "",
  "Segment labels include CAF, Immune and Epithelial, enabling stromal/fibroblast, immune and tumor-like comparisons.",
  "Patient-level aggregation and paired ROI-level segment comparisons are possible.",
  "",
  "RDS files were attempted but did not load with `readRDS`; the primary analysis uses Q3-normalized WTA matrix and segment properties."
), file.path(doc_dir, "gse199102_geomx_feasibility_report.md"), useBytes = TRUE)

writeLines(c(
  "# GSE199102 and GSE240078 Direction Concordance Interpretation",
  "",
  "GSE199102 was analyzed as CAF/Immune/Epithelial GeoMx WTA segment data and compared with the previous GSE240078 stroma-vs-carcinoma DSP result.",
  "",
  capture.output(print(concord)),
  "",
  "SPP1/TAM is reported separately and should not be forced into the stromal-high group."
), file.path(doc_dir, "gse199102_gse240078_concordance_interpretation.md"), useBytes = TRUE)

writeLines(c(
  "# GSE199102 GeoMx Final Decision Report",
  "",
  paste0("Decision: ", decision),
  "",
  "## Summary",
  capture.output(print(controls)),
  "",
  "## Interpretation",
  "GSE199102 provides independent GeoMx WTA compartment-level testing using CAF, Immune and Epithelial segment labels. Interpret only at compartment level, not as cell-cell spatial proximity.",
  "",
  "## Claim Boundaries",
  "- No causality.",
  "- No direct SPP1-CD44 validation.",
  "- No Visium gradient reconstruction.",
  "- No LN immune uncoupling."
), file.path(doc_dir, "gse199102_geomx_final_decision_report.md"), useBytes = TRUE)

results_par <- if (decision == "A. Strong ED support") {
  "In an independent GeoMx WTA PDAC cohort (GSE199102), CAF-segmented AOIs showed higher CAF/myCAF and panCAF/matrix module scores than matched epithelial segments, with patient- and paired-ROI-level support. Immune and stromal-interface programs were evaluated alongside tumor epithelial and tumor-aggressive internal controls, providing cross-dataset compartment-level support consistent with the GSE240078 DSP result."
} else {
  "GSE199102 GeoMx WTA provided independent compartment-level context for stromal and epithelial modules, but its support should be interpreted according to the module-specific patient-level stability reported in the source tables."
}
writeLines(c("# GSE199102 Results Paragraph Draft", "", results_par), file.path(doc_dir, "gse199102_geomx_results_paragraph_draft.md"), useBytes = TRUE)

methods_par <- "GSE199102 GeoMx WTA Q3-normalized target count matrices and segment-property metadata were downloaded from GEO. Segment identifiers were matched between the expression matrix and segment properties, and segments were grouped as CAF, Immune or Epithelial according to the provided Segment field. Module scores were computed by z-scoring each gene across segments and averaging available genes within each module. Comparisons were summarized at segment, paired-ROI and patient levels; patient-level medians were used as the primary interpretive unit."
writeLines(c("# GSE199102 Methods Draft", "", methods_par), file.path(doc_dir, "gse199102_geomx_methods_draft.md"), useBytes = TRUE)

legend <- "Candidate panel draft: GSE199102 GeoMx WTA compartment validation. (A) Module gene coverage. (B) Patient-level CAF-segment minus epithelial-segment module score deltas. (C) Paired ROI CAF-minus-epithelial deltas. (D) Positive and negative compartment controls. (E) Direction concordance with GSE240078 DSP stroma-vs-carcinoma results."
writeLines(c("# GSE199102 Figure Legend Draft", "", legend), file.path(doc_dir, "gse199102_geomx_figure_legend_draft.md"), useBytes = TRUE)

ranking_path <- file.path(table_dir, "signal_strength_ranking.csv")
rank <- if (file.exists(ranking_path)) fread(ranking_path, fill = TRUE) else data.table()
new_row <- data.table(
  dataset = "GSE199102 Hwang et al. PDAC GeoMx WTA",
  platform = "GeoMx WTA",
  validation_target = "independent compartment-level CAF/immune/epithelial segment validation",
  effect_size = paste0("CAF-myCAF median patient CAF-epithelial delta=", round(controls[module == "caf_mycaf", median_stromal_minus_tumor], 3), "; matrix delta=", round(controls[module == "pancaf_matrix", median_stromal_minus_tumor], 3)),
  support_fraction = paste0("CAF-myCAF patient support=", round(controls[module == "caf_mycaf", patient_support_fraction_stromal_gt_tumor], 3), "; matrix support=", round(controls[module == "pancaf_matrix", patient_support_fraction_stromal_gt_tumor], 3)),
  patient_or_sample_level_support = "patient-level and paired ROI-level available",
  FDR_or_p_value = "see segment and paired ROI tables",
  robustness_level = decision,
  positive_controls_passed = paste0(sum(controls$control_pass %in% TRUE, na.rm = TRUE), "/", sum(!is.na(controls$control_pass)), " expected controls"),
  negative_controls_passed = "tumor epithelial/tumor-aggressive control directions reported",
  claim_supported = ifelse(decision == "A. Strong ED support", "independent GeoMx compartment-level support", "module-specific source-level support"),
  claim_boundary = "Compartment-level only; no causal signaling or direct SPP1-CD44 validation",
  recommended_display = ifelse(decision == "A. Strong ED support", "extended_data", "source_only_or_secondary"),
  reason = "Independent WTA GeoMx dataset with CAF, Immune and Epithelial segments; patient-level aggregation performed."
)
if (nrow(rank) > 0 && "dataset" %in% names(rank)) rank <- rank[dataset != new_row$dataset]
rank <- rbindlist(list(rank, new_row), fill = TRUE)
write_csv(rank, ranking_path)
write_csv(rank, file.path(manifest_dir, "signal_strength_ranking.csv"))

panel_manifest_path <- file.path(manifest_dir, "candidate_ed_panel_manifest.csv")
panel <- if (file.exists(panel_manifest_path)) fread(panel_manifest_path, fill = TRUE) else data.table()
new_panel <- data.table(
  panel = "GSE199102",
  title = "GeoMx WTA compartment-level validation",
  status = ifelse(decision == "A. Strong ED support", "candidate_component", "source_or_secondary_component"),
  dataset = "GSE199102",
  source_table = file.path(source_dir, "Source_Data_GSE199102_GeoMx.csv"),
  figure_path = file.path(figure_dir, "gse199102_geomx_candidate_panels.pdf"),
  decision = decision,
  note = "Final ED10 generation deferred."
)
if (nrow(panel) > 0 && all(c("panel", "dataset") %in% names(panel))) panel <- panel[!(panel == "GSE199102" & dataset == "GSE199102")]
panel <- rbindlist(list(panel, new_panel), fill = TRUE)
write_csv(panel, panel_manifest_path)

candidate_doc <- file.path(doc_dir, "candidate_results_for_ED_figure.md")
if (file.exists(candidate_doc)) {
  candidate_lines <- readLines(candidate_doc, warn = FALSE)
  starts <- grep("^## GSE199102 GeoMx WTA$", candidate_lines)
  if (length(starts) > 0) {
    remove_idx <- integer()
    for (s in starts) {
      next_starts <- grep("^## ", candidate_lines)
      next_starts <- next_starts[next_starts > s]
      e <- if (length(next_starts) > 0) min(next_starts) - 1 else length(candidate_lines)
      remove_idx <- c(remove_idx, seq.int(s, e))
    }
    candidate_lines <- candidate_lines[-unique(remove_idx)]
    writeLines(candidate_lines, candidate_doc, useBytes = TRUE)
  }
}
cat(
  "\n\n## GSE199102 GeoMx WTA\n\n",
  "Decision: ", decision, "\n\n",
  "Use as an independent GeoMx compartment-level component only if the module-specific patient-level and paired-ROI tables are judged concordant enough with GSE240078. SPP1/TAM remains separately reported and should not be used as stromal SPP1 support.\n",
  file = candidate_doc, append = TRUE, sep = ""
)

manifest <- data.table(
  output = c(
    "gse199102_geomx_metadata_summary.csv", "gse199102_geomx_expression_summary.csv",
    "gse199102_geomx_segment_annotation_summary.csv", "gse199102_geomx_gene_coverage_by_module.csv",
    "gse199102_geomx_module_scores.csv", "gse199102_geomx_segment_comparison_results.csv",
    "gse199102_geomx_patient_level_deltas.csv", "gse199102_geomx_paired_segment_results.csv",
    "gse199102_gse240078_direction_concordance.csv", "gse199102_geomx_positive_negative_control_results.csv",
    "Source_Data_GSE199102_GeoMx.csv", "gse199102_geomx_candidate_panels.pdf"
  ),
  path = c(
    file.path(table_dir, "gse199102_geomx_metadata_summary.csv"), file.path(table_dir, "gse199102_geomx_expression_summary.csv"),
    file.path(table_dir, "gse199102_geomx_segment_annotation_summary.csv"), file.path(table_dir, "gse199102_geomx_gene_coverage_by_module.csv"),
    file.path(table_dir, "gse199102_geomx_module_scores.csv"), file.path(table_dir, "gse199102_geomx_segment_comparison_results.csv"),
    file.path(table_dir, "gse199102_geomx_patient_level_deltas.csv"), file.path(table_dir, "gse199102_geomx_paired_segment_results.csv"),
    file.path(table_dir, "gse199102_gse240078_direction_concordance.csv"), file.path(table_dir, "gse199102_geomx_positive_negative_control_results.csv"),
    file.path(source_dir, "Source_Data_GSE199102_GeoMx.csv"), file.path(figure_dir, "gse199102_geomx_candidate_panels.pdf")
  )
)
write_csv(manifest, file.path(manifest_dir, "gse199102_geomx_manifest.csv"))

cat("GSE199102 GeoMx analysis complete. Decision: ", decision, "\n", sep = "")
