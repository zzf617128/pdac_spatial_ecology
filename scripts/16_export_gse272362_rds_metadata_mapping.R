#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "E:/PDAC_TLS/pdac_spatial_ecology"
rds_path <- file.path(project_root, "data", "external", "GSE272362_zenodo", "PDAC_Updated.rds")
out_path <- file.path(project_root, "metadata", "gse272362_rds_sample_mapping.csv")

object <- readRDS(rds_path)
meta <- attributes(object)[["meta.data"]]
cols <- intersect(c("orig.ident", "Sample_ID2", "SlideName", "AreaCode", "Origin", "patient"), colnames(meta))
mapping <- unique(meta[, cols, drop = FALSE])
mapping <- mapping[order(mapping[["orig.ident"]]), , drop = FALSE]
write.csv(mapping, out_path, row.names = FALSE, fileEncoding = "UTF-8")
print(mapping)
cat("Wrote", out_path, "\n")

