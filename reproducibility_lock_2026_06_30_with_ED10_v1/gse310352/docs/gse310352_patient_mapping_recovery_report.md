# GSE310352 Patient/Sample Mapping Recovery Report

## Scope
This pass attempted to recover patient ID, specimen ID, sample ID, treatment group and tissue block mapping for the eight GSE310352 CosMx slides without adding new biology analyses and without modifying manuscript or ED10 files.

## Summary Finding
Patient/specimen/block mapping could not be recovered from the available public GEO records or the processed CosMx files inspected locally. The only high-confidence mapping is slide title to GEO sample accession, plus treatment/tissue descriptors shared across all eight samples.

## Recovered High-Confidence Fields
- GEO sample IDs: GSM9294399-slide1, GSM9294400-slide2, GSM9294401-slide3, GSM9294402-slide4, GSM9294404-slide5, GSM9294405-slide5b2, GSM9294406-slide5b3, GSM9294403-slide6.
- Treatment group: treatment naive for all GSE310352 CosMx samples.
- Tissue/source: PDAC tumor resection, FFPE 5 um section, Nanostring CosMx TAP 1k universal RNA panel.

## Unrecovered Fields
- Patient ID: not recovered for any slide.
- Specimen ID: not recovered for any slide beyond the GEO sample/slide accession itself.
- Tissue block: not recovered for any slide.
- Cross-modal key linking CosMx slide names to Visium/other subseries patient identifiers: not recovered.

## Sources Checked
- Local/official GEO metadata files: GSE310352_series_matrix.txt.gz, GSE310352_family.soft.gz, GSE310352_family.xml.tgz, GSE310353_family.soft.gz, GSE282302_family.soft.gz, GSE310388_family.soft.gz.
- GSE310352 family SOFT and series matrix list the eight slide titles and GSM accessions and state tissue/treatment, but do not include patient, specimen, block, donor, case, or anonymized subject fields.
- GSE310352 MINiML gives the same slide/GSM/tissue/treatment/supplementary-file information and no patient/specimen/block fields.
- GSE310353 SuperSeries SOFT confirms GSE310352 is a CosMx subseries alongside GSE282302 and GSE310388, but does not provide a key connecting CosMx slide names to patient/specimen IDs.
- GSE282302 Spatial Transcriptomics SOFT uses ROI/sample names such as C1_D10_ROI1_s1 and treatment descriptors, but no public key was found linking these names to GSE310352 slide1-slide6/slide5b2/slide5b3.
- GSE310388 RNA-seq SOFT describes cell-line normoxia/hypoxia samples and is not informative for CosMx patient mapping.
- Processed local GSE310352 files inspected: GSM9294399_slide1_cell_by_gene_counts.csv.gz; GSM9294399_slide1_cell_metadata.csv.gz; GSM9294400_slide2_cell_by_gene_counts.csv.gz; GSM9294400_slide2_cell_metadata.csv.gz; GSM9294401_slide3_cell_by_gene_counts.csv.gz; GSM9294401_slide3_cell_metadata.csv.gz; GSM9294402_slide4_cell_by_gene_counts.csv.gz; GSM9294402_slide4_cell_metadata.csv.gz; GSM9294403_slide6_cell_by_gene_counts.csv.gz; GSM9294403_slide6_cell_metadata.csv.gz; GSM9294404_slide5_cell_by_gene_counts.csv.gz; GSM9294404_slide5_cell_metadata.csv.gz; GSM9294405_slide5b2_cell_by_gene_counts.csv.gz; GSM9294405_slide5b2_cell_metadata.csv.gz; GSM9294406_slide5b3_cell_by_gene_counts.csv.gz; GSM9294406_slide5b3_cell_metadata.csv.gz.
- Processed cell metadata headers contain fov, cell_ID, Area, AspectRatio, local/global coordinates, width/height, and IF intensity fields only.
- GEO supplementary FTP directories list per-slide counts, cell metadata, raw_core.tar.gz, and Seurat RDS files, but no standalone sample metadata/clinical mapping table.

## Raw Archive / Seurat Attempts
- Raw core and Seurat downloads were attempted for the smallest slide5b3 files, but the FTP transfers ended with unexpected EOF and the incomplete files were removed. These incomplete downloads were not used as evidence.
- GEO data_processing states raw_core archives contain CosMx tx_file.csv, polygons.csv, metadata_file.csv and fov_positions_file.csv. This description suggests coordinate/segmentation content, not a patient/specimen key; nevertheless, this remains a residual limitation unless full raw archives can be downloaded later.

## Slide5 Family Note
The names slide5, slide5b2 and slide5b3 are compatible with a possible slide/block family, but public GEO metadata does not define what b2/b3 mean. This filename pattern is therefore not sufficient to merge them into one patient or specimen.

## Recommendation For Manuscript Wording
Continue to describe GSE310352 CosMx support as slide-level/FOV-level rather than patient-level or specimen-level. Avoid wording such as patient-level replication, independent patients, per-patient consistency, or specimen-level support for this dataset unless an explicit external sample key is obtained.

## Output
- Mapping attempt table: `E:/PDAC_TLS/pdac_spatial_ecology/results/orthogonal_validation_strong_search/tables/gse310352_slide_patient_mapping_attempt.csv`.
