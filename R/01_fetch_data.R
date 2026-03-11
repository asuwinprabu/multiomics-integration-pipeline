# ==============================================================================
# Script 01: TCGA Data Acquisition (Real Clinical Data)
# Purpose: Downloads Clinical and RNA-seq data for the specified cancer cohort
# Author: Multi-Omics Pipeline
# ==============================================================================

# Specify the Bioconductor packages we need
packages <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "readr")

# Ensure required packages are installed
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# 1. Initialization and Configuration
# ------------------------------------------------------------------------------
# Set the project and cohort
project_id <- "TCGA-SKCM" # Skin Cutaneous Melanoma
data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE)

cat(sprintf("\n[INFO] Starting TCGA data download pipeline for %s...\n", project_id))

# ------------------------------------------------------------------------------
# 2. Downloading Clinical Data (JSON API)
# ------------------------------------------------------------------------------
cat("[INFO] Querying Clinical Data via JSON endpoint...\n")
# Using the GDCquery_clinic function bypasses the buggy file download process
patient_data <- GDCquery_clinic(project = project_id, type = "clinical")

if (!is.null(patient_data) && nrow(patient_data) > 0) {
    # Match the expected barcode format for downstream scripts (TCGA-XX-XXXX)
    if (!"bcr_patient_barcode" %in% names(patient_data) && "submitter_id" %in% names(patient_data)) {
        patient_data <- patient_data %>% rename(bcr_patient_barcode = submitter_id)
    }
    write_csv(patient_data, file.path(data_dir, "tcga_skcm_clinical.csv"))
    cat(sprintf("[SUCCESS] Saved clinical metadata with %d patients.\n", nrow(patient_data)))
} else {
    stop("Failed to retrieve clinical data from the GDC API. The server may be down.")
}

# ------------------------------------------------------------------------------
# 3. Downloading RNA-seq Expression Data (Counts)
# ------------------------------------------------------------------------------
cat("[INFO] Querying RNA-seq Gene Expression quantification (STAR - Counts)...\n")
# Note: For testing/demonstration, we restrict to a few samples to save time and bandwidth.
query_rnaseq <- GDCquery(
  project = project_id,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Metastatic")
)

# [TESTING MODE] - only taking first 10 samples to prevent massive download right now
cat("[WARNING] Running in testing mode: Subsetting to first 10 barcodes to save download time.\n")
query_rnaseq$results[[1]] <- head(query_rnaseq$results[[1]], 10)

# Force the API download method to prevent timeout failures on standard client downloads
tryCatch({
  GDCdownload(query_rnaseq, directory = file.path(data_dir, "GDCdata"), method = "api")
}, error = function(e) {
  cat("[WARN] API download failed, trying client method...\n")
  GDCdownload(query_rnaseq, directory = file.path(data_dir, "GDCdata"), method = "client")
})


cat("[INFO] Preparing SummarizedExperiment object...\n")
se_rna <- GDCprepare(query_rnaseq, directory = file.path(data_dir, "GDCdata"))

# Save the raw R object (SummarizedExperiment contains counts and sample metadata)
saveRDS(se_rna, file.path(data_dir, "tcga_skcm_rnaseq_se.rds"))
cat(sprintf("[SUCCESS] Saved RNA-seq SummarizedExperiment with %d genes and %d samples.\n", 
            nrow(se_rna), ncol(se_rna)))

cat("\n[INFO] Data acquisition complete. Ready for Feature Engineering.\n")
