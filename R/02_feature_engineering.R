# ==============================================================================
# Script 02: Feature Engineering (Real TCGA Data)
# Purpose: Calculate TMB from mutation data and extract immune gene signatures
# Author: Multi-Omics Pipeline
# ==============================================================================

packages <- c("maftools", "TCGAbiolinks", "dplyr", "tidyr", "SummarizedExperiment")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
library(maftools)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)

# ------------------------------------------------------------------------------
# 1. Initialization
# ------------------------------------------------------------------------------
data_dir <- "data"
project_id <- "TCGA-SKCM"

cat("[INFO] Starting Feature Engineering...\n")

# ------------------------------------------------------------------------------
# 2. Calculate Tumor Mutational Burden (TMB)
# ------------------------------------------------------------------------------
cat("[INFO] Querying Mutation Data (MAF)...\n")
# Download MAF (Mutation Annotation Format)
query_maf <- GDCquery(
    project = project_id, 
    data.category = "Simple Nucleotide Variation", 
    access = "open", 
    data.type = "Masked Somatic Mutation", 
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
tryCatch({
  GDCdownload(query_maf, directory = file.path(data_dir, "GDCdata"), method = "api")
}, error = function(e) {
  GDCdownload(query_maf, directory = file.path(data_dir, "GDCdata"), method = "client")
})
maf_data <- GDCprepare(query_maf, directory = file.path(data_dir, "GDCdata"))

# Create a maftools object
skcm_maf <- read.maf(maf = maf_data, verbose = FALSE)

# Calculate TMB (Assuming ~38MB exome size for TCGA WES)
cat("[INFO] Calculating TMB...\n")
tmb_results <- tmb(maf = skcm_maf, captureSize = 38, logScale = FALSE)

# Format the output dataframe
tmb_df <- tmb_results %>% 
  select(Tumor_Sample_Barcode, total_perMB) %>%
  rename(Sample = Tumor_Sample_Barcode, TMB = total_perMB) %>%
  # TCGA barcodes are long, extract the patient portion (first 12 chars: TCGA-XX-XXXX)
  mutate(Patient = substr(Sample, 1, 12))

# ------------------------------------------------------------------------------
# 3. Extract RNA-seq Signatures (e.g. IFNg, CD8 T effector)
# ------------------------------------------------------------------------------
cat("[INFO] Extracting Gene Expression Signatures...\n")

# Load the SummarizedExperiment object saved in Script 01
se_file <- file.path(data_dir, "tcga_skcm_rnaseq_se.rds")
if (file.exists(se_file)) {
    se_rna <- readRDS(se_file)
    
    # Extract counts matrix
    counts <- assay(se_rna)
    
    # Simple log2(CPM+1) normalization
    cpm <- apply(counts, 2, function(x) (x/sum(x))*1000000)
    log_cpm <- log2(cpm + 1)
    
    # Define a simplified Interferon-Gamma (IFNg) signature
    # (Adapted from Ayers et al., JCI, 2017)
    ifng_genes <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")
    
    # The rownames of TCGA SE objects are usually Ensembl IDs. 
    gene_info <- rowData(se_rna)
    
    # Match symbols
    matched_genes <- rownames(se_rna)[gene_info$gene_name %in% ifng_genes]
    
    if (length(matched_genes) > 0) {
        # Calculate signature score (mean of constituent genes)
        ifng_score <- colMeans(log_cpm[matched_genes, , drop=FALSE])
        
        sig_df <- data.frame(
            Sample = names(ifng_score),
            IFNg_Score = ifng_score
        ) %>%
        mutate(Patient = substr(Sample, 1, 12))
        
        # Merge TMB and Signatures
        feature_matrix <- inner_join(tmb_df, sig_df, by="Patient")
        
        dir.create(file.path(data_dir, "processed"), showWarnings = FALSE)
        write_csv(feature_matrix, file.path(data_dir, "processed", "feature_matrix.csv"))
        cat(sprintf("[SUCCESS] Saved Feature Matrix with %d patients.\n", nrow(feature_matrix)))
    } else {
        warning("Could not map signature genes to rownames.")
    }
} else {
    warning("RNA-seq data not found. Please run 01_fetch_data.R first.")
}

cat("[INFO] Feature Engineering Complete.\n")
