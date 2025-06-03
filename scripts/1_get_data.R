##########################################################################################
### 1_get_data.R
### 
### This script reads the bulk RNA sequencing and microarray datasets and filters them
### to only include genes present in one common gene mapping list. It transforms the
### microarray data using 2^(***) to match the scale of the bulk RNA sequencing data values
### – this is used for InstaPrism deconvolution. It also transforms the bulk RNA sequencing
### data using log10(***+1) to match the scale of the microarray data – this is used for
### clustering. All of these matrices are saved in a uniform format containing on sample
### ID (rows) and genes (columns); file names are appended with either “asImported” or
### “transformed.” It also saves any metadata information (ex. prior clustering) about
### the samples in a separate file for reference.
##########################################################################################

# These are large files that take up a lot of memory - clear R's memory to make room
rm(list = ls())
gc()

# Load the renv library from the lockfile
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")  # Install renv if not already installed
}
library(renv)
renv::init()
renv::restore()

# Load Required Libraries
library(data.table)
library(dplyr)
library(Biobase)  # for ExpressionSet access
library(curatedOvarianData)
library(here)
library(readxl)

# Set base directories
# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
input_data <- file.path(here(), "input_data")
output_data <- file.path(here(), "output_data/bulk_datasets")

# Ensure output_data directory exists
dir.create(output_data, recursive = TRUE, showWarnings = FALSE)

# Load Gene List
MAD_genes <- data.frame(fread(file.path(input_data,"GlobalMAD_genelist.csv")))
colnames(MAD_genes) <- "hgnc_symbol"

##########################################################
# 1) Helper reading/merging functions 
##########################################################

# read Schildkraut data from file
read_format_expr <- function(in_file, metadata_table){
  # 1) Read file with fill=TRUE to handle irregular columns
  rnaseq_expr_df <- fread(in_file,
                          header = TRUE,
                          data.table = FALSE,
                          fill   = TRUE)
  
  # 2) Remove any trailing empty columns
  empty_cols <- sapply(rnaseq_expr_df, function(x) all(is.na(x)))
  if (any(empty_cols)) {
    rnaseq_expr_df <- rnaseq_expr_df[, !empty_cols, drop = FALSE]
  }
  
  # 3) The first column is gene IDs
  gene_ids <- rnaseq_expr_df[, 1]
  rnaseq_expr_df <- rnaseq_expr_df[, -1, drop = FALSE]
  
  # 4) Transpose
  rnaseq_expr_df <- data.frame(t(rnaseq_expr_df))
  colnames(rnaseq_expr_df) <- gene_ids
  
  # 5) Create sample ID column
  rnaseq_expr_df$ID <- gsub("^Sample_", "", rownames(rnaseq_expr_df))
  
  # 6) Merge with metadata by ID
  full_df <- merge(metadata_table, rnaseq_expr_df, by = "ID")
  
  # Return the combined data.frame as [[1]] and gene IDs as [[2]]
  return(list(full_df, gene_ids))
}

# read curatedOvarianData data from an in_df object
read_format_cOD_expr <- function(in_df, metadata_table){
  rnaseq_expr_df <- as.data.frame(in_df)  # ensure data.frame
  gene_ids <- rownames(rnaseq_expr_df)
  sample_ids <- colnames(rnaseq_expr_df)
  
  # Transpose
  rnaseq_expr_df <- data.frame(t(rnaseq_expr_df))
  colnames(rnaseq_expr_df) <- gene_ids
  
  # Create sample ID col
  rnaseq_expr_df$ID <- sample_ids
  
  # Merge with metadata
  full_df <- merge(metadata_table, rnaseq_expr_df, by = "ID")
  
  return(list(full_df, gene_ids))
}

##########################################################
# 2) Filter to MAD_genes & write out as CSV plus metadata
#    We'll define a function that saves both the 'as
#    imported' version and the 'transformed' version.
##########################################################

save_dual_versions <- function(expr_merged, dataset_name, transform_type = c("rnaseq", "microarray")) {
  expr_merged <- as.data.frame(expr_merged)
  
  # 1) Save as-imported (raw/log2) version
  save_filtered(expr_merged, dataset_name, suffix = "asImported")
  
  # 2) Create a “transformed” version
  #    For bulk RNA seq (SchildkrautB, SchildkrautW, TCGA_bulk), log10(...+1)
  #    For microarray (TCGA_microarray, Tothill, Yoshihara), do 2^(...)
  transform_type <- match.arg(transform_type)
  if (transform_type == "rnaseq") {
    expr_transformed <- transform_rnaseq(expr_merged)
    
    # 3) Save the transformed version
    save_filtered(expr_transformed, dataset_name, suffix = "transformed")
  } 
  if (transform_type == "microarray") {
    expr_transformed <- transform_microarray(expr_merged)
  
  # 3) Save the transformed version
  save_filtered(expr_transformed, dataset_name, suffix = "transformed")
}}

# This sub-function saves ONLY ID + gene columns in the expression CSV,
# and everything else in a separate metadata file.
save_filtered <- function(expr_merged, dataset_name, suffix) {
  # 1) Gene columns = those in the MAD_genes list
  valid_genes <- intersect(colnames(expr_merged), MAD_genes$hgnc_symbol)
  
  # 2) The final expression CSV has only: ID + (valid_genes)
  final_order <- c("ID", valid_genes)
  final_order <- intersect(final_order, colnames(expr_merged))
  expr_filtered <- expr_merged[, final_order, drop = FALSE]

  # Write expression CSV
  expr_output_file <- file.path(output_data,
    paste0(dataset_name, "_filtered_", suffix, ".csv")
  )
  write.csv(expr_filtered, expr_output_file, row.names = FALSE)
  
  # 3) Metadata file: everything else (excluding ID + valid_genes)
  meta_cols <- setdiff(colnames(expr_merged), final_order)
  metadata_only <- expr_merged[, meta_cols, drop = FALSE]
  
  meta_output_file <- file.path(output_data,
    paste0(dataset_name, "_metadata_", suffix, ".csv")
  )
  write.csv(metadata_only, meta_output_file, row.names = FALSE)
}

##########################################
# Helper transformations for each data type
##########################################

# 1) For RNA-seq => log10(... + 1)
transform_rnaseq <- function(expr_df) {
  df <- expr_df
  # We'll define 'ID' as non-gene. 
  # Every other column that is also in 'MAD_genes' is numeric. 
  # So we can just do:
  numeric_cols <- intersect(colnames(df), MAD_genes$hgnc_symbol)
  
  for (cc in numeric_cols) {
    df[[cc]] <- as.numeric(df[[cc]])
    df[[cc]] <- log10(df[[cc]] + 1)
  }
  return(df)
}

# 2) For microarray => 2^(...) to revert from log2
transform_microarray <- function(expr_df) {
  df <- expr_df
  # Again, only exponentiate gene columns
  numeric_cols <- intersect(colnames(df), MAD_genes$hgnc_symbol)
  
  for (cc in numeric_cols) {
    df[[cc]] <- as.numeric(df[[cc]])
    df[[cc]] <- 2^(df[[cc]])
  }
  return(df)
}

##########################################################
# 3) Load microarray metadata
##########################################################

clust_file <- file.path(input_data, "FullClusterMembership.csv")

if (file.exists(clust_file)) {
  clust_df <- data.frame(fread(clust_file))
  colnames(clust_df)[1] <- "ID"
} else {
  stop("Could not find FullClusterMembership.csv for microarray metadata.")
}

##########################################################
# 4) Process, transform, and save datasets
##########################################################

## A) SchildkrautB (RNA-seq)
message("Processing SchildkrautB (RNA-seq)")
schildB_metadata <- fread(file.path(input_data, "main_AA_metadata_table.tsv"))
schildB_res <- read_format_expr(
  in_file = file.path(input_data, "salmon_raw_counts_for_way_pipeline.tsv"),
  metadata_table = schildB_metadata
)
save_dual_versions(
  expr_merged = schildB_res[[1]],
  dataset_name = "SchildkrautB",
  transform_type = "rnaseq"
)

## B) SchildkrautW (RNA-seq)
message("Processing SchildkrautW (RNA-seq)")
schildW_metadata <- fread(file.path(input_data, "main_white_metadata_table.tsv"))
schildW_res <- read_format_expr(
  in_file = file.path(input_data, "salmon_raw_counts_for_way_pipeline_whites.tsv"),
  metadata_table = schildW_metadata
)
save_dual_versions(
  expr_merged = schildW_res[[1]],
  dataset_name = "SchildkrautW",
  transform_type = "rnaseq"
)

## C) TCGA_bulk (RNA-seq) and TCGA_microarray (Microarray)
## We will do these together to isolate only the samples in common between datasets

message("Processing TCGA_bulk (RNA-seq)")

  # Read TCGA bulk RNAseq
  tcga_bulk_dta <- fread(file.path(input_data,"EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"),
                                   data.table = FALSE)
  # Set the first column (genes) as rownames
  rownames(tcga_bulk_dta) <- tcga_bulk_dta[,1]
  tcga_bulk_dta <- tcga_bulk_dta[,-1]
  # These sample IDs (column names) do not exactly match the format of the sample IDs in the supplemental info table
  # Trim them down to match the TCGA-XX-XXXX patient barcode format of the supplemental info table
  # Ex. "TCGA-OR-A5J3-01A-11R-A29S-07" --> "TCGA-OR-A5J3"
  colnames(tcga_bulk_dta) <- sapply(strsplit(colnames(tcga_bulk_dta), "-"),
                           function(x) paste(x[1:3], collapse = "-"))
  # Read TCGA bulk supplemental sample information
  tcga_bulk_info <- readxl::read_excel(file.path(input_data,"TCGA-CDR-SupplementalTableS1.xlsx"),
                                                   sheet = "TCGA-CDR", na = "#N/A")
  # Use the supplemental info table to find the HGSOC sample IDs
  hgsoc_sample_indices <- tcga_bulk_info$type == "OV" &
                          tcga_bulk_info$histological_type == "Serous Cystadenocarcinoma" &
                          tcga_bulk_info$histological_grade %in% c("G3", "G4") 
  hgsoc_sample_IDs <- tcga_bulk_info$bcr_patient_barcode[hgsoc_sample_indices]
  # Filter TCGA bulk to only include HGSOC
  tcga_bulk_dta <- tcga_bulk_dta[,colnames(tcga_bulk_dta) %in% hgsoc_sample_IDs]
  # Change the hyphens in the TCGA_bulk sample IDs to periods to match the TCGA microarray sample ID formatting
  colnames(tcga_bulk_dta) <- gsub("-",".",colnames(tcga_bulk_dta))
  # The genes (row names) are denoted as "hgnc_symbol|entrez_ID"
  # Some have "?" for the hgnc_symbol, we will remove those genes
  # We will also remove any duplicated hgnc genes
  gene_names <- strsplit(rownames(tcga_bulk_dta), split = "\\|")
  gene_names <- sapply(gene_names, "[[", 1)
  rows_to_remove <- union(which(gene_names == "?"),
                          which(duplicated(gene_names)))
  tcga_bulk_dta <- tcga_bulk_dta[-rows_to_remove,]
  rownames(tcga_bulk_dta) <- gene_names[-rows_to_remove]
  # Set any NA values to zero
  tcga_bulk_dta[is.na(tcga_bulk_dta)] <- 0
  # Set any negative values to zero
  tcga_bulk_dta[tcga_bulk_dta < 0] <- 0

message("Processing TCGA_microarray (microarray)")
  
  # Read TCGA microarray
  data("TCGA_eset", package = "curatedOvarianData")
  tcga_microarray_dta <- exprs(TCGA_eset)
  # Filter TCGA microarray to only include HGSOC
  pheno <- pData(TCGA_eset)
  hgsoc_samples <- pheno$histological_type=="ser" & pheno$sample_type=="tumor" & pheno$summarygrade=="high"
  hgsoc_samples[is.na(hgsoc_samples)] <- FALSE
  tcga_microarray_dta <- tcga_microarray_dta[, hgsoc_samples]
  
  # Find samples in common between the two datasets
  tcga_common_samples <- intersect(colnames(tcga_bulk_dta), colnames(tcga_microarray_dta))
  # Filter datasets to only include common samples
  tcga_bulk_dta <- tcga_bulk_dta[, tcga_common_samples]
  tcga_microarray_dta <- tcga_microarray_dta[, tcga_common_samples]
  
  # Get metadata
  tcga_metadata <- subset(clust_df, Dataset == "TCGA")
  
  # Read TCGA bulk
  res_tcga_bulk <- read_format_cOD_expr(tcga_bulk_dta, tcga_metadata)
  save_dual_versions(
    expr_merged = res_tcga_bulk[[1]],
    dataset_name = "TCGA_bulk",
    transform_type = "rnaseq"
  )
  
  # Read TCGA microarray
  res_tcga_microarray <- read_format_cOD_expr(tcga_microarray_dta, tcga_metadata)
  save_dual_versions(
    expr_merged = res_tcga_microarray[[1]],
    dataset_name = "TCGA_microarray",
    transform_type = "microarray"
  )
  
## D) Tothill (Microarray)
message("Processing Tothill (microarray)")
data("GSE9891_eset", package = "curatedOvarianData")
tothill_dta <- exprs(GSE9891_eset)
# Filter to only include HGSOC
pheno <- pData(GSE9891_eset)
hgsoc_samples <- pheno$histological_type=="ser" & pheno$sample_type=="tumor" & pheno$summarygrade=="high"
hgsoc_samples[is.na(hgsoc_samples)] <- FALSE
tothill_dta <- tothill_dta[, hgsoc_samples]
# Get metadata
tothill_metadata <- subset(clust_df, Dataset == "Tothill")

res_tothill <- read_format_cOD_expr(tothill_dta, tothill_metadata)
save_dual_versions(
  expr_merged = res_tothill[[1]],
  dataset_name = "Tothill",
  transform_type = "microarray"
)

## E) Yoshihara (Microarray)
message("Processing Yoshihara (microarray)")
data("GSE32062.GPL6480_eset", package = "curatedOvarianData")
yoshi_dta <- exprs(GSE32062.GPL6480_eset)
# Filter to only include HGSOC
pheno <- pData(GSE32062.GPL6480_eset)
hgsoc_samples <- pheno$histological_type=="ser" & pheno$sample_type=="tumor" & pheno$summarygrade=="high"
hgsoc_samples[is.na(hgsoc_samples)] <- FALSE
yoshi_dta <- yoshi_dta[, hgsoc_samples]
# Get metadata
yoshi_metadata <- subset(clust_df, Dataset == "Yoshihara")

res_yoshi <- read_format_cOD_expr(yoshi_dta, yoshi_metadata)
save_dual_versions(
  expr_merged = res_yoshi[[1]],
  dataset_name = "Yoshihara",
  transform_type = "microarray"
)

message("All datasets transformed and saved as both asImported & transformed versions!")