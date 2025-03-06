##########################################################################################
### 2_get_clustering.R
### 
### This script performs k-means clustering (k=2,3,4), NMF clustering (k=2,3,4), and 
### consensusOV subtyping for each bulk dataset. For k-means clustering, it uses log10(...+1)
### transformed data (for RNAseq), and raw log2 data (for microarray). For NMF clustering
### and consensusOV subtyping, it uses raw counts (for RNAseq), and 2^(...) scaled
### "pseudocounts" data (for microarray). It saves the results in one csv file per dataset,
### where each row corresponds to one sample. It also saves a csv file containing the
### results for all datasets.
##########################################################################################

##########################################################
# 1) Set up the environment
##########################################################

# Clear R's memory
rm(list = ls())
gc()

# Load environment
renv::load()

# Load Required Libraries
library(data.table)
library(dplyr)
library(caret)
library(cluster)
library(NMF)
library(consensusOV)
library(here)

# Set base directories
# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
input_data <- file.path(here(), "input_data")
output_data <- file.path(here(), "output_data")
clustering_results_dir <- file.path(output_data, "bulk_data_clustering_and_subtypes")

# Ensure clustering_results_dir directory exists
dir.create(clustering_results_dir, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(5)

# Load gene_map, needed for real Entrez IDs for consensusOV
gene_map_file <- file.path(input_data, "ensembl_hgnc_entrez.tsv")
if (!file.exists(gene_map_file)) {
  stop("gene_map_file not found!")
}
gene_map <- data.frame(fread(gene_map_file)) 
# We expect columns like: ensembl_gene_id, hgnc_symbol, entrezgene_id, etc.

##########################################################
# 1) Helper reading, clustering, and subtyping functions 
##########################################################

# read_dataset_counts():
# Reads a CSV (with ID + gene columns) that has *rows=samples* after we set rownames=ID.
# This function reads raw counts (for RNAseq data - asImported), and
# 2^(...) scaled "pseudocounts" data (for microarray - transformed).
read_dataset_counts <- function(dataset_name) {
  if (dataset_name == "SchildkrautB" | dataset_name == "SchildkrautW" | dataset_name == "TCGA_bulk"){
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_asImported.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)  # row=sample, col=gene
  } else {
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_transformed.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)
  }
}

# read_dataset_log():
# Reads a CSV (with ID + gene columns) that has *rows=samples* after we set rownames=ID.
# This function reads log10(...+1) transformed data (for RNAseq data - transformed), and
# raw log2 data (for microarray - asImported).
read_dataset_log <- function(dataset_name) {
  if (dataset_name == "SchildkrautB" | dataset_name == "SchildkrautW" | dataset_name == "TCGA_bulk"){
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_transformed.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)  # row=sample, col=gene
  } else {
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_asImported.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)
  }
}

# run_kmeans():
# We will run the kmeans clustering on raw log2/log10(...+1) scaled data
# The matrix is (samples × genes).
run_kmeans <- function(samp_x_gene) {
  out_df <- data.frame(ID = rownames(samp_x_gene))
  for (k in c(2, 3, 4)) {
    km <- kmeans(samp_x_gene, centers=k, nstart=25)
    out_df[[paste0("ClusterK", k, "_kmeans")]] <- km$cluster
  }
  return(out_df)
}

# run_nmf():
# We will run the NMF clustering on raw counts/2^(...) scaled data
# For NMF, we need row=genes, col=samples => transpose
run_nmf <- function(samp_x_gene) {
  gene_x_samp <- t(samp_x_gene)  # row=genes, col=samples
  # ensure positivity
  minval <- min(gene_x_samp, na.rm=TRUE)
  if (minval < 0) {
    gene_x_samp <- gene_x_samp - minval + 1e-3
  }
  out_df <- data.frame(ID=rownames(samp_x_gene)) 
  for (k in c(2,3,4)) {
    nmf_res <- nmf(gene_x_samp, rank=k, nrun=10, .options='v')
    clust <- predict(nmf_res)
    # reorder to match rownames(samp_x_gene)
    sampleIDs <- colnames(gene_x_samp)  # same as rownames(samp_x_gene)
    out_df[[paste0("ClusterK", k, "_NMF")]] <- clust[sampleIDs]
  }
  return(out_df)
}

# run_consensusOV():
# We will run cOV subtyping on raw counts/2^(...) scaled data
# cOV expects row=genes, col=samples, plus *actual Entrez IDs* for geneIDs.
# We'll do the same transpose, then map gene symbols -> entrez IDs
run_consensusOV <- function(samp_x_gene) {
  # transpose => row=genes, col=samples
  gene_x_samp <- t(samp_x_gene)
  
  # rownames(gene_x_samp) are gene symbols. We must find their entrez IDs in gene_map
  # We'll keep only genes that appear in gene_map. 
  # Because cOV needs a numeric vector of geneIDs in the same order as the rows.
  all_genes <- rownames(gene_x_samp)
  # subset gene_map to these symbols
  map_sub <- subset(gene_map, hgnc_symbol %in% all_genes & !is.na(entrezgene_id))
  # remove duplicates if any
  map_sub <- unique(map_sub[, c("hgnc_symbol","entrezgene_id")])
  # reorder map_sub to match rownames(gene_x_samp)
  map_sub <- map_sub[match(all_genes, map_sub$hgnc_symbol), ]
  
  # Some genes might be missing from map_sub => NA in entrezgene_id => exclude them
  keep_idx <- which(!is.na(map_sub$entrezgene_id))
  # filter gene_x_samp to only those rows
  gene_x_samp_filt <- gene_x_samp[keep_idx, , drop=FALSE]
  gene_ids_W <- map_sub$entrezgene_id[keep_idx]
  
  if (nrow(gene_x_samp_filt) < 2) {
    # too few genes => can't run cOV meaningfully
    # return all NA
    return(data.frame(
      ID=rownames(samp_x_gene),
      consensusOV=rep(NA_character_, nrow(samp_x_gene))
    ))
  }
  
  # Now run get.subtypes
  subtypes_res <- get.subtypes(gene_x_samp_filt, gene_ids_W, method='consensus')
  
  # Generate a dataframe that aligns sample ID with consensus subtype
  out_df <- data.frame(ID = rownames(subtypes_res$rf.probs),
                       consensusOVsubtype = subtypes_res$consensusOV.subtypes)
  
  # Ensure that no samples were accidentally filtered out
  if (nrow(out_df) != nrow(samp_x_gene)) {
    warning("consensusOV did not generate a result for all samples in ", ds)
  }
  
  return(out_df)
}

# combine_clusterings():
# merges multiple data frames with the same 'ID' column.
combine_clusterings <- function(...) {
  dfs <- list(...)
  out <- Reduce(function(x, y) merge(x, y, by='ID', all=TRUE), dfs)
  return(out)
}

##################################################
# 3) Main Script
##################################################

dataset_list <- c("SchildkrautB", "SchildkrautW", "TCGA_bulk",
                  "TCGA_microarray", "Tothill", "Yoshihara")

all_clustering_res <- list()

for (ds in dataset_list) {
  message("Clustering dataset: ", ds)
  
  # 1) Load counts and log transformed (row=sample, col=gene)
  expr_df_counts <- read_dataset_counts(ds)
  expr_df_log <- read_dataset_log(ds)
  
  # 2) K-means
  kmeans_out <- run_kmeans(expr_df_log)
  
  # 3) NMF
  nmf_out <- run_nmf(expr_df_counts)
  
  # 4) consensusOV
  cov_out <- run_consensusOV(expr_df_counts)
  
  # 5) Combine
  combined <- combine_clusterings(kmeans_out, nmf_out, cov_out)
  combined$Dataset <- ds
  
  # 6) Write per-dataset
  outf <- file.path(clustering_results_dir, paste0(ds, "_clustering_labels.csv"))
  write.csv(combined, outf, row.names=FALSE)
  
  all_clustering_res[[ds]] <- combined
}

final_df <- dplyr::bind_rows(all_clustering_res)
final_out <- file.path(clustering_results_dir, "all_datasets_clustering_labels.csv")
write.csv(final_df, final_out, row.names=FALSE)

message("Done! Clustering results written to ", final_out)
