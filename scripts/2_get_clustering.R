###############################################
# 2_get_clustering.R
# Perform KMeans, NMF, and consensusOV clusterings
# for each dataset. Saves results in a CSV,
# matching original approach (rows=samples for k-means).
###############################################

# 1) Set up environment
# ======================
# Set the working directory
setwd("/Users/grace/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/greene_lab_projects/hgsoc_adipocyte_deconvolution")

# Load Required Libraries
library(data.table)
library(dplyr)
library(caret)
library(cluster)
library(NMF)        # for NMF-based clustering
library(consensusOV)
library(ggplot2)

# Set Base Directories
proj_dir <- file.path(getwd(), "prior_data")
bulk_data_dir <- file.path(getwd(), "data/bulks")
clustering_results_dir <- file.path(bulk_data_dir, "subtype_clusters")
dir.create(clustering_results_dir, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(5)

############################################################
# 1) Load gene_map, needed for real Entrez IDs for cOV
############################################################
gene_map_file <- file.path(proj_dir, "reference_data/ensembl_hgnc_entrez.tsv")
if (!file.exists(gene_map_file)) {
  stop("gene_map_file not found!")
}
gene_map <- data.frame(fread(gene_map_file)) 
# We expect columns like: ensembl_gene_id, hgnc_symbol, entrezgene_id, etc.

############################################################
# 2) Read data
############################################################

# We want the log10 transformed version of these datasets as that is what was used in Davidson 2024:
# (https://github.com/greenelab/hgsc_characterization/blob/master/figure_notebooks/compare_centroids.Rmd)
# The Schildkraut RNAseq transformed files and the Mayo microarray asImported file
# are already in this log10 format so we can use them directly.
# We will import the TCGA, Tothill, and Yoshihara microarray transformed files (in "counts"),
# then perform log10(...+1) transformation on them.
dataset_list <- c("SchildkrautB", "SchildkrautW", "TCGA", "Mayo", "Tothill", "Yoshihara")

for (ds in dataset_list) {
  if (ds == "Mayo") {
    assign(paste0(ds,"_bulk"),
           fread(file = (file.path(bulk_data_dir,paste0(ds, "_filtered_asImported.csv"))), header = TRUE, data.table = FALSE))
  } else {
    assign(paste0(ds,"_bulk"),
           fread(file = (file.path(bulk_data_dir,paste0(ds, "_filtered_transformed.csv"))), header = TRUE, data.table = FALSE))
  }
}

# Set column 1 (sample IDs) as rownames
for (ds in dataset_list) {
  assign(paste0(ds,"_bulk"),
         data.frame(get(paste0(ds,"_bulk"))[,-1], row.names=get(paste0(ds,"_bulk"))[,1]))
}

# Use log10(...+1) to transform the TCGA, Tothill, and Yoshihara 
log10transform <- function(df){
  transformed <- log10(df + 1)
  return(transformed)
}

TCGA_bulk <- log10transform(TCGA_bulk)
Tothill_bulk <- log10transform(Tothill_bulk)
Yoshihara_bulk <- log10transform(Yoshihara_bulk)
# Now all 6 bulk datasets are in a log10 transformed format
# Rows = samples, columns = genes

# Isolate only the columns (genes) in common between the datasets
all_genes <- list(colnames(SchildkrautB_bulk), colnames(SchildkrautW_bulk), colnames(TCGA_bulk),
                  colnames(Mayo_bulk), colnames(Tothill_bulk), colnames(Yoshihara_bulk))
common_genes <- Reduce(intersect, all_genes)

# Subset each dataset to only include common genes
for (ds in dataset_list) {
  assign(paste0(ds,"_subset"),
         get(paste0(ds,"_bulk"))[,common_genes])
}

############################################################
# 3) Helper Functions
############################################################

# run_kmeans():
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
# 4) Main Script
##################################################

dataset_list <- c("SchildkrautB", "SchildkrautW", "TCGA", "Mayo", "Tothill", "Yoshihara")

all_clustering_res <- list()

for (ds in dataset_list) {
  message("Clustering dataset: ", ds)
  
  # 1) Load => row=sample, col=gene
  expr_df <- get(paste0(ds,"_subset"))
  
  # 2) K-means
  kmeans_out <- run_kmeans(expr_df)
  
  # 3) NMF
  nmf_out <- run_nmf(expr_df)
  
  # 4) consensusOV
  cov_out <- run_consensusOV(expr_df)
  
  # 5) Combine
  combined <- combine_clusterings(kmeans_out, nmf_out, cov_out)
  combined$Dataset <- ds
  
  # 6) Write per-dataset
  outf <- file.path(clustering_results_dir, paste0(ds, "_clustering_labels.csv"))
  write.csv(combined, outf, row.names=FALSE)
  
  all_clustering_res[[ds]] <- combined
}

final_df <- dplyr::bind_rows(all_clustering_res)
final_out <- file.path(clustering_results_dir, "AllDatasets_ClusteringLabels.csv")
write.csv(final_df, final_out, row.names=FALSE)

message("Done! Clustering results written to ", final_out)

