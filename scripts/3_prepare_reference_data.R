##########################################################################################
### 3_prepare_reference_data.R
### 
### This script creates a single-cell/single-nucleus reference matrix for use as input into
### InstaPrism. It requires single cell RNA sequencing data of HGSOC
### (from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517) and cell type labels
### for those cells (from https://github.com/greenelab/deconvolution_pilot/tree/main/data/cell_labels).
### It also requires single nucleus RNA sequencing data of adipocytes
### (from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171). Please refer to the
### “Input and output data preparation and organization” subsection of README for specific
### files necessary. It reads in the HGSOC single cell and adipocyte single nucleus RNA
### sequencing data and performs some pre-processing on the adipocyte data (removing duplicate
### samples, removing non-adipocyte nuclei, removing samples with many mitochondrial gene
### reads, and using Seurat to remove low-quality nuclei, empty droplets, and nuclei
### doublets/multiplets). It then combines the single cell and single nucleus data to
### generate an expression matrix where each row corresponds to a gene (GeneCards symbols)
### and each column corresponds to a sample (each assigned a unique numerical ID). It also
### generates a cell type file which serves as a key, and associates each sample ID to
### its cell type. 
##########################################################################################

# These are very large files that take up a lot of memory - clear R's memory to make room
rm(list = ls())
gc()

# Load environment
renv::load()

# Load required libraries
library(dplyr)
library(Matrix)
library(data.table)
library(scuttle)
library(here)
library(Seurat)

# Set base directories
# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
input_data <- file.path(here(), "input_data")
output_data <- file.path(here(), "output_data")
reference_data_dir <- file.path(output_data,"sc_sn_reference_data")

# Ensure reference_data_dir directory exists
dir.create(reference_data_dir, recursive = TRUE, showWarnings = FALSE)

##########################################################
# 1) Read the scRNAseq reference datasets and associate
#    each cell with a pre-identified cell state
##########################################################

# The eight different datasets have numerous identifying numbers/names, for simplicity we will use rep1-rep8.
# Corresponding identifiers can be verified at (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517)
# We will create a data frame associating the different identifiers for reference.
identifiers <-
  data.frame(rep = c("rep1", "rep2", "rep3", "rep4", "rep5", "rep6", "rep7", "rep8"),
             gsm = c("GSM6720925", "GSM6720926", "GSM6720927", "GSM6720928", "GSM6720929", "GSM6720930", "GSM6720931", "GSM6720932"),
             num = c("2251", "2267", "2283", "2293", "2380", "2428", "2467", "2497"))

# Create functions for reading the files
read_matrix <- function(n) {
  matrix <- Matrix::readMM(file.path(input_data, paste0(identifiers[n,2],"_single_cell_matrix_", identifiers[n,3], ".mtx")))
  return(matrix)
}
read_barcodes <- function(n){
  df <- fread(file = (file.path(input_data, paste0(identifiers[n,2], "_single_cell_barcodes_", identifiers[n,3], ".tsv"))),
              header = FALSE,
              data.table = FALSE)
  return(df)
}
read_features <- function(n){
  df <- fread(file = (file.path(input_data, paste0(identifiers[n,2], "_single_cell_features_", identifiers[n,3], ".tsv"))),
              header = FALSE,
              data.table = FALSE)
  return(df)
}
read_labels <- function(n){
  df <- read.delim(file = (file.path(input_data, paste0(identifiers[n,3], "_labels.txt"))),
                   header = TRUE)
  return(df)
}

# Reading the files
for (n in c(seq.int(1,8))) {
  print(paste0("Reading scRNAseq files, rep ",n))
  assign(paste0("rep",n,"_sc_matrix"), read_matrix(n))
  assign(paste0("rep",n,"_sc_barcodes"), read_barcodes(n))
  assign(paste0("rep",n,"_sc_features"), read_features(n))
  assign(paste0("rep",n,"_sc_labels"), read_labels(n))
}

# Now we have read the files, including the matrices (they are sparse matrices of class dgTMatrix)
# The "barcodes" files are the column names for the matrices (cells), and the
# "features" files are the row names for the matrices (genes).
# Let's make sure that the numbers of rows and columns align.

for (n in c(seq.int(1,8))) {
  if (nrow(get(paste0("rep",n,"_sc_matrix"))) != nrow(get(paste0("rep",n,"_sc_features"))))
  {
    stop(paste0("Warning! Number of rows in matrix and number of rows in features do not align for rep ", n))
  }
  if (ncol(get(paste0("rep",n,"_sc_matrix"))) != nrow(get(paste0("rep",n,"_sc_barcodes"))))
  {
    stop(paste0("Warning! Number of columns in matrix and number of rows in barcodes do not align for rep ", n))
  }
}

# The cells contained in the "labels" files don't perfectly overlap with the cells contained in
# the matrices/"barcodes" files. Let's isolate the overlap.
for (n in c(seq.int(1,8))) {
  assign(paste0("rep",n,"_sc_overlap_cells"),
         intersect((get(paste0("rep",n,"_sc_barcodes")))[,1], (get(paste0("rep",n,"_sc_labels")))[,1]))
}

# Find the row numbers of overlap cells in "barcodes"
get_row_numbers <- function(df, v){
  return(which(df[,1] %in% v))
}

for (n in c(seq.int(1,8))) {
  assign(paste0("rep",n,"_sc_overlap_rows"),
         get_row_numbers(get(paste0("rep",n,"_sc_barcodes")), get(paste0("rep",n,"_sc_overlap_cells"))))
}

# Now we need to isolate only the columns of the matrices that match with these overlapping cells
for (n in c(seq.int(1,8))) {
  assign(paste0("rep",n,"_sc_matrix_subset"),
         get(paste0("rep",n,"_sc_matrix"))[,get(paste0("rep",n,"_sc_overlap_rows"))])
}

# Removing the large matrices from R's memory to save space
rm(rep1_sc_matrix, rep2_sc_matrix, rep3_sc_matrix, rep4_sc_matrix, rep5_sc_matrix, rep6_sc_matrix, rep7_sc_matrix, rep8_sc_matrix)
gc()


##########################################################
# 2) For the scRNA seq data, filter out any duplicate samples
#    within reps then create a large matrix with rows = genes
#    and columns = samples. Also create a key that associates
#    each sample ID with its cell type.
##########################################################

# First we will create dataframes that associate the overlapping cell barcodes with their cell type labels
for (n in c(seq.int(1,8))) {
  assign(paste0("rep",n,"_sc_overlap_cells_labeled"),
         filter(get(paste0("rep",n,"_sc_labels")),
                get(paste0("rep",n,"_sc_labels"))[,1] %in% get(paste0("rep",n,"_sc_overlap_cells"))))
}

# Removing the sc_labels from R's memory to save space
rm(rep1_sc_labels, rep2_sc_labels, rep3_sc_labels, rep4_sc_labels, rep5_sc_labels, rep6_sc_labels, rep7_sc_labels, rep8_sc_labels)
gc()

# Now let's remove any duplicate samples in both the _sc_overlap_cells_labeled cell type keys
# and the _sc_matrix_subset expression matrices.
for (n in c(1:8)){
  cell_types <- get(paste0("rep",n,"_sc_overlap_cells_labeled"))
  matrix <- get(paste0("rep",n,"_sc_matrix_subset"))
  
  # Ensuring that the number of samples is the same
  if(nrow(cell_types) != ncol(matrix)){
    stop(paste0("Warning! Number of samples in scRNAseq data rep ",n," do not match up between the cell type labels and the expression matrix."))
  }
  
  # Identifying duplicate sample barcode indices
  indices <- which(duplicated(cell_types$Barcode))
  print(paste0("Identified ",length(indices)," duplicate samples in scRNAseq data rep ",n,", of ",nrow(cell_types)," total samples."))
  if (length(indices) == 0){
    next
  }

  # Removing any duplicate samples
  cell_types_filt <- cell_types[-indices,]
  matrix_filt <- matrix[,-indices]

  # Ensuring that the remaining number of samples matchs up between filtered cell type key and filtered matrix
  if(nrow(cell_types_filt) != ncol(matrix_filt)){
    stop("Warning! Number of samples does not match up between cell type labels and expression matrix for scRNAseq data rep ",n," after removal of duplicate samples.")
  }
  
  # Assign back to original object names
  assign(paste0("rep",n,"_sc_overlap_cells_labeled"),
         cell_types_filt)
  assign(paste0("rep",n,"_sc_matrix_subset"),
         matrix_filt)
}

# Then we will create combined key that aligns cell barcodes with cell type
all_sc_overlap_cells_labeled <- rbind(rep1_sc_overlap_cells_labeled,
                                      rep2_sc_overlap_cells_labeled,
                                      rep3_sc_overlap_cells_labeled,
                                      rep4_sc_overlap_cells_labeled,
                                      rep5_sc_overlap_cells_labeled,
                                      rep6_sc_overlap_cells_labeled,
                                      rep7_sc_overlap_cells_labeled,
                                      rep8_sc_overlap_cells_labeled)

# Removing the overlap cells objects from R's memory to save space
rm(rep1_sc_overlap_cells_labeled,
   rep2_sc_overlap_cells_labeled,
   rep3_sc_overlap_cells_labeled,
   rep4_sc_overlap_cells_labeled,
   rep5_sc_overlap_cells_labeled,
   rep6_sc_overlap_cells_labeled,
   rep7_sc_overlap_cells_labeled,
   rep8_sc_overlap_cells_labeled)
gc()

# Let's ensure that the row names (genes), contained in the "features" objects, are
# identical between reps so that the matrices can be accurately bound together horizontally.
for (n in c(seq.int(2,8))) {
  if(!identical(get(paste0("rep",n,"_sc_features")),rep1_sc_features)){
    stop("Warning! Row names (genes) are not identical between scRNAseq reps.")
  }
}

# Now we can create a big sparse matrix that combines the information in the matrices of all 8 reps
all_sc_expr <- cbind(rep1_sc_matrix_subset,
                     rep2_sc_matrix_subset,
                     rep3_sc_matrix_subset,
                     rep4_sc_matrix_subset,
                     rep5_sc_matrix_subset,
                     rep6_sc_matrix_subset,
                     rep7_sc_matrix_subset,
                     rep8_sc_matrix_subset)

colnames(all_sc_expr) <- all_sc_overlap_cells_labeled[,1]
#using cell barcodes for column names
rownames(all_sc_expr) <- rep1_sc_features[,2]

# ────────────────────────────────────────────────────────────────
# ► REMOVE cells labelled "Unknown1" or "Unknown2"
unwanted <- c("Unknown1", "Unknown2")
keep      <- !(all_sc_overlap_cells_labeled$cellType %in% unwanted)

all_sc_overlap_cells_labeled <- all_sc_overlap_cells_labeled[keep, ]
all_sc_expr <- all_sc_expr[, keep]    # drop the same columns in the matrix
# ────────────────────────────────────────────────────────────────

#using GeneCards symbols for gene names as that is what is used in the adipocyte snRNAseq data
# Removing the subsetted matrices from R's memory to save space
rm(rep1_sc_matrix_subset,
   rep2_sc_matrix_subset,
   rep3_sc_matrix_subset,
   rep4_sc_matrix_subset,
   rep5_sc_matrix_subset,
   rep6_sc_matrix_subset,
   rep7_sc_matrix_subset,
   rep8_sc_matrix_subset)
gc()

##########################################################
# 3) Read the adipose tissue snRNAseq data, subset to only
#    include genes in common with scRNAseq data
##########################################################

# Read the files
adipose_file_names <- c("GSM5359325_Hs_OAT_01-1.dge.tsv",
                        "GSM5359326_Hs_OAT_01-2.dge.tsv",
                        "GSM5359327_Hs_OAT_253-1.dge.tsv",
                        "GSM5359328_Hs_OAT_254-1.dge.tsv",
                        "GSM5359329_Hs_OAT_255-1.dge.tsv",
                        "GSM5359330_Hs_OAT_256-1.dge.tsv",
                        "GSM5359331_Hs_SAT_01-1.dge.tsv",
                        "GSM5359332_Hs_SAT_02-1.dge.tsv",
                        "GSM5359333_Hs_SAT_04-1.dge.tsv",
                        "GSM5359334_Hs_SAT_253-1.dge.tsv",
                        "GSM5359335_Hs_SAT_254-1.dge.tsv",
                        "GSM5359336_Hs_SAT_255-1.dge.tsv",
                        "GSM5359337_Hs_SAT_256-1.dge.tsv",
                        "GSM5820679_Hs_OAT_09-1.dge.tsv",
                        "GSM5820680_Hs_OAT_10-1.dge.tsv",
                        "GSM5820684_Hs_SAT_09-1.dge.tsv",
                        "GSM5820685_Hs_SAT_10-1.dge.tsv",
                        "GSM5820686_Hs_SAT_11-1.dge.tsv")

# Use readSparseCounts from "scuttle" package to read these large files
for (file in adipose_file_names) {
  print(paste0("Reading ",file," (adipose",which(adipose_file_names == file),")"))
  assign(paste0("adipose",(which(adipose_file_names == file))),
         readSparseCounts(file.path(input_data, file), sep="\t", quote=NULL, comment.char="", row.names=TRUE,
                          col.names=TRUE, ignore.row=0L, skip.row=0L, ignore.col=0L, skip.col=0L, chunk=1000L))
}

# Let's append the sample name to the column names (cell IDs), to match the formatting of the metadata file
for (n in c(seq.int(1,length(adipose_file_names)))) {
  matrix <- get(paste0("adipose",n))
  new_colnames <- colnames(matrix)
  sample_name <- sub("GSM\\d+_([^_]+_[^_]+_[^_]+).dge.tsv", "\\1", adipose_file_names[n])
  new_colnames <- paste(sample_name, new_colnames, sep="_")
  colnames(matrix) <- new_colnames
  assign(paste0("adipose",n), matrix)
}

# Make a list of the sn adipose objects
adipose_object_list <- mget(paste0("adipose", 1:length(adipose_file_names)))

# Let's find the row names (genes) that are shared in common between adipose files
adipose_genes <- lapply(adipose_object_list, rownames)
adipose_common_genes <- Reduce(intersect, adipose_genes)

# Now we will find the genes (as GeneCards symbols) in common between the the single cell and the adipose data.
all_overlap_genes <- intersect(rep1_sc_features[,2], adipose_common_genes)
print(paste0("Identified ",length(all_overlap_genes)," total genes in common between sc and sn RNAseq data."))

# Now we will subset each adipose matrix to only include common genes
for (n in c(seq.int(1,length(adipose_file_names)))) {
  assign(paste0("adipose",n,"_subset"),
         get(paste0("adipose",n))[all_overlap_genes,])
}

# Removing the large matrices from R's memory to save space
rm(list = ls(pattern = "^sn_adipose\\d+$"))
gc()

##########################################################
# 4) Preprocessing adipose snRNAseq data: remove duplicates,
#    remove non-adipocytes nuclei, remove nuclei with high
#    mitochondrial gene read levels, remove low-quality
#    cells, empty droplets, and cell doublets/multiplets
#    using Seurat. Then, combine all processed samples
#    together into one large gene expression matrix.
##########################################################

# Remove any duplicate samples (columns) from the expression matrices
for (n in c(seq.int(1,length(adipose_file_names)))) {
  matrix <- get(paste0("adipose",n,"_subset"))
  indices <- which(duplicated(colnames(matrix)))
  print(paste0("Identified ",length(indices)," duplicate samples in adipose",n,", of ",ncol(matrix)," total samples."))
  if (length(indices) == 0){
    next
  }
  filt_matrix <- matrix[,-indices]
  assign(paste0("adipose",n,"_subset"), filt_matrix)
}

# Removing the large objects from R's memory to save space
rm(matrix, indices)
gc()

# The cell types in these tissue have been labeled by Emont et al. (Nature 2022). 
# Use these labels to remove all non-adipocyte nuclei.

# Read cell metadata
# This is a table where each row corresponds to one nucleus (sample)
adipose_metadata <- fread(file.path(input_data,"GSE176171_cell_metadata.tsv"), data.table = FALSE)

# Isolate the nucleus (sample) IDs of the adipocytes
adipocytes_IDs <- adipose_metadata[adipose_metadata$cell_type__custom == "adipocyte",]
adipocytes_IDs <- adipocytes_IDs$cell_id

# Subset each expression matrix to only include the adipocytes
for (n in c(seq.int(1,length(adipose_file_names)))) {
  matrix <- get(paste0("adipose",n,"_subset"))
  adipocytes <- intersect(colnames(matrix), adipocytes_IDs)
  print(paste0("Identified ",length(adipocytes)," adipocytes in adipose",n," of ",ncol(matrix)," total samples."))
  assign(paste0("adipose",n,"_subset"),
         matrix[, colnames(matrix) %in% adipocytes])
}

# Removing the large objects from R's memory to save space
rm(matrix, adipocytes)
gc()

# As these are nuclei, there should not be any mitochondrial gene RNAs present.
# However, nearly all of these samples have some amount of mitochondrial gene RNA.
# Therefore we will filter out any samples with >10 total mitochondrial gene reads,
# as well as any sample with >5% of all reads coming from mitochondrial genes.

# Removing samples with >50 total mitochondrial gene reads
for (n in c(seq.int(1,length(adipose_file_names)))) {
  # Find mitochondrial genes (gene names starting with "MT-")
  mito_gene_indices <- grep("MT-", rownames(get(paste0("adipose",n,"_subset"))))
  mito_genes <- get(paste0("adipose",n,"_subset"))[mito_gene_indices,]
  
  # Select the samples (columns) that have >10 reads of any mitochondrial genes
  mito_samples <- which(colSums(mito_genes) > 10)
  print(paste0("Identified ",length(mito_samples)," samples with >50 mitochondrial gene reads in adipose",
               n,"_subset, of ",ncol(mito_genes)," total samples."))
  
  if(length(mito_samples) != 0){
    # Remove the identified columns from the matrix
    assign(paste0("adipose",n,"_subset"),
           get(paste0("adipose",n,"_subset"))[,-mito_samples])
    }
}

# Removing samples with >% of all reads coming from mitochondrial genes
for (n in c(seq.int(1,length(adipose_file_names)))) {
  # Find mitochondrial genes (gene names starting with "MT-")
  mito_gene_indices <- grep("MT-", rownames(get(paste0("adipose",n,"_subset"))))
  
  # Find the sum of all mitochondrial reads and the sum of all gene reads,
  # then divide to find the proportion
  mito_genes_sum <- as.vector(colSums(get(paste0("adipose",n,"_subset"))[mito_gene_indices,]))
  all_genes_sum <- as.vector(colSums(get(paste0("adipose",n,"_subset"))))
  proportions <- mito_genes_sum / all_genes_sum
  # Select the samples which have >2% of all reads coming from mitochondrial genes
  indices <- which(proportions > 0.02)
  print(paste0("Identified ",length(indices)," samples with >2% of all reads coming from mitochondrial genes in adipose",
               n,"_subset, of ",length(mito_genes_sum)," total samples."))
  
  if(length(indices) != 0){
    # Remove the identified columns from the matrix
    assign(paste0("adipose",n,"_subset"),
           get(paste0("adipose",n,"_subset"))[,-indices])
  }
}

# Now we will use the Seurat package to remove low-quality nuclei, empty droplets, and nuclei doublets/multiplets.
# We will filter out samples that express fewer than 200 unique genes (likely low quality), 
# as well as samples that express greater than 2500 unique genes (likely doublet/multiplet).
for (n in seq_along(adipose_file_names)) {
  mat <- get(paste0("adipose",n,"_subset"))
  if (ncol(mat) == 0) {
    message("Skipping adipose", n, ": no cells remain after mito filtering")
    next
  }
  obj <- CreateSeuratObject(counts = mat)
  low_genes <- which(obj$nFeature_RNA < 200)
  high_genes <- which(obj$nFeature_RNA > 2500)
  print(paste0("Identified ",length(low_genes)," samples with <200 unique genes expressed and ",length(high_genes),
               " samples with >4000 unique genes expressed in adipose",n,"_subset, of ",ncol(obj)," total samples."))
  remove <- c(low_genes, high_genes)
  if(length(remove) != 0){
    assign(paste0("adipose",n,"_subset"),
           get(paste0("adipose",n,"_subset"))[,-remove])
  }
}

# Let's ensure that the row names (genes), are identical between adipose
# samples so that the matrices can be accurately bound together horizontally.
for (n in c(seq.int(2,length(adipose_file_names)))) {
  if(!identical(rownames(get(paste0("adipose",n,"_subset"))),rownames(adipose1_subset))){
    stop("Warning! Row names (genes) are not identical between adipose samples.")
  }
}

# Finally we can combine these into one large adipose matrix
# Create a list of objects
adipose_subset_object_list <- mget(paste0("adipose",1:length(adipose_file_names),"_subset"))

# Extract expression matrices from each object
get_expression_matrices <- lapply(adipose_subset_object_list, function(obj) {
  return(obj) 
})

# Bind the matrices together horizontally
all_adipose_expr <- do.call(cbind, get_expression_matrices)

# Removing the subsetted matrices from R's memory to save space
rm(adipose_subset_object_list)
rm(list = ls(pattern = "^sn_adipose\\d+_subset$"))
gc()

##########################################################
# 5) Combine the scRNAseq data and the adipose snRNAseq data
##########################################################

# Subsetting the scRNAseq data to only contain overlapping genes with the snRNAseq data
all_sc_expr_overlap <- (all_sc_expr)[all_overlap_genes, ]

# Removing all_sc_expr from R's memory to save space
rm(all_sc_expr)
gc()

# Creating a large matrix with all of the scRNAseq and adipose cells as columns and all genes in common as rows.
all_expr_final <- cbind(all_sc_expr_overlap,
                        all_adipose_expr)

# Removing all_sc_expr_overlap from R's memory to save space
rm(all_sc_expr_overlap)
gc()

# Adding the adipocyte data to the key
adipocyte_cells_labeled <- data.frame(Barcode = colnames(all_adipose_expr),
                                      cellType = "Adipocytes",
                                      stringsAsFactors = FALSE)

# Removing all_adipose_expr from R's memory to save space
rm(all_adipose_expr)
gc()

# Combining the adipocyte cell type key with the scRNA seq cell type key
all_cell_type_final <- rbind(all_sc_overlap_cells_labeled, adipocyte_cells_labeled)

# Removing all_sc_overlap_cells_labeled, adipocyte_cells_labeled from R's memory to save space
rm(all_sc_overlap_cells_labeled, adipocyte_cells_labeled)
gc()

# Replacing all cell barcodes with a unique numeric ID number
all_cell_type_final <- all_cell_type_final %>%
  mutate(cell_id = 1:nrow(all_cell_type_final))

colnames(all_expr_final) <- 1:ncol(all_expr_final)

##########################################################
# 6) Write final files
##########################################################

writeMM(all_expr_final, file.path(reference_data_dir, "reference_expr_final.csv"))
write.csv(rownames(all_expr_final), file.path(reference_data_dir, "reference_expr_final_genes.csv"))
write.csv(all_cell_type_final, file.path(reference_data_dir, "reference_cell_types_final.csv"), row.names=FALSE)
