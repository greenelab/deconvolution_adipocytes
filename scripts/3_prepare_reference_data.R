
##############################
### The following script creates reference data from single cell and single nucleus RNA seq
### for use with InstaPrism. 
##############################

# This script requires data that can be downloaded from
# (https://github.com/greenelab/deconvolution_pilot/tree/main/data/cell_labels) and
# (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517)
# Please see README for specific files necessary

# Set the working directory
setwd("/Users/grace/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/greene_lab_projects/hgsoc_adipocyte_deconvolution")

# Set directory for location of reference data
reference_data_dir <- file.path(getwd(), "reference_data")

# Load required libraries
library(dplyr)

##############################
### 1) Read the scRNAseq reference datasets and associate each cell with a pre-identified cell state
##############################

# The eight different datasets have numerous identifying numbers/names, for simplicity we will use rep1-rep8.
# Corresponding identifiers can be verified at (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517)
# We will create a data frame associating the different identifiers for reference.
identifiers <-
  data.frame(rep = c("rep1", "rep2", "rep3", "rep4", "rep5", "rep6", "rep7", "rep8"),
             gsm = c("GSM6720925", "GSM6720926", "GSM6720927", "GSM6720928", "GSM6720929", "GSM6720930", "GSM6720931", "GSM6720932"),
             num = c("2251", "2267", "2283", "2293", "2380", "2428", "2467", "2497"))

# Create functions for reading the files
read_matrix <- function(n) {
  matrix <- Matrix::readMM(file.path(reference_data_dir, paste0(identifiers[n,2],"_single_cell_matrix_", identifiers[n,3], ".mtx")))
  return(matrix)
}
read_barcodes <- function(n){
  df <- read.table(file = (file.path(reference_data_dir, paste0(identifiers[n,2], "_single_cell_barcodes_", identifiers[n,3], ".tsv"))),
                   header = FALSE)
  return(df)
}
read_features <- function(n){
  df <- read.table(file = (file.path(reference_data_dir, paste0(identifiers[n,2], "_single_cell_features_", identifiers[n,3], ".tsv"))),
                   header = FALSE)
  return(df)
}
read_labels <- function(n){
  df <- read.delim(file = (file.path(reference_data_dir, paste0(identifiers[n,3], "_labels.txt"))),
                   header = TRUE)
  return(df)
}

# Reading the files
for (n in c(seq.int(1,8))) {
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

##############################
### 2) Create the information required to input into InstaPrism
##############################

# First we will create dataframes that associate the overlapping cell barcodes with their cell type labels
for (n in c(seq.int(1,8))) {
  assign(paste0("rep",n,"_sc_overlap_cells_labeled"),
         filter(get(paste0("rep",n,"_sc_labels")),
                get(paste0("rep",n,"_sc_labels"))[,1] %in% get(paste0("rep",n,"_sc_overlap_cells"))))
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

# Let's ensure that the row names (genes), contained in the "features" objects, are
# identical between reps so that the matrices can be accurately bound together horizontally.
for (n in c(seq.int(2,8))) {
  if(!identical(get(paste0("rep",n,"_sc_features")),rep1_sc_features)){
    stop("Warning! Row names (genes) are not identical between reps.")
  }
}

# Now we can create a data frame that combines the information in the matrices of all 8 reps
for (n in c(seq.int(1,8))) {
  assign(paste0("rep",n,"_sc_df"),
         as.data.frame(as.matrix(get(paste0("rep",n,"_sc_matrix_subset")))))
}

all_sc_expr <- cbind(rep1_sc_df,
                     rep2_sc_df,
                     rep3_sc_df,
                     rep4_sc_df,
                     rep5_sc_df,
                     rep6_sc_df,
                     rep7_sc_df,
                     rep8_sc_df)

##############################
### 3) Read the adipocyte snRNAseq data
##############################


##############################
### 4) Combine the scRNAseq data and the adipocyte snRNAseq data
##############################


#only use genes in common








