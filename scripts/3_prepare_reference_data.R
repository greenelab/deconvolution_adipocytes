
##############################
### The following script creates reference data from single cell and single nucleus RNA seq
### for use with InstaPrism. 
##############################

# This script requires data that can be downloaded from
# (https://github.com/greenelab/deconvolution_pilot/tree/main/data/cell_labels),
# (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517),
# and (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171)
# Please see README for specific files necessary

# These are very large files that take up a lot of memory - clear R's memory to make room
rm(list = ls())
gc()

# # Set the working directory
# setwd("projects/gakatsu@xsede.org/hgsoc_adipocyte_deconvolution")

# Set directory for location of reference data
reference_data_dir <- file.path("/projects/gakatsu@xsede.org/hgsoc_adipocyte_deconvolution/reference_data")

# Load required libraries
library(dplyr)
library(Matrix)
library(data.table)

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
  df <- fread(file = (file.path(reference_data_dir, paste0(identifiers[n,2], "_single_cell_barcodes_", identifiers[n,3], ".tsv"))),
              header = FALSE,
              data.table = FALSE)
  return(df)
}
read_features <- function(n){
  df <- fread(file = (file.path(reference_data_dir, paste0(identifiers[n,2], "_single_cell_features_", identifiers[n,3], ".tsv"))),
              header = FALSE,
              data.table = FALSE)
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

# Removing the large matrices from R's memory to save space
rm(rep1_sc_matrix, rep2_sc_matrix, rep3_sc_matrix, rep4_sc_matrix, rep5_sc_matrix, rep6_sc_matrix, rep7_sc_matrix, rep8_sc_matrix)
gc()

##############################
### 2) Create the information required to input into InstaPrism
##############################

# First we will create dataframes that associate the overlapping cell barcodes with their cell type labels
for (n in c(seq.int(1,8))) {
  assign(paste0("rep",n,"_sc_overlap_cells_labeled"),
         filter(get(paste0("rep",n,"_sc_labels")),
                get(paste0("rep",n,"_sc_labels"))[,1] %in% get(paste0("rep",n,"_sc_overlap_cells"))))
}

# Removing the sc_labels from R's memory to save space
rm(rep1_sc_labels, rep2_sc_labels, rep3_sc_labels, rep4_sc_labels, rep5_sc_labels, rep6_sc_labels, rep7_sc_labels, rep8_sc_labels)
gc()

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
    stop("Warning! Row names (genes) are not identical between reps.")
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

##############################
### 3) Read the adipose tissue snRNAseq data and combine the files together
##############################

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
                        "GSM5359337_Hs_SAT_256-1.dge.tsv")

# Here is a function to read these large tsv files as sparse matrices (dgCMatrix)
## SOURCE: "readSparseCounts" package by Aaron Lun (https://rdrr.io/bioc/scuttle/man/readSparseCounts.html)
# Package was not available for direct installation (not compatible with R version)
readSparseCounts <- function(file, sep="\t", quote=NULL, comment.char="", row.names=TRUE, col.names=TRUE,
                             ignore.row=0L, skip.row=0L, ignore.col=0L, skip.col=0L, chunk=1000L) {
  if (is.character(file)) {
    fhandle <- file(file, open='r')
    on.exit(close(fhandle))
  } else if (is(file, "connection")) {
    fhandle <- file
    if (!isOpen(fhandle, "read")) {
      stop("'file' should be a connection in read mode")
    }
  } else {
    stop("'file' should be a connection or a character string")
  }

  # Scanning through rows.
  if (ignore.row) {
    readLines(fhandle, n=ignore.row)
  }
  if (col.names) {
    cell.names <- scan(fhandle, sep=sep, nlines=1, what=character(), quote=quote, comment.char=comment.char, quiet=TRUE)
  } else {
    cell.names <- NULL
  }
  if (skip.row) {
    readLines(fhandle, n=skip.row)
  }

  # Figuring out how to extract the columns.
  first <- scan(fhandle, sep=sep, quote=quote, what=character(), comment.char=comment.char, nlines=1L, quiet=TRUE)

  nentries <- length(first)
  what <- vector("list", nentries)
  if (row.names) {
    row.name.col <- ignore.col + 1L
    what[[row.name.col]] <- "character"
    ignore.col <- row.name.col
  }

  skip.col <- skip.col + ignore.col
  ncells <- nentries - skip.col
  cell.cols <- skip.col + seq_len(ncells)
  what[cell.cols] <- 0

  # Processing the first element.
  gene.names <- NULL
  if (row.names) {
    gene.names <- first[[row.name.col]]
  }

  output <- list(as(rbind(as.double(first[cell.cols])), "dgCMatrix"))
  it <- 2L

  # Reading it in, chunk by chunk (see behavior of nmax= when what= is a list).
  repeat {
    current <- scan(fhandle, what=what, sep=sep, quote=quote, comment.char=comment.char, nmax=chunk, quiet=TRUE)
    if (row.names) {
      gene.names <- c(gene.names, current[[row.name.col]])
    }
    output[[it]] <- as(do.call(cbind, current[cell.cols]), "dgCMatrix")
    it <- it + 1L
    if (chunk<0 || length(current[[1]]) < chunk) {
      break
    }
  }
  output <- do.call(rbind, output)

  # Adding row and column names, if available.
  if (row.names) {
    rownames(output) <- gene.names
  }
  if (col.names) {
    colnames(output) <- tail(cell.names, ncells)
  }
  return(output)
}

# Use readSparseCounts to read the files
for (file in adipose_file_names) {
  print(paste0("reading ",file))
  assign(paste0("adipose",(which(adipose_file_names == file))),
         readSparseCounts(file.path(reference_data_dir, file), sep="\t", quote=NULL, comment.char="", row.names=TRUE,
                          col.names=TRUE, ignore.row=0L, skip.row=0L, ignore.col=0L, skip.col=0L, chunk=1000L))
}

# Let's find the row names (genes) that are shared in common between files so that the
# adipose matrices can be accurately bound together horizontally.
adipose_genes <- list(rownames(adipose1), rownames(adipose2), rownames(adipose3), rownames(adipose4), rownames(adipose5),
                      rownames(adipose6), rownames(adipose7), rownames(adipose8), rownames(adipose9),
                      rownames(adipose10), rownames(adipose11), rownames(adipose12), rownames(adipose13))

adipose_common_genes <- Reduce(intersect, adipose_genes)

# Now we will find the row names (genes) in common between the the single cell and the adipose data.
all_overlap_genes <- intersect(rep1_sc_features[,2], adipose_common_genes)

# Now we will subset each adipose matrix to only include common genes
for (n in c(seq.int(1,13))) {
  assign(paste0("adipose",n,"_subset"),
         get(paste0("adipose",n))[all_overlap_genes, ])
}

# Removing the large matrices from R's memory to save space
rm(adipose1, adipose2, adipose3, adipose4, adipose5, adipose6, adipose7, adipose8, adipose9, adipose10, adipose11, adipose12, adipose13)
gc()

# Let's ensure that the row names (genes), are now identical between adipose
# samples so that the matrices can be accurately bound together horizontally.
for (n in c(seq.int(2,13))) {
  if(!identical(rownames(get(paste0("adipose",n,"_subset"))),rownames(adipose1_subset))){
    stop("Warning! Row names (genes) are not identical between adipose samples.")
  }
}

# Finally we can combine these into one large adipose matrix
all_adipose_expr <- cbind(adipose1_subset, adipose2_subset, adipose3_subset, adipose4_subset, adipose5_subset,
                          adipose6_subset, adipose7_subset, adipose8_subset, adipose9_subset,
                          adipose10_subset, adipose11_subset, adipose12_subset, adipose13_subset)

# Removing the subsetted matrices from R's memory to save space
rm(adipose1_subset, adipose2_subset, adipose3_subset, adipose4_subset, adipose5_subset,
   adipose6_subset, adipose7_subset, adipose8_subset, adipose9_subset,
   adipose10_subset, adipose11_subset, adipose12_subset, adipose13_subset)
gc()

##############################
### 4) Combine the scRNAseq data and the adipose snRNAseq data
##############################

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

# Ensuring that the cell barcodes match up between all_expr_final and all_cell_type_final
if(!identical(colnames(all_expr_final),all_cell_type_final[,1])){
  stop("Warning! Cell barcodes do not match up between all_expr_final and all_cell_type_final.")
}

# Replacing the cell barcodes with a unique ID number
all_cell_type_final <- all_cell_type_final %>%
  mutate(cell_id = 1:nrow(all_cell_type_final))

colnames(all_expr_final) <- 1:ncol(all_expr_final)

# Writing final files
writeMM(all_expr_final, file.path(reference_data_dir, "000combined_reference_expr_final.csv"), row.names=TRUE, col.names=TRUE)
write.csv(colnames(all_expr_final), file.path(reference_data_dir, "000combined_reference_expr_final_colnames.csv"))
write.csv(rownames(all_expr_final), file.path(reference_data_dir, "000combined_reference_expr_final_rownames.csv"))
write.csv(all_cell_type_final, file.path(reference_data_dir, "000combined_reference_cell_types_final.csv"), row.names=FALSE)


