
###############
# The following script uses InsaPrism for deconvolution of the bulks in BayesPrism framework
# using the following tutorial for reference (https://humengying0907.github.io/InstaPrism_tutorial.html)
##############

# These are large files that take up a lot of memory - clear R's memory to make room
rm(list = ls())
gc()

# Install InstaPrism 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("Biobase", quietly = TRUE))
  BiocManager::install("Biobase")
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("humengying0907/InstaPrism")

# Load required libraries
library(InstaPrism)
library(data.table)
library(Matrix)

# Set seed for reproducibility
set.seed(13)

# Define import and export directories
import_dir <- file.path("/projects/gakatsu@xsede.org/hgsoc_adipocyte_deconvolution")
export_dir <- file.path("/projects/gakatsu@xsede.org/hgsoc_adipocyte_deconvolution/instaprism_output")

# Read files
scExpr_all <- Matrix::readMM(file.path(import_dir,"reference_data/000combined_reference_expr_final.csv"))
colnames(scExpr_all) <- (fread(file.path(import_dir,"reference_data/000combined_reference_expr_final_colnames.csv"), header = TRUE, data.table = FALSE))$"x"
rownames(scExpr_all) <- (fread(file.path(import_dir,"reference_data/000combined_reference_expr_final_rownames.csv"), header = TRUE, data.table = FALSE))$"x"

cell_type_labels <- fread(file.path(import_dir,"reference_data/000combined_reference_cell_types_final.csv"), data.table = FALSE)
cell_type_labels <- cell_type_labels[,"cellType"]

# Function to select a random subset of each cell type to decrease the size of the single cell/nucleus reference matrix
select_cell_subsets <- function(cell_type_labels, n){
  cell_types <- unique(cell_type_labels) # Getting a list of cell types present
  all_selected_cells <- c()
  for(cell in cell_types){
    all <- grep(cell, cell_type_labels) # Get vector indices for cell type
    if(length(all) < n){
      all_selected_cells <- c(all_selected_cells, all)
      print(paste0("Fewer than ",n," cells of type ",cell," so will not subset. Total cell count ",length(all)))
      next
    }
    selected_cells <- sample(all, n) # Select a random subset of these vector inidices
    all_selected_cells <- c(all_selected_cells, selected_cells)
  }
  return(all_selected_cells) # Indices for selected cells
}

# Using the function defined above to subset the genes x cells matrix and the cell type labels vector
selected_cells <- select_cell_subsets(cell_type_labels, n = 500)
scExpr_subset <- scExpr_all[, selected_cells]
cell_type_labels_subset <- cell_type_labels[selected_cells]

# Removing large objects from R's memory to save room 
rm(scExpr_all, cell_type_labels)
gc()

# Create reference object for input into InstaPrism
refPhi_obj = refPrepare(sc_Expr = scExpr_subset, cell.type.labels = cell_type_labels_subset, cell.state.labels = cell_type_labels_subset)

# List of datasets
dataset_list <- c("SchildkrautB", "SchildkrautW", "TCGA", "Mayo", "Tothill", "Yoshihara")

# Run InstaPrism for each dataset
for (ds in dataset_list) {
  print(paste0("Running InstaPrism on ", ds))
  # InstaPrism requires non-log-transformed bulk data, so we will read the pre-transformation Schildkraut data (before doing log10(...+1))
  # The other datasets are microarray and we must read the post-transformation (pseudo "counts") version for those (after doing 2^(...) or 10^(...))
  if (ds == "SchildkrautB" | ds == "SchildkrautW") {
    bulk_expr <- fread(file = (file.path(paste0(import_dir, "/data/bulks/", ds, "_filtered_asImported.csv"))), header = TRUE, data.table = FALSE)
  } else {
    bulk_expr <- fread(file = (file.path(paste0(import_dir, "/data/bulks/", ds, "_filtered_transformed.csv"))), header = TRUE, data.table = FALSE)
  }
  bulk_expr_rownames <- data.frame(bulk_expr[,-1], row.names=bulk_expr[,1]) # Setting column 1 (sample IDs) as the rownames
  assign(paste0("bulk_expr_", ds), t(bulk_expr_rownames)) # Transposing the rows and columns
  assign(paste0("instaprism_output_", ds),
         InstaPrism(bulk_Expr = get(paste0("bulk_expr_", ds)),
                    refPhi_cs = refPhi_obj))
}

# Writing final files
for (ds in dataset_list) {
  write.csv(t((get(paste0("instaprism_output_", ds))@Post.ini.ct@theta)),
            file.path(export_dir, paste0("instaprism_output_", ds)),
            row.names=TRUE)
}

