##########################################################################################
### 4_run_instaprism.R
### 
### This script runs InstaPrism, an R package that performs deconvolution with a similar
### but faster method to BayesPrism. First, it loads in the reference single cell plus
### single nucleus RNA sequencing reference dataset created in the previous script. Since
### the reference dataset is so large, it randomly selects 500 of each cell type to use.
### It then generates and saves two reference objects for input into InstaPrism, one with
### adipocytes and one without adipocytes. It runs InstaPrism twice on each of the six
### bulk datasets, both with and without the adipocytes in the reference data. Of note,
### Instaprism requires non-log-transformed bulk data, so it is performed on the original,
### non-transformed data for the bulk RNA sequencing datasets and on the 2^(…) transformed
### data for the microarray datasets.
##########################################################################################

# These are large files that take up a lot of memory - clear R's memory to make room
rm(list = ls())
gc()

# Load environment
renv::load()

# Load required libraries
library(InstaPrism)
library(data.table)
library(Matrix)
library(here)

# Set base directories
# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
output_data <- file.path(here(), "output_data")
instaprism_folder <- file.path(output_data, "instaprism")

# Ensure instaprism_folder directory exists
dir.create(instaprism_folder, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(5)

##########################################################
# 0) Genes to remove (for adipocyte reference only)
##########################################################
genes_to_remove <- fread(
  file.path(here(), "input_data", "intersect_3ds.csv"),
  header = FALSE,        # your file has no header row
  data.table = FALSE
)[, 2]                   # second column holds the symbols
genes_to_remove <- unique(na.omit(genes_to_remove))

##########################################################
# 1) Read sn+sc reference files
##########################################################

# Read expression matrix
scExpr_all <- Matrix::readMM(file.path(output_data,"sc_sn_reference_data/reference_expr_final.csv"))
rownames(scExpr_all) <- (fread(file.path(output_data,"sc_sn_reference_data/reference_expr_final_genes.csv"),
                               header = TRUE, data.table = FALSE))$"x"
colnames(scExpr_all) <- (1:ncol(scExpr_all))

# Read cell type label (key)
cell_type_labels <- fread(file.path(output_data,"sc_sn_reference_data/reference_cell_types_final.csv"),
                          data.table = FALSE)
cell_type_labels <- cell_type_labels[,"cellType"]

#  Remove the unwanted genes *only for the adipocyte run* intersection genes based on previous work
# scExpr_all <- scExpr_all[
#   !rownames(scExpr_all) %in% genes_to_remove, ]

##########################################################
# 2) Select 500 cells of each cell type from the large 
#    reference dataset
##########################################################

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
scExpr_subset_all <- scExpr_all[, selected_cells]
cell_type_labels_subset_all <- cell_type_labels[selected_cells]

# Removing large objects from R's memory to save room 
rm(scExpr_all, cell_type_labels)
gc()

# Now create a version with no adipocytes
adipos <- grep("Adipocytes", cell_type_labels_subset_all)
scExpr_subset_no_adipos <- scExpr_subset_all[,-adipos]
cell_type_labels_subset_no_adipos <- cell_type_labels_subset_all[-adipos]

##########################################################
# 3) Create and save reference objects for input
#    into InstaPrism
##########################################################

# Create reference objects for input into InstaPrism
refPhi_obj_all = refPrepare(sc_Expr = scExpr_subset_all,
                        cell.type.labels = cell_type_labels_subset_all,
                        cell.state.labels = cell_type_labels_subset_all)
refPhi_obj_no_adipos = refPrepare(sc_Expr = scExpr_subset_no_adipos,
                            cell.type.labels = cell_type_labels_subset_no_adipos,
                            cell.state.labels = cell_type_labels_subset_no_adipos)

# Write these reference objects as files
dir.create(file.path(instaprism_folder,"instaprism_reference_objects"), recursive = TRUE, showWarnings = FALSE)
write.csv(slot(refPhi_obj_all,"phi.cs"),
        file.path(instaprism_folder, "instaprism_reference_objects/instaprism_reference_with_adipos"))
write.csv(slot(refPhi_obj_no_adipos,"phi.cs"),
        file.path(instaprism_folder, "instaprism_reference_objects/instaprism_reference_no_adipos"))

##########################################################
# 4) Read bulk dataset files
##########################################################

# List of datasets
dataset_list <- c("SchildkrautB", "SchildkrautW", "TCGA_bulk", "TCGA_microarray", "Tothill", "Yoshihara")

# Reading the bulk datasets
# Instaprism requires non-log-transformed bulk data, so it is performed on the original, non-transformed data
# for the bulk RNA sequencing datasets and on the 2^(…) transformed data for the microarray datasets.
for (ds in dataset_list){
  if (ds == "SchildkrautB" | ds == "SchildkrautW" | ds == "TCGA_bulk"){
    bulk_expr <- fread(file.path(output_data, paste0("bulk_datasets/",ds,"_filtered_asImported.csv")),
          header = TRUE, data.table = FALSE)
  } else {
    bulk_expr <- fread(file.path(output_data, paste0("bulk_datasets/",ds,"_filtered_transformed.csv")),
                       header = TRUE, data.table = FALSE)
  }
  # Setting column 1 (sample IDs) as the rownames
  bulk_expr <- data.frame(bulk_expr[,-1], row.names=bulk_expr[,1]) 
  # Assigning object and transposing so that rows = genes and columns = samples
  assign(paste0("bulk_expr_",ds), t(bulk_expr)) 
}

##########################################################
# 5) Run InstaPrism twice on each bulk dataset, both with
#    and without adipocytes
##########################################################

# Run Instaprism with adipocytes
for (ds in dataset_list){
  print(paste0("Running InstaPrism (with adipocytes) on ", ds))
  instaprism_output <- InstaPrism(bulk_Expr = get(paste0("bulk_expr_", ds)),
                                  refPhi_cs = refPhi_obj_all, input_type = "refPhi_cs")
  assign(paste0("instaprism_output_", ds, "_with_adipocytes"),
         instaprism_output)
}

# Run Instaprism without adipocytes
for (ds in dataset_list){
  print(paste0("Running InstaPrism (no adipocytes) on ", ds))
  instaprism_output <- InstaPrism(bulk_Expr = get(paste0("bulk_expr_", ds)),
                                  refPhi_cs = refPhi_obj_no_adipos, input_type = "refPhi_cs")
  assign(paste0("instaprism_output_", ds, "_no_adipocytes"),
         instaprism_output)
}

# Creating directory to store Instaprism outputs
dir.create(file.path(instaprism_folder,"instaprism_outputs"), recursive = TRUE, showWarnings = FALSE)

# Writing final files (with adipocytes)
for (ds in dataset_list) {
  write.csv(t((get(paste0("instaprism_output_", ds, "_with_adipocytes"))@Post.ini.ct@theta)),
            file.path(instaprism_folder, paste0("instaprism_outputs/instaprism_output_", ds,"_with_adipocytes.csv")),
            row.names=TRUE)
}

# Writing final files (no adipocytes)
for (ds in dataset_list) {
  write.csv(t((get(paste0("instaprism_output_", ds, "_no_adipocytes"))@Post.ini.ct@theta)),
            file.path(instaprism_folder, paste0("instaprism_outputs/instaprism_output_", ds,"_no_adipocytes.csv")),
            row.names=TRUE)
}
