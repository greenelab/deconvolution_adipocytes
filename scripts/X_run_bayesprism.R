#The following script opens files created in Python (deconvolution_differences/scripts/prepare_deconvolution.py)
# and uses InsaPrism for deconvolution of the bulks in BayesPrism framework, using references with SC and SN, wihth DEG filtered and not.

# Import necessary libraries
library(InstaPrism)
library(BayesPrism)

data_type = "ADP"

# Define the combinations of gene_type and reference_type
combinations <- expand.grid(reference_type = c("sc", "sn"), gene_type = c("filtered", "notfiltered"))

# Loop through each combination
for (i in 1:nrow(combinations)) {
  # Extract the current combination
  gene_type <- combinations$gene_type[i]
  reference_type <- combinations$reference_type[i]
  
  ###################
  #check if running in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    # Get the directory where the R script resides in RStudio
    script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    # If not in RStudio, use a fallback method (e.g., set it manually)
    script_dir <- ""
  }
  
  # Define the relative path based on the current combination
  import_path <- file.path(script_dir, "..", "data/deconvolution/", data_type, "/")
  
  # Combine the script directory and relative path to get the full path
  export_path <- paste0(relative_path, "results/", data_type, "/")
  
  # Define fil names based on the current combination
  mixture_file <- paste0(import_path, "pseudobulks.csv")
  signal_file <- paste0(import_path,  data_type, "_", reference_type,"_", gene_type, "_signal.csv")
  cell_state_file <- paste0(import_path ,  data_type, "_", reference_type,"_", gene_type, "_cell_state.csv")
  export_file_prop <- paste0(export_path,  data_type, "_", reference_type,"_", gene_type, "_BayesPrism_proportionsresults.csv")
  export_file_ref <- paste0(export_path, data_type, "_", reference_type,"_", gene_type, "_BayesPrism_usedref.csv")
  
  # read, inspect, mixture file as RDS
  mixture_file_rds <- mixture_file
  mixture_data <- read.csv(mixture_file, stringsAsFactors = FALSE, row.names = 1)
  cat("Mixture Data:\n")
  print(head(mixture_data[,1:3]))
  
  # read, inspect, signal files as RDS
  signal_files_rds <- signal_file
  signal_data <- read.csv(signal_file, stringsAsFactors = FALSE, row.names = 1)
  cat("Signal Data:\n")
  print(head(signal_data[,1:3]))
  
  # read, inspect, and save the cell state files as CSV
  cell_state <- read.csv(cell_state_file)
  cat("Cell State Data:\n")
  print(head(cell_state))
  
  ########################################################
  ###### Roughly following IntaPrism tutorial: https://humengying0907.github.io/InstaPrism_tutorial.html
  bulk_Expr <- mixture_data
  sc_Expr <- signal_data
  cell_type_labels <- t(cell_state[1])
  cell_state_labels <- t(cell_state[2])
  
  InstaPrism.res = InstaPrism(input_type = 'raw',sc_Expr = sc_Expr,bulk_Expr = bulk_Expr,
                              cell.type.labels = cell_type_labels,cell.state.labels = cell_state_labels)
  
  ##Export files: proportions estimated and references used
  cell_frac = InstaPrism.res@Post.ini.cs@theta
  write.table(cell_frac, export_file_prop, sep="\t", quote=F,  row.names = TRUE, col.names = TRUE,)
  cell_ref = InstaPrism.res@initial.reference@phi.cs
  write.table(cell_ref, export_file_ref, sep="\t", quote=F, row.names = TRUE, col.names = TRUE,)
  
  # Printing the value of num_missing and combination to keep track
  cat(paste("Processing with genes = ", gene_type, "and reference type = ", reference_type), "\n")
}

