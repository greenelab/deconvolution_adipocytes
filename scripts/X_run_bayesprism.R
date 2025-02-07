# Set up CAN repository (needed in SLUR M)
options(repose = c(CRANE = "https://cloud.r-project.org"))

if(!requireNamespace("scran", quietly = TRUE)){
  BiocManager::install("scran")
}
library(scran)

# Load necessary libraries
library(InstaPrism)
library(BayesPrism)

# Get data type from SLURM argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No data type provided. Run script with: Rscript BayesPrism.R <DATA_TYPE>")
}
data_type <- args[1]

# Define the combinations of gene_type and reference_type
combinations <- expand.grid(reference_type = c("sc", "sn"), gene_type = c("filtered", "notfiltered"))

# Loop through each combination
for (i in 1:nrow(combinations)) {
  # Extract the current combination
  gene_type <- combinations$gene_type[i]
  reference_type <- combinations$reference_type[i]
  
  # Detect script directory
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    script_dir <- getwd()  # Default to current working directory if not in RStudio
  }
  
  # Define paths
  import_path <- file.path(script_dir, "..", "data/deconvolution", data_type, "/")
  export_path <- file.path(script_dir, "..", "results/", data_type, "/")
  
  # Define file names
  mixture_file <- paste0(import_path, "pseudobulks.csv")
  signal_file <- paste0(import_path, reference_type, "_", gene_type, "_signal.csv")
  cell_state_file <- paste0(import_path, reference_type, "_", gene_type, "_cell_state.csv")
  export_file_prop <- paste0(export_path, reference_type, "_", gene_type, "_BayesPrism_proportionsresults_new.csv")
  export_file_ref <- paste0(export_path, reference_type, "_", gene_type, "_BayesPrism_usedref_new.csv")
  
  # Print paths to debug
  cat("Import path:", import_path, "\n")
  cat("Export path:", export_path, "\n")
  cat("Mixture file:", mixture_file, "\n")
  cat("Signal file:", signal_file, "\n")
  cat("Cell state file:", cell_state_file, "\n")
  
  # Read mixture (pseudobulks)
  mixture_data <- read.csv(mixture_file, stringsAsFactors = FALSE, row.names = 1)
  cat("Mixture Data (before filtering):\n")
  print(head(mixture_data[, 1:3]))
  
  # Read signal (reference)
  signal_data <- read.csv(signal_file, stringsAsFactors = FALSE, row.names = 1)
  cat("Signal Data:\n")
  print(head(signal_data[, 1:3]))
  
  # Read cell state
  cell_state <- read.csv(cell_state_file)
  cat("Cell State Data:\n")
  print(head(cell_state))
  
  # Ensure genes in pseudobulks match the reference
  common_genes <- intersect(rownames(mixture_data), rownames(signal_data))
  
  if (length(common_genes) < 500) {
    stop(paste("ERROR: Too few matching genes (", length(common_genes), ") between mixture and reference!"))
  }
  
  # Filter and reorder mixture_data to match reference genes
  mixture_data <- mixture_data[common_genes, ]
  signal_data <- signal_data[common_genes, ]
  
  cat("Mixture Data (after filtering):\n")
  print(head(mixture_data[, 1:3]))
  
  ########################################################
  bulk_Expr <- mixture_data
  cell_type_labels <- t(cell_state[1])
  cell_state_labels <- t(cell_state[2])
  
  scExpr <- SingleCellExperiment(assays = list(counts = as.matrix(signal_data)))
  
  refPhi_obj = refPrepare(
    sc_Expr = assay(scExpr, "counts"),
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels
  )
  results = InstaPrism(bulk_Expr = bulk_Expr, refPhi_cs = refPhi_obj)
  
  # Export results
  cell_frac <- results@Post.ini.cs@theta
  write.table(cell_frac, export_file_prop, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  cell_ref <- results@initial.reference@phi.cs
  write.table(cell_ref, export_file_ref, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # Log progress
  cat(paste("Processing completed for genes = ", gene_type, " and reference type = ", reference_type, "\n"))
}

