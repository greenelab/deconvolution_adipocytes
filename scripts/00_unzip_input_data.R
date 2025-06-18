##########################################################################################
### 00_unzip_input_data.R
### 
### This script unzips any gzip and zip files that may exist in the input_data directory.
##########################################################################################

# Load here library
library(here)

# Ensure the R.utils package is installed for gunzip functionality
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}
library(R.utils)

# Set base directory for location of .zip and .gz files
# We will also be saving the unzipped files here 
input_data <- file.path(here(), "input_data")

##########################################################
# 1) Function to check for and unzip .zip and .gz files in directory
##########################################################

unzip_compressed_files <- function(directory_path) {
  
  # Check if the provided directory exists
  if (!dir.exists(directory_path)) {
    stop("Error: The specified directory does not exist: ", directory_path)
  }
  
  message(paste("Checking directory:", directory_path, "for compressed files..."))
  
  # List all files in the directory
  all_files <- list.files(directory_path, full.names = TRUE)
  
  # Initialize counters for unzipped files
  zip_unzipped_count <- 0
  gz_unzipped_count <- 0
  
  # --- Handle .zip files ---
  zip_files <- all_files[grep("\\.zip$", all_files, ignore.case = TRUE)]
  
  if (length(zip_files) > 0) {
    message(paste("Found", length(zip_files), ".zip file(s). Attempting to unzip..."))
    for (zip_file in zip_files) {
      tryCatch({
        # Unzip the file into the same directory
        unzip(zip_file, exdir = directory_path)
        message(paste("Successfully unzipped:", basename(zip_file)))
        zip_unzipped_count <- zip_unzipped_count + 1
      }, error = function(e) {
        warning(paste("Failed to unzip", basename(zip_file), ":", e$message))
      })
    }
  } else {
    message("No .zip files found.")
  }
  
  # --- Handle .gz files ---
  gz_files <- all_files[grep("\\.gz$", all_files, ignore.case = TRUE)]
  
  if (length(gz_files) > 0) {
    message(paste("Found", length(gz_files), ".gz file(s). Attempting to gunzip..."))
    for (gz_file in gz_files) {
      tryCatch({
        # gunzip the file. The output file name will be the original name without .gz
        # e.g., "data.txt.gz" will become "data.txt"
        R.utils::gunzip(gz_file, remove = TRUE, overwrite = TRUE)
        message(paste("Successfully gunzipped:", basename(gz_file)))
        gz_unzipped_count <- gz_unzipped_count + 1
      }, error = function(e) {
        warning(paste("Failed to gunzip", basename(gz_file), ":", e$message))
      })
    }
  } else {
    message("No .gz files found.")
  }
  
  message(paste("\nSummary of unzipping operations:"))
  message(paste("- .zip files unzipped:", zip_unzipped_count))
  message(paste("- .gz files unzipped:", gz_unzipped_count))
  
  if (zip_unzipped_count == 0 && gz_unzipped_count == 0) {
    message("No compressed files were found or unzipped.")
  } else {
    message("Unzipping process complete.")
  }
}

##########################################################
# 2) Run function
##########################################################

unzip_compressed_files(input_data)