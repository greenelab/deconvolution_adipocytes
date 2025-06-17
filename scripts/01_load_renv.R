##########################################################################################
### 01_load_renv.R
### 
### This script loads the renv from the lockfile.
##########################################################################################


# Load the renv library from the lockfile
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", version = "1.1.0")  # Install renv if not already installed
}
renv::restore()