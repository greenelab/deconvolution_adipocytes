#!/bin/bash
#SBATCH --job-name=hgsoc_nf
#SBATCH --account=amc-general
#SBATCH --output=log_hgsoc_nf.log
#SBATCH --error=errorhgsoc_nf.err
#SBATCH --time=08:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1 

set -euo pipefail

# --------------------------------------------------
# 0) Resolve project root (directory of this script)
# --------------------------------------------------
PRJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "${PRJ_DIR}"

# --------------------------------------------------
# 1) Load Conda and activate env_hgsoc
# --------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV_YML="${PRJ_DIR}/environments/env_hgsoc.yml"

# create the env once; reuse afterwards
if ! conda env list | grep -q '^env_hgsoc '; then
    echo "••• Creating Conda environment env_hgsoc"
    conda env create -f "${ENV_YML}"
fi
conda activate env_hgsoc          # puts R, Nextflow, compilers on PATH

# --------------------------------------------------
# 2) Install R packages that Conda cannot provide
# --------------------------------------------------
Rscript - <<'RSCRIPT'
pkgs_cran  <- c("InstaPrism")                    # GitHub-only
pkgs_bioc  <- c()                               # add names if Bioc, but absent on bioconda

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")

for (p in pkgs_bioc)
  if (!requireNamespace(p, quietly = TRUE))
      BiocManager::install(p, ask = FALSE, update = FALSE)

# example GitHub install (replace with the repo that actually hosts InstaPrism)
if (!requireNamespace("InstaPrism", quietly = TRUE))
    remotes::install_github("JinmiaoChenLab/InstaPrism")
RSCRIPT

# --------------------------------------------------
# 3) Run the pipeline
# --------------------------------------------------
echo "••• Launching Nextflow"
nextflow run main.nf -profile slurm -resume

echo "••• Pipeline finished 🎉"
