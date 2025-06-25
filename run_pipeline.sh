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

set -eo pipefail

# --------------------------------------------------
# 0) Resolve project root (directory of this script)
# --------------------------------------------------
# Check if SLURM_SUBMIT_DIR is set (indicates Slurm environment)
if [ -n "${SLURM_SUBMIT_DIR}" ]; then
    # Running on HPC via Slurm
    PRJ_DIR="${SLURM_SUBMIT_DIR}"
    RUN_MODE="HPC"
    echo "Running on HPC (Slurm). Project directory: ${PRJ_DIR}"
else
    # Running locally
    # Get the directory where the script itself is located when run locally
    PRJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    RUN_MODE="LOCAL"
    echo "Running locally. Project directory: ${PRJ_DIR}"
fi

# Change to the project directory
cd "${PRJ_DIR}"
echo "Current working directory: $(pwd)"

# --------------------------------------------------
# 0.5) Specifying paths specific to Alpine
# --------------------------------------------------
# export paths specific to Alpine
if [ "${RUN_MODE}" == "HPC" ]; then
    export PATH=/usr/include:$PATH
    export CPATH=/usr/include/:$CPATH
    export C_INCLUDE_PATH=/usr/include/:$C_INCLUDE_PATH
fi

# --------------------------------------------------
# 1) Load miniforge, activate Conda environment env_hgsoc
# --------------------------------------------------
module load miniforge

ENV_YML="${PRJ_DIR}/env_hgsoc.yml"

# create the env once; reuse afterwards
 if ! conda env list | grep -q '^env_hgsoc '; then
    echo "••• Creating Conda environment env_hgsoc"
    conda env create -f "${ENV_YML}"
 fi

echo "••• Activating Conda environment env_hgsoc"
conda activate env_hgsoc

# --------------------------------------------------
# 2) Ensure that Nextflow is installed
# --------------------------------------------------
if ! command -v nextflow &> /dev/null
then
    echo "Error: Nextflow is not installed or not found in your PATH."
    echo "Please install Nextflow and ensure it's accessible in your system's PATH."
    echo "You can typically install it with: curl -s https://get.nextflow.io | bash"
    echo "Then move the 'nextflow' executable to a directory in your PATH (e.g., ~/bin or /usr/local/bin)."
    exit 1 # Exit the script with an error code
fi

# --------------------------------------------------
# OPTIONAL: unzip .gz and .zip files in input_data
# --------------------------------------------------
# python scripts/00_unzip_input_data.py

# --------------------------------------------------
# 3) Run the pipeline
# --------------------------------------------------
echo "••• Launching Nextflow"

# Define Nextflow's work directory based on run mode
NEXTFLOW_WORK_DIR="${PRJ_DIR}/nextflow_work"


if [ "${RUN_MODE}" == "HPC" ]; then
    nextflow run main.nf -profile slurm -resume -w "${NEXTFLOW_WORK_DIR}"
else
    nextflow run main.nf -profile local -resume -w "${NEXTFLOW_WORK_DIR}"
fi

echo "••• Pipeline finished 🎉"