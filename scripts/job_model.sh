#!/bin/bash
#SBATCH --job-name="miRNA diagnostic analysis"
#SBATCH --output=/mnt/users/reich/rockerprojects/bestageing2022/logs/job_model_%j.txt
#SBATCH --error=/mnt/users/reich/rockerprojects/bestageing2022/logs/job_model_%j.txt
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --mail-user=q050cr@gmail.com
#SBATCH --mail-type=END,FAIL

# This script runs an R script with 3 arguments. The arguments are passed 
# to the R script and also incorporated into the output filename.

# Usage: sbatch job_model.sh arg1 arg2 arg3
# ex: sbatch job_model.sh TRUE TRUE 10
# arg1: use miRetrieveBiomarker; pass (TRUE, FALSE)
# arg2: use random_selection of miRNAs; pass (TRUE, FALSE) 
# arg3: no of repeats of cross-val, pass integer

echo "== Diagnostic analysis of research miRNAs  =="
echo "Job started at $(date)"
echo "Using miRetrieve Biomarkers selection: $1"
echo "Using random selection: $2"
echo "Using 5-fold CV with repeats: $3"

SINGULARITY_IMAGE="/mnt/users/reich/programs/christophs_custom_rocker_24.sif"
R_SCRIPT_PATH="/mnt/users/reich/rockerprojects/bestageing2022/scripts/model.R"

# Set Singularity to use a different temporary directory
export SINGULARITY_TMPDIR="/mnt/users/reich/tmp"
export R_LIBS="/mnt/users/reich/programs/R43/lib:$R_LIBS"

singularity exec ${SINGULARITY_IMAGE} Rscript ${R_SCRIPT_PATH} $1 $2 $3

echo "Job ended at $(date)"
echo "== End of Job =="
