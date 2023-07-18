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


echo "== Diagnostic analysis of research miRNAs  =="
echo "Job started at $(date)"


SINGULARITY_IMAGE="/mnt/users/reich/programs/christophs_custom_rocker_24.sif"
R_SCRIPT_PATH="/mnt/users/reich/rockerprojects/bestageing2022/scripts/model.R"

export R_LIBS="/mnt/users/reich/programs/R43/lib:$R_LIBS"

singularity exec ${SINGULARITY_IMAGE} Rscript ${R_SCRIPT_PATH}

echo "Job ended at $(date)"
echo "== End of Job =="
