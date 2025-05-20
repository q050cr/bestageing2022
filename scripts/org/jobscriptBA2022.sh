#!/bin/bash
#SBATCH --job-name="miRNA diagnostic analysis"
#SBATCH --output=R_job_output.txt
#SBATCH --error=R_job_error.txt
#SBATCH --mail-user=q050cr@gmail.com
#SBATCH --mail-type=END,FAIL

echo "== Diagnostic analysis of research miRNAs  =="

Rscript scripts/jobRunScript.R

echo "== End of Job =="