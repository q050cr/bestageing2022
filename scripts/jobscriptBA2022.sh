#!/bin/bash
#SBATCH --job-name="miRNA diagnostic analysis"


echo "== Diagnostic analysis of research miRNAs  =="

Rscript scripts/jobRunScript.R

echo "== End of Job =="