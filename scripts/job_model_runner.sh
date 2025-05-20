#!/bin/bash
# filepath: /mnt/users/reich/rockerprojects/bestageing2022/scripts/job_model_runner.sh
#SBATCH --job-name="BestAgeing miRNA Analysis"
#SBATCH --output=/mnt/users/reich/rockerprojects/bestageing2022/logs/job_model_%j.txt
#SBATCH --error=/mnt/users/reich/rockerprojects/bestageing2022/logs/err/job_model_err_%j.txt
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --mail-user=q050cr@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -w benjamin

# Usage: sbatch job_model_runner.sh [options]
#
# Options:
#   --script SCRIPT_PATH      R script to run (required)
#   --disease DISEASE         Disease to analyze (default: "MI")
#   --model MODEL_TYPE        Model type: glmnet, rf, xgb (default: "glmnet")
#   --mirna SELECTION         miRNA selection: all, significant, literature, random (default: "all")
#   --matched                 Use matched samples (flag)
#   --cv-repeats N            Number of CV repeats (default: 10)
#   --seed N                  Random seed (default: 42)
#   --features N              Number of features to use (default: 20)
#
# Example:
#   sbatch job_model_runner.sh --script scripts/02_machine_learning_models.R --disease MI --model rf --mirna significant --matched --cv-repeats 10

# Default values
SCRIPT_PATH=""
DISEASE="CAD"
MODEL_TYPE="glmnet"
MIRNA_SELECTION="all"
MATCHED=false
CV_REPEATS=10
SEED=42
NUM_FEATURES=20

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --script)
      SCRIPT_PATH="$2"
      shift
      shift
      ;;
    --disease)
      DISEASE="$2"
      shift
      shift
      ;;
    --model)
      MODEL_TYPE="$2"
      shift
      shift
      ;;
    --mirna)
      MIRNA_SELECTION="$2"
      shift
      shift
      ;;
    --matched)
      MATCHED=true
      shift
      ;;
    --cv-repeats)
      CV_REPEATS="$2"
      shift
      shift
      ;;
    --seed)
      SEED="$2"
      shift
      shift
      ;;
    --features)
      NUM_FEATURES="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

# Check if script path is provided
if [ -z "$SCRIPT_PATH" ]; then
  echo "Error: Script path is required"
  echo "Usage: sbatch job_model_runner.sh --script SCRIPT_PATH [options]"
  exit 1
fi

# Check if script exists
if [ ! -f "$SCRIPT_PATH" ]; then
  echo "Error: Script not found: $SCRIPT_PATH"
  exit 1
fi

# Create unique job ID for logging
JOB_ID=${SLURM_JOB_ID:-$(date +%Y%m%d%H%M%S)}
LOG_FILE="/mnt/users/reich/rockerprojects/bestageing2022/logs/job_model_${JOB_ID}.log"

# Print job information
echo "========== Job Information ==========" | tee -a "$LOG_FILE"
echo "Job ID: $JOB_ID" | tee -a "$LOG_FILE"
echo "Script: $SCRIPT_PATH" | tee -a "$LOG_FILE"
echo "Disease: $DISEASE" | tee -a "$LOG_FILE"
echo "Model Type: $MODEL_TYPE" | tee -a "$LOG_FILE"
echo "miRNA Selection: $MIRNA_SELECTION" | tee -a "$LOG_FILE"
echo "Matched Samples: $MATCHED" | tee -a "$LOG_FILE"
echo "CV Repeats: $CV_REPEATS" | tee -a "$LOG_FILE"
echo "Seed: $SEED" | tee -a "$LOG_FILE"
echo "Number of Features: $NUM_FEATURES" | tee -a "$LOG_FILE"
echo "======================================" | tee -a "$LOG_FILE"

# Set up R environment
module load r/4.3.0

# Convert matched boolean to R logical
if [ "$MATCHED" = true ]; then
  MATCHED_R="TRUE"
else
  MATCHED_R="FALSE"
fi

# Run R script with parameters
echo "Running R script with parameters..." | tee -a "$LOG_FILE"
Rscript -e "
  # Set parameters
  DISEASE <- '$DISEASE'
  MODEL_TYPE <- '$MODEL_TYPE'
  MIRNA_SELECTION <- '$MIRNA_SELECTION'
  MATCHED <- $MATCHED_R
  CV_REPEATS <- $CV_REPEATS
  SEED <- $SEED
  NUM_FEATURES <- $NUM_FEATURES
  
  # Source the script
  source('$SCRIPT_PATH')
" 2>&1 | tee -a "$LOG_FILE"

# Check exit status
STATUS=$?
if [ $STATUS -eq 0 ]; then
  echo "Job completed successfully" | tee -a "$LOG_FILE"
  mkdir -p /mnt/users/reich/rockerprojects/bestageing2022/logs/success
  cp "$LOG_FILE" "/mnt/users/reich/rockerprojects/bestageing2022/logs/success/$(basename "$LOG_FILE")"
else
  echo "Job failed with status: $STATUS" | tee -a "$LOG_FILE"
  mkdir -p /mnt/users/reich/rockerprojects/bestageing2022/logs/err
  cp "$LOG_FILE" "/mnt/users/reich/rockerprojects/bestageing2022/logs/err/$(basename "$LOG_FILE")"
fi

exit $STATUS
