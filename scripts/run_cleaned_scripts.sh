#!/bin/bash
# filepath: /mnt/users/reich/rockerprojects/bestageing2022/scripts/run_cleaned_scripts.sh
# This script runs the cleaned versions of BestAgeing2022 analysis scripts
# Author: Christoph Reich
# Date: 2024-10-18

# Usage information
show_usage() {
  echo "Usage: $0 [options]"
  echo
  echo "Options:"
  echo "  -h, --help                 Display this help message"
  echo "  -s, --script SCRIPT        Specify script to run (de, ml, survival)"
  echo "  -d, --disease DISEASE      Disease to analyze (dcm, acs, cad, hfref)"
  echo "  -m, --mirna SELECTION      miRNA selection method (all, significant, literature, random)"
  echo "  --matched                  Use matched samples"
  echo "  -r, --repeats N            Number of CV repeats for ML models (default: 10)"
  echo "  --model MODEL              Model type for ML (glmnet, rf, xgb, etc.)"
  echo "  --seed N                   Random seed (default: 42)"
  echo
  echo "Examples:"
  echo "  $0 -s de -d dcm            # Run differential expression analysis for DCM"
  echo "  $0 -s ml -d cad --matched  # Run ML model for CAD with matched samples"
  echo "  $0 -s survival -d hfref    # Run survival analysis for HFREF"
}

# Default values
SCRIPT=""
DISEASE=""
SELECTION="all"
MATCHED=false
REPEATS=10
MODEL="glmnet"
SEED=42

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -h|--help)
      show_usage
      exit 0
      ;;
    -s|--script)
      SCRIPT="$2"
      shift
      shift
      ;;
    -d|--disease)
      DISEASE="$2"
      shift
      shift
      ;;
    -m|--mirna)
      SELECTION="$2"
      shift
      shift
      ;;
    --matched)
      MATCHED=true
      shift
      ;;
    -r|--repeats)
      REPEATS="$2"
      shift
      shift
      ;;
    --model)
      MODEL="$2"
      shift
      shift
      ;;
    --seed)
      SEED="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown parameter: $1"
      show_usage
      exit 1
      ;;
  esac
done

# Check required parameters
if [ -z "$SCRIPT" ]; then
  echo "Error: Script type is required"
  show_usage
  exit 1
fi

if [ -z "$DISEASE" ]; then
  echo "Error: Disease is required"
  show_usage
  exit 1
fi

# Set up R environment
module load r/4.3.0

# Create log directory if it doesn't exist
LOG_DIR="/mnt/users/reich/rockerprojects/bestageing2022/logs"
mkdir -p "$LOG_DIR"
mkdir -p "$LOG_DIR/success"
mkdir -p "$LOG_DIR/err"

# Generate timestamp for logging
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/run_${SCRIPT}_${DISEASE}_${TIMESTAMP}.log"

# Print job information
echo "========== Job Information ==========" | tee -a "$LOG_FILE"
echo "Script: $SCRIPT" | tee -a "$LOG_FILE"
echo "Disease: $DISEASE" | tee -a "$LOG_FILE"
echo "miRNA Selection: $SELECTION" | tee -a "$LOG_FILE"
echo "Matched Samples: $MATCHED" | tee -a "$LOG_FILE"
if [ "$SCRIPT" = "ml" ]; then
  echo "Model Type: $MODEL" | tee -a "$LOG_FILE"
  echo "CV Repeats: $REPEATS" | tee -a "$LOG_FILE"
fi
echo "Seed: $SEED" | tee -a "$LOG_FILE"
echo "======================================" | tee -a "$LOG_FILE"

# Convert matched boolean to R logical
if [ "$MATCHED" = true ]; then
  MATCHED_R="TRUE"
else
  MATCHED_R="FALSE"
fi

# Set script path based on script type
case $SCRIPT in
  "de")
    R_SCRIPT="/mnt/users/reich/rockerprojects/bestageing2022/scripts/001c_model_de_analysis_DET_MATRIX_cleaned.R"
    ;;
  "ml")
    R_SCRIPT="/mnt/users/reich/rockerprojects/bestageing2022/scripts/003c_ml_model_matchIt_cleaned.R"
    ;;
  "survival")
    R_SCRIPT="/mnt/users/reich/rockerprojects/bestageing2022/scripts/03_survival_analysis.R"
    ;;
  *)
    echo "Error: Invalid script type. Must be 'de', 'ml', or 'survival'" | tee -a "$LOG_FILE"
    exit 1
    ;;
esac

# Check if script exists
if [ ! -f "$R_SCRIPT" ]; then
  echo "Error: Script not found: $R_SCRIPT" | tee -a "$LOG_FILE"
  exit 1
fi

# Run R script with parameters
echo "Running R script: $R_SCRIPT" | tee -a "$LOG_FILE"
echo "Command started at: $(date)" | tee -a "$LOG_FILE"

# Set up R command based on script type
if [ "$SCRIPT" = "ml" ]; then
  # ML script takes different parameters
  Rscript -e "
    # Set up disease and analysis
    all_combis <- data.frame(
      diseases = c('$DISEASE'),
      analysis = c('selected')
    )
    
    # Command line arguments for script
    args <- c($MATCHED_R, FALSE, $REPEATS)
    
    # Source the script
    source('$R_SCRIPT')
  " 2>&1 | tee -a "$LOG_FILE"
else
  # DE or survival scripts
  Rscript -e "
    # Set parameters
    DISEASE <- '$DISEASE'
    MIRNA_SELECTION <- '$SELECTION'
    MATCHED <- $MATCHED_R
    SEED <- $SEED
    
    # For ML models
    if ('$SCRIPT' == 'ml') {
      MODEL_TYPE <- '$MODEL'
      CV_REPEATS <- $REPEATS
    }
    
    # Source the script
    source('$R_SCRIPT')
  " 2>&1 | tee -a "$LOG_FILE"
fi

# Check exit status
STATUS=$?
if [ $STATUS -eq 0 ]; then
  echo "Job completed successfully" | tee -a "$LOG_FILE"
  cp "$LOG_FILE" "$LOG_DIR/success/$(basename "$LOG_FILE")"
else
  echo "Job failed with status: $STATUS" | tee -a "$LOG_FILE"
  cp "$LOG_FILE" "$LOG_DIR/err/$(basename "$LOG_FILE")"
fi

echo "Command finished at: $(date)" | tee -a "$LOG_FILE"

exit $STATUS
