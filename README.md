# The European BestAgeing Study on microRNA Biomarkers

**Background:** In the present work, we present the outcome of a multinational, prospective biomarker validation study of the BestAgeing consortium. We systematically evaluated the diagnostic potential of miRNA signatures using ML approaches for a range of cardiovascular diseases (CVD), including acute coronary syndrome (ACS), chronic coronary artery disease (CAD), dilated cardiomyopathy (DCM), and ischemic cardiomyopathy (ICM) in a uniform, standardized fashion. We also aimed to investigate the impact on patient survival and disease severity.

**Conclusion:** The European BestAgeing miRNA study reveals the potential of several miRNA biomarkers for disease diagnosis and prognostication. Especially signatures for CAD and heart failure detection might improve clinical decision making for a larger patient group.

## Project Structure

This repository contains code and data for the BestAgeing2022 project, which analyzes miRNA expression patterns in relation to aging and disease.

### Directory Structure

```
bestageing2022/
├── data/               # Raw and processed data
├── data-literature/    # Literature references and resources
├── figures/            # Generated figures
├── logs/               # Log files from job executions
│   ├── err/            # Error logs
│   ├── success/        # Success logs
├── output/             # Analysis outputs
│   ├── de_results/     # Differential expression results
│   ├── models/         # Trained models
│   ├── plots/          # Generated plots
│   └── tables/         # Generated tables
├── scripts/            # Analysis scripts
│   ├── config/         # Configuration files
│   ├── helper/         # Helper functions
│   ├── figures_create/ # Figure generation scripts
│   └── tables_create/  # Table generation scripts
└── scripts_2024/       # Scripts for 2024 analyses
```

## Script Organization

The scripts have been reorganized for clarity and maintainability:

1. **Configuration Files** (`scripts/config/`): Centralized configuration for paths and settings.
2. **Helper Functions** (`scripts/helper/`): Reusable functions for data processing and visualization.
3. **Analysis Scripts** (`scripts/`): Main analysis scripts, numbered by workflow order.
4. **Job Scripts** (`scripts/job_model_runner.sh`): Scripts for running analyses on the compute cluster.

## Main Analysis Scripts

- `01_differential_expression.R`: Differential expression analysis of miRNAs.
- `02_machine_learning_models.R`: Machine learning model training and evaluation.
- `03_survival_analysis.R`: Survival analysis using miRNA expression data.

## How to Run

### Interactive Mode

To run the analysis scripts interactively:

```R
# In R
source("scripts/01_differential_expression.R")
```

You can modify the parameters at the top of each script to customize the analysis.

### Cluster Mode

To run the analysis on the compute cluster:

```bash
# Submit a differential expression analysis job
sbatch scripts/job_model_runner.sh --script scripts/01_differential_expression.R --disease MI --matched

# Submit a machine learning model job
sbatch scripts/job_model_runner.sh --script scripts/02_machine_learning_models.R --disease MI --model rf --mirna significant --matched --cv-repeats 10
```

## Dependencies

This project requires the following R packages:

- tidyverse
- limma
- edgeR
- caret
- randomForest
- glmnet
- xgboost
- pROC
- mlr3
- mlr3verse
- survival
- survminer
