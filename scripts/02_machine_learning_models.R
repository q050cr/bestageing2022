# Machine Learning Model Training and Evaluation
#
# This script trains and evaluates ML models on miRNA data
# Author: Christoph Reich
# Date: 2025-05-20 (updated)

# ----- Setup -----
# Load configuration and utilities
source("scripts/config/config.R")
source("scripts/helper/utils.R")

# Get paths
paths <- get_project_paths()

# Load required libraries
load_project_libraries(c(
  "tidyverse",
  "caret",
  "randomForest",
  "glmnet",
  "xgboost",
  "pROC",
  "mlr3",
  "mlr3verse"
))

# ----- Parameters -----
# Analysis parameters - edit these as needed
DISEASE <- "CAD" # Choose from available diseases
MATCHED <- TRUE # Use matched samples?
CV_REPEATS <- 10 # Number of cross-validation repeats
SEED <- 42 # Random seed for reproducibility
MODEL_TYPE <- "glmnet" # Model type: "glmnet", "rf", "xgb", "ensemble"
MIRNA_SELECTION <- "all" # How to select miRNAs: "all", "significant", "literature", "random"
NUM_FEATURES <- 20 # Number of features to select (if not using "all")

# Set seed for reproducibility
set.seed(SEED)

# ----- Data Import -----
# Import data
clinical_data <- import_project_data("clinical")
mirna_data <- import_project_data("mirna")
metadata <- import_project_data("metadata")

# If using significant miRNAs from DE analysis, load those results
if (MIRNA_SELECTION == "significant") {
  de_results_path <- list.files(
    file.path(paths$output_path, "de_results"),
    pattern = paste0(".*de_results_", DISEASE, ".*\\.csv$"),
    full.names = TRUE
  )[1]

  if (!is.na(de_results_path)) {
    de_results <- read.csv(de_results_path)
  } else {
    stop("No DE results found. Run the differential expression analysis first.")
  }
}

# ----- Data Preparation -----
# Prepare data for analysis
prepare_ml_data <- function(clinical_data, mirna_data, metadata, disease,
                            mirna_selection = "all", num_features = 20,
                            de_results = NULL, matched = FALSE) {
  # Join clinical data with miRNA data
  merged_data <- inner_join(clinical_data, mirna_data, by = "patient_id")

  # Filter by disease
  disease_data <- merged_data %>%
    filter(Disease == disease | Disease == "Control") %>%
    mutate(Disease = factor(Disease, levels = c("Control", disease)))

  # Match samples if required
  if (matched) {
    # Create matched groups based on age, sex, etc.
    matched_formula <- as.formula(paste("Disease ~", paste(c("Age", "Sex", "BMI"), collapse = " + ")))
    matched_data <- matchIt::matchIt(
      matched_formula,
      data = disease_data,
      method = "nearest",
      ratio = 1
    )

    # Extract matched samples
    disease_data <- matchIt::match.data(matched_data)
  }

  # Select miRNAs based on selection method
  mirna_cols <- colnames(disease_data)[grep("^hsa-", colnames(disease_data))]

  selected_mirnas <- switch(mirna_selection,
    "all" = mirna_cols,
    "significant" = {
      if (is.null(de_results)) {
        stop("DE results required for 'significant' miRNA selection")
      }
      sig_mirnas <- de_results %>%
        filter(adj.P.Val < 0.05) %>%
        pull(miRNA)
      intersect(sig_mirnas, mirna_cols)
    },
    "literature" = {
      lit_mirnas <- read.csv(file.path(paths$data_path, "researchMiRNAAccession.csv"))
      lit_mirnas_vec <- lit_mirnas %>% pull(miRNA)
      intersect(lit_mirnas_vec, mirna_cols)
    },
    "random" = sample(mirna_cols, min(num_features, length(mirna_cols)))
  )

  # If not "all", limit to specified number
  if (mirna_selection != "all") {
    selected_mirnas <- selected_mirnas[1:min(num_features, length(selected_mirnas))]
  }

  # Create feature matrix
  features <- disease_data %>%
    select(Disease, all_of(selected_mirnas))

  return(list(
    data = features,
    selected_mirnas = selected_mirnas
  ))
}

# Prepare ML data
ml_data <- prepare_ml_data(
  clinical_data = clinical_data,
  mirna_data = mirna_data,
  metadata = metadata,
  disease = DISEASE,
  mirna_selection = MIRNA_SELECTION,
  num_features = NUM_FEATURES,
  de_results = if (exists("de_results")) de_results else NULL,
  matched = MATCHED
)

# ----- Model Training and Evaluation -----
# Train and evaluate model
train_evaluate_model <- function(data, model_type = "glmnet", cv_repeats = 10) {
  # Prepare training control
  train_control <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = cv_repeats,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )

  # Set up model parameters
  model_params <- switch(model_type,
    "glmnet" = expand.grid(
      alpha = seq(0, 1, by = 0.2),
      lambda = 10^seq(-5, 0, length.out = 10)
    ),
    "rf" = expand.grid(
      mtry = seq(2, min(20, floor(sqrt(ncol(data) - 1) * 3)), by = 2)
    ),
    "xgb" = expand.grid(
      nrounds = c(50, 100, 150),
      max_depth = c(3, 6, 9),
      eta = c(0.01, 0.1, 0.3),
      gamma = c(0, 0.1),
      colsample_bytree = c(0.6, 0.8, 1.0),
      min_child_weight = c(1, 3, 5),
      subsample = c(0.6, 0.8, 1.0)
    ),
    stop("Unsupported model type")
  )

  # Train model
  model <- train(
    Disease ~ .,
    data = data,
    method = model_type,
    trControl = train_control,
    tuneGrid = model_params,
    metric = "ROC"
  )

  # Calculate performance metrics
  predictions <- predict(model, data, type = "prob")
  roc_curve <- roc(data$Disease, predictions[, 2])
  auc <- auc(roc_curve)

  # Create performance metrics dataframe
  perf_metrics <- data.frame(
    Model = model_type,
    AUC = auc,
    Accuracy = max(model$results$Accuracy),
    Sensitivity = max(model$results$Sens),
    Specificity = max(model$results$Spec)
  )

  # Variable importance
  if (model_type %in% c("glmnet", "rf")) {
    var_imp <- varImp(model)
  } else {
    var_imp <- NULL
  }

  return(list(
    model = model,
    predictions = predictions,
    roc_curve = roc_curve,
    performance = perf_metrics,
    importance = var_imp
  ))
}

# Train and evaluate model
model_results <- train_evaluate_model(
  data = ml_data$data,
  model_type = MODEL_TYPE,
  cv_repeats = CV_REPEATS
)

# ----- Visualization -----
# Create ROC curve plot
create_roc_plot <- function(roc_curve, title = "ROC Curve") {
  plot_data <- data.frame(
    specificity = 1 - roc_curve$specificities,
    sensitivity = roc_curve$sensitivities
  )

  roc_plot <- ggplot(plot_data, aes(x = specificity, y = sensitivity)) +
    geom_line(color = "blue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(
      title = title,
      subtitle = paste("AUC =", round(auc(roc_curve), 3)),
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )

  return(roc_plot)
}

# Create variable importance plot
create_importance_plot <- function(importance, n_features = 10) {
  if (is.null(importance)) {
    return(NULL)
  }

  imp_df <- importance$importance %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    arrange(desc(Overall)) %>%
    head(n_features)

  imp_plot <- ggplot(imp_df, aes(x = reorder(Feature, Overall), y = Overall)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = "Variable Importance",
      x = NULL,
      y = "Importance"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  return(imp_plot)
}

# Create ROC plot
roc_plot <- create_roc_plot(
  roc_curve = model_results$roc_curve,
  title = paste("ROC Curve -", DISEASE, "-", MODEL_TYPE)
)

# Create importance plot
imp_plot <- create_importance_plot(
  importance = model_results$importance,
  n_features = 10
)

# ----- Save Results -----
if (SAVE_FILES) {
  output_dir <- file.path(paths$output_path, "models")
  perf_dir <- file.path(paths$output_path, "performance_summary_df")
  plot_dir <- file.path(paths$output_path, "plots")

  # Ensure directories exist
  ensure_dir_exists(output_dir)
  ensure_dir_exists(perf_dir)
  ensure_dir_exists(plot_dir)

  # Save model
  save_project_object(
    model_results$model,
    name = paste0("model_", DISEASE, "_", MODEL_TYPE, "_", MIRNA_SELECTION),
    type = ifelse(MATCHED, "matched", "all"),
    extension = "rds",
    dir = output_dir
  )

  # Save performance metrics
  save_project_object(
    model_results$performance,
    name = paste0("performance_", DISEASE, "_", MODEL_TYPE, "_", MIRNA_SELECTION),
    type = ifelse(MATCHED, "matched", "all"),
    extension = "csv",
    dir = perf_dir
  )

  # Save plots
  ggsave(
    file.path(plot_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_roc_", DISEASE, "_", MODEL_TYPE, ".png")),
    roc_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  if (!is.null(imp_plot)) {
    ggsave(
      file.path(plot_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_imp_", DISEASE, "_", MODEL_TYPE, ".png")),
      imp_plot,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}

# Print performance metrics
print(model_results$performance)

# ----- Session Info -----
# Print session information for reproducibility
print(sessionInfo())
