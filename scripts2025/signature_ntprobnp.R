## ML signatures, ntproBNP analysis

# ----- Setup -----
data_path_bestageing2022 <- "/mnt/nas185/reich/rockerprojects/bestageing2022"
data_path_BestAgeing <- "/mnt/nas185/reich/BestAgeing"


# Load required libraries
library(tidyverse)
library(limma)
library(edgeR)
library(ggplot2)
library(readxl)
library(writexl)
library(pROC)
library(readxl)
library(janitor)
library(glue)
library(ggthemes)
library(pROC)
library(dplyr) # Ensure dplyr is loaded for case_when
library(stringr) # For str_extract function
# Load required libraries for ML predictions
library(workflows)
library(tune)
library(tidymodels)
library(colino)

source(
  file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R")
)

# load data
dcm_cleaned <- readRDS("revision2025/data/dcm_cleaned_data.rds")
hfref_cleaned <- readRDS("revision2025/data/hfref_cleaned_data.rds")


# ML Model Prediction Analysis Function
add_ml_predictions <- function(data, disease_name) {
  cat("\n=== ADDING ML MODEL PREDICTIONS FOR", toupper(disease_name), "===\n")

  # Define models to test
  MODEL <- c("logistic_reg_norm", "logistic_reg_simple", "RF", "XGB")

  # Initialize data with predictions
  data_with_predictions <- data

  successful_models <- c()

  for (model in MODEL) {
    cat("Loading model:", model, "\n")

    # Try loading full model first, then selected model
    model_paths <- c(
      glue(
        "{data_path_bestageing2022}/output/test_set_results/test_results_202311/{disease_name}/004c_TEST_RESULTS_{model}_analysis_full_m20.rds"
      ),
      glue(
        "{data_path_bestageing2022}/output/test_set_results/test_results_202311/{disease_name}/selected_analysis/004c_TEST_RESULTS_{model}_analysis_selected.rds"
      ),
      glue(
        "{data_path_bestageing2022}/output/test_set_results/test_results_202311/{disease_name}/selected_analysis/004c_20240121_TEST_RESULTS_{model}_analysis_selected_MATCHED.rds"
      )
    )

    loaded_models <- c()

    for (i in seq_along(model_paths)) {
      if (file.exists(model_paths[i])) {
        tryCatch(
          {
            cat("  Loading from:", model_paths[i], "\n")
            test_results_tmp <- readRDS(model_paths[i])

            # Extract workflow
            extract_workflow_loop <- extract_workflow(test_results_tmp)

            # Prepare data for prediction (keep only miRNA columns)
            prediction_data <- data %>%
              select(pat_id, age, sex, starts_with("hsa_"))

            # Make predictions
            predictions <- predict(
              extract_workflow_loop,
              new_data = prediction_data,
              type = "prob"
            )

            # Create column names
            model_type <- case_when(
              i == 1 ~ "full",
              str_detect(model_paths[i], "MATCHED") ~ "selected_matched",
              TRUE ~ "selected"
            )
            colname_pred_control <- glue(
              ".pred_control_{model}_{model_type}_{disease_name}"
            )
            colname_pred_disease <- glue(
              ".pred_{disease_name}_{model}_{model_type}_{disease_name}"
            )

            # Add predictions to data
            temp_dat <- data %>%
              select(pat_id) %>%
              bind_cols(predictions) %>%
              rename(
                !!colname_pred_control := !!colnames(predictions)[1],
                !!colname_pred_disease := !!colnames(predictions)[2]
              )

            data_with_predictions <- data_with_predictions %>%
              left_join(temp_dat, by = "pat_id", relationship = "many-to-many")

            loaded_models <- c(loaded_models, paste0(model, "_", model_type))

            cat(
              "  Successfully added predictions for",
              model,
              "(",
              model_type,
              ")\n"
            )
          },
          error = function(e) {
            cat("  Error loading", model_paths[i], ":", e$message, "\n")
          }
        )
      } else {
        cat("  File not found:", model_paths[i], "\n")
      }
    }

    successful_models <- c(successful_models, loaded_models)

    if (length(loaded_models) == 0) {
      cat("  Failed to load any model for", model, "\n")
    }
  }

  cat(
    "Successfully loaded",
    length(successful_models),
    "models:",
    paste(successful_models, collapse = ", "),
    "\n"
  )

  return(list(
    data = data_with_predictions,
    successful_models = successful_models
  ))
}

# ML Model Performance Analysis Function
analyze_ml_performance <- function(data, disease_name, successful_models) {
  cat("\n=== ML MODEL PERFORMANCE ANALYSIS FOR", toupper(disease_name), "===\n")

  if (length(successful_models) == 0) {
    cat("No successful models to analyze\n")
    return(tibble())
  }

  # Get prediction columns
  pred_cols <- names(data)[grepl(
    paste0(".pred_", disease_name, "_"),
    names(data)
  )]

  if (length(pred_cols) == 0) {
    cat("No prediction columns found\n")
    return(tibble())
  }

  results <- tibble()

  for (pred_col in pred_cols) {
    # Extract model info from column name
    model_info <- str_extract(pred_col, "(?<=.pred_)[^_]+_[^_]+_[^_]+")

    # Prepare data for analysis
    analysis_data <- data %>%
      filter(
        !is.na(.data[[pred_col]]),
        !is.na(ntprobnp),
        !is.na(age),
        !is.na(sex)
      ) %>%
      mutate(disease_binary = ifelse(disease == disease_name, 1, 0))

    if (nrow(analysis_data) < 10) {
      cat("Insufficient data for", pred_col, "\n")
      next
    }

    # Model 1: ML signature alone + age + sex
    f1 <- as.formula(paste("disease_binary ~", pred_col, "+ age + sex"))
    glm1 <- tryCatch(
      {
        glm(f1, data = analysis_data, family = binomial())
      },
      error = function(e) NULL
    )

    # Model 2: NT-proBNP + age + sex (baseline)
    f2 <- as.formula("disease_binary ~ ntprobnp + age + sex")
    glm2 <- tryCatch(
      {
        glm(f2, data = analysis_data, family = binomial())
      },
      error = function(e) NULL
    )

    # Model 3: Combined ML signature + NT-proBNP + age + sex
    f3 <- as.formula(paste(
      "disease_binary ~",
      pred_col,
      "+ ntprobnp + age + sex"
    ))
    glm3 <- tryCatch(
      {
        glm(f3, data = analysis_data, family = binomial())
      },
      error = function(e) NULL
    )

    if (!is.null(glm1) && !is.null(glm2) && !is.null(glm3)) {
      # Calculate AUCs
      auc1 <- as.numeric(
        roc(
          analysis_data$disease_binary,
          predict(glm1, type = "response"),
          quiet = TRUE
        )$auc
      )
      auc2 <- as.numeric(
        roc(
          analysis_data$disease_binary,
          predict(glm2, type = "response"),
          quiet = TRUE
        )$auc
      )
      auc3 <- as.numeric(
        roc(
          analysis_data$disease_binary,
          predict(glm3, type = "response"),
          quiet = TRUE
        )$auc
      )

      # Get p-values
      p_ml <- tryCatch(coef(summary(glm1))[2, 4], error = function(e) NA)
      p_ntprobnp <- tryCatch(coef(summary(glm2))[2, 4], error = function(e) NA)
      p_combined_ml <- tryCatch(coef(summary(glm3))[2, 4], error = function(e) {
        NA
      })
      p_combined_ntprobnp <- tryCatch(
        coef(summary(glm3))[3, 4],
        error = function(e) NA
      )

      # Likelihood ratio test for incremental benefit
      lrt_p <- tryCatch(
        {
          anova(glm2, glm3, test = "Chisq")$`Pr(>Chi)`[2]
        },
        error = function(e) NA
      )

      results <- bind_rows(
        results,
        tibble(
          model = model_info,
          n_samples = nrow(analysis_data),
          auc_ml_alone = auc1,
          auc_ntprobnp_baseline = auc2,
          auc_combined = auc3,
          incremental_benefit = auc3 - auc2,
          p_ml_signature = p_ml,
          p_ntprobnp_baseline = p_ntprobnp,
          p_combined_ml = p_combined_ml,
          p_combined_ntprobnp = p_combined_ntprobnp,
          p_incremental = lrt_p
        )
      )

      cat(sprintf("  %s:\n", model_info))
      cat(sprintf(
        "    ML signature + age + sex:     AUC=%.3f (p=%.3f)\n",
        auc1,
        p_ml
      ))
      cat(sprintf(
        "    NT-proBNP + age + sex:        AUC=%.3f (p=%.3f)\n",
        auc2,
        p_ntprobnp
      ))
      cat(sprintf("    Combined model:               AUC=%.3f\n", auc3))
      cat(sprintf(
        "    Incremental benefit:          Δ=%.3f (p=%.3f)\n",
        auc3 - auc2,
        lrt_p
      ))
      cat(sprintf(
        "    Sample size:                  n=%d\n\n",
        nrow(analysis_data)
      ))
    }
  }

  return(results)
}

# ML Gray Zone Analysis Function
ml_gray_zone_analysis <- function(data, disease_name, successful_models) {
  cat("\n=== ML GRAY ZONE ANALYSIS FOR", toupper(disease_name), "===\n")

  if (length(successful_models) == 0) {
    cat("No successful models to analyze\n")
    return(tibble())
  }

  # Use data with actual NT-proBNP measurements only
  data_original <- data %>%
    filter(!is.na(ntprobnp))

  # Remove imputed control values using the stored imputation value
  imputation_value <- attr(data, "control_imputation_value")
  if (!is.null(imputation_value)) {
    data_original <- data_original %>%
      filter(!(disease == "control" & abs(ntprobnp - imputation_value) < 0.01))
  }

  # Define gray zone (clinical threshold: 125-450 pg/mL)
  gray_zone_data <- data_original %>%
    filter(ntprobnp >= 125, ntprobnp <= 450)

  if (nrow(gray_zone_data) < 20) {
    cat("Insufficient data in gray zone for", disease_name, "\n")
    return(tibble())
  }

  cat(sprintf("Gray zone samples (125-450 pg/mL): %d\n", nrow(gray_zone_data)))

  # Get prediction columns
  pred_cols <- names(gray_zone_data)[grepl(
    paste0(".pred_", disease_name, "_"),
    names(gray_zone_data)
  )]

  results <- tibble()

  for (pred_col in pred_cols) {
    # Extract model info
    model_info <- str_extract(pred_col, "(?<=.pred_)[^_]+_[^_]+_[^_]+")

    # Prepare data
    analysis_data <- gray_zone_data %>%
      filter(!is.na(.data[[pred_col]]), !is.na(age), !is.na(sex)) %>%
      mutate(disease_binary = ifelse(disease == disease_name, 1, 0))

    if (nrow(analysis_data) < 10) {
      next
    }

    # Models in gray zone
    f1 <- as.formula(paste("disease_binary ~", pred_col, "+ age + sex"))
    f2 <- as.formula("disease_binary ~ ntprobnp + age + sex")
    f3 <- as.formula(paste(
      "disease_binary ~",
      pred_col,
      "+ ntprobnp + age + sex"
    ))

    glm1 <- tryCatch(
      glm(f1, data = analysis_data, family = binomial()),
      error = function(e) NULL
    )
    glm2 <- tryCatch(
      glm(f2, data = analysis_data, family = binomial()),
      error = function(e) NULL
    )
    glm3 <- tryCatch(
      glm(f3, data = analysis_data, family = binomial()),
      error = function(e) NULL
    )

    if (!is.null(glm1) && !is.null(glm2) && !is.null(glm3)) {
      auc1 <- as.numeric(
        roc(
          analysis_data$disease_binary,
          predict(glm1, type = "response"),
          quiet = TRUE
        )$auc
      )
      auc2 <- as.numeric(
        roc(
          analysis_data$disease_binary,
          predict(glm2, type = "response"),
          quiet = TRUE
        )$auc
      )
      auc3 <- as.numeric(
        roc(
          analysis_data$disease_binary,
          predict(glm3, type = "response"),
          quiet = TRUE
        )$auc
      )

      lrt_p <- tryCatch(
        anova(glm2, glm3, test = "Chisq")$`Pr(>Chi)`[2],
        error = function(e) NA
      )

      results <- bind_rows(
        results,
        tibble(
          model = model_info,
          n_gray_zone = nrow(analysis_data),
          auc_ml_gray = auc1,
          auc_ntprobnp_gray = auc2,
          auc_combined_gray = auc3,
          incremental_benefit_gray = auc3 - auc2,
          p_incremental_gray = lrt_p
        )
      )

      cat(sprintf("  %s Gray Zone:\n", model_info))
      cat(sprintf("    ML signature:     AUC=%.3f\n", auc1))
      cat(sprintf("    NT-proBNP:        AUC=%.3f\n", auc2))
      cat(sprintf("    Combined:         AUC=%.3f\n", auc3))
      cat(sprintf(
        "    Incremental:      Δ=%.3f (p=%.3f)\n",
        auc3 - auc2,
        lrt_p
      ))
    }
  }

  return(results)
}

# Add ML predictions to cleaned data
cat("\n", rep("=", 70), "\n")
cat("MACHINE LEARNING MODEL PREDICTIONS ANALYSIS\n")
cat(rep("=", 70), "\n")

# Add predictions for DCM
dcm_ml_result <- add_ml_predictions(dcm_cleaned$glm_data, "dcm")
dcm_data_with_ml <- dcm_ml_result$data
dcm_successful_models <- dcm_ml_result$successful_models

# Add predictions for HFrEF
hfref_ml_result <- add_ml_predictions(hfref_cleaned$glm_data, "hfref")
hfref_data_with_ml <- hfref_ml_result$data
hfref_successful_models <- hfref_ml_result$successful_models

# Analyze ML model performance
dcm_ml_performance <- analyze_ml_performance(
  dcm_data_with_ml,
  "dcm",
  dcm_successful_models
)
hfref_ml_performance <- analyze_ml_performance(
  hfref_data_with_ml,
  "hfref",
  hfref_successful_models
)

# ML Gray zone analysis
dcm_ml_gray <- ml_gray_zone_analysis(
  dcm_data_with_ml,
  "dcm",
  dcm_successful_models
)
hfref_ml_gray <- ml_gray_zone_analysis(
  hfref_data_with_ml,
  "hfref",
  hfref_successful_models
)

# Save ML analysis results
if (nrow(dcm_ml_performance) > 0) {
  write_csv(dcm_ml_performance, "output/dcm_ml_performance_analysis.csv")
  cat(
    "DCM ML performance results saved to output/dcm_ml_performance_analysis.csv\n"
  )
}

if (nrow(hfref_ml_performance) > 0) {
  write_csv(hfref_ml_performance, "output/hfref_ml_performance_analysis.csv")
  cat(
    "HFrEF ML performance results saved to output/hfref_ml_performance_analysis.csv\n"
  )
}

if (nrow(dcm_ml_gray) > 0) {
  write_csv(dcm_ml_gray, "output/dcm_ml_gray_zone_analysis.csv")
}

if (nrow(hfref_ml_gray) > 0) {
  write_csv(hfref_ml_gray, "output/hfref_ml_gray_zone_analysis.csv")
}

# Summary of ML model performance
cat("\n=== ML MODEL PERFORMANCE SUMMARY ===\n")

if (nrow(dcm_ml_performance) > 0) {
  cat("\nDCM - Best performing ML models:\n")
  dcm_ml_best <- dcm_ml_performance %>%
    arrange(desc(auc_combined)) %>%
    select(
      model,
      auc_ml_alone,
      auc_ntprobnp_baseline,
      auc_combined,
      incremental_benefit
    ) %>%
    head(3)
  print(dcm_ml_best)
}

if (nrow(hfref_ml_performance) > 0) {
  cat("\nHFrEF - Best performing ML models:\n")
  hfref_ml_best <- hfref_ml_performance %>%
    arrange(desc(auc_combined)) %>%
    select(
      model,
      auc_ml_alone,
      auc_ntprobnp_baseline,
      auc_combined,
      incremental_benefit
    ) %>%
    head(3)
  print(hfref_ml_best)
}
