# Survival Analysis Using miRNA Expression
#
# This script performs survival analysis on miRNA data
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
  "survival",
  "survminer",
  "ggplot2",
  "glmnet",
  "readxl",
  "writexl"
))

# ----- Parameters -----
# Analysis parameters - edit these as needed
DISEASE <- "MI" # Choose from available diseases
MATCHED <- TRUE # Use matched samples?
SURVIVAL_ENDPOINT <- "OS" # Overall survival (OS) or Event-free survival (EFS)
TIME_UNIT <- "months" # Unit for time-to-event
MIRNA_SELECTION <- "significant" # How to select miRNAs: "all", "significant", "literature", "random"
NUM_FEATURES <- 10 # Number of features to use (if not "all")

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
    warning("No DE results found. Using all miRNAs instead.")
    MIRNA_SELECTION <- "all"
  }
}

# ----- Data Preparation -----
# Prepare data for survival analysis
prepare_survival_data <- function(clinical_data, mirna_data, metadata, disease,
                                  mirna_selection = "all", num_features = 10,
                                  de_results = NULL, matched = FALSE,
                                  survival_endpoint = "OS") {
  # Join clinical data with miRNA data
  merged_data <- inner_join(clinical_data, mirna_data, by = "patient_id")

  # Filter by disease and ensure survival data is available
  disease_data <- merged_data %>%
    filter(Disease == disease) %>%
    filter(!is.na(survival_time) & !is.na(survival_status))

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
        warning("DE results required for 'significant' miRNA selection. Using all miRNAs.")
        mirna_cols
      } else {
        sig_mirnas <- de_results %>%
          filter(adj.P.Val < 0.05) %>%
          pull(miRNA)
        intersect(sig_mirnas, mirna_cols)
      }
    },
    "literature" = {
      lit_mirnas <- read.csv(file.path(paths$data_path, "researchMiRNAAccession.csv"))
      lit_mirnas_vec <- lit_mirnas %>% pull(miRNA)
      intersect(lit_mirnas_vec, mirna_cols)
    },
    "random" = sample(mirna_cols, min(num_features, length(mirna_cols)))
  )

  # If not "all", limit to specified number
  if (mirna_selection != "all" && length(selected_mirnas) > num_features) {
    selected_mirnas <- selected_mirnas[1:num_features]
  }

  # Create feature matrix
  surv_data <- disease_data %>%
    select(patient_id, survival_time, survival_status, Age, Sex, BMI, all_of(selected_mirnas)) %>%
    mutate(
      status = as.numeric(survival_status == "Dead"),
      time = survival_time
    )

  # Ensure all miRNA values are numeric
  surv_data <- surv_data %>%
    mutate(across(all_of(selected_mirnas), as.numeric))

  return(list(
    data = surv_data,
    selected_mirnas = selected_mirnas
  ))
}

# Prepare survival data
surv_data <- prepare_survival_data(
  clinical_data = clinical_data,
  mirna_data = mirna_data,
  metadata = metadata,
  disease = DISEASE,
  mirna_selection = MIRNA_SELECTION,
  num_features = NUM_FEATURES,
  de_results = if (exists("de_results")) de_results else NULL,
  matched = MATCHED,
  survival_endpoint = SURVIVAL_ENDPOINT
)

# ----- Survival Analysis -----
# Perform basic survival analysis
perform_survival_analysis <- function(data, selected_mirnas) {
  # Create a survival object
  surv_obj <- Surv(data$time, data$status)

  # Basic Kaplan-Meier curve
  fit <- survfit(surv_obj ~ 1, data = data)

  # Initialize results list
  results <- list(
    overall_fit = fit,
    individual_fits = list(),
    cox_models = list(),
    risk_score_model = NULL
  )

  # Analyze each miRNA individually
  for (mirna in selected_mirnas) {
    # Create a binary variable based on median expression
    median_expr <- median(data[[mirna]], na.rm = TRUE)
    data[[paste0(mirna, "_high")]] <- as.factor(data[[mirna]] > median_expr)

    # Fit Kaplan-Meier
    mirna_fit <- survfit(surv_obj ~ get(paste0(mirna, "_high")), data = data)
    results$individual_fits[[mirna]] <- mirna_fit

    # Fit Cox proportional hazards model
    cox_formula <- as.formula(paste("surv_obj ~", mirna))
    cox_model <- coxph(cox_formula, data = data)
    results$cox_models[[mirna]] <- cox_model
  }

  # Fit multivariate Cox model
  if (length(selected_mirnas) > 1) {
    multi_formula <- as.formula(paste("surv_obj ~", paste(selected_mirnas, collapse = " + ")))
    multi_cox <- coxph(multi_formula, data = data)
    results$multivariate_cox <- multi_cox

    # Create a risk score model
    risk_score <- predict(multi_cox, newdata = data, type = "risk")
    data$risk_score <- risk_score
    data$high_risk <- as.factor(risk_score > median(risk_score))

    # Fit Kaplan-Meier by risk group
    risk_fit <- survfit(surv_obj ~ high_risk, data = data)
    results$risk_score_fit <- risk_fit
    results$risk_score_data <- data
  }

  return(results)
}

# Run survival analysis
surv_results <- perform_survival_analysis(
  data = surv_data$data,
  selected_mirnas = surv_data$selected_mirnas
)

# ----- Visualization -----
# Create Kaplan-Meier plot
create_km_plot <- function(fit, title = "Kaplan-Meier Curve", time_unit = "months") {
  km_plot <- ggsurvplot(
    fit,
    data = surv_data$data,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = paste("Time (", time_unit, ")", sep = ""),
    ggtheme = theme_bw(),
    risk.table.y.text = FALSE,
    title = title
  )

  return(km_plot)
}

# Create forest plot for hazard ratios
create_forest_plot <- function(cox_models, title = "Hazard Ratios") {
  # Extract hazard ratios and CIs
  hr_data <- lapply(names(cox_models), function(mirna) {
    model <- cox_models[[mirna]]
    summary_model <- summary(model)

    data.frame(
      mirna = mirna,
      hr = summary_model$conf.int[1, "exp(coef)"],
      lower = summary_model$conf.int[1, "lower .95"],
      upper = summary_model$conf.int[1, "upper .95"],
      p.value = summary_model$coefficients[1, "Pr(>|z|)"]
    )
  }) %>% bind_rows()

  # Sort by hazard ratio
  hr_data <- hr_data %>%
    arrange(desc(hr))

  # Create forest plot
  forest_plot <- ggplot(hr_data, aes(x = hr, y = reorder(mirna, hr), xmin = lower, xmax = upper)) +
    geom_point(size = 3, aes(color = p.value < 0.05)) +
    geom_errorbarh(height = 0.2, aes(color = p.value < 0.05)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("gray", "red"), guide = "none") +
    labs(
      title = title,
      y = NULL,
      x = "Hazard Ratio (95% CI)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(hjust = 0)
    )

  return(list(
    plot = forest_plot,
    data = hr_data
  ))
}

# Create plots
km_overall <- create_km_plot(
  fit = surv_results$overall_fit,
  title = paste("Overall Survival -", DISEASE),
  time_unit = TIME_UNIT
)

if (length(surv_results$cox_models) > 0) {
  forest_hr <- create_forest_plot(
    cox_models = surv_results$cox_models,
    title = paste("Hazard Ratios -", DISEASE)
  )
}

if (!is.null(surv_results$risk_score_fit)) {
  km_risk <- create_km_plot(
    fit = surv_results$risk_score_fit,
    title = paste("Survival by Risk Group -", DISEASE),
    time_unit = TIME_UNIT
  )
}

# ----- Save Results -----
if (SAVE_FILES) {
  output_dir <- file.path(paths$output_path, "survival_analysis")
  ensure_dir_exists(output_dir)

  # Save data
  save_project_object(
    surv_data$data,
    name = paste0("survival_data_", DISEASE),
    type = ifelse(MATCHED, "matched", "all"),
    extension = "rds",
    dir = output_dir
  )

  # Save results
  save_project_object(
    surv_results,
    name = paste0("survival_results_", DISEASE),
    type = ifelse(MATCHED, "matched", "all"),
    extension = "rds",
    dir = output_dir
  )

  # Save plots
  if (exists("km_overall")) {
    png(
      file.path(output_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_km_overall_", DISEASE, ".png")),
      width = 10, height = 8, units = "in", res = 300
    )
    print(km_overall)
    dev.off()
  }

  if (exists("forest_hr")) {
    ggsave(
      file.path(output_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_forest_hr_", DISEASE, ".png")),
      forest_hr$plot,
      width = 10,
      height = 8,
      dpi = 300
    )

    write.csv(
      forest_hr$data,
      file.path(output_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_hazard_ratios_", DISEASE, ".csv")),
      row.names = FALSE
    )
  }

  if (exists("km_risk")) {
    png(
      file.path(output_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_km_risk_", DISEASE, ".png")),
      width = 10, height = 8, units = "in", res = 300
    )
    print(km_risk)
    dev.off()
  }
}

# ----- Summary -----
# Print summary information
cat("===== Survival Analysis Summary =====\n")
cat("Disease:", DISEASE, "\n")
cat("Number of samples:", nrow(surv_data$data), "\n")
cat("Number of events:", sum(surv_data$data$status), "\n")
cat("Number of miRNAs used:", length(surv_data$selected_mirnas), "\n")
cat("miRNA selection method:", MIRNA_SELECTION, "\n")

# Print median survival
median_surv <- surv_results$overall_fit$table["median"]
cat("Median survival time:", ifelse(is.na(median_surv), "Not reached", paste(median_surv, TIME_UNIT)), "\n")

# Print top significant miRNAs
if (length(surv_results$cox_models) > 0) {
  top_mirnas <- lapply(names(surv_results$cox_models), function(mirna) {
    model <- surv_results$cox_models[[mirna]]
    summary_model <- summary(model)

    data.frame(
      mirna = mirna,
      hr = summary_model$conf.int[1, "exp(coef)"],
      p.value = summary_model$coefficients[1, "Pr(>|z|)"]
    )
  }) %>%
    bind_rows() %>%
    filter(p.value < 0.05) %>%
    arrange(p.value) %>%
    head(5)

  if (nrow(top_mirnas) > 0) {
    cat("\nTop significant miRNAs:\n")
    print(top_mirnas)
  } else {
    cat("\nNo significant miRNAs found.\n")
  }
}

# ----- Session Info -----
# Print session information for reproducibility
print(sessionInfo())
