


library(dplyr)
library(yardstick)  # for metric calculation
library(vip)        # for permutation importance
library(boot)
library(purrr)
library(progress)

# Assuming you have a fitted model `final_model` and data `test_data`
original_data <- results$dat_test
original_data <- original_data %>%
  mutate(original_pred = predict(final_model, original_data)$.pred)
original_rsquared <- yardstick::rsq_vec(truth = original_data |> pull({target}), estimate = original_data$original_pred)

# Function to calculate drop in R-squared after permutation
calculate_importance <- function(data, feature) {
  # Permute the feature
  data_permuted <- data
  data_permuted[[feature]] <- sample(data[[feature]])
  
  # Add predictions for permuted data
  prediction <- predict(final_model, data_permuted)$.pred
  
  # Calculate R-squared for permuted data
  permuted_rsquared <- rsq_vec(truth = data_permuted |> pull({target}), estimate = prediction)
  
  # Calculate the drop in performance
  drop_in_performance <- original_rsquared - permuted_rsquared
  return(drop_in_performance)
}

# Bootstrapping function to return importance scores directly
bootstrap_importance <- function(data, indices, feature) {
  sample_data <- data[indices, ]
  return(calculate_importance(sample_data, feature))
}


# Calculate importance for all features ----------
set.seed(123)
importance_scores <- map_dbl(names(results$dat_test)[-c(1:2)], ~ calculate_importance(original_data, .x), .progress = "Calculating Importance Scores")
importance_scores_df <- tibble(Feature = names(results$dat_test)[-c(1:2)], Importance = importance_scores)

## bootstrap importance CI for features with importance !=0 -------
importance_scores_df <- importance_scores_df |> 
  as_tibble() |> 
  filter(Importance!=0)

# Perform bootstrapping to calculate importance scores and confidence intervals
set.seed(123)
n_bootstraps <- 20
non_zero_features <- importance_scores_df |> pull(Feature)

results_bootstrapped_importance <- map_dfr(non_zero_features, ~ {
  boot_res <- boot(original_data, bootstrap_importance, R = n_bootstraps, feature = .x)
  boot_scores <- boot_res$t
  mean_importance <- mean(boot_scores)
  
  # Attempt to calculate BCa CI, fallback to Percentile on error
  ci_result <- tryCatch({
    boot.ci(boot_res, type = "bca")
  }, error = function(e) {
    message(paste("BCa CI calculation failed for feature:", .x, "- Falling back to percentile CI."))
    boot.ci(boot_res, type = "perc")
  })
  
  # Check if BCa or Percentile CIs have been used
  if (!is.null(ci_result$bca)) {
    ci_low <- ci_result$bca[4]
    ci_high <- ci_result$bca[5]
  } else if (!is.null(ci_result$percent)) {
    ci_low <- ci_result$percent[4]
    ci_high <- ci_result$percent[5]
  } else {
    ci_low <- NA  # In case both methods fail
    ci_high <- NA
  }
  
  tibble(
    Feature = .x,
    MeanImportance = mean_importance,
    Low = ci_low,
    High = ci_high
  )
}, .progress = "Calculating Bootstrapped Importances and CIs")


# Plot the results
ggplot(
  results_bootstrapped_importance |> mutate(Feature = convert_mir_name_V(Feature)), 
  aes(x = reorder(Feature, MeanImportance), y = MeanImportance)
) +
  geom_col(fill = ggthemes::gdocs_pal()(1)) +
  geom_errorbar(aes(ymin = Low, ymax = High), width = 0.2) +
  coord_flip() +
  labs(
    y = "Mean Bootstrapped Importance", 
    x = "Feature", 
    title = "", # "Feature Importance with 95% Confidence Intervals from Bootstrapping"
  ) +
  theme_minimal(base_size = 16, base_family = 'Arial') +
  theme(legend.position = "none", axis.title.y = element_blank())


