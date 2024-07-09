
require(tidymodels)
require(finetune)
library(glue)
library(dplyr)
library(caret)
library(pROC)
library(yardstick)
library(rsample)  # For bootstrapping
library(broom)    # For augment()
library(boot)
library(purrr)
library(ggplot2)
library(ggthemes)

convert_mir_name <- function(name) {
  # Replace 'mir' with 'miR'
  name <- gsub("mir", "miR", name)
  
  # Replace underscores with hyphens
  name <- gsub("_", "-", name)
  
  return(name)
}
convert_mir_name_V <- Vectorize(convert_mir_name)

# Regression function ------
run_regression_analysis <- function(
    modeldata, 
    .outcome_var, 
    .predictors,
    .id_col = NULL,
    .exclude_cols = NULL,  # must be kept after feature selection!
    split_proportion, 
    n_feature_select, 
    no_folds, 
    no_repeats,
    grid_size
) {
  ###
  #* Performs regression analysis 
  #* Data split/ Preprocessing/ Tuning
  #* returns tuning results 
  ###
  
  set.seed(123)
  
  # Prepare the data: Ensure the outcome variable is treated as numeric
  modeldata <- modeldata %>%
    mutate(!!.outcome_var := as.numeric(.data[[.outcome_var]])) |> 
    select(all_of(c(.id_col, .outcome_var)), all_of(.predictors)) |> 
    # drop missings in outcome
    drop_na(!!.outcome_var)
  
  if(nrow(modeldata) < 10){
    stop(call. = FALSE, "Too many missings in outcome variable!")
  }
  
  ### Split the data ----
  dat_split <- initial_split(modeldata, prop = split_proportion)
  dat_train <- training(dat_split)
  dat_test <- testing(dat_split)
  folds <- vfold_cv(dat_train, v = no_folds, repeats = no_repeats)
  
  # Prepare a recipe for preprocessing the data
  
  # * Tidyeval ------- https://stackoverflow.com/questions/69742572/create-recipes-and-passing-column-names-dynamically
  # outcome_col_var_expr      <- rlang::enquo(.outcome_var)
  # pred_col_var_expr      <- rlang::enquo(.predictors)
  
  # .outcome_var  <- rlang::sym(names(modeldata)[[1]])
  
  ### recipes ------
  f <- as.formula(paste(.outcome_var, " ~ ."))
  
  base_recipe <- recipe(formula = f , data = dat_train) |> 
    update_role(!!.id_col, new_role = "ID") |> 
    step_zv(all_predictors()) |> 
    step_impute_mean(all_numeric_predictors()) |> 
    step_impute_mode(all_nominal_predictors()) |> 
    step_corr(all_numeric_predictors(), threshold = 0.9)
  
  # Normalized recipe with specific feature selection
  normalized_rec <- base_recipe %>%
    step_normalize(all_numeric_predictors()) |> 
    step_dummy(all_nominal_predictors()) |> 
    colino::step_select_forests(
      all_predictors(), 
      -all_of( contains(.exclude_cols)), 
      outcome = .outcome_var, 
      top_p = n_feature_select
    )
  
  models <- list(
    "LM" = linear_reg(penalty = tune(), mixture = tune()) %>% set_engine("glmnet") %>% set_mode("regression"),
    "RF" = rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>% set_engine("ranger") %>% set_mode("regression")
  )
  
  workflows <- workflow_set(
    preproc = list(normalized = normalized_rec),
    models = models
  )
  
  # Finalize mtry based on the preprocessed data
  # Assuming 'dat_train' is available at this point and has been processed by `normalized_rec`
  # We need to apply the preprocessing steps to calculate the actual number of predictors
  prepped_data <- prep(normalized_rec, training = dat_train) %>% 
    bake(new_data = NULL)
  
  # Number of predictors could be counted like so:
  num_predictors <- ncol(prepped_data) - 2  # remove outcome variable + id var
  
  # Grid definition and tuning setup
  grid_list <- lapply(models, function(model_spec) {
    # Extract parameters and finalize mtry based on number of predictors
    params <- extract_parameter_set_dials(model_spec)
    if ("mtry" %in% params$name) {
      # Typically, a good starting point is to 
      ## a) try square root of the total number of predictors for classification tasks `ceiling(sqrt(num_predictors))` 
      ## b) one-third for regression tasks.
      params <- params %>%
        update(mtry = mtry(range = c(1, num_predictors))) |> 
        update(
          trees = trees(range = c(10, 500)),  # default c(1, 2000)   # CODE: trees() %>% range_get()
          min_n = min_n(range = c(2, 40)) # default c(2, 40)  control complexity  -> increasing min_child_weight, forces algorithm to create simpler trees
        )
    }
    
    params |> 
      grid_latin_hypercube(size = grid_size)
  })
  
  
  for (i in seq_along(grid_list)) {
    # Construct the wflow_id using the pattern from preprocessing name and model name
    wflow_id <- paste("normalized", names(grid_list)[i], sep = "_")
    workflows <- option_add(workflows, grid = grid_list[[i]], id = wflow_id)
  }
  
  ### Tune the workflows ---------
  race_ctrl <- control_race(verbose = FALSE, allow_par = TRUE, save_pred = TRUE, burn_in = 3, num_ties = 10, alpha = 0.05, randomize = TRUE, parallel_over = "everything", save_workflow = FALSE)
  race_results <- workflows %>%
    workflow_map(
      "tune_race_anova",
      resamples = folds,
      metrics = metric_set(rmse, rsq),
      control = race_ctrl,
      seed = 123,
      verbose = TRUE
    )
  
  return(
    list(
      race_results = race_results, 
      modeldata = modeldata, 
      dat_split = dat_split, dat_train = dat_train, dat_test = dat_test, 
      recipe = normalized_rec
    )
  )
}


# evaluate ----------------------------------------------------------------

## evaluate plotting helpers  --------

create_prediction_error_plot <- function(plot_data, trait, wflow_id) {
  ggplot(plot_data, aes(x = actual, y = predicted)) +
    geom_point(alpha = 0.5, shape=16) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = glue("Predictions vs Actual for {trait} - {wflow_id}"),
         x = "Actual", y = "Predicted") +
    theme_minimal(base_size = 16, base_family = 'Arial')
}

create_residual_plot <- function(plot_data) {
  ggplot(plot_data, aes(x = predicted, y = actual - predicted)) +
    geom_point(shape=16) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Residuals vs. Fitted", x = "Fitted Values", y = "Residuals") +
    theme_minimal(base_size = 16, base_family = 'Arial')
}

create_qq_plot <- function(plot_data) {
  ggplot(plot_data, aes(sample = actual - predicted)) +
    stat_qq(alpha=0.5, shape=16) +
    stat_qq_line() +
    labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal(base_size = 16, base_family = 'Arial')
}

## evaluate process model function  --------

# Define a functionfor bootstrapping to calculate metrics on each resample
bootstrap_metrics <- function(split, data, trait_col) {
  # Resampled data
  resampled_data <- analysis(split)
  
  # Calculate metrics for train and test combined
  combined_metrics <- resampled_data %>%
    yardstick::metrics(trait_col, .pred)
  
  # Return the metrics
  return(combined_metrics)
}

#' Process Model Results
#'
#' This function takes the results of a model evaluation, processes it to extract
#' and finalize the best model based on hyperparameters, and calculates the model's
#' performance metrics on test data. It also generates diagnostic plots to evaluate
#' the model's predictions.
#'
#' @param results A list containing the results from multiple model evaluations,
#'                typically read from an RDS file, which includes data like model
#'                tuning results, data split, test data, and best hyperparameter settings.
#' @param trait A string representing the specific trait being analyzed, which
#'              influences the plots' labeling and subset of predictions to be
#'              analyzed.
#' @param model_idx An integer index specifying which model's results to process
#'                  from the `results$race_results` data.
#'
#' @return A list containing the model's performance metrics (R-squared, RMSE)
#'         and a list of ggplot objects for the prediction error plot, residual
#'         plot, and QQ plot both for ONLY test set and combined train/test set
#'
process_model <- function(results, trait, model_idx) {
  wflow_id <- results$race_results[[1]][model_idx]
  
  # select best hyperparams
  best_results_tmp <- 
    results$race_results %>% 
    extract_workflow_set_result(wflow_id) %>% 
    select_best(metric = "rsq")
  
  # fit to training data, and calculate on test data only
  test_results_tmp <- 
    results$race_results %>% 
    extract_workflow(wflow_id) %>% 
    finalize_workflow(best_results_tmp) %>% 
    last_fit(split = results$dat_split)
  
  # Calculate metrics
  last_fit_metrics <- collect_metrics(test_results_tmp)
  
  # Store metrics
  model_results <- list(
    R_squared = last_fit_metrics %>% filter(.metric == "rsq") %>% pull(.estimate),
    RMSE = last_fit_metrics %>% filter(.metric == "rmse") %>% pull(.estimate)
  )
  
  # Generate plots
  predictions <- collect_predictions(test_results_tmp)
  plot_data <- tibble(actual = predictions[[trait]], predicted = predictions$.pred)
  
  model_results$plots <- list(
    prediction_error = create_prediction_error_plot(plot_data, trait, wflow_id),
    residual_plot = create_residual_plot(plot_data),
    qq_plot = create_qq_plot(plot_data)
  )
  
  # predict on train and test data together https://stackoverflow.com/questions/68124804/tidymodels-get-predictions-and-metrics-on-training-data-using-workflow-recipe
  final_model <- test_results_tmp$.workflow[[1]]
  
  all_predictions <- bind_rows(
    augment(final_model, new_data = results$dat_train) %>% 
      mutate(type = "train"),
    augment(final_model, new_data = results$dat_test) %>% 
      mutate(type = "test")
  )
  
  train_test_metrics_grouped <- all_predictions %>%
    group_by(type) %>%
    yardstick::metrics(trait, .pred)
  # add train/test combined
  train_test_metrics_grouped <- train_test_metrics_grouped |> 
    rbind(all_predictions %>%
            yardstick::metrics(trait, .pred) |> 
            mutate(type="traintest")
    ) |> 
    arrange(.metric, desc(type))
  
  ### bootstrap CIs --------
  set.seed(123)  
  bootstrap_samples <- 200
  bootstraps_traintest <- rsample::bootstraps(all_predictions, times = bootstrap_samples)
  bootstraps_train <- rsample::bootstraps(all_predictions |> filter(type=="train"), times = bootstrap_samples)
  bootstraps_test <- rsample::bootstraps(all_predictions|> filter(type=="test"), times = bootstrap_samples)
  
  # Apply function to each bootstrap sample
  bootstraps_traintest <- bootstraps_traintest %>%
    mutate(boot_results = map(
      splits, 
      bootstrap_metrics, # function from above, returns tibble with metrics
      data = all_predictions, trait_col = trait # fn args
    )
    ) %>%
    unnest(boot_results)
  
  bootstraps_train <- bootstraps_train %>%
    mutate(boot_results = map(splits, bootstrap_metrics, data = all_predictions |> filter(type=="train"), trait_col = trait)) %>%
    unnest(boot_results)
  bootstraps_test <- bootstraps_test %>%
    mutate(boot_results = map(splits, bootstrap_metrics, data = all_predictions|> filter(type=="test"), trait_col = trait)) %>%
    unnest(boot_results)
  
  # Compute the confidence intervals for each metric
  ci_traintest <- bootstraps_traintest %>%
    group_by(.metric) %>%
    summarise(
      lower = quantile(.estimate, probs = 0.025),
      upper = quantile(.estimate, probs = 0.975)
    ) |> 
    mutate(type = "traintest")
  ci_train <- bootstraps_train %>%
    group_by(.metric) %>%
    summarise(
      lower = quantile(.estimate, probs = 0.025),
      upper = quantile(.estimate, probs = 0.975)
    ) |> 
    mutate(type = "train")
  ci_test <- bootstraps_test %>%
    group_by(.metric) %>%
    summarise(
      lower = quantile(.estimate, probs = 0.025),
      upper = quantile(.estimate, probs = 0.975)
    ) |> 
    mutate(type = "test")
  ci_combined <- ci_traintest |> rbind(ci_train, ci_test)
  
  train_test_metrics_grouped <- train_test_metrics_grouped |> 
    left_join(ci_combined, by=join_by(type, .metric))
  
  model_results$train_test_metrics_grouped <- train_test_metrics_grouped
  
  # Generate traintest plots
  plot_data_traintest <- tibble(actual = all_predictions[[trait]], predicted = all_predictions$.pred, type=all_predictions$type)
  prediction_error_traintest <- create_prediction_error_plot(plot_data_traintest, trait, wflow_id) +
    geom_point(aes(color = type), alpha = 0.5) +
    ggthemes::scale_color_gdocs() +
    theme(legend.title = element_blank())
  residual_plot_traintest <- create_residual_plot(plot_data_traintest) +
    geom_point(aes(color = type), alpha = 0.5) +
    ggthemes::scale_color_gdocs() +
    theme(legend.title = element_blank())
  
  model_results$plots_traintest <- list(
    prediction_error = prediction_error_traintest,
    residual_plot = residual_plot_traintest,
    qq_plot = create_qq_plot(plot_data_traintest)
  )
  
  return(model_results)
}


## calculate feature importance manually with bootstrapping function  --------

calculate_feature_importance <- function(results, trait, model_idx, n_bootstraps) {
  # 1) get model
  wflow_id <- results$race_results[[1]][model_idx]
  best_results_tmp <- 
    results$race_results %>% 
    extract_workflow_set_result(wflow_id) %>% 
    select_best(metric = "rsq")
  test_results_tmp <- 
    results$race_results %>% 
    extract_workflow(wflow_id) %>% 
    finalize_workflow(best_results_tmp) %>% 
    last_fit(split = results$dat_split)
  final_model <- test_results_tmp$.workflow[[1]]
  
  
  # 2) define data
  original_data <- results$dat_train
  
  # Assume `final_model` is a fitted model available globally or passed as a parameter
  # If it should be part of `results`, you'd access it similarly
  # final_model <- results$final_model
  
  # Adding predictions to the dataset
  original_data <- original_data %>%
    mutate(original_pred = predict(final_model, original_data)$.pred)
  
  # 3) Calculating the original R-squared
  original_rsquared <- yardstick::rsq_vec(
    truth = original_data %>% pull({{ trait }}),
    estimate = original_data$original_pred
  )
  
  # 4) Function to calculate the drop in R-squared after permutation
  calculate_importance <- function(data, feature) {
    data_permuted <- data
    data_permuted[[feature]] <- sample(data[[feature]])
    prediction <- predict(final_model, data_permuted)$.pred
    permuted_rsquared <- rsq_vec(truth = data_permuted %>% pull({{ trait }}), estimate = prediction)
    drop_in_performance <- original_rsquared - permuted_rsquared
    return(drop_in_performance)
  }
  
  # Bootstrapping function to return importance scores directly
  bootstrap_importance <- function(data, indices, feature) {
    sample_data <- data[indices, ]
    return(calculate_importance(sample_data, feature))
  }
  
  # 5) Calculate importance for all features
  set.seed(123)
  importance_scores <- map_dbl(
    names(results$dat_train)[-c(1:2)], 
    ~ calculate_importance(original_data, .x), 
    .progress = "Calculating Importance Scores"
  )
  
  importance_scores_df <- tibble(
    Feature = names(results$dat_train)[-c(1:2)], 
    Importance = importance_scores
  ) %>% filter(Importance != 0)
  
  # 6) Perform bootstrapping only for != 0 features to calculate mean importance scores and confidence intervals
  set.seed(123)
  n_bootstraps <- n_bootstraps
  non_zero_features <- importance_scores_df %>% pull(Feature)
  
  results_bootstrapped_importance <- map_dfr(non_zero_features, ~ {
    boot_res <- boot(original_data, bootstrap_importance, R = n_bootstraps, feature = .x)
    boot_scores <- boot_res$t
    mean_importance <- mean(boot_scores)
    
    ci_result <- tryCatch({
      boot.ci(boot_res, type = "bca")
    }, error = function(e) {
      message(paste("BCa CI calculation failed for feature:", .x, "- Falling back to percentile CI."))
      boot.ci(boot_res, type = "perc")
    })
    
    if (!is.null(ci_result$bca)) {
      ci_low <- ci_result$bca[4]
      ci_high <- ci_result$bca[5]
    } else if (!is.null(ci_result$percent)) {
      ci_low <- ci_result$percent[4]
      ci_high <- ci_result$percent[5]
    } else {
      ci_low <- NA
      ci_high <- NA
    }
    
    tibble(
      Feature = .x,
      MeanImportance = mean_importance,
      Low = ci_low,
      High = ci_high
    )
  }, .progress = "Calculating Bootstrapped Importances and CIs")
  
  results_bootstrapped_importance <- results_bootstrapped_importance |> arrange(desc(MeanImportance))
  
  # Plot the results
  importance_plot <- ggplot(
    results_bootstrapped_importance |> mutate(Feature=convert_mir_name_V(Feature)), 
    aes(x = reorder(Feature, MeanImportance), y = MeanImportance)
  ) +
    geom_col(fill = gdocs_pal()(1)) +
    geom_errorbar(aes(ymin = Low, ymax = High), width = 0.2) +
    coord_flip() +
    labs(
      y = "Mean Bootstrapped Importance", 
      x = "Feature", 
      title = "", # "Feature Importance with 95% Confidence Intervals from Bootstrapping"
    ) +
    theme_minimal(base_size = 16, base_family = 'Arial') +
    theme(legend.position = "none", axis.title.y = element_blank())
  
  return(
    list(
      results_bootstrapped_importance=results_bootstrapped_importance,
      plot_bootstrapped_vip=importance_plot
      )
    )
}

