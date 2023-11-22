

# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
SAVE.files <- TRUE

no_folds = 5
no_repeats = 10

# Define library and data paths based on system
if (system_name == "MacBook-Pro-CR-2065.local" | grepl("laptop-zim.uni-heidelberg.de", system_name)) {
  lib_path <- .libPaths()[1]
  data_path_bestageing2022 <- "/Volumes/T7CR/data/bestageing2022"
  data_path_BestAgeing <- "/Volumes/T7CR/data/BestAgeing"
  if(mount_filesystem == TRUE) {
    data_path_bestageing2022 <- "/Users/christophreich/Desktop/mount/rockerprojects/bestageing2022"  # mount -t nfs 10.55.1.185:/data/users/reich/ ~/Desktop/mount/
    data_path_BestAgeing <- "/Users/christophreich/Desktop/mount/BestAgeing"
  }
} else {  # assuming cluster
  .libPaths("/mnt/users/reich/programs/R43/lib")
  lib_path <- "/mnt/users/reich/programs/R43/lib" 
  data_path_bestageing2022 <- "/mnt/users/reich/rockerprojects/bestageing2022"
  data_path_BestAgeing <- "/mnt/users/reich/BestAgeing"
}


# libraries -------------------------------------------------------------------
require(tidyverse)
require(probably)
require(patchwork)
require(probably)
library(DALEXtra)
library(flextable)
library(officer)
library(glue)
require(colino, lib.loc = lib_path)
require(tidymodels)

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_calibration_plot.R"))
source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))
source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_vip_plot.R"))


# for word export
sect_properties <- prop_section(
  page_size = page_size(
    orient = "portrait",  # "portrait" "landscape"
    width = 8.3, height = 14
  ),
  type = "continuous",
  page_margins = page_mar(),
  #section_columns = section_columns(widths = c(4.75, 4.75))
)

diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="full")


# run analysis on test set for each disease ....................................
for (i in 1:nrow(all_combis)){
  # need to get data split right to finalize model
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_{all_combis$diseases[i]}_data01.rds")
  if(!file.exists(path2dataprocessed)) {
    next
  }
  
  data01 <- readRDS(file = path2dataprocessed)
  no.mirnas <- ncol(data01) - 4
  text_disease <- stringr::str_to_upper(all_combis$diseases[i])
  set.seed(123) # get same indices 
  modeldat <- data01
  
  # make control first factor for all analyses ;)
  modeldat <- modeldat %>% 
    mutate(disease = factor(disease, levels = c("control", all_combis$diseases[i])))
  
  dat_split <- rsample::initial_split(modeldat, strata = disease)
  dat_train <- training(dat_split)
  dat_test <- testing(dat_split)
  folds <- vfold_cv(dat_train, strata = disease, v = no_folds, repeats = no_repeats)
  
  # load race results ---------------------------------------------------------
  # more logit
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_randomMIR_FALSE_more_logit.rds"))
  
  # place to store
  best_results <- list()  # hyperparams selected
  test_results <- list()  # last fit
  
  aucs <- c()
  ci.auc.texts <- c()
  accuracies <- c()
  sensitiv <- c()
  specific <- c()
  ppv <- c()
  npv <- c()
  precision <- c()
  recall <- c()
  f1 <- c()
  
  accuracies_sens09 <- c()
  sensitiv_sens09 <- c()
  specific_sens09 <- c()
  ppv_sens09 <- c()
  npv_sens09 <- c()
  precision_sens09 <- c()
  recall_sens09 <- c()
  f1_sens09 <- c()
  
  accuracies_youden <- c()
  sensitiv_youden <- c()
  specific_youden <- c()
  ppv_youden <- c()
  npv_youden <- c()
  precision_youden <- c()
  recall_youden <- c()
  f1_youden <- c()
  
  conf_mat_list <- list()
  predictions_loop <- list()
  caret_conf_mat <- list()
  caret_conf_mat_sens09 <- list()
  caret_conf_mat_youden <- list()
  coefficients_loop <- list()
  
  # TEST SET -----------------------------------------------------------------
  
  for (models in seq_along(1:nrow(race_results)) ) {
    ## retrieve ID
    wflow_id <- race_results[[1]][models]
    # select best hyperparams
    best_results[[models]] <- 
      race_results %>% 
      extract_workflow_set_result(wflow_id) %>% 
      select_best(metric = "roc_auc")
    # fit to training data, and calculate on test data
    test_results[[models]] <- 
      race_results %>% 
      extract_workflow(wflow_id) %>% 
      finalize_workflow(best_results[[models]]) %>% 
      last_fit(split = dat_split)
    
    aucs <- append(aucs, collect_metrics(test_results[[models]])[[3]][2])
    accuracies <- append(accuracies, collect_metrics(test_results[[models]])[[3]][1])
    
    predictions_loop[[models]] <- test_results[[models]] %>%
      collect_predictions() 
    
    # CONF MATRIX - original model cutoff ------------------------------------
    caret_conf_mat[[models]] <-  caret::confusionMatrix(predictions_loop[[models]]$.pred_class,                
                                                        predictions_loop[[models]]$disease, 
                                                        positive = all_combis$diseases[i])
    # print(caret_conf_mat[[models]])
    #accuracies <- append(accuracies, caret_conf_mat[[models]]$overall[[1]])
    sensitiv <- append(sensitiv, caret_conf_mat[[models]][[4]][1])
    specific <- append(specific, caret_conf_mat[[models]][[4]][2])
    ppv <- append(ppv, caret_conf_mat[[models]][[4]][3])
    npv <- append(npv, caret_conf_mat[[models]][[4]][4])
    precision <- append(precision, caret_conf_mat[[models]][[4]][5])
    recall <- append(recall, caret_conf_mat[[models]][[4]][6])
    f1 <- append(f1, caret_conf_mat[[models]][[4]][7])
    
    pred_string <- paste0(".pred_", all_combis$diseases[i])
    
    ## YOUDEN INDEX -----------------------------------------------------------
    my_roc <- pROC::roc(predictions_loop[[models]]$disease, predictions_loop[[models]][[pred_string]])
    youden.threshold <- coords(my_roc, x="best", best.method="youden", ret = "threshold")[[1]]
    
    ## GOAL: Sensitivity = 90% ------------------------------------------------
    target_sensitivity <- 0.90 
    threshold_90 <- my_roc$thresholds[which.min(abs(my_roc$sensitivities - 0.9))] 
    
    ci.auc <- ci.auc(my_roc)
    ci.auc.text <- paste("AUC CI: ", round(ci.auc[1],2), "-",round(ci.auc[3],2), " DeLong" , sep = "")
    ci.auc.texts <- append(ci.auc.texts, paste("AUC CI: ", round(ci.auc[1],2), "-",round(ci.auc[3],2), " DeLong" , sep = ""))
    
    # threshold for spec & sens -----------------------------------------------
    predictions_loop[[models]] <- predictions_loop[[models]] %>% 
      mutate(
        # youden
        .pred_class_youden = ifelse(.data[[pred_string]] < youden.threshold, "control", all_combis$diseases[i]),
        .pred_class_youden = factor(.pred_class_youden, levels = levels(predictions_loop[[models]]$disease)),
        # sens 90%
        .pred_class_sens09 = ifelse(.data[[pred_string]] > threshold_90, all_combis$diseases[i], "control"),
        .pred_class_sens09 = factor(.pred_class_sens09, levels = levels(predictions_loop[[models]]$disease))
      )
    
    # CONF MATRIX - sens optimized ----
    caret_conf_mat_sens09[[models]] <-  caret::confusionMatrix(predictions_loop[[models]]$.pred_class_sens09,
                                                               predictions_loop[[models]]$disease, 
                                                               positive = all_combis$diseases[i])
    #print(caret_conf_mat_sens09[[models]])
    accuracies_sens09 <- append(accuracies_sens09, caret_conf_mat_sens09[[models]]$overall[[1]])
    sensitiv_sens09 <- append(sensitiv_sens09, caret_conf_mat_sens09[[models]][[4]][1])
    specific_sens09 <- append(specific_sens09, caret_conf_mat_sens09[[models]][[4]][2])
    ppv_sens09 <- append(ppv_sens09, caret_conf_mat_sens09[[models]][[4]][3])
    npv_sens09 <- append(npv_sens09, caret_conf_mat_sens09[[models]][[4]][4])
    precision_sens09 <- append(precision_sens09, caret_conf_mat_sens09[[models]][[4]][5])
    recall_sens09 <- append(recall_sens09, caret_conf_mat_sens09[[models]][[4]][6])
    f1_sens09 <- append(f1_sens09, caret_conf_mat_sens09[[models]][[4]][7])
    
    # CARET CONF MATRIX - youden optimized ----
    caret_conf_mat_youden[[models]] <-  caret::confusionMatrix(predictions_loop[[models]]$.pred_class_youden, 
                                                               predictions_loop[[models]]$disease, 
                                                               positive = all_combis$diseases[i])
    # print(caret_conf_mat_youden[[models]])
    accuracies_youden <- append(accuracies_youden, caret_conf_mat_youden[[models]]$overall[[1]])
    sensitiv_youden <- append(sensitiv_youden, caret_conf_mat_youden[[models]][[4]][1])
    specific_youden <- append(specific_youden, caret_conf_mat_youden[[models]][[4]][2])
    ppv_youden <- append(ppv_youden, caret_conf_mat_youden[[models]][[4]][3])
    npv_youden <- append(npv_youden, caret_conf_mat_youden[[models]][[4]][4])
    precision_youden <- append(precision_youden, caret_conf_mat_youden[[models]][[4]][5])
    recall_youden <- append(recall_youden, caret_conf_mat_youden[[models]][[4]][6])
    f1_youden <- append(f1_youden, caret_conf_mat_youden[[models]][[4]][7])
    

    
    # TEST ROC CURVE
    ci_auc_plusminus <- aucs[models] - ci.auc[1]
    subtitle <- paste(wflow_id, " - Test set", sep = "")
    roc_plot_test <- predictions_loop[[models]] %>%
      ## factor order usually alphabetically! (acs < control)
      mutate(disease = factor(disease, levels=c("control", all_combis$diseases[i]))) %>%  
      roc_curve(disease, !!sym(pred_string), event_level="second") %>%
      ggplot(aes(1 - specificity, sensitivity)) +
      geom_line(linewidth = 1.5, color = ggthemes_data$few$colors$Dark[2,2][[1]]) +
      geom_path(show.legend = FALSE) + # black line in the middle
      geom_abline(lty = 2, color = ggthemes_data$few$colors$Dark[1,2][[1]], size = 1.2, alpha=0.5) +
      coord_equal() +
      labs(title = NULL, subtitle = NULL) +
      annotate("text", x = 0.4, y=0.2, label = glue("AUC={round(aucs[models],2)}±{round(ci_auc_plusminus,2)}"), size=4, hjust = 0) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6)) +
      my_base_theme()
    print(roc_plot_test)
    
    if(SAVE.files == TRUE) {
      filename_plot_ROC_wflow_id <- 
        glue("{data_path_bestageing2022}/output/plots/roc_curves/{all_combis$diseases[i]}/004c_ROC_{wflow_id}_analysis_full_analysis_m20.svg")
      ggsave(filename = filename_plot_ROC_wflow_id, 
             plot = roc_plot_test, 
             width = 3, height = 3, 
             units = "in"
      )
    }
    
    # CALIBRATION ------------------------------------------------------------
    calibration_plot <- predictions_loop[[models]] %>% 
      cal_plot_breaks(disease, .data[[pred_string]], event_level = "second", num_breaks = 8,
                      include_points = FALSE) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
      my_base_theme()
    #print(calibration_plot)
    
    custom_cal_plots <- custom_calibration_plot(predicted_probabilities = predictions_loop[[models]][[pred_string]], 
                             observed_outcomes = predictions_loop[[models]]$disease, 
                             outcome=all_combis$diseases[i],
                             recalibration=FALSE)
    
    # adjust_legend
    custom_cal_plots$ridges_probs_plot <- custom_cal_plots$ridges_probs_plot+
      scale_fill_manual(name = NULL, labels = toupper, values = thematic::okabe_ito(6))
    
    
    ## Patchwork plot -------------------------------------------------------
    patchwork_roc_calib <- (calibration_plot | (roc_plot_test / custom_cal_plots$ridges_probs_plot)) +  # 2nd brackets needed for labels
      plot_annotation(tag_levels = 'A')
    patchwork_roc_calib <- patchwork_roc_calib + plot_layout(widths = c(2, 1))
    # patchwork_roc_calib
    
    # 4 plots
    # layout <- (calibration_plot + (roc_plot_test / custom_cal_plots$ridges_probs_plot) + calibration_plot) +
    #   plot_annotation(tag_levels = 'A')
    
    filename_plot_patchwork_ROC_calib <- 
      glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_patchwork_ROC_calib_{wflow_id}_analysis_full_analysis_m20")
    
    ggsave(filename = glue("{filename_plot_patchwork_ROC_calib}.svg"), 
           plot = patchwork_roc_calib, 
           width = 12, height = 8, 
           units = "in"
    )
    saveRDS(object = patchwork_roc_calib, file = glue("{filename_plot_patchwork_ROC_calib}.rds"))  # store for later patchwork
    
    # save each component to modify later ;)
    # 1) calib
    saveRDS(object = calibration_plot, file = glue("{filename_plot_patchwork_ROC_calib}_01_calib_only.rds"))
    # 2) roc
    saveRDS(object = roc_plot_test, file = glue("{filename_plot_patchwork_ROC_calib}_02_roc_only.rds"))
    # 3) ridges probs
    saveRDS(object = custom_cal_plots$ridges_probs_plot, file = glue("{filename_plot_patchwork_ROC_calib}_03_ridges_probs_only.rds"))
    
    # next model
  }
  
  # summarize all models --------------------------------------
  performance.summary.table.sens09 <- tibble(
    model = race_results$wflow_id,
    AUC = round(aucs,2), 
    accuracy = round(accuracies_sens09,2),
    sensitivity = round(sensitiv_sens09, 2),
    specificity = round(specific_sens09, 2),
    ppv = round(ppv_sens09,2),
    npv= round(npv_sens09,2), 
    precision= round(precision_sens09,2), 
    recall= round(recall_sens09,2), 
    f1.score =round(f1_sens09,2)
  ) %>% 
    arrange(desc(AUC))
  
  filename_df_performance_summary_sens09 <- 
    glue("{data_path_bestageing2022}/output/performance_summary_df/{all_combis$diseases[i]}/004c_performance_summary_table_analysis_full_analysis_m20.rds")
  saveRDS(object = performance.summary.table.sens09, file = filename_df_performance_summary_sens09)
  
  # next disease
  message(glue("---------Run finished for disease: {toupper(all_combis$diseases[i])}---------------"))
  
}



# Feature Importance ---------------------------------------------------------

for (i in 1:nrow(all_combis)){
  # need to get data split right to finalize model
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_{all_combis$diseases[i]}_data01.rds")
  if(!file.exists(path2dataprocessed)) {
    next
  }
  
  data01 <- readRDS(file = path2dataprocessed)
  no.mirnas <- ncol(data01) - 4
  text_disease <- stringr::str_to_upper(all_combis$diseases[i])
  set.seed(123) # get same indices 
  modeldat <- data01
  
  # make control first factor for all analyses ;)
  modeldat <- modeldat %>% 
    mutate(disease = factor(disease, levels = c("control", all_combis$diseases[i])))
  
  dat_split <- rsample::initial_split(modeldat, strata = disease)
  dat_train <- training(dat_split)
  dat_test <- testing(dat_split)
  
  folds <- vfold_cv(dat_train, strata = disease, v = no_folds, repeats = no_repeats)
  
  n_feature_select <- 20  # how many mirna select in training?
  normalized_rec <- 
    recipe(disease ~ ., data = dat_train) %>%
    ### https://recipes.tidymodels.org/articles/Ordering.html
    #  To make sure we don’t get any unexpected results, it’s best to use 
    #  the following ordering of high-level transformations:
    #     Skewness Transformations - step_YeoJohnson()
    #     Centering, Scaling, or Normalization on Numeric Predictors
    #     Dummy Variables for Categorical Data
    update_role(pat_id, new_role="ID") %>% 
    step_zv(all_predictors()) %>%  
    step_impute_mean(all_numeric_predictors()) %>% 
    step_impute_mode(all_nominal_predictors(), -disease) %>% 
    step_corr(all_numeric_predictors(), threshold = 0.9) %>% 
    step_YeoJohnson() %>% 
    step_normalize(all_numeric_predictors()) %>%  
    step_dummy(all_nominal_predictors(),-disease) %>% 
    step_select_forests(all_predictors(), -c(age, sex_Male), outcome = "disease", top_p = n_feature_select)
  # B
  poly_rec <- 
    recipe(disease ~ ., data = dat_train) %>%
    update_role(pat_id, new_role="ID") %>% 
    step_zv(all_predictors()) %>%  
    step_impute_mean(all_numeric_predictors()) %>% 
    step_impute_mode(all_nominal_predictors(), -disease) %>% 
    step_corr(all_numeric_predictors(), threshold = 0.9) %>% 
    step_normalize(all_numeric_predictors()) %>%  
    step_poly(all_numeric_predictors()) %>% 
    step_dummy(all_nominal_predictors(),-disease) %>% 
    step_select_forests(all_predictors(), -c(starts_with("age"), starts_with("sex")), outcome = "disease", top_p = n_feature_select)
  #step_interact( ~all_predictors():all_predictors())
  
  # C
  simple_rec <- 
    recipe(disease ~ ., data = dat_train) %>%
    update_role(pat_id, new_role="ID") %>% 
    # ZERO VARIANCE
    step_zv(all_predictors()) %>%  
    # IMPUTE
    step_impute_mean(all_numeric_predictors()) %>% 
    step_impute_mode(all_nominal_predictors(), -disease) %>% 
    # DECORRELATE
    step_corr(all_numeric_predictors(), threshold = 0.9) %>% 
    step_dummy(all_nominal_predictors(),-disease) %>% 
    step_select_forests(all_predictors(), -c(age, sex_Male), outcome = "disease", top_p = n_feature_select)
  
  # more logit
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_randomMIR_FALSE_more_logit.rds"))
  
  ### filter first, vip takes long, and we do not like all models ;) ------------
  race_results <- race_results %>% 
    filter(!wflow_id %in% c("full_quad_logistic_reg")) %>% 
    filter(!stringr::str_detect(wflow_id, "KNN"))
  
  # common selected vars from models
  commmon_vars_list <- list()
  
  # now vip for all models
  
  best_results <- list()  # hyperparams selected
  test_results <- list()  # last fit
  
  for (models in seq_along(1:nrow(race_results)) ) {
    ## retrieve ID
    wflow_id <- race_results[[1]][models]  # wflow_id column
    # select best hyperparams
    best_results[[models]] <- 
      race_results %>% 
      extract_workflow_set_result(wflow_id) %>% 
      select_best(metric = "roc_auc")
    # fit to training data, and calculate on test data
    test_results[[models]] <- 
      race_results %>% 
      extract_workflow(wflow_id) %>% 
      finalize_workflow(best_results[[models]]) %>% 
      last_fit(split = dat_split)
    
    last_fit_metrics <- collect_metrics(test_results[[models]])
    # extract model
    last_fit_model <- extract_workflow(test_results[[models]])
    
    ## Global Explanations --------------------------------------------------------
    # We compute variable importance by permutating features. If shuffling a column causes a large degradation in model performance, it is important and vica versa. This approach is model agnostic. 
    
    model_vars <- simple_rec %>% summary()
    
    if(str_detect(wflow_id, "full_quad")) {
      model_vars <- poly_rec %>% summary()
    }
    
    index <- #model_vars$variable != "id" &    # needed for global importance since in recipe
      model_vars$variable != "disease" #& (model_vars$role == "predictor" | model_vars$role == "outcome")
    # modify data
    predictor_vars <- model_vars[["variable"]][index]
    vip_train <- dat_train %>% select(all_of(predictor_vars))  # although we used feature selection, all vars must be present.. takes time
    
    # create explainers | top model
    explainer1 <- 
      explain_tidymodels(
        last_fit_model, 
        data = vip_train, #eval(as.symbol(recipe_models$recipe[models])) %>% prep() %>% bake(new_data = NULL, all_predictors(), all_outcomes()),  #vip_train, 
        y = as.numeric(dat_train$disease), 
        label = wflow_id,
        verbose=TRUE
      )
    
    # variable importance takes time
    vip1 <- model_parts(explainer = explainer1, 
                        loss_function = loss_default("classification")    #  # loss_default(explainer_xgb$model_info$type) "One minus AUC"
    )
    
    saveRDS(vip1, file = glue("{data_path_bestageing2022}/output/feature_importance/{all_combis$diseases[i]}/004c_full_{wflow_id}_vip1.rds"))
    
    # cave used feature selection, but cannot supply less vars to explain tidymodels
    predictor_vars_tmp <- extract_recipe(last_fit_model) %>% 
      summary() %>% 
      filter(str_detect(variable, "hsa_mir")) %>% 
      pull(variable) %>% 
      c("pat_id", "age", "sex")
    
    vip1_tmp <- vip1 %>% 
      filter(variable %in% predictor_vars_tmp | variable == "_full_model_" | variable == "_baseline_")
    
    metric_name <- attr(vip1_tmp, "loss_name")
    
    # 
    features_importance_table <- vip1_tmp %>% 
      group_by(variable) %>% 
      filter(variable != "_baseline_" & variable != "pat_id") %>% 
      summarize(dropout_loss = mean(dropout_loss)) %>% 
      arrange((dropout_loss))
    
    saveRDS(features_importance_table, file = glue("{data_path_bestageing2022}/output/feature_importance/{all_combis$diseases[i]}/004c_full_{wflow_id}_features_importance_table.rds"))
    
    common_vars_new <- setNames(list(features_importance_table$variable), wflow_id)
    
    commmon_vars_list <- c(commmon_vars_list, common_vars_new)
    
    variable_imp_plot <- ggplot_imp(vip1_tmp)
    variable_imp_plot
    
    vip_plot_filename <- glue("{data_path_bestageing2022}/output/plots/feature_importance/{all_combis$diseases[i]}/004c_{wflow_id}_variable_imp_plot")
    ggsave(
      filename = glue("{vip_plot_filename}.svg"), plot = variable_imp_plot, 
      width = 6, height = 6, 
      units = "in"  # default
    )
    
    saveRDS(object = variable_imp_plot, file = glue("{vip_plot_filename}.rds"))
    
    message(glue("wflow_id: {wflow_id} for {toupper(all_combis$diseases[i])} done.\ncont..."))
    # next model
  }
  
  ## Venn | Common Vars in final models ----------------------------------------------
  # selected features calculating intersections
  # remove some vars we are not interested
  commmon_vars_list_modified <- commmon_vars_list[names(commmon_vars_list) %in% c("XGB", "RF", "SVM_linear", "SVM_poly") | str_detect(names(commmon_vars_list), pattern= "logistic_reg")]
  
  commmon_vars_list_modified <- lapply(commmon_vars_list_modified, function(x) {
    x[!x %in% c("age", "_full_model_", "sex")]
  })
  
  list_names <- names(commmon_vars_list_modified)
  
  # Function to find intersection for a given combination of lists
  get_intersection <- function(combination) {
    lists_to_intersect <- commmon_vars_list_modified[combination]
    Reduce(intersect, lists_to_intersect)
  }
  
  # Generate all combinations of list names and find intersections
  list_of_intersections <- lapply(1:length(list_names), function(i) {
    combn(list_names, i, function(combination) {
      setNames(list(get_intersection(combination)), paste(combination, collapse = ":"))
    }, simplify = FALSE)
  }) %>% unlist(recursive = FALSE)
  
  # Filter out empty lists and check for colons (we are only interested in intersections)
  filtered_list <- Filter(function(x) length(x[[1]]) > 0, list_of_intersections)
  filtered_list <- Filter(function(x) any(str_detect(names(x), pattern = ":")), filtered_list)
  
  # Flatten the list and convert to data frame
  flattened_list <- do.call(c, filtered_list)  # concatenate function
  df_intersections <- data.frame(
    ListName = names(flattened_list),
    Elements = I(flattened_list)  # I() is used to prevent the list from being simplified into individual columns
  )
  df_intersections$Elements <- sapply(df_intersections$Elements, function(x) paste(x, collapse = ", "))
  df_intersections$Count <- sapply(flattened_list, length)
  df_intersections <- df_intersections %>% 
    as_tibble() %>% 
    arrange(desc(Count))
  df_intersections
  saveRDS(object = df_intersections, file = glue("{data_path_bestageing2022}/output/feature_importance/{all_combis$diseases[i]}/004c_full_m20_selected_vars_intersection.rds"))
  
  # docx
  file_path <- glue("{data_path_bestageing2022}/output/feature_importance/{all_combis$diseases[i]}/004c_full_m20_selected_vars_intersection.docx")
  df_intersections_flextable <- flextable(df_intersections) %>% 
    colformat_double(
      big.mark = ",", digits = 3, na_str = "N/A"
    ) %>% 
    flextable::set_table_properties(layout = "autofit") %>%
    fontsize(size = 8, part = "all") %>% 
    flextable::font(fontname = "Times New Roman", part = "all") %>%
    autofit() %>% 
    save_as_docx(path=file_path, pr_section = sect_properties)
  
  # # previous approach # only until 9 intersection elements, but these also would be huge lists
  # intersections <- gplots::venn(commmon_vars_list[1:3])
  # list_of_intersections <- attr(intersections, "intersections")
  # 
  # # Displaying the intersections
  # intersections_diseases <- tibble(
  #   Intersection = character(),
  #   Variables = character(), 
  #   count_intersections = integer()
  # )
  # for (name in seq_along(list_of_intersections)) {
  #   intersection_name <- names(list_of_intersections[name])
  #   variables_joined <- paste(gsub("_", "-", list_of_intersections[[name]]), collapse = ", ")
  #   no_variables <- length(list_of_intersections[[name]])
  #   
  #   cat(intersection_name, ":\n", variables_joined, "\n\n")
  #   
  #   intersections_diseases <- intersections_diseases %>% 
  #     add_row(Intersection = intersection_name, Variables = variables_joined, count_intersections = no_variables)
  # }
  

  #write.csv2(x = intersections_diseases, file=glue("{data_path_bestageing2022}/output/feature_importance/{all_combis$diseases[i]}/004c_full_m20_selected_vars_intersection.csv"))
  
  # intersections_diseases %>% 
  #   kableExtra::kable(digits = 3, caption = "Intersection of variables in the top 10 VIPs of each model",) %>%
  #   kableExtra::kable_styling(font_size=12) %>% 
  #   kableExtra::kable_classic(full_width = T) %>% 
  #   kableExtra::save_kable()

  
  # next disease
  message(glue("-------------Run finished for disease: {toupper(all_combis$diseases[i])}-------------------"))
}


