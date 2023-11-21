

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_calibration_plot.R"))
source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))

save.FILE=TRUE


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="full")

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
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_randomMIR_FALSE.rds"))
  
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
    
    
    ## factor order usually alphabetically! (acs < control)
    predictions_loop[[models]] <- predictions_loop[[models]] %>% 
      mutate(.pred_class_youden = factor(.pred_class_youden, levels=c("control", params$disease)),
             disease = factor(disease, levels=c("control", params$disease)),
      )   
    
    gbm::calibrate.plot(y = predictions_loop[[models]]$disease, p = predictions_loop[[models]][[pred_string]], 
                        main=paste(wflow_id, "-Calibration Plot | ", sep=""), 
                        line.par = list(col = "black"), shade.col="grey", xlim=c(0,0.6))
  }
  
}


