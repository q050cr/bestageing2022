




# ml performance summary plot

library(dplyr)
library(tidyr)
library(glue)
library(ggplot2)
library(patchwork)

# see also /Volumes/T7CR/data/observeACS/scripts/figures_create/performance_ml_models.R
# Assuming you have these data frames:
# performance_table_HFREF, performance_table_ACS, performance_table_Disease3, performance_table_Disease4

mount_filesystem <- TRUE

if (system_name == "MacBook-Pro-CR-2065.local" | stringr::str_detect(string = system_name, "laptop-zim.uni-heidelberg.de")) {
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

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_performance_summary_barplot.R"))
source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))


# A) selected analysis -------------------------------------------------------

diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="full")

for (i in 1:nrow(all_combis) ) {
  
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_20231122_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_SELECTED_miRetrieve_TRUE_randomMIR_FALSE_more_logit.rds"))
  # not included in calculating vip
  race_results <- race_results %>% 
    filter(!wflow_id %in% c("full_quad_logistic_reg")) %>% 
    filter(!stringr::str_detect(wflow_id, "KNN"))
  
  for (models in seq_along(1:nrow(race_results)) ) {
    wflow_id <- race_results[[1]][models]

    # load the 3 components created in script 004c_{filename}.R
    filename_plot_patchwork_ROC_calib <- 
      glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_patchwork_ROC_calib_{wflow_id}_analysis_selected_analysis")
    # 1) calib
    calibration_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_01_calib_only.rds"))
    # 2) roc
    roc_plot_test <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_02_roc_only.rds"))
    # 3) ridges probs
    ridges_probs_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_03_ridges_probs_only.rds")) 
    ridges_probs_plot <- ridges_probs_plot +theme(legend.position = "bottom")
    # load vip plot created later in script 004c_ml_test_set_full_analysis.R
    variable_imp_plot <- readRDS(glue("{data_path_bestageing2022}/output/plots/feature_importance/{all_combis$diseases[i]}/004c_{wflow_id}_variable_imp_plot_selectedmiRetrieve.rds"))
    
    # create performance plot
    # load data created in script  004c_ml_test_set_full_analysis.R
    performance.summary.table.sens09 <- readRDS(glue("{data_path_bestageing2022}/output/performance_summary_df/{all_combis$diseases[i]}/004c_performance_summary_table_analysis_selected_analysis.rds")) %>% 
      # filter for model in loop
      filter(model %in% wflow_id)
    
    performance_barplot <- my_performance_plot(mydata = performance.summary.table.sens09)
    
    # Arrange the plots
    combined_plot <- (roc_plot_test / performance_barplot) | (calibration_plot / ridges_probs_plot) |  variable_imp_plot
    combined_plot <- combined_plot + plot_layout(ncol = 3, nrow = 1, heights = c(1, 1), widths = c(1,1,2)) + plot_annotation(tag_levels = 'A')
    combined_plot
    
    filename_plot_patchwork_ml_summary_plot <- 
      glue("{data_path_bestageing2022}/output/plots/ml_summary_plot_patchwork/{all_combis$diseases[i]}/004c_patchwork_ml_summary_plot_{wflow_id}_analysis_selected")
    
    ggsave(filename = glue("{filename_plot_patchwork_ml_summary_plot}.svg"), 
           plot = combined_plot, 
           width = 14*1.25, height = 8*1.25, 
           units = "in"
    )
    
  }
  # next disease
  message(glue("-------------Run finished for disease: {toupper(all_combis$diseases[i])}-------------------"))
}


# B) full analysis not matched -----------------------------------------------


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="full")

for (i in 1:nrow(all_combis) ) {
  
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_randomMIR_FALSE_more_logit.rds"))
  
  # not included in calculating vip
  race_results <- race_results %>% 
    filter(!wflow_id %in% c("full_quad_logistic_reg")) %>% 
    filter(!stringr::str_detect(wflow_id, "KNN"))
  
  for (models in seq_along(1:nrow(race_results)) ) {
    wflow_id <- race_results[[1]][models]
    patchwork_roc_calib <- readRDS(glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_patchwork_ROC_calib_{wflow_id}_analysis_full_analysis_m20.rds"))
    
    # load the 3 components created in script 004c_ml_test_set_full_analysis.R
    filename_plot_patchwork_ROC_calib <- 
      glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_patchwork_ROC_calib_{wflow_id}_analysis_full_analysis_m20")
    # 1) calib
    calibration_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_01_calib_only.rds"))
    # 2) roc
    roc_plot_test <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_02_roc_only.rds"))
    # 3) ridges probs
    ridges_probs_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_03_ridges_probs_only.rds")) 
    ridges_probs_plot <- ridges_probs_plot +theme(legend.position = "bottom")
    # load vip plot created later in script 004c_ml_test_set_full_analysis.R
    variable_imp_plot <- readRDS(glue("{data_path_bestageing2022}/output/plots/feature_importance/{all_combis$diseases[i]}/004c_{wflow_id}_variable_imp_plot.rds"))
    
    # create performance plot
    # load data created in script  004c_ml_test_set_full_analysis.R
    performance.summary.table.sens09 <- readRDS(glue("{data_path_bestageing2022}/output/performance_summary_df/{all_combis$diseases[i]}/004c_performance_summary_table_analysis_full_analysis_m20.rds")) %>% 
      # filter for model in loop
      filter(model %in% wflow_id)
    
    performance_barplot <- my_performance_plot(mydata = performance.summary.table.sens09)
    
    # Arrange the plots
    combined_plot <- (roc_plot_test / performance_barplot) | (calibration_plot / ridges_probs_plot) |  variable_imp_plot
    combined_plot <- combined_plot + plot_layout(ncol = 3, nrow = 1, heights = c(1, 1), widths = c(1,1,2)) + plot_annotation(tag_levels = 'A')
    combined_plot
    
    filename_plot_patchwork_ml_summary_plot <- 
      glue("{data_path_bestageing2022}/output/plots/ml_summary_plot_patchwork/{all_combis$diseases[i]}/004c_patchwork_ml_summary_plot_{wflow_id}_analysis_full_analysis_m20")
    
    ggsave(filename = glue("{filename_plot_patchwork_ml_summary_plot}.svg"), 
           plot = combined_plot, 
           width = 14, height = 8, 
           units = "in"
    )
    
  }
  # next disease
  message(glue("-------------Run finished for disease: {toupper(all_combis$diseases[i])}-------------------"))
}



# C) full analysis matched -----------------------------------------------

for (i in 1:nrow(all_combis) ) {
  
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_randomMIR_FALSE_more_logit_matchIt.rds"))

  # not included in calculating vip
  race_results <- race_results %>% 
    filter(!wflow_id %in% c("full_quad_logistic_reg")) %>% 
    filter(!stringr::str_detect(wflow_id, "KNN"))
  
  for (models in seq_along(1:nrow(race_results)) ) {
    wflow_id <- race_results[[1]][models]
    patchwork_roc_calib <- readRDS(glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_patchwork_ROC_calib_{wflow_id}_analysis_full_analysis_m20_matchIt.rds"))
    
    # load the 3 components created in script 004c_ml_test_set_full_analysis.R
    filename_plot_patchwork_ROC_calib <- 
      glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_patchwork_ROC_calib_{wflow_id}_analysis_full_analysis_m20_matchIt")
    # 1) calib
    calibration_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_01_calib_only.rds"))
    # 2) roc
    roc_plot_test <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_02_roc_only.rds"))
    # 3) ridges probs
    ridges_probs_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_03_ridges_probs_only.rds")) 
    ridges_probs_plot <- ridges_probs_plot +theme(legend.position = "bottom")
    # load vip plot created later in script 004c_ml_test_set_full_analysis.R
    variable_imp_plot <- readRDS(glue("{data_path_bestageing2022}/output/plots/feature_importance/{all_combis$diseases[i]}/004c_{wflow_id}_variable_imp_plot_matchIt.rds"))
    
    # create performance plot
    # load data created in script  004c_ml_test_set_full_analysis.R
    performance.summary.table.sens09 <- readRDS(glue("{data_path_bestageing2022}/output/performance_summary_df/{all_combis$diseases[i]}/004c_performance_summary_table_analysis_full_analysis_m20_matchIt.rds")) %>% 
      # filter for model in loop
      filter(model %in% wflow_id)
    
    performance_barplot <- my_performance_plot(mydata = performance.summary.table.sens09)
    
    # Arrange the plots
    combined_plot <- (roc_plot_test / performance_barplot) | (calibration_plot / ridges_probs_plot) |  variable_imp_plot
    combined_plot <- combined_plot + plot_layout(ncol = 3, nrow = 1, heights = c(1, 1), widths = c(1,1,2)) + plot_annotation(tag_levels = 'A')
    combined_plot
    
    filename_plot_patchwork_ml_summary_plot <- 
      glue("{data_path_bestageing2022}/output/plots/ml_summary_plot_patchwork/{all_combis$diseases[i]}/004c_20240121_patchwork_ml_summary_plot_{wflow_id}_analysis_full_analysis_m20_matchIt")
    
    ggsave(filename = glue("{filename_plot_patchwork_ml_summary_plot}.svg"), 
           plot = combined_plot, 
           width = 14, height = 8, 
           units = "in"
    )
    
  }
  # next disease
  message(glue("-------------Run finished for disease matched analysis: {toupper(all_combis$diseases[i])}-------------------"))
}


# D) selected analysis matched -------------------------------------------------------

diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="selected")

for (i in 1:nrow(all_combis) ) {
  # updated 20240125
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_20240125_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_SELECTED_miRetrieve_TRUE_randomMIR_FALSE_more_logit_matchIt.rds"))
  # not included in calculating vip
  race_results <- race_results %>% 
    filter(!wflow_id %in% c("full_quad_logistic_reg")) %>% 
    filter(!stringr::str_detect(wflow_id, "KNN"))
  
  for (models in seq_along(1:nrow(race_results)) ) {
    wflow_id <- race_results[[1]][models]
    
    # load the 3 components created in script 004c_{filename}.R
    filename_plot_patchwork_ROC_calib <- 
      glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_20240121_patchwork_ROC_calib_{wflow_id}_analysis_selected_analysis_matchIt")
    # 1) calib
    calibration_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_01_calib_only.rds"))
    # 2) roc
    roc_plot_test <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_02_roc_only.rds"))
    # 3) ridges probs
    ridges_probs_plot <- readRDS(file = glue("{filename_plot_patchwork_ROC_calib}_03_ridges_probs_only.rds")) 
    ridges_probs_plot <- ridges_probs_plot +theme(legend.position = "bottom")
    # load vip plot created later in script 004c_ml_test_set_full_analysis.R
    variable_imp_plot <- readRDS(glue("{data_path_bestageing2022}/output/plots/feature_importance/{all_combis$diseases[i]}/004c_20240121_{wflow_id}_variable_imp_plot_selectedmiRetrieve_matchIt.rds"))
    
    # create performance plot
    # load data created in script  004c_ml_test_set_full_analysis.R
    performance.summary.table.sens09 <- readRDS(glue("{data_path_bestageing2022}/output/performance_summary_df/{all_combis$diseases[i]}/004c_20240121_performance_summary_table_analysis_selected_analysis_matchIt.rds")) %>% 
      # filter for model in loop
      filter(model %in% wflow_id)
    
    performance_barplot <- my_performance_plot(mydata = performance.summary.table.sens09)
    
    # Arrange the plots
    combined_plot <- (roc_plot_test / performance_barplot) | (calibration_plot / ridges_probs_plot) |  variable_imp_plot
    combined_plot <- combined_plot + plot_layout(ncol = 3, nrow = 1, heights = c(1, 1), widths = c(1,1,2)) + plot_annotation(tag_levels = 'A')
    combined_plot
    
    filename_plot_patchwork_ml_summary_plot <- 
      glue("{data_path_bestageing2022}/output/plots/ml_summary_plot_patchwork/{all_combis$diseases[i]}/004c_20240121_patchwork_ml_summary_plot_{wflow_id}_analysis_selected_matchIt")
    
    ggsave(filename = glue("{filename_plot_patchwork_ml_summary_plot}.svg"), 
           plot = combined_plot, 
           width = 14*1.25, height = 8*1.25, 
           units = "in"
    )
    
  }
  # next disease
  message(glue("-------------Run finished for disease: {toupper(all_combis$diseases[i])}-------------------"))
}


