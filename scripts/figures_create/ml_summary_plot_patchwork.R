




# ml performance summary plot

library(dplyr)
library(tidyr)
library(glue)
library(ggplot2)
library(patchwork)

# see also /Volumes/T7CR/data/observeACS/scripts/figures_create/performance_ml_models.R
# Assuming you have these data frames:
# performance_table_HFREF, performance_table_ACS, performance_table_Disease3, performance_table_Disease4

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


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="full")

for (i in 1:nrow(all_combis) ) {
  
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_randomMIR_FALSE_more_logit.rds"))
  
  for (models in seq_along(1:nrow(race_results)) ) {
    wflow_id <- race_results[[1]][models]
    patchwork_roc_calib <- readRDS(glue("{data_path_bestageing2022}/output/plots/ml_roc_calib_patchwork/{all_combis$diseases[i]}/004c_patchwork_ROC_calib_{wflow_id}_analysis_full_analysis_m20.rds"))
    
    variable_imp_plot <- readRDS(glue("{data_path_bestageing2022}/output/plots/feature_importance/{all_combis$diseases[i]}/004c_{wflow_id}_variable_imp_plot.rds"))
    
    # add plot to existing patchwork plot
    new_combined_plot <- (patchwork_roc_calib | variable_imp_plot) +
      plot_layout(ncol=3, guides = "collect") +  # guides = "collect" is important to have only one legend 
      plot_annotation(tag_levels = 'A') +
      theme(legend.position = "bottom")
    new_combined_plot
  }
  
  patchwork_roc_calibtest <- patchwork_roc_calib + theme(legend.position = "bottom")
  patchwork_roc_calibtest
  
}




ggsave(filename = glue("{data_path_bestageing2022}/output/plots/IFL_poster2023/performance_sum_plot/{Sys.Date()}_performance_summary_barplot.svg"), 
       plot = performance_barplot, 
       width = 6, height = 3, 
       units = "in"
)
