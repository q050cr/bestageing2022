





diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="full")

for (i in 1:nrow(all_combis)){
  race_results <- readRDS(glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_10_folds_5_{toupper(all_combis$diseases[i])}_analysis_randomMIR_FALSE.rds"))
  
  
  
  
  
}