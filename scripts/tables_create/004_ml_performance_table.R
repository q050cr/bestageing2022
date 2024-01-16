

# performance bar plot

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(gt)

# see also /Volumes/T7CR/data/observeACS/scripts/figures_create/performance_ml_models.R
# Assuming you have these data frames:
# performance_table_HFREF, performance_table_ACS, performance_table_Disease3, performance_table_Disease4

system_name <- Sys.info()["nodename"]
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



# load data ---------------------------------------------------------------

# load data created in script "main.Rmd"
performance_table_ACS_XGB_full <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/acs/models_full/004c_XGB_epiR_crosstabs.csv"))
performance_table_CAD_XGB_full <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/cad/models_full/004c_XGB_epiR_crosstabs.csv"))
performance_table_DCM_XGB_full <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/dcm/models_full/004c_XGB_epiR_crosstabs.csv"))
performance_table_HFREF_XGB_full <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/hfref/models_full/004c_XGB_epiR_crosstabs.csv"))

performance_table_ACS_XGB_selected <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/acs/models_selected/004c_XGB_epiR_crosstabs.csv"))
performance_table_CAD_XGB_selected <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/cad/models_selected/004c_XGB_epiR_crosstabs.csv"))
performance_table_DCM_XGB_selected <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/dcm/models_selected/004c_XGB_epiR_crosstabs.csv"))
performance_table_HFREF_XGB_selected <- read.csv2(file = glue("{data_path_bestageing2022}/output/performance_summary_df/hfref/models_selected/004c_XGB_epiR_crosstabs.csv"))


# create for only full model -------------------------------------------------

# Step 1: Combine the performance tables
all_performance <- bind_rows(
  performance_table_ACS_XGB_full %>% mutate(disease = "ACS"),
  performance_table_CAD_XGB_full %>% mutate(disease = "CAD"),
  performance_table_DCM_XGB_full %>% mutate(disease = "DCM"),
  performance_table_HFREF_XGB_full %>% mutate(disease = "HFREF"),
) %>%
  select(-X) |> 
  filter(statistic %in% c("sp", "se", "pv.pos", "pv.neg", "auc")) |> 
  mutate(statistic = case_when(
    statistic == "sp"     ~ "Specificity",
    statistic == "se"     ~ "Sensitivity",
    statistic == "pv.pos" ~ "PPV",
    statistic == "pv.neg" ~ "NPV",
    statistic == "auc"    ~ "AUC",
    TRUE                  ~ statistic
  )) 

# Create the GT table
gt_table <- all_performance %>%
  group_by(disease) %>%
  gt(groupname_col = "disease") %>%
  #tab_row_group(groupvar = disease) %>%
  cols_label( statistic = "Statistic", est = "Estimate", lower = "Lower CI", upper = "Upper CI") |> 
  fmt_number(columns = c(est, lower, upper), decimals = 3)

gt_table


# combine full and selected model stats -----------------------------------

all_performance_full <- bind_rows(
  performance_table_ACS_XGB_full %>% mutate(disease = "ACS"),
  performance_table_CAD_XGB_full %>% mutate(disease = "CAD"),
  performance_table_DCM_XGB_full %>% mutate(disease = "DCM"),
  performance_table_HFREF_XGB_full %>% mutate(disease = "HFREF"),
) %>%
  select(-X) |> 
  filter(statistic %in% c("sp", "se", "pv.pos", "pv.neg", "auc")) |> 
  mutate(statistic = case_when(
    statistic == "sp"     ~ "Specificity",
    statistic == "se"     ~ "Sensitivity",
    statistic == "pv.pos" ~ "PPV",
    statistic == "pv.neg" ~ "NPV",
    statistic == "auc"    ~ "AUC",
    TRUE                  ~ statistic
  )) |> 
  rename(est_full = est, lower_full = lower , upper_full = upper)

all_performance_selected <- bind_rows(
  performance_table_ACS_XGB_selected %>% mutate(disease = "ACS"),
  performance_table_CAD_XGB_selected %>% mutate(disease = "CAD"),
  performance_table_DCM_XGB_selected %>% mutate(disease = "DCM"),
  performance_table_HFREF_XGB_selected %>% mutate(disease = "HFREF"),
) %>%
  select(-X) |> 
  filter(statistic %in% c("sp", "se", "pv.pos", "pv.neg", "auc")) |> 
  mutate(statistic = case_when(
    statistic == "sp"     ~ "Specificity",
    statistic == "se"     ~ "Sensitivity",
    statistic == "pv.pos" ~ "PPV",
    statistic == "pv.neg" ~ "NPV",
    statistic == "auc"    ~ "AUC",
    TRUE                  ~ statistic
  )) |> 
  rename(est_selected = est, lower_selected = lower , upper_selected = upper)

all_performance_combined <- left_join(all_performance_full, all_performance_selected, by = c("disease", "statistic"))

# Create the GT table
gt_table_combined <- all_performance_combined %>%
  group_by(disease) %>%
  gt(groupname_col = "disease") %>%
  #tab_row_group(groupvar = disease) %>%
  cols_label( 
    statistic = "Statistic", est_full = "Estimate", lower_full = "Lower CI", upper_full = "Upper CI",
    est_selected = "Estimate", lower_selected = "Lower CI", upper_selected = "Upper CI"
  ) |> 
  tab_spanner(label = "Full model", columns = c("est_full", "lower_full", "upper_full")) |>
  tab_spanner(label = "Selected model", columns = c("est_selected", "lower_selected", "upper_selected")) |>
  fmt_number(columns = c("est_full", "lower_full", "upper_full", "est_selected", "lower_selected", "upper_selected"), decimals = 3)

gt_table_combined

# merge cols
gt_table_combined <- all_performance_combined %>%
  group_by(disease) %>%
  gt(groupname_col = "disease") %>%
  fmt_number(columns = c("est_full", "lower_full", "upper_full", "est_selected", "lower_selected", "upper_selected"), decimals = 3) |> 
  cols_merge(
    columns = c("est_full", "lower_full", "upper_full"),
    pattern = "{1} ({2}; {3})",
  ) |> 
  cols_merge(
    columns = c("est_selected", "lower_selected", "upper_selected"),
    pattern = "{1} ({2}; {3})",
  ) |> 
  cols_label( 
    statistic = "Statistic", 
    est_full = "Feature selection", 
    est_selected = "A-priori selected features"
  )

gt_table_combined

gt_table_combined %>% 
  gtsave(glue("{data_path_bestageing2022}/output/tables/ml_performance_table/table_ml_performance_summary.docx"))
