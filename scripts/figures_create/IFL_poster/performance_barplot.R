

# IFL performance bar plot

library(dplyr)
library(tidyr)
library(ggplot2)

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

# load data created in script "main.Rmd"
performance_table_ACS <- readRDS(file = glue("{data_path_bestageing2022}/output/plots/IFL_poster2023/performance_sum_plot/data/performance_table_ACS.rds"))
performance_table_CAD <- readRDS(file = glue("{data_path_bestageing2022}/output/plots/IFL_poster2023/performance_sum_plot/data/performance_table_CAD.rds"))
performance_table_DCM <- readRDS(file = glue("{data_path_bestageing2022}/output/plots/IFL_poster2023/performance_sum_plot/data/performance_table_DCM.rds"))
performance_table_HFREF <- readRDS(file = glue("{data_path_bestageing2022}/output/plots/IFL_poster2023/performance_sum_plot/data/performance_table_HFREF.rds"))


# Step 1: Combine the performance tables
all_performance <- bind_rows(
  performance_table_ACS %>% mutate(disease = "ACS"),
  performance_table_CAD %>% mutate(disease = "CAD"),
  performance_table_DCM %>% mutate(disease = "DCM"),
  performance_table_HFREF %>% mutate(disease = "HFREF"),
)

# Step 2: Filter top 3 models by AUC for each disease
top_models <- all_performance %>%
  group_by(disease) %>%
  arrange(desc(AUC)) %>%
  slice_head(n = 3) %>%
  ungroup()

# Step 3: Reshape the data to long format
long_data <- top_models %>%
  select(disease, model, AUC, accuracy, f1.score) %>%
  pivot_longer(cols = c(AUC, accuracy, f1.score), names_to = "metric", values_to = "value") %>% 
  mutate(metric = ifelse(metric == "f1.score", "F1-score", metric),
         metric = ifelse(metric == "accuracy", "Accuracy", metric)) %>% 
  mutate(disease = factor(disease, levels = rev(c("ACS", "CAD", "DCM", "HFREF"))))


# Step 4: Plot with ggplot2
performance_barplot <- ggplot(long_data, aes(x = disease, y = value, fill = metric,
                                             #group = model)
)) +
  geom_bar(stat = "identity", position = "dodge", alpha=0.5) +
  #facet_wrap(~ metric, scales = "free_y") + # Use facet_wrap just for y-axis variations
  theme_minimal() +
  ylab("Value") +
  xlab("Disease") +
  labs(fill = "Model") +
  labs(fill=NULL)+
  xlab(NULL)+
  ylab(NULL)+
  theme(legend.position = "top")+
  coord_flip()+
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) + # Set the desired breaks for the y-axis
  my_base_theme()
performance_barplot

ggsave(filename = glue("{data_path_bestageing2022}/output/plots/IFL_poster2023/performance_sum_plot/{Sys.Date()}_performance_summary_barplot.svg"), 
       plot = performance_barplot, 
       width = 6, height = 3, 
       units = "in"
)


# paper plot selected --------------------------------------------------------------

performance_table_ACS_selected <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/acs/004c_performance_summary_table_analysis_selected_analysis.rds"))
performance_table_CAD_selected <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/cad/004c_performance_summary_table_analysis_selected_analysis.rds"))
performance_table_DCM_selected <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/dcm/004c_performance_summary_table_analysis_selected_analysis.rds"))
performance_table_HFREF_selected <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/hfref/004c_performance_summary_table_analysis_selected_analysis.rds"))

performance_table_ACS_selected_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/acs/004c_performance_summary_table_analysis_selected_analysis_matchIt.rds"))
performance_table_CAD_selected_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/cad/004c_performance_summary_table_analysis_selected_analysis_matchIt.rds"))
performance_table_DCM_selected_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/dcm/004c_performance_summary_table_analysis_selected_analysis_matchIt.rds"))
performance_table_HFREF_selected_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/hfref/004c_performance_summary_table_analysis_selected_analysis_matchIt.rds"))
# had to remove metafeatures from matching
performance_table_ACS_selected_matched2 <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/acs/004c_20240121_performance_summary_table_analysis_selected_analysis_matchIt.rds"))
performance_table_CAD_selected_matched2 <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/cad/004c_20240121_performance_summary_table_analysis_selected_analysis_matchIt.rds"))
performance_table_DCM_selected_matched2 <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/dcm/004c_20240121_performance_summary_table_analysis_selected_analysis_matchIt.rds"))
performance_table_HFREF_selected_matched2 <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/hfref/004c_20240121_performance_summary_table_analysis_selected_analysis_matchIt.rds"))


# Step 1: Combine the performance tables, select models
top_models <- bind_rows(
  performance_table_ACS_selected %>% mutate(disease = "ACS") %>% filter(model=="XGB"),
  performance_table_CAD_selected_matched2 %>% mutate(disease = "CAD") %>% filter(model == "logistic_reg_simple"),
  performance_table_DCM_selected %>% mutate(disease = "DCM") %>% filter(model=="XGB"),
  performance_table_HFREF_selected_matched2 %>% mutate(disease = "HFREF") %>% filter(model == "RF"),
)
top_models

# Step 3: Reshape the data to long format
long_data <- top_models %>%
  select(disease, model, sensitivity, specificity, AUC) %>%
  pivot_longer(cols = c(sensitivity, specificity, AUC), names_to = "metric", values_to = "value") %>% 
  mutate(
    #metric = ifelse(metric == "f1.score", "F1-score", metric),
    metric = ifelse(metric == "accuracy", "Accuracy", metric),
    metric = ifelse(metric == "sensitivity", "Sensitivity", metric),
    metric = ifelse(metric == "specificity", "Specificity", metric)) %>% 
  mutate(disease = factor(disease, levels = rev(c("ACS", "CAD", "DCM", "HFREF")))) %>% 
  mutate(metric = factor(metric, levels = (c("AUC", "Sensitivity", "Specificity"))))

# Step 4: Plot with ggplot2
performance_barplot_selected <- ggplot(
  long_data, aes(x = disease, y = value, fill = metric,
                 #group = model)
  )) +
  geom_bar(stat = "identity", position = "dodge", alpha=0.8) +
  #facet_wrap(~ metric, scales = "free_y") + # Use facet_wrap just for y-axis variations
  theme_minimal() +
  ylab("Value") +
  xlab("Disease") +
  labs(fill = "Model") +
  labs(fill=NULL)+
  xlab(NULL)+
  ylab(NULL)+
  theme(legend.position = "top")+
  coord_flip()+
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) + # Set the desired breaks for the y-axis
  my_base_theme()
performance_barplot_selected

ggsave(filename = glue("{data_path_bestageing2022}/output/plots/ml_summary_plot_patchwork/performance_summary_barplot_selected.svg"), 
       plot = performance_barplot_selected, 
       width = 6, height = 3, 
       units = "in"
)

# paper plot m=20 feature selection --------------------------------------------------------------

performance_table_ACS_full <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/acs/004c_performance_summary_table_analysis_full_analysis_m20.rds"))
performance_table_CAD_full <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/cad/004c_performance_summary_table_analysis_full_analysis_m20.rds"))
performance_table_DCM_full <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/dcm/004c_performance_summary_table_analysis_full_analysis_m20.rds"))
performance_table_HFREF_full <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/hfref/004c_performance_summary_table_analysis_full_analysis_m20.rds"))

performance_table_ACS_full_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/acs/004c_performance_summary_table_analysis_full_analysis_m20_matchIt.rds"))
performance_table_CAD_full_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/cad/004c_performance_summary_table_analysis_full_analysis_m20_matchIt.rds"))
performance_table_DCM_full_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/dcm/004c_performance_summary_table_analysis_full_analysis_m20_matchIt.rds"))
performance_table_HFREF_full_matched <- readRDS(file = glue("{data_path_bestageing2022}/output/performance_summary_df/hfref/004c_performance_summary_table_analysis_full_analysis_m20_matchIt.rds"))


# Step 1: Combine the performance tables, select models
top_models_full <- bind_rows(
  performance_table_ACS_full %>% mutate(disease = "ACS") %>% filter(model=="XGB"),
  performance_table_CAD_full_matched %>% mutate(disease = "CAD") %>% filter(model == "XGB"),
  performance_table_DCM_full %>% mutate(disease = "DCM") %>% filter(model=="logistic_reg_simple"),
  performance_table_HFREF_full_matched %>% mutate(disease = "HFREF") %>% filter(model == "logistic_reg_simple"),
)
top_models_full

# Step 3: Reshape the data to long format
long_data <- top_models_full %>%
  select(disease, model, sensitivity, specificity, AUC) %>%
  pivot_longer(cols = c(sensitivity, specificity, AUC), names_to = "metric", values_to = "value") %>% 
  mutate(
    #metric = ifelse(metric == "f1.score", "F1-score", metric),
    metric = ifelse(metric == "accuracy", "Accuracy", metric),
    metric = ifelse(metric == "sensitivity", "Sensitivity", metric),
    metric = ifelse(metric == "specificity", "Specificity", metric)) %>% 
  mutate(disease = factor(disease, levels = rev(c("ACS", "CAD", "DCM", "HFREF")))) %>% 
  mutate(metric = factor(metric, levels = (c("AUC", "Sensitivity", "Specificity"))))

# Step 4: Plot with ggplot2
performance_barplot_full <- ggplot(
  long_data, aes(x = disease, y = value, fill = metric,
                 #group = model)
  )) +
  geom_bar(stat = "identity", position = "dodge", alpha=0.8) +
  #facet_wrap(~ metric, scales = "free_y") + # Use facet_wrap just for y-axis variations
  theme_minimal() +
  ylab("Value") +
  xlab("Disease") +
  labs(fill = "Model") +
  labs(fill=NULL)+
  xlab(NULL)+
  ylab(NULL)+
  theme(legend.position = "top")+
  coord_flip()+
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) + # Set the desired breaks for the y-axis
  my_base_theme()
performance_barplot_full

ggsave(filename = glue("{data_path_bestageing2022}/output/plots/ml_summary_plot_patchwork/performance_summary_barplot_full.svg"), 
       plot = performance_barplot_full, 
       width = 6, height = 3, 
       units = "in"
)

