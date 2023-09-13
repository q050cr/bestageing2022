

# age analysis

# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
SAVE.files <- TRUE

# Define library and data paths based on system
if (system_name == "MacBook-Pro-CR-2065.local") {
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

# dependencies ---------------------------------------------------------------
library(readxl, lib.loc = lib_path)
library(janitor, lib.loc = lib_path)
library(glue, lib.loc = lib_path)
library(gt, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(tidyr, lib.loc = lib_path)
library(stringr, lib.loc = lib_path)
library(purrr, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(broom, lib.loc = lib_path)
library(ggplot2, lib.loc = lib_path)
library(fmsb, lib.loc = lib_path)  # radar plot
library(RColorBrewer, lib.loc = lib_path)
library(ggdist, lib.loc = lib_path)
library(gghalves, lib.loc = lib_path)
library(ggrepel, lib.loc = lib_path)
library(rstatix, lib.loc = lib_path)
library(ggthemes, lib.loc = lib_path)
library(ggpubr, lib.loc = lib_path)
library(pROC, lib.loc = lib_path)
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(janitor::make_clean_names)

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) 


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="full")


# load data ---------------------------------------------------------------

# MIRNA DAT
model_data1 <- clean_names(readRDS(file = glue('{data_path_BestAgeing}/data_new/model_data1.RDS')))  # has also multiclass col + diagnoses
load(file = glue('{data_path_BestAgeing}/data/mirnas.rda'))  # "UKL-HD" n=765
load(file = glue('{data_path_BestAgeing}/data/data.rda'))  # "UKL-HD" n=731
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

# load-mirnas-from_research
# create vector of described mirnas
load(glue("{data_path_BestAgeing}/data_research/fromR/researchMiRNAAccession.rda"))
# check if all mirnas are named the same
researchMiRNAAccession$miRNAName_v21 <-  make_clean_names(researchMiRNAAccession$miRNAName_v21) %>% 
  str_replace(pattern = "mi_r", replacement = "mir")
# researchMiRNAAccession$miRNAName_v21[(!researchMiRNAAccession$miRNAName_v21 %in% colnames(all_mirnas))] 
## "hsa_mir_106a_5p" --> should be --> "hsa_mir_106b_5p"
## all_mirnas[,(str_detect(string = colnames(all_mirnas), pattern = "106"))]
## researchMiRNAAccession[which(!researchMiRNAAccession$miRNAName_v21 %in% colnames(all_mirnas)),]
researchMiRNAAccession$miRNAName_v21[researchMiRNAAccession$miRNAName_v21 == "hsa_mir_106a_5p"] <- "hsa_mir_106b_5p"

###
# load-metadat --
## DIAGNOSES DAT
load(glue("{data_path_BestAgeing}/data/diagnoses_df.rda"))

## SURVIVAL DAT
survival_dat <- clean_names(readRDS(glue("{data_path_bestageing2022}/data/202211908_XMELD_abfrage_best_ageing.rds"))) # %>% 
# original path "/mnt/users/reich/XMeldPortal_neu/meldeportal-tools-meldeportalclient-9.3/Rout/202211908_XMELD_abfrage_best_ageing.rds"

## metadata from DB
# https://www.bestageing.org/Pages/Login.aspx?ReturnUrl=%2f&AspxAutoDetectCookieSupport=1
load(glue("{data_path_BestAgeing}/data/clean_all_meta.rda"))  # created in "scripts/_prepare_metadata.R"
clean_all_meta <- clean_all_meta %>% 
  mutate(age = ifelse(age < 18, NA, age))  # wrong age remove
# cath data? "hkdb"

## load all original metadat xlsx files again to make sure that also overlapped 
#patients (e.g. dcm+cad) are in each group
control_ids <- read_excel(glue("{data_path_BestAgeing}/data/pheno_controls.xlsx")) %>% 
  dplyr::pull(BestAgeingCode)

# "UKL-HD-00318" both in Control and CAD dataset, looked it up (HK Nr 1289-2015): KHK ohne hg Stenosen, LV gut --> assign to CAD only
control_ids <- control_ids[control_ids != "UKL-HD-00318"]


# Stratified Analysis ---------------------------------------------------------
# Perform analysis separately for control and disease groups. This could reveal differing relationships between age and biomarker expression in each group.
mirnas_age_analysis <- all_mirnas %>% 
  left_join(model_data1 %>% select(pat_id, multiclass, control, acs, cad, dcm, ref), by=c("pat_id"="pat_id")) %>% 
  left_join(clean_all_meta %>% select(patID, age, sex),  by=c("pat_id"="patID")) %>% 
  relocate(age, sex, multiclass, control, acs, cad, dcm, ref, .after = pat_id)

mirnas_age_analysis <- mirnas_age_analysis %>% 
  filter(multiclass != "pef")

# nested analysis
##### Using maps instead idea ~/maps_instead_of_for_loops.R  (but here we have two iterables: group and miRNA vector)

miRNA_vector <- names(mirnas_age_analysis)[10:length(names(mirnas_age_analysis))]

mirnas_age_analysis$multiclass <- relevel(mirnas_age_analysis$multiclass, ref = "control")

nested_data <- mirnas_age_analysis |> 
  group_by(multiclass) |> 
  nest() %>% 
  arrange(multiclass)

nested_data %>% 
  mutate(
    lin_mod = map(
      .x = data,
      function(x) lm(data = x, hsa_let_7a_3p ~ age)
    ),
    coefficients = map(lin_mod, coefficients),
    slope = map_dbl(coefficients, \(x) x[2]),
    slope_short = map_dbl(coefficients, 2)
  )


nested_data <- nested_data %>% 
  mutate(results_lm_cor = list(NULL), 
         results_interaction_term = list(NULL), 
         plots_lm_scatter = list(NULL), 
         plots_interaction_scatter = list(NULL)
         )

for(mynest in 1:nrow(nested_data)) {
  result_list <- miRNA_vector %>% 
    # x ... miRNA from vector
    map(.f = function(x) summary(lm(data = nested_data$data[[mynest]], as.formula(paste0(x, "~ I(age/10)")) )), .progress = TRUE) %>% 
    set_names(nm = miRNA_vector)
  
  result_list_cor_coef <- miRNA_vector %>% 
    #map(.f = function(x) cor(nested_data$data[[mynest]]$age, nested_data$data[[mynest]][[x]], use = "complete" ), .progress = TRUE) %>% 
    map(.f = function(x) cor.test(nested_data$data[[mynest]]$age, nested_data$data[[mynest]][[x]], use = "complete" ), .progress = TRUE) %>%   # more info
    set_names(nm = miRNA_vector)
  
  results_tibble_tmp <- tibble(miRNA = names(result_list), lin_mod = result_list, cor_analysis = result_list_cor_coef) %>% 
    mutate(raw_pval_age = map_dbl(.x=lin_mod, ~ .x$coefficients[2, 4] ,.progress = TRUE),
           padj_age = map_dbl(raw_pval_age, function(x) p.adjust(x, method = "BH", n = length(miRNA_vector))),
           OR_10y_age = map_dbl(.x=lin_mod, ~ .x$coefficients[2, 1], .progress = TRUE),
           raw_pval_age_cor = map_dbl(.x=cor_analysis, ~ .x$p.value, .progress = TRUE),
           padj_age_cor = map_dbl(raw_pval_age_cor, function(x) p.adjust(x, method = "BH", n = length(miRNA_vector))),
           cor_pearson_age = map_dbl(.x=cor_analysis, ~ .x$estimate, .progress = TRUE),
           )
  nested_data$results_lm_cor[[mynest]] <- results_tibble_tmp
}

# outline: get significant miRNAs, make scatter plots like in BM publication. make strata 40-50y, 50-60y ... -> spider plot for mean expression of each miRNA



# Interaction Term each biomarker~I(Age/10):Diseases ------------------------------------------------------------
# In a regression model, include an interaction term between age and disease status. 
# This would allow you to see if the effect of age on biomarker expression differs by disease status.

# In the output, look for the coefficient for the interaction term (age:disease). 
# If it's significant, it suggests that the relationship between age and biomarker expression differs by group!

for(mydisease in seq_along(names(mirnas_age_analysis)[6:9]) ){  # checked that variable names are ordered!
  colname_disease <- names(mirnas_age_analysis)[6:9][mydisease]
  mirnas_interaction_age_analysis <- mirnas_age_analysis %>% 
    filter(control == 1 | !!rlang::sym(colname_disease) == 1)
  
  result_interaction_list <- miRNA_vector %>% 
    map(.f = function(x) summary(lm(data = mirnas_interaction_age_analysis, as.formula(paste0(x, "~ I(age/10) * ", colname_disease)) )), .progress = TRUE) %>% 
    set_names(nm = miRNA_vector)
    
  results_interaction_tibble_tmp <- tibble(miRNA = names(result_list), lin_mod_interaction = result_interaction_list) %>% 
    mutate(raw_pval_ageXdisease = map_dbl(.x=lin_mod_interaction, ~ .x$coefficients[4, 4] ,.progress = TRUE),
           padj_ageXdisease = map_dbl(raw_pval_ageXdisease, function(x) p.adjust(x, method = "BH", n = length(miRNA_vector)))
    )
  
  nested_data$results_interaction_term[[mydisease+1]] <- results_interaction_tibble_tmp  # need the "+1" for index, because we are comparing control vs diseases (control stays empty)
}

# interaction term interpretation always difficult.. 
# "The interaction term between age and group is significant for microRNA X. This suggests that the association between age and 
# the expression of microRNA X is not the same across all disease groups. For example, as age increases, microRNA X expression might 
# increase in Group A, decrease in Group B, and stay the same in the control group (Group C). It indicates that the dynamics of 
# microRNA X expression as a function of age differs depending on the disease group. Therefore, both age and disease type must 
# be considered together when studying the expression of microRNA X."

## plots -------------------------------------------------------------------
labels_interaction_plot <- c("Control", "ACS", "CAD", "DCM", "HFrEF")

for(mynest in 1:nrow(nested_data)){
  # A) lm scatter plot
  mirnas_age_analysis_lm_significant <- nested_data$results_lm_cor[[mynest]] %>% 
    filter(padj_age < 0.05) %>% 
    pull(miRNA)
  
  if(length(mirnas_age_analysis_lm_significant) > 0) {  # could also use trycatch, but throws error Error in `pivot_longer()`:! `cols` must select at least one column.
    df_long <- nested_data$data[[mynest]] %>% 
      select(age, all_of(mirnas_age_analysis_lm_significant)) %>% 
      pivot_longer(cols = -age, names_to = "miRNA", values_to = "log2expression") %>% 
      mutate(miRNA = str_replace_all(miRNA, "_", "-"))
    
    lm_scatter_plot <- ggplot(df_long, aes(x = age, y = log2expression)) + 
      geom_point(alpha=0.3, shape=16) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      facet_wrap(~ miRNA, scales = "free", nrow=5) +
      labs(x = "Age", y = expression(paste("log"[2], " expression"))) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6)) +
      my_base_theme()
    # store plot
    nested_data$plots_lm_scatter[[mynest]] <- lm_scatter_plot
    
    height_lm_scatter_plot <- ifelse(as.character(nested_data$multiclass[mynest]) != "control", 10, 50)  # 10*10  .. does also not look good
    
    if (SAVE.files ==TRUE) {
      path_lm_scatter_plot <- glue("{data_path_bestageing2022}/output/plots/age_analysis/stratified_lm/{Sys.Date()}_{nested_data$multiclass[mynest]}_lm_scatter_significant_wrap.svg")
      ggsave(filename = path_lm_scatter_plot, plot = lm_scatter_plot, 
             width = 10, height = height_lm_scatter_plot, limitsize = FALSE,
             units = "in"  # default
      )
    }
    
    # radar plots -----------------------------------------------------------
    df_radar <- nested_data$data[[mynest]] %>% 
      select(age, all_of(mirnas_age_analysis_lm_significant)) %>% 
      drop_na(age) %>% 
      mutate(age_group = cut(age, breaks = c(0, 50, 70, Inf), labels = c("50-", "50-70", "70+")))
    
    # Calculate the mean of each miRNA for each age group
    df_mean <- df_radar %>% 
      group_by(age_group) %>% 
      summarise(across(all_of(mirnas_age_analysis_lm_significant), mean, na.rm = TRUE))
    # Create a dataframe suitable for radarchart()
    df_mean <- df_mean %>% as.tibble()
    ## The first row needs to be the max values for the chart
    #max_values <- apply(df_mean[-1], 2, max)
    #df_radar_prepared_plot <- rbind(max_values, df_mean)
    
    # dotchart -----------------------------------------------------------
    # https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
    df_dotchart_age_bins <- df_mean %>% 
      pivot_longer(cols = -age_group, names_to = "miRNA", values_to = "log2expression") %>% 
      mutate(miRNA = str_replace_all(miRNA, "_", "-"))
    agegroup_dotchart_plot <- ggdotchart(data = df_dotchart_age_bins, x = "miRNA", y = "log2expression", group="age_group", color="age_group", 
               add = "segment", position = position_dodge(0.3), 
               dot.size = 3,
               sorting = "descending",
               rotate = TRUE) +
      labs(color="Age Group")+
      xlab(NULL) +
      ylab("Mean log 2 expression")+
      scale_y_log10() +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_color_manual(values = thematic::okabe_ito(6)) +
      #theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
      my_base_theme()
    
    if (SAVE.files ==TRUE) {
      path_agegroup_dotchart_plot <- glue("{data_path_bestageing2022}/output/plots/age_analysis/dotchart_agegroups/{Sys.Date()}_{nested_data$multiclass[mynest]}_agegroup_dotchart_plot_significant.svg")
      ggsave(filename = path_agegroup_dotchart_plot, plot = agegroup_dotchart_plot, 
             width = 10, height = 14, limitsize = FALSE,
             units = "in"  # default
      )
    }
    
  }

  # end of lm scatter plot analysis
  
  # B) interaction plot
  if(nested_data$multiclass[mynest] != "control") {
    # run interaction plots
    colname_disease <- as.character(nested_data$multiclass[mynest])  # convert factor variable to string
    mirnas_interaction_age_analysis <- mirnas_age_analysis %>% 
      filter(control == 1 | !!rlang::sym(colname_disease) == 1)
    
    # get significant miRNAs
    mirnas_ageXdisease_int_analysis_lm__significant <- nested_data$results_interaction_term[[mynest]] %>% 
      filter(padj_ageXdisease < 0.05) %>% 
      pull(miRNA)
    
    if(length(mirnas_ageXdisease_int_analysis_lm__significant) > 0) {
      df_long_interaction <- mirnas_interaction_age_analysis %>% 
        select(age, all_of(colname_disease), all_of(mirnas_ageXdisease_int_analysis_lm__significant)) %>% 
        pivot_longer(cols = -c(age, all_of(colname_disease)), names_to = "miRNA", values_to = "log2expression") %>% 
        mutate(!!colname_disease := recode_factor(!!sym(colname_disease), "0" = "Control", "1" = labels_interaction_plot[mynest])) %>% 
        mutate(miRNA = str_replace_all(miRNA, "_", "-"))
      
      lm_interaction_scatter_plot <- ggplot(df_long_interaction, aes(x=age, y=log2expression, color=!!sym(colname_disease))) + 
        geom_point(alpha=0.3, shape=16) +
        geom_smooth(method="lm", se=TRUE) + 
        facet_wrap(~ miRNA, scales = "free") +
        theme_minimal() +
        labs(title = NULL, # "Interaction of age and disease on biomarker expression",
             x = "Age",
             y = expression(paste("log"[2], " expression")),
             color = NULL) +
        theme_minimal(base_size = 16, base_family = 'Arial')+
        scale_color_manual(values = thematic::okabe_ito(6)) +
        my_base_theme()
      #theme(legend.position="bottom")
      # store plot
      nested_data$plots_interaction_scatter[[mynest]] <- lm_interaction_scatter_plot
      if (SAVE.files ==TRUE) {
        path_lm_interaction_scatter_plot <- glue("{data_path_bestageing2022}/output/plots/age_analysis/interaction_ageXdisease/{Sys.Date()}_{nested_data$multiclass[mynest]}_lm_scatter_significant_wrap.svg")
        ggsave(filename = path_lm_interaction_scatter_plot, plot = lm_interaction_scatter_plot, 
               width = 10, height = 10, 
               units = "in"  # default
        )
      }
    }
    # end of interaction analysis
  }
  # end of loop
}




# Mediation Analysis ------------------------------------------------------------
#  If you believe that age affects disease status, and disease status in turn affects 
#  biomarker expression (or vice versa), a mediation analysis could be used to quantify the direct and indirect effects of age on biomarker expression.





