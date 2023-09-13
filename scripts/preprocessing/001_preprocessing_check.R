

# check miRNA preprocessing again 2023-09-05 


# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
SAVE.files <- TRUE
runTests <- TRUE

# Define library and data paths based on system
if (system_name == "MacBook-Pro-CR-2065.local" | system_name == "dhcp172-619.laptop-zim.uni-heidelberg.de") {
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
library(AgiMicroRna)
require(readxl, lib.loc = lib_path)
require(janitor, lib.loc = lib_path)
require(e1071, lib.loc = lib_path)
require(reshape2, lib.loc = lib_path)
require(glue, lib.loc = lib_path)
require(gt, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(tidyr, lib.loc = lib_path)
require(stringr, lib.loc = lib_path)
require(purrr, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(ggplot2, lib.loc = lib_path)
require(gplots, lib.loc = lib_path)  # for Venn plot/ data intersections
require(gridExtra, lib.loc = lib_path)  # for adding marginal density plots
require(cowplot, lib.loc = lib_path)  # for adding marginal density plots
require(RColorBrewer, lib.loc = lib_path)
require(ggdist, lib.loc = lib_path)
require(gghalves, lib.loc = lib_path)
require(ggrepel, lib.loc = lib_path)
require(ggvenn, lib.loc = lib_path)
require(rstatix, lib.loc = lib_path)
require(ggthemes, lib.loc = lib_path)
require(ggpubr, lib.loc = lib_path)
require(pROC, lib.loc = lib_path)
require(reshape2, lib.loc = lib_path)
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(janitor::make_clean_names)
conflicted::conflicts_prefer(pROC::roc)

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

###
# We apply both parametric t-tests and nonparametric U-tests ----------------

###
# START LOOP for specified combinations -------------------------------------
###

p2_list <- list()

for (i in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])) {  # only consider all miRNAs for DE (not research miRNAs only)
  ## reassign disease since only full analysis here
  disease <- all_combis[all_combis[["analysis"]] == "full", ]$diseases[i]
  
  ## DATA FIRST
  # create parameter specific {disease}_ids
  filename <- paste0(data_path_BestAgeing , "/data/pheno_", disease, ".xlsx")
  disease_vector <- paste0(disease, "_ids")
  # assign vector to string variable above    (CAVE: to evaluate -> `eval(as.symbol(X))`)
  assign(x = disease_vector,
         value = read_excel(filename) %>% 
           dplyr::pull(BestAgeingCode))
  # df with pat_id and disease identifier
  disease_ident_df <- tibble(pat_id = c(eval(as.symbol(disease_vector)),control_ids), 
                             disease = c(rep(disease, length(eval(as.symbol(disease_vector)))),
                                         rep("control", length(control_ids)))
  )
  
  # prepare-combine-data
  data01 <- all_mirnas %>% 
    # filter control patients and from {disease}_ids
    filter(pat_id %in% !!sym(disease_vector) |
             pat_id %in% control_ids
    ) %>% 
    # bind age and gender from metadata
    left_join(clean_all_meta %>% select(patID, age, sex), 
              by = c("pat_id"= "patID") ) %>% 
    left_join(disease_ident_df, by = "pat_id") %>% 
    relocate(c(disease, age, sex), .after = pat_id) %>% 
    mutate(disease = factor(disease),
           sex = factor(sex) ) %>% 
    select(-age, -sex, -pat_id) %>% 
    as_tibble()
  
  # check log2 transformation --------------------------------------------------
  # there should not be any negative or zero numbers!
  
  # Filter out only numeric columns
  numeric_data <- data01[sapply(data01, is.numeric)]
  
  # Count negative values for each numeric column
  negative_counts <- sapply(numeric_data, function(col) sum(col < 0, na.rm = TRUE))
  # Count zero values for each numeric column
  zero_counts <- sapply(numeric_data, function(col) sum(col == 0, na.rm = TRUE))
  
  # Print the results
  glue("Negative counts for each column: {sum(negative_counts != 0)}")
  glue("Zero counts for each column: {sum(zero_counts != 0)}")
  # this looks good
  
  ## summary stats: skewness and kurtosis --------------------------------------
  #  If the values for most features are close to 0, it suggests the data might be log-transformed and/or normalized
  
  skewness_values <- apply(data01[, -1], 2, skewness)
  kurtosis_values <- apply(data01[, -1], 2, kurtosis)
  sum_stats_skew_kurt <- tibble(miRNA = names(data01[, -1]), skewness_miRNA = skewness_values, kurtosis_miRNA = kurtosis_values)
  data_long <- pivot_longer(sum_stats_skew_kurt, cols = !miRNA, names_to = "statistic", values_to = "value")
  
  # Plot histogram with skewness and kurtosis
  histo_skew_kurt <- ggplot(data_long, aes(x = value, fill = statistic)) +
    geom_histogram(bins = 30, position = 'identity', alpha = 0.5) +
    labs(title = "Distribution of Skewness and Kurtosis", 
         x = "Value", 
         y = "Count") +
    theme_minimal() +
    facet_wrap(~ statistic, scales = "free_x")  # Separate plots for skewness and kurtosis
  print(histo_skew_kurt)  # features with high kurtosis and skewness are prevalent!
  
  set.seed(1234)  # randomly pick miRNAs
  sampled_columns <- sample(names(data01)[-1], 20)
  df_sampled <- data01[, c("disease", sampled_columns)]
  df_long_sampled <- melt(df_sampled, id.vars="disease") %>% 
    mutate(variable = sub("hsa_", "", variable))
  
  # Plot histograms to check log2 transformation
  histo_random_distribution <-  ggplot(df_long_sampled, aes(x=value)) + 
    geom_histogram(bins=30, fill="skyblue", alpha=0.7) + 
    facet_wrap(~variable, scales="free") + 
    theme_minimal() +
    labs(title="Histogram of Sampled Features")
  print(histo_random_distribution)
  
  # Q-Q plots
  qqplot_dist <- ggplot(df_long_sampled, aes(sample=value)) + 
    stat_qq(alpha=0.3, shape=16, size=1, color="skyblue") + 
    stat_qq_line() + # aes(color=variable)) + 
    theme_minimal() + 
    labs(title="Q-Q Plot of Sampled Features")+
    facet_wrap(~variable, scales="free")
  print(qqplot_dist)
  
  # density plot all features
  # df_melt <- melt(data01, id.vars="disease")
  # ggplot(df_melt, aes(x=value, fill=variable)) + 
  #   geom_density(alpha=0.5) +
  #   theme_minimal() +
  #   labs(title="Density Plot of All Features")
  
  # apply qq-normalization (again?) -----------------------------------------
  require(preprocessCore, lib.loc = lib_path)

  # Normalize within each group, cave transposing needed
  group_levels <- unique(data01$disease)
  for (group in group_levels) {
    group_indices <- which(data01$disease == group)
    data01[group_indices, -1] <- t(normalize.quantiles(t(data01[group_indices, -1])))
  }
  
  # Normalize the entire dataset
  data01[,-1] <- t(normalize.quantiles(t(data01[,-1])))
  
  # plots after qq-normalization again
    ## summary stats: skewness and kurtosis --------------------------------------
  #  If the values for most features are close to 0, it suggests the data might be log-transformed and/or normalized
  
  skewness_values <- apply(data01[, -1], 2, skewness)
  kurtosis_values <- apply(data01[, -1], 2, kurtosis)
  sum_stats_skew_kurt <- tibble(miRNA = names(data01[, -1]), skewness_miRNA = skewness_values, kurtosis_miRNA = kurtosis_values)
  data_long <- pivot_longer(sum_stats_skew_kurt, cols = !miRNA, names_to = "statistic", values_to = "value")
  
  # Plot histogram with skewness and kurtosis
  histo_skew_kurt <- ggplot(data_long, aes(x = value, fill = statistic)) +
    geom_histogram(bins = 30, position = 'identity', alpha = 0.5) +
    labs(title = "Distribution of Skewness and Kurtosis", 
         x = "Value", 
         y = "Count") +
    theme_minimal() +
    facet_wrap(~ statistic, scales = "free_x")  # Separate plots for skewness and kurtosis
  print(histo_skew_kurt)  # features with high kurtosis and skewness are prevalent!
  
  set.seed(1234)  # randomly pick miRNAs
  sampled_columns <- sample(names(data01)[-1], 20)
  df_sampled <- data01[, c("disease", sampled_columns)]
  df_long_sampled <- melt(df_sampled, id.vars="disease") %>% 
    mutate(variable = sub("hsa_", "", variable))
  
  # Plot histograms to check log2 transformation
  histo_random_distribution <-  ggplot(df_long_sampled, aes(x=value)) + 
    geom_histogram(bins=30, fill="skyblue", alpha=0.7) + 
    facet_wrap(~variable, scales="free") + 
    theme_minimal() +
    labs(title="Histogram of Sampled Features")
  print(histo_random_distribution)
  
  # Q-Q plots
  qqplot_dist <- ggplot(df_long_sampled, aes(sample=value)) + 
    stat_qq(alpha=0.3, shape=16, size=1, color="skyblue") + 
    stat_qq_line() + # aes(color=variable)) + 
    theme_minimal() + 
    labs(title="Q-Q Plot of Sampled Features")+
    facet_wrap(~variable, scales="free")
  print(qqplot_dist)
  
}
