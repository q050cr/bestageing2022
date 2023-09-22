
### INFO ----------------------------------------------------------------------
# Differential miRNA Expression Analysis Sensitivity Analysis

# update in 001d script: compared to 001c this script considers only a subgroup analysis that removes batch center (only possible if both disease and controls were sampled from same center)
#                         therefore, it only considers ACS vs controls (UKHD and Uppsala university) n =205; DCM vs controls (only UKHD n=105); CAD vs controls (only UKHD n =338)

# this script is sourced from `scripts/render_param_reports.R`
# selection provided by `all_combis$diseases` and `all_combis$analysis`

# script creates plots: "fig01vogel2013", "fig02vogel2013"


# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
SAVE.files <- TRUE
runTests <- TRUE
runqqnorm <- FALSE  # already done before
filterDetMatrix <- TRUE  # new in this 001c script
runCOMBAT <- FALSE  # problem perfect correlation recruitment and disease
runLimmaRemoveBatchEffect <- TRUE
runSVA <- TRUE  # adds sva if computed to dataframe and adjusts for it when testing
runRobustRegression <- TRUE

# Define library and data paths based on system
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

# dependencies ---------------------------------------------------------------
library(readxl, lib.loc = lib_path)
library(janitor, lib.loc = lib_path)
library(glue, lib.loc = lib_path)
# NEW
require(preprocessCore, lib.loc = lib_path)  # preprocess qq normalization again
library(sva) # , lib.loc = lib_path)
library(limma) # , lib.loc = lib_path)
library(robustbase, lib.loc = lib_path)
library(robust, lib.loc = lib_path)

library(gt, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(tidyr, lib.loc = lib_path)
library(stringr, lib.loc = lib_path)
library(purrr, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(tibble, lib.loc = lib_path)
library(ggplot2, lib.loc = lib_path)
library(gplots, lib.loc = lib_path)  # for Venn plot/ data intersections
library(gridExtra, lib.loc = lib_path)  # for adding marginal density plots
library(cowplot, lib.loc = lib_path)  # for adding marginal density plots
library(RColorBrewer, lib.loc = lib_path)
library(ggpubr, lib.loc = lib_path)
library(patchwork, lib.loc = lib_path)
library(ggridges, lib.loc = lib_path)
library(ggdist, lib.loc = lib_path)
library(gghalves, lib.loc = lib_path)
library(ggrepel, lib.loc = lib_path)
library(ggvenn, lib.loc = lib_path)
library(rstatix, lib.loc = lib_path)
library(ggthemes, lib.loc = lib_path)

library(pROC, lib.loc = lib_path)
library(reshape2, lib.loc = lib_path)
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(janitor::make_clean_names)
conflicted::conflicts_prefer(pROC::roc)

# functions ---------------------------------------------------------------
mode_function <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

remove_digits <- function(vec) {
  # Use regex to replace patterns like '-DIGITS' with an empty string
  return(sub("-\\d+$", "", vec))
}

# Define a function for the permutation test
permute_test <- function(data, column_name, nperm = 1000) {
  # Split data by disease status
  group1 <- data[data$disease == 'control', column_name]
  group2 <- data[data$disease != 'control', column_name]  # Assuming all other are disease
  
  # Compute observed t-statistic
  observed_t <- abs(t.test(group1, group2)$statistic)
  
  # INIT vector with length "nperm" to store permuted t-values
  permuted_t <- numeric(nperm)
  
  # Loop for permutations
  for(i in 1:nperm) {
    shuffled_disease <- sample(data$disease)  # Shuffle disease status!!
    perm_group1 <- data[shuffled_disease == 'control', column_name]
    perm_group2 <- data[shuffled_disease != 'control', column_name]
    permuted_t[i] <- abs(t.test(perm_group1, perm_group2)$statistic)
  }
  
  # Compute p-value as proportion of times permuted t-values exceed observed t-value
  p_value <- mean(permuted_t >= observed_t)
  
  return(list(observed_t = observed_t, p_value = p_value))
}

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))


# load data ---------------------------------------------------------------

diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) 


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="full")

# MIRNA DAT
model_data1 <- clean_names(readRDS(file = glue('{data_path_BestAgeing}/data_new/model_data1.RDS')))  # has also multiclass col + diagnoses
load(file = glue('{data_path_BestAgeing}/data/mirnas.rda'))  # "UKL-HD" n=765
load(file = glue('{data_path_BestAgeing}/data/data.rda'))  # "UKL-HD" n=731
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

## seq metadat 09-2023
# batch info
hbdx_metadat <- clean_names(readxl::read_excel(path = glue('{data_path_bestageing2022}/data/kahraman2023/230907_annotation_chrstoph_reich.xlsx')))  # has also multiclass col + diagnoses
#RIN
rin_mean <- mean(as.numeric(hbdx_metadat$rin), na.rm=TRUE)
rin_sd <- sd(as.numeric(hbdx_metadat$rin), na.rm=TRUE)

# detection matrix
det_mat_all_mirnas <- read.table("/Volumes/T7CR/data/bestageing2022/data/kahraman2023/det_mat_all_mirnas.txt") %>% 
  t() %>% 
  as.data.frame() %>% 
  clean_names()
colnames(det_mat_all_mirnas) <- str_replace(string = colnames(det_mat_all_mirnas), pattern = "mi_r", replacement = "mir")
rownames(det_mat_all_mirnas) <- gsub("\\.", "-", rownames(det_mat_all_mirnas))
# convert to tibble
det_mat_all_mirnas <- det_mat_all_mirnas %>% 
  rownames_to_column(var = "pat_id") %>% 
  as_tibble()


# hbdx_metadat %>% select(slide_id, slide_array_id)
unique_slide_ids <- length(unique(hbdx_metadat$slide_id))
arrays_per_slide <- nrow(hbdx_metadat) / unique_slide_ids
sum_duplicated_samples <- sum(duplicated(hbdx_metadat$customer_id)); index_duplicated_samples <- duplicated(hbdx_metadat$customer_id)
duplicate_ids <- hbdx_metadat[index_duplicated_samples, ] %>% pull(customer_id)
# check
checked_duplicated_hbdx <- hbdx_metadat %>% 
  filter(customer_id %in% duplicate_ids) %>% 
  select(customer_id, date_analysis, sample_qc, array_qc) %>% arrange(customer_id)
# remove
hbdx_metadat <- hbdx_metadat %>% 
  distinct(customer_id, .keep_all = TRUE)

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
# We apply both parametric t-tests and nonparametric U-tests


# START LOOP for specified combinations -------------------------------------

p2_list <- list()

for (i in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])) {  # only consider all miRNAs for DE (not research miRNAs only)
  # reassign disease since only full analysis here
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
    as_tibble()
  
  if(runqqnorm == TRUE) {  # already done
    ## apply qq-normalization ---------------------------------------------------
    # Normalize within each group, cave transposing needed
    group_levels <- unique(data01$disease)
    for (group in group_levels) {
      group_indices <- which(data01$disease == group)
      data01[group_indices, -c(1:4)] <- t(normalize.quantiles(t(data01[group_indices, -c(1:4)])))
    }
    
    # Normalize the entire dataset
    data01[,-c(1:4)] <- t(normalize.quantiles(t(data01[,-c(1:4)])))
  }
  
  
  ## FILTERING DET MATRIX ----------------------------------------------------
  if (filterDetMatrix == TRUE){
    filter_dat <- data01 %>% select(-c(age,sex))
    
    det_mat_all_mirnas_tmp <- det_mat_all_mirnas %>% 
      filter(pat_id %in% filter_dat$pat_id) %>% 
      left_join(filter_dat %>% select(pat_id, disease), by=c("pat_id")) %>% 
      relocate(disease, .after = pat_id)
    
    print(dim(filter_dat) == dim(det_mat_all_mirnas_tmp))
    
    # EDA - Plotting a density of average log2 expressions 
    data_long <- filter_dat %>% 
      pivot_longer(cols = -c(disease, pat_id), names_to = "miRNA", values_to = "expression")
    
    # Stratified density plot
    density_plot_filter_eda_stratified <- ggplot(data_long, aes(x=expression, fill=disease)) + 
      geom_density(alpha=0.5) + 
      labs(x="Log2 Expression") +  # title="Density plot of log2 expressions", 
      facet_wrap(~ disease, nrow=2) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
      theme(axis.text.y = element_blank())
    #print(density_plot_filter_eda_stratified)
    
    density_ridges_plot_filter_eda_stratified <- ggplot(data_long, aes(x=expression, y=disease, fill=disease)) + 
      geom_density_ridges(alpha=0.5) +  # , jittered_points = TRUE, position = position_points_jitter(width = 0.05)) + # too much data involved if points plotted
      labs(x="Log2 Expression", y="Density") +  # title="Density Ridge plot of log2 expressions", 
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
      theme(axis.text.y = element_blank())
    #print(density_ridges_plot_filter_eda_stratified)
    
    if (SAVE.files ==TRUE) {
      density_ridges_plot_filter_eda_stratified_path <- glue("{data_path_bestageing2022}/output/plots/preprocessing/filtering_miRNA/001d_{disease}_density_ridges.svg")
      ggsave(filename = density_ridges_plot_filter_eda_stratified_path, plot = density_ridges_plot_filter_eda_stratified, 
             width = 8, height = 5, 
             units = "in"  # default
      )
    }
    
    threshold <- 0.9  # more conservative threshold
    
    # Filter for control samples
    lowExprControl <- colMeans(det_mat_all_mirnas_tmp[det_mat_all_mirnas_tmp$disease == "control", -(1:2)]) <= threshold
    
    # Filter for disease samples
    lowExprDisease <- colMeans(det_mat_all_mirnas_tmp[det_mat_all_mirnas_tmp$disease == disease, -(1:2)]) <= threshold
    
    # Identify miRNAs that are lowly expressed in both conditions
    toFilter <- lowExprControl & lowExprDisease
    
    # Filter out those miRNAs from the dataset
    filteredData <- filter_dat[, c(TRUE, TRUE, !toFilter)] ## TRUE TRUE corresponds to keep pat_id and disease
    
    data01 <- data01 %>% 
      select(1:4) %>% 
      dplyr::bind_cols(filteredData[, -(1:2)])
  }
  
  
  # COMBAT ------------------------------------------------------------------
  if (runCOMBAT == TRUE){
    
    # batch 1 center
    batch_center <- remove_digits(data01$pat_id)
    table(batch_center, data01$disease)
    # batch_center control hfref  # problem with COMBAT! perfect correlation of recruitment status with disease variable
    # AMC          0    39
    # GUF          0    37
    # INSERM       0    81
    # NSC          0   156
    # SERMAS       0    16
    # UCSC         0     5
    # UKL-HD     149     0
    # UNIPD      180     0
    # UOA          0    56
    # UU         519     0
    
    # batch 2 biochip id
    # batch 2
    batch_slide_id <- data01 %>% 
      left_join(hbdx_metadat %>% select(customer_id, slide_id),
                by = c("pat_id"="customer_id")) %>% 
      pull(slide_id)
    
    exprs_data <- as.matrix(data01[,-c(1:4)])
    rownames(exprs_data) <- data01$pat_id
    
    # pheno data
    pheno_data <- data01 %>% 
      select(disease, age, sex) %>% 
      as.data.frame()
    rownames(pheno_data) <- data01$pat_id
    
    # mean/ mode imputation necessary
    # lapply(pheno_data, function(x) sum(is.na(x)))
    pheno_data$age[is.na(pheno_data$age)] <- mean(pheno_data$age, na.rm = TRUE)
    #pheno_data$disease[is.na(pheno_data$disease)] <- mode_function(pheno_data$disease)
    pheno_data$sex[is.na(pheno_data$sex)] <- mode_function(pheno_data$sex)
    
    ## model matrix 
    modcombat <- model.matrix(~1, data=pheno_data) # had thisdesign , but perfect correlation   (~ disease + sex + age, data=pheno_data)
    #null_model <- model.matrix(~ sex + age, data=pheno_data)
    
    combat_edata = ComBat(
      dat=t(exprs_data), batch=batch_center, mod=modcombat, 
      par.prior=TRUE, prior.plots=FALSE, 
      #mean.only = FALSE
    ) 
    
    combat_edata_t <- combat_edata %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "pat_id") %>% 
      as_tibble()
    
    pValuesComBat = f.pvalue(combat_edata,design,null_model)
    qValuesComBat = p.adjust(pValuesComBat,method="BH")
    sum(qValuesComBat<0.05)
    pValuesComBat = f.pvalue(t(exprs_data),design,null_model)
    qValuesComBat = p.adjust(pValuesComBat,method="BH")
    sum(qValuesComBat<0.05)
  }
  
  # LIMMA BATCH EFFECT ------------------------------------------------------
  
  if(runLimmaRemoveBatchEffect == TRUE) {
    batch_center <- remove_digits(data01$pat_id)
    table(batch_center, data01$disease)
    
    
    data01 <- data01 %>% 
      # batch 1
      mutate(batch_center = batch_center) %>% 
      # NEW NEW NEW only keep centers that have both controls and disease sampled NEW NEW NEW !!!!
      group_by(batch_center) %>% 
      filter(all(c("control", {{ disease }} ) %in% disease) ) %>% 
      ungroup()
    
    if(nrow(data01) == 0) {
      warning(glue("CUSTOM ERROR: For {toupper(disease)}: all samples (control/ {disease}) sampled from different centers"))
      next
    }
    
    # add metadata
    exprs_metadat <- data01 %>% 
      # batch 2
      left_join(hbdx_metadat %>% select(customer_id, slide_id, sample_qc, array_qc),
                by = c("pat_id"="customer_id")) %>% 
      relocate(batch_center, slide_id, .after = pat_id) %>% 
      select(-c(sample_qc, array_qc))
    
    exprs_data <- as.matrix(exprs_metadat[,-c(1:6)])
    
    # pheno data
    pheno_data <- data01 %>% 
      mutate(batch_center = batch_center) %>% 
      group_by(batch_center) %>% 
      filter(all(c("control", {{ disease }} ) %in% disease) ) %>% 
      ungroup() %>% 
      select(disease, age, sex) %>% 
      as.data.frame()
    rownames(pheno_data) <- exprs_metadat$pat_id
    pheno_data$age[is.na(pheno_data$age)] <- mean(pheno_data$age, na.rm = TRUE)
    pheno_data$sex[is.na(pheno_data$sex)] <- mode_function(pheno_data$sex)
    
    # update
    batch_slide_id <- exprs_metadat %>% pull(slide_id)
    batch_center <- exprs_metadat %>% pull(batch_center)
    
    data01 <- data01 %>% 
      select(-batch_center)
    
    # PCA after batch removal -------------------------------------------------
    runPCA_batch <- TRUE
    # PCA
    if(runPCA_batch == TRUE) {
      # pca
      pca_res <- prcomp(exprs_data, center = TRUE, scale. = TRUE)
      
      # Get the first two principal components
      pca_df <- as_tibble(as.data.frame(pca_res$x[,1:5]))
      
      exprs_metadat <- exprs_metadat %>% 
        bind_cols(pca_df)
      
      # PC-plot
      # # Custom labels for the legend
      # custom_labels <- c(
      #   "acs" = "ACS",
      #   "cad" = "CAD",
      #   "control" = "Control",
      #   "dcm" = "DCM",
      #   "ref" = "HFrEF"
      # )
      
      palette <- brewer.pal(10, "Set3")
      if(length(unique(batch_center)) < 6)  {
        palette <- thematic::okabe_ito(6)   # only works up to 6 centers!
      }
      # takes some time
      ggplot(exprs_metadat, aes(PC1, PC2, color = batch_center)) +
        geom_point(alpha=0.9, aes(shape=disease)) +
        #coord_fixed(ratio=1)+
        labs(title = NULL,  # "PCA Plot Colored by Disease",
             x = "PC1",
             y = "PC2") +
        theme_minimal(base_size = 16, base_family = 'Arial')+
        scale_fill_manual(values = palette) +
        scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
        scale_shape_discrete(name=NULL, labels = c("Control", toupper(disease))) +
        my_base_theme() -> plot_PC1_PC2
      # plot_PC1_PC2
      
      ggplot(exprs_metadat, aes(PC1, PC3, color = batch_center)) +
        geom_point(alpha=0.9, aes(shape=disease)) +
        labs(title = NULL,  # "PCA Plot Colored by Disease",
             x = "PC1",
             y = "PC3") +
        theme_minimal(base_size = 16, base_family = 'Arial')+
        scale_fill_manual(values = palette) +
        scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
        scale_shape_discrete(name=NULL, labels = c("Control", toupper(disease))) +
        my_base_theme() -> plot_PC1_PC3
      #plot_PC1_PC3
      
      if (SAVE.files ==TRUE) {
        pca12_beforebatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001d_{disease}_pca12_before_batch.svg")
        ggsave(filename = pca12_beforebatch, plot = plot_PC1_PC2, 
               width = 8, height = 6, 
               units = "in"  # default
        )
        pca13_beforebatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001d_{disease}_pca13_before_batch.svg")
        ggsave(filename = pca13_beforebatch, plot = plot_PC1_PC3, 
               width = 6, height = 6, 
               units = "in"  # default
        )
      }
    }
    
    # run limma batch remove https://evayiwenwang.github.io/Managing_batch_effects/adjust.html#accounting-for-batch-effects
    # modlimma <- model.matrix(~1, data=pheno_data) 
    
    # model
    modlimma <- model.matrix( ~ disease, data=pheno_data)
    if (length(unique(batch_center)) > 1) {
      # batch1
      adjustedMatrix <- removeBatchEffect(t(exprs_data), batch=batch_center, design = modlimma)  # collinearity issues with center batch:  one-to-one mapping with any condition
      # batch2
      adjustedMatrix <- removeBatchEffect(adjustedMatrix, batch=batch_slide_id, design = modlimma)  
    } else {
      # batch1
      adjustedMatrix <- removeBatchEffect(t(exprs_data), batch=batch_slide_id, design = modlimma) 
    }
    
    limma_edata <-  adjustedMatrix %>% t() %>%  as.data.frame() %>% as_tibble() %>% 
      mutate(pat_id = exprs_metadat$pat_id) %>% 
      select(pat_id, everything())
    
    # inspect random boxplots 
    set.seed(123)
    random_index <- sample(x = 1:nrow(exprs_data), size = 50, replace = FALSE) %>% sort()  # samples, not miRNAs
    
    # fast boxplot inspect
    #par(mfrow=c(1,2))
    #boxplot(as.data.frame(t(exprs_data))[, random_index],main="Original")
    #boxplot(as.data.frame(adjustedMatrix)[, random_index],main="Batch corrected")
    
    # ggboxplot 
    df_original <- as.data.frame(t(exprs_data)) %>% 
      select(all_of(random_index)) %>% 
      tibble::rownames_to_column(var = "miRNA") %>% 
      pivot_longer(cols = -miRNA, names_to = "Sample", values_to = "Expression")
    
    df_corrected <- as.data.frame(adjustedMatrix) %>% 
      select(all_of(random_index)) %>% 
      tibble::rownames_to_column(var = "miRNA") %>% 
      pivot_longer(cols = -miRNA, names_to = "Sample", values_to = "Expression")
    
    # Find the global y limits to make them same for both plots
    y_limits <- c(min(c(df_original$Expression, df_corrected$Expression)), 
                  max(c(df_original$Expression, df_corrected$Expression)))
    
    # Plot original data
    p1 <- ggplot(df_original, aes(x = Sample, y = Expression)) +
      geom_boxplot(aes(group = Sample), alpha=0.5) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      labs(title = NULL, x=NULL) +
      ylim(y_limits) +
      scale_color_manual(values = palette, name=NULL) +
      theme(axis.text.x = element_blank()) +
      my_base_theme() 
    # Plot batch-corrected data
    p2 <- ggplot(df_corrected, aes(x = Sample, y = Expression)) +
      geom_boxplot(aes(group = Sample), alpha=0.5) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      labs(title = NULL, x=NULL) +
      ylim(y_limits) +
      scale_color_manual(values = palette, name=NULL) +
      theme(axis.text.x = element_blank(),
            axis.title.y=element_blank(), axis.text.y=element_blank()  # remove y labels from 2nd plot
      ) +
      my_base_theme() 
    
    # Combine plots
    combined_plot <- p1 | p2
    if (SAVE.files ==TRUE) {
      boxplot_batch_correct <- glue("{data_path_bestageing2022}/output/plots/preprocessing/boxplot_batch_corr/001_d{disease}_boxplot_before_after_batch.svg")
      ggsave(filename = boxplot_batch_correct, plot = combined_plot, 
             width = 8, height = 6, 
             units = "in"  # default
      )
    }
    # pValuesLimma = f.pvalue(adjustedMatrix,design,null_model)
    # qValuesLimma = p.adjust(pValuesLimma,method="BH")
    # sum(qValuesLimma<0.05)
    
    runPCA_batch2 <- TRUE
    # PCA
    if(runPCA_batch2 == TRUE) {
      # pca
      pca_res <- prcomp(limma_edata[-1], center = TRUE, scale. = TRUE)
      
      # Get the first two principal components
      pca_df <- as_tibble(as.data.frame(pca_res$x[,1:5]))
      
      # add metadata
      exprs_metadat <- limma_edata %>% 
        bind_cols(pca_df) %>% 
        mutate(batch_center = batch_center, 
               disease = exprs_metadat$disease)
      
      # takes some time
      ggplot(exprs_metadat, aes(PC1, PC2, color = batch_center)) +
        geom_point(alpha=0.9, aes(shape=disease)) +
        labs(title = NULL,  # "PCA Plot Colored by Disease",
             x = "PC1",
             y = "PC2") +
        theme_minimal(base_size = 16, base_family = 'Arial')+
        scale_fill_manual(values = palette) +
        scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
        scale_shape_discrete(name=NULL, labels = c("Control", toupper(disease))) +
        my_base_theme() -> plot_PC1_PC2_corrected
      # plot_PC1_PC2_corrected
      
      ggplot(exprs_metadat, aes(PC1, PC3, color = batch_center)) +
        geom_point(alpha=0.9, aes(shape=disease)) +
        labs(title = NULL,  # "PCA Plot Colored by Disease",
             x = "PC1",
             y = "PC3") +
        theme_minimal(base_size = 16, base_family = 'Arial')+
        scale_fill_manual(values = palette) +
        scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
        scale_shape_discrete(name=NULL, labels = c("Control", toupper(disease))) +
        my_base_theme() -> plot_PC1_PC3_corrected
      # plot_PC1_PC3_corrected
      
      if (SAVE.files ==TRUE) {
        pca12_afterbatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001d_{disease}_pca12_after_batch.svg")
        ggsave(filename = pca12_afterbatch, plot = plot_PC1_PC2_corrected, 
               width = 8, height = 6, 
               units = "in"  # default
        )
        pca13_afterbatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001d_{disease}_pca13_after_batch.svg")
        ggsave(filename = pca13_afterbatch, plot = plot_PC1_PC3_corrected, 
               width = 6, height = 6, 
               units = "in"  # default
        )
      }
    }
  }
  
  if(runLimmaRemoveBatchEffect == TRUE) {
    data01 <- data01 %>% select(1:4) %>% 
      left_join(limma_edata, by = join_by(pat_id))
  }
  
  # SVA ---------------------------------------------------------------------
  if (runSVA == TRUE){
    exprs_data <- as.matrix(data01[,-c(1:4)])
    rownames(exprs_data) <- data01$pat_id
    pheno_data <- data01 %>% select(disease, age, sex) %>% as.data.frame()
    rownames(pheno_data) <- data01$pat_id
    
    # mean/ mode imputation necessary
    sum(is.na(pheno_data))
    pheno_data$age[is.na(pheno_data$age)] <- mean(pheno_data$age, na.rm = TRUE)
    
    #pheno_data$disease[is.na(pheno_data$disease)] <- mode_function(pheno_data$disease)
    pheno_data$sex[is.na(pheno_data$sex)] <- mode_function(pheno_data$sex)
    
    design <- model.matrix(~ disease + sex + age, data=pheno_data)
    null_model <- model.matrix(~ sex + age, data=pheno_data)
    
    #dim(exprs_data)
    #dim(pheno_data)
    #dim(design)
    #dim(null_model)
    
    n.sv <- num.sv(t(exprs_data), design, method="leek")
    if (n.sv != 0 & exists(x = "n.sv")) {
      sva_obj <- sva(t(exprs_data), design, null_model, n.sv=n.sv)
      
      sva_var_df <- sva_obj$sv
      colnames(sva_var_df) <- paste0("sva_", seq(n.sv))
      
      # update design matrix
      sva_design <- cbind(design, sva_var_df)
      sva_null_model <- cbind(null_model, sva_var_df)
      
      # append sva to adjust later
      data01 <- cbind(data01, sva_var_df)  # new cols 
    }
  }
  
  # if (SAVE.files ==TRUE) {
  #   filename.data01 <- glue("{data_path_bestageing2022}/output/de_results/{disease}/data01.rds")
  #   saveRDS(object = data01, file = filename.data01)
  # }
  
  ## SUBSET of literature miRNAs
  # selected_mirna_dat <- paste0("data01_", length(researchMiRNAAccession$miRNAName_v21), "mirnas")  
  # assign(x = selected_mirna_dat, 
  #        value = data01 %>% 
  #          select(pat_id, disease, age, sex, researchMiRNAAccession$miRNAName_v21))
  
  if (SAVE.files ==TRUE) {
    path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001d_{disease}_data01.rds")
    saveRDS(object = data01, file = path2dataprocessed)
  }
  
  
  ## RUN TESTS & log reg model  ----------------------------------------------------
  # create space to store results for ALL miRNAs
  all_filtered_mirnas <- data01 %>% select(-c(disease, age, sex))  # bring to same str as "all_filtered_mirnas" to work with existing code
  
  if (n.sv != 0 & exists(x = "n.sv")) {
    all_filtered_mirnas <- all_filtered_mirnas %>% 
      select(!all_of(colnames(sva_var_df)))
  }
  
  pval.t.test<-rep(NA,ncol(all_filtered_mirnas)-1)
  pval.u.test<-rep(NA,ncol(all_filtered_mirnas)-1)
  pval.glm <- rep(NA,ncol(all_filtered_mirnas)-1)
  pval.glm_sva  <- rep(NA,ncol(all_filtered_mirnas)-1)  # adjusting for sv1, sv2, ... sv_n
  pval.t.test.permute <- rep(NA,ncol(all_filtered_mirnas)-1)  # robust testing
  pval.glm_rob <- rep(NA,ncol(all_filtered_mirnas)-1)  # robust regression to remove outliers
  pval.glm_pca  <- rep(NA,ncol(all_filtered_mirnas)-1)  # including PC1 + PC2 + PC3
  pval.glm_rob_pca <- rep(NA,ncol(all_filtered_mirnas)-1)  # robust regression to remove outliers including PCs
  #average.difference <- rep(NA,ncol(all_filtered_mirnas)-1)
  log2FoldChange <- rep(NA,ncol(all_filtered_mirnas)-1)
  median.cont <- rep(NA,ncol(all_filtered_mirnas)-1)
  median.case <- rep(NA,ncol(all_filtered_mirnas)-1)
  mean.cont <- rep(NA,ncol(all_filtered_mirnas)-1)
  mean.case <- rep(NA,ncol(all_filtered_mirnas)-1)
  empse.cont <- rep(NA,ncol(all_filtered_mirnas)-1)
  empse.case <- rep(NA,ncol(all_filtered_mirnas)-1)
  aucs <- rep(NA,ncol(all_filtered_mirnas)-1)
  auc_glm <- rep(NA,ncol(all_filtered_mirnas)-1)
  auc_glm_sva  <- rep(NA,ncol(all_filtered_mirnas)-1)
  auc_glm_rob  <- rep(NA,ncol(all_filtered_mirnas)-1)
  auc_glm_pca <- rep(NA,ncol(all_filtered_mirnas)-1)
  auc_glm_rob_pca  <- rep(NA,ncol(all_filtered_mirnas)-1)
  
  name.mir <- rep(NA,ncol(all_filtered_mirnas)-1)
  
  # indexing
  cont.index <- data01$disease == "control"
  case.index <- data01$disease == disease
  # run tests
  
  
  if (runTests == TRUE){  #takes time
    total <- ncol(all_filtered_mirnas)-1
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for(miRNA in 1:(ncol(all_filtered_mirnas)-1)) {
      # update index since first colnames are [1] "pat_id"  "disease"  "age"   "sex"  "hsa_let_7a_3p" 
      miRNA_col <- miRNA+4
      cont <- data01[cont.index, miRNA_col] %>% as_vector()
      case <- data01[case.index, miRNA_col] %>% as_vector()
      # median for log-median expression plot (Figure 1 Vogel2013)
      mean.cont[miRNA] <- mean(cont)
      mean.case[miRNA] <- mean(case)
      median.cont[miRNA] <- median(cont)
      median.case[miRNA] <- median(case)
      empse.cont[miRNA] <- sd(cont)/ sqrt(length(cont))
      empse.case[miRNA] <- sd(case)/ sqrt(length(case))
      aucs[miRNA] <- suppressMessages(pROC::roc(controls=cont, cases=case)$auc[[1]])
      name.mir[miRNA] <- names(data01[miRNA_col])
      
      # average difference and logfold
      mean.control <- mean(cont)
      mean.case <- mean(case)
      log2FoldChange[miRNA] <- mean.case - mean.control
      ##log2FoldChange[miRNA] <- (mean.case - mean.control)/mean.control   # we are already on the log2 scale (https://support.bioconductor.org/p/117881/)
      # pvals
      pval.t.test[miRNA] <-t.test(cont,case)$p.value 
      pval.u.test[miRNA] <- wilcox.test(as.numeric(cont), as.numeric(case), exact = FALSE)$p.value
      # glm
      name_miRNA <- colnames(data01)[miRNA_col]
      f <- as.formula(paste("disease ~ ", name_miRNA, " + sex + age"))
      logreg <- glm(formula = f, data = data01, family = binomial(link = "logit") )
      predictions <- predict(logreg, newdata = data01[ , c("age", "sex", name_miRNA)], type = "response")
      auc_glm[miRNA] <- suppressMessages(roc(data01$disease, predictions)$auc[[1]])
      
      pval.glm[miRNA] <- coef(summary(logreg))[2,4]
      
      # 2023-09-05 SVA adjusted. again 
      if (n.sv != 0 & exists(x = "n.sv")) {
        f2 <- as.formula(paste("disease ~ ", name_miRNA, " + sex + age +", paste0(colnames(sva_var_df), collapse=" + ")))
        logreg2 <- glm(formula = f2, data = data01, family = binomial(link = "logit") )
        predictions2 <- predict(logreg2, newdata = data01[ , c("age", "sex", name_miRNA, paste0(colnames(sva_var_df)))], type = "response")
        auc_glm_sva[miRNA] <- suppressMessages(roc(data01$disease, predictions2)$auc[[1]])
        pval.glm_sva[miRNA] <- coef(summary(logreg2))[2,4]
      }
      
      # 2023-09-05 Permutation Test
      result <- permute_test(data01, name_miRNA, nperm = 500)
      pval.t.test.permute[miRNA] <- result$p_value
      
      
      # 2023-09-08 robust regression against outliers bc still a lot of pval inflation
      if (runRobustRegression == TRUE) {
        logreg_rob <- robustbase::glmrob(formula = f, data = data01, family = binomial(link = "logit") )
        predictions_rob <- predict(logreg_rob, newdata = data01[ , c("age", "sex", name_miRNA)],  type = "response")
        auc_glm_rob[miRNA] <- suppressMessages(roc(data01$disease, predictions_rob)$auc[[1]])
        pval.glm_rob[miRNA] <- coef(summary(logreg_rob))[2,4]
      }
      
      # including PC1, PC2, PC3 in glm
      if(runPCA_batch2 == TRUE) {
        combined_data <- cbind(data01, pca_df[, c("PC1", "PC2", "PC3")])
        formula_glm <- as.formula(paste("disease ~ ", name_miRNA, " + sex + age + PC1 + PC2 + PC3"))
        # fit least squares glm
        logreg_ls <- glm(formula = formula_glm, data = combined_data, family = binomial(link = "logit"))
        predictions_ls <- predict(logreg_ls, newdata = combined_data[ , c("age", "sex", name_miRNA, "PC1", "PC2", "PC3")], type = "response")
        
        # fit robust regression
        logreg_rob <- glmrob(formula = formula_glm, data = combined_data, family = binomial(link = "logit"))
        predictions_rob <- predict(logreg_rob, newdata = combined_data[ , c("age", "sex", name_miRNA, "PC1", "PC2", "PC3")], type = "response")
        # store results
        pval.glm_pca[miRNA] <- coef(summary(logreg_ls))[2,4]
        auc_glm_pca[miRNA] <- suppressMessages(roc(combined_data$disease, predictions_ls)$auc[[1]])
        pval.glm_rob_pca[miRNA] <- coef(summary(logreg_rob))[2,4]
        auc_glm_rob_pca[miRNA] <- suppressMessages(roc(combined_data$disease, predictions_rob)$auc[[1]])
      }
      
      # status bar
      setTxtProgressBar(pb, miRNA)
    }
    
    # GATHER Results -----------------------------------------------------------
    # we conducted 2549 t-tests and 2549 glm-models for each gene
    de.results <- tibble(miRNA = colnames(all_filtered_mirnas)[-1],  # all miRNAs without patID
                         # average.difference = average.difference,
                         log2FoldChange = log2FoldChange,
                         pval.t.test = pval.t.test,
                         pval.u.test = pval.u.test,
                         pval.glm = pval.glm,
                         pval.glm.sva = pval.glm_sva,
                         pval.t.test.permute =pval.t.test.permute, 
                         pval.glm_rob =pval.glm_rob,
                         aucs_univariate = aucs,
                         aucs_glm = auc_glm,
                         aucs_glm_sva = auc_glm_sva,
                         pval.glm_pca = pval.glm_pca,
                         auc_glm_pca = auc_glm_pca,
                         pval.glm_rob_pca =pval.glm_rob_pca,
                         auc_glm_rob_pca=auc_glm_rob_pca)
    
    # for figure 1 vogel 2013 
    adjusting_method <- "BH" # "holm" # "BH"
    results_logmedians <- tibble(miR =name.mir, auc=aucs, aucs_glm = auc_glm, aucs_glm_sva = auc_glm_sva,
                                 pval.t.test = pval.t.test, pval.u.test = pval.u.test, pval.glm = pval.glm, pval.glm_sva = pval.glm_sva, pval.t.test.permute=pval.t.test.permute,
                                 pval.glm_rob = pval.glm_rob, pval.glm_rob_pca =pval.glm_rob_pca,
                                 pval.glm_pca = pval.glm_pca, 
                                 auc_glm_pca = auc_glm_pca,
                                 auc_glm_rob_pca=auc_glm_rob_pca,
                                 logmedian.cont = median.cont, logmedian.case = median.case,
                                 logmean.cont = mean.cont, logmean.case = mean.case, empse.case, empse.cont) %>% 
      mutate(auc= ifelse(auc<0.5, 1-auc, auc)) %>% 
      # changed on 2023-08-02
      mutate(padj = p.adjust(pval.t.test, method = adjusting_method, n = length(name.mir)),  # inflation with "BH", use Bonferroni-Holm
             padj.u.test = p.adjust(pval.u.test, method = adjusting_method, n = length(name.mir)), 
             padj.glm = p.adjust(pval.glm, method = adjusting_method, n = length(name.mir)),
             padj.glm_sva = p.adjust(pval.glm_sva, method = adjusting_method, n = length(name.mir)),
             padj.t.test.permute = p.adjust(pval.t.test.permute, method = adjusting_method, n = length(name.mir)),
             padj.glm_pca = p.adjust(pval.glm_pca, method = adjusting_method, n = length(name.mir)),
             padj.glm_rob = p.adjust(pval.glm_rob, method = adjusting_method, n = length(name.mir)),
             padj.glm_rob_pca = p.adjust(pval.glm_rob_pca, method = adjusting_method, n = length(name.mir)),
      ) %>% 
      mutate(sign_indicator = ifelse(padj < 0.05, "p.adj≤0.05", "n.s."),
             sign_indicator_sva = ifelse(padj.glm_sva < 0.05, "p.adj≤0.05", "n.s."))
    
    if (SAVE.files ==TRUE) {
      filename.de.tibble <- glue("{data_path_bestageing2022}/output/de_results/{disease}/001d_de_results_batch_corrected.rds")
      saveRDS(object = de.results, file = filename.de.tibble)
      filename.logmedians <- glue("{data_path_bestageing2022}/output/de_results/{disease}/001d_results_logmedians_batch_corrected.rds")
      saveRDS(object = results_logmedians, file = filename.logmedians)
    }
    
    # inspect results
    
    # qqplot ----------------------------------------------------------
    pvalues <- results_logmedians$pval.t.test
    pvalues <- results_logmedians$pval.glm_pca
    pvalues <- results_logmedians$pval.glm_rob_pca
    pvalues <- results_logmedians$padj     # simple t-test
    pvalues <- results_logmedians$padj.glm
    pvalues <- results_logmedians$padj.glm_pca
    pvalues <- results_logmedians$padj.t.test.permute
    pvalues <- results_logmedians$padj.glm_rob
    pvalues <- results_logmedians$padj.glm_rob_pca
    
    qqplot_data <- data.frame(observed = -log10(sort(pvalues)),
                              expected = -log10(ppoints(length(pvalues))))
    
    # genomic inflation not for microarray data (GWAS!)
    qqplot_pvalues <- ggplot(qqplot_data, aes(x = expected, y = observed)) +
      geom_point(alpha=0.5, shape=16) +
      #coord_equal()+
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(x = "-log10(Expected P-values)",
           y = "-log10(Observed P-values)",
           title = NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6)) +
      my_base_theme()
    #qqplot_pvalues
    
    
    # qqplot with 6 adjusted pval qqplots
    # Creating a data frame of p-values in long format
    pvalues_long <- results_logmedians %>%
      select(padj, padj.glm, padj.glm_pca, padj.t.test.permute, padj.glm_rob, padj.glm_rob_pca) %>%
      pivot_longer(cols = everything(), names_to = "method", values_to = "pval")
    
    # Human-friendly labels
    method_labels <- c(
      padj = "T-Test",
      padj.glm = "GLM",
      padj.glm_pca = "GLM with PCA",
      padj.t.test.permute = "T-Test Permute",
      padj.glm_rob = "Robust GLM",
      padj.glm_rob_pca = "Robust GLM with PCA"
    )
    
    pvalues_long$method <- recode(pvalues_long$method, !!!method_labels)
    
    # Creating the QQ plot data
    qqplot_data <- pvalues_long %>%
      group_by(method) %>%
      mutate(observed = -log10(sort(pval)),
             expected = -log10(ppoints(length(pval))))
    
    # Making the QQ plot
    qqplot_pvalues <- ggplot(qqplot_data, aes(x = expected, y = observed, color = method)) +
      geom_point(alpha=0.5, shape=16) +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(x = "-log10(Expected P-values)",
           y = "-log10(Observed P-values)",
           title = NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      theme(legend.position = "bottom") +
      scale_color_manual(values = thematic::okabe_ito(6), name=NULL) +
      my_base_theme()
    
    # Counting the number of significant results for each method
    table_data <- pvalues_long %>%
      group_by(method) %>%
      summarise(significant = sum(pval < 0.05)) %>% 
      arrange(significant) %>% 
      rename(Method = method, 
             Significant=significant)
    
    # Create the table as a ggplot object
    table_plot <- ggtexttable(
      table_data, 
      rows = NULL, 
      theme = ttheme("mBlue")
    )
    
    # Combine the QQ plot and table
    combined_plot <- ggarrange(
      qqplot_pvalues,
      table_plot,
      ncol = 2, 
      #heights = c(1, 1),
      widths = c(2,1),
      align = "hv"
    )
    
    combined_plot
    
    if (SAVE.files ==TRUE) {
      qqplot_pvalues_path <- glue("{data_path_bestageing2022}/output/plots/qqplot/de/001d_{disease}_de_qqplot.svg")
      ggsave(filename = qqplot_pvalues_path, plot = combined_plot, 
             width = 10, height = 8, 
             units = "in"  # default
      )
    }
    
    if (n.sv != 0 & exists(x = "n.sv")) {
      pvalues <- results_logmedians$pval.glm_sva
      
      qqplot_data <- data.frame(observed = -log10(sort(pvalues)),
                                expected = -log10(ppoints(length(pvalues))))
      
      qqplot_pvalues_sva <- ggplot(qqplot_data, aes(x = expected, y = observed)) +
        geom_point(alpha=0.5, shape=16) +
        #coord_equal()+
        geom_abline(intercept = 0, slope = 1, color = "red") +
        labs(x = "-log10(Expected P-values)",
             y = "-log10(Observed P-values)",
             title = NULL) +
        theme_minimal(base_size = 16, base_family = 'Arial')+
        scale_fill_manual(values = thematic::okabe_ito(6)) +
        my_base_theme()
      qqplot_pvalues_sva
      
      if (SAVE.files ==TRUE) {
        qqplot_pvalues_path_sva <- glue("{data_path_bestageing2022}/output/plots/qqplot/de/001d_{disease}_de_qqplot_sva_included.svg")
        ggsave(filename = qqplot_pvalues_path_sva, plot = qqplot_pvalues_sva, 
               width = 6, height = 6, 
               units = "in"  # default
        )
      }
    }
  }
  
  print(paste0("|||-----------------------Run finished for disease: ", toupper(disease), " -----------------------|||"))
}





rm(disease)  # errors in plot otherwise because colname and varname


# Figures Vogel -----------------------------------------------------------

# changed from "rockerprojects/bestageing2022/data-literature/miRetrieve/top50mirnas_all_diseases.rds"
miRetrieve_alldiseases <- readRDS(glue("{data_path_bestageing2022}/data-literature/miRetrieve/2023-07-27_top50mirnas_all_diseases_pmids_gpt.rds")) # created in "scripts/miRetrieve/miRetrieve_topmirnas_all_diseases.R"
length(unique(miRetrieve_alldiseases$Accession))  ## 112

# UPDATE!! so that for each disease can be looked at all miRNAs above biomarker threshold/ not as before only top 50 miRNAS
miRetrieve_alldiseases <- readRDS(glue("{data_path_bestageing2022}/data-literature/miRetrieve/2023-09-20_ALL_bm_mirnas_all_diseases.rds"))
length(unique(miRetrieve_alldiseases$Accession))  ## 254

p2_listnew <- list()
for(mydisease in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])){
  
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001d_{all_combis$diseases[mydisease]}_data01.rds")
  
  if(!file.exists(path2dataprocessed)) {
    next
  }
  
  data01 <- readRDS(file = path2dataprocessed)
  de_results <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/001d_de_results_batch_corrected.rds"))
  results_logmedians <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/001d_results_logmedians_batch_corrected.rds"))
  
  # figure 01 vogel plot ------------------------------------------------------
  p1.1 <- ggplot(results_logmedians, aes(x = logmedian.cont, y = logmedian.case)) +
    geom_point(alpha=0.3, shape=16) +
    #geom_count(color="black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    labs(x = glue("controls\n[log-median expression]"), y = glue("{toupper(all_combis$diseases[mydisease])} patients\n[log-median expression]")) +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  p1.2 <- ggplot() +
    geom_histogram(data = results_logmedians, aes(x = padj.glm), fill = "white", color = "black", # bins = 500,
                   breaks = seq(0, 1, by = 0.05) ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.05), expand = c(0, 0)) +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
    labs(x = glue("GLM P-value\n({adjusting_method} adjusted)"), y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  p1.2
  p1.3 <- ggplot() +
    geom_histogram(data = results_logmedians, aes(x = auc), fill = "white", color = "black", binwidth = 0.01) +
    scale_x_continuous(breaks = seq(0.5, 0.7, by = 0.1), limits = c(0.45, 0.8), expand = c(0, 0)) +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
    labs(x = "Univariate AUC", y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  baseline_glm_model <- glm(formula = disease~age+sex, data = data01, family = binomial(link = "logit") )
  predictions <- predict(baseline_glm_model, newdata = data01[ , c("age", "sex")], type = "response")
  auc_glm_baseline <- suppressMessages(roc(data01$disease, predictions)$auc[[1]])
  p1.3b <- ggplot() +  # adjusted glm for age and sex
    geom_histogram(data = results_logmedians, aes(x = aucs_glm), fill = "white", color = "black", binwidth = 0.01) +
    scale_x_continuous(breaks = seq(0.6, 0.9, by = 0.1), limits = c(0.55, 0.85), expand = c(0, 0)) +
    geom_vline(xintercept = auc_glm_baseline, color = "red", linetype = "dashed", size = 1) +
    labs(x = "Univariate AUC", y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  # Arrange the plots into two rows
  plot_row1 <- p1.1 + labs(title = "A")
  plot_row2 <- ggarrange(p1.2 + labs(title = "B"), p1.3 + labs(title = "C"), nrow = 1, ncol = 2, widths = c(1, 1))
  plot_row2b <- ggarrange(p1.2 + labs(title = "B"), p1.3b + labs(title = "C"), nrow = 1, ncol = 2, widths = c(1, 1))
  # Combine the plots into a single plot with two rows
  fig1vogel2013 <- ggarrange(plot_row1, plot_row2, nrow = 2, heights = c(2, 1))
  fig1vogel2013b <- ggarrange(plot_row1, plot_row2b, nrow = 2, heights = c(2, 1))
  print(fig1vogel2013)
  
  if (SAVE.files ==TRUE) {
    fig1vogel2013_path <- glue("{data_path_bestageing2022}/output/plots/fig01vogel2013/001d_{toupper(all_combis$diseases[mydisease])}_linear_correlation.svg")
    fig1vogel2013_path_b <- glue("{data_path_bestageing2022}/output/plots/fig01vogel2013/001d_{toupper(all_combis$diseases[mydisease])}_linear_correlation_b.svg")
    ggsave(filename = fig1vogel2013_path, plot = fig1vogel2013, 
           width = 8, height = 8, 
           units = "in"  # default
    )
    ggsave(filename = fig1vogel2013_path_b, plot = fig1vogel2013b, 
           width = 8, height = 8, 
           units = "in"  # default
    )
  }
  
  ## figure 02 vogel 2013. ----------------------------
  # Step 1: select miRNAs from miRetrieve and disease
  disease_miRetrieve_mirnas <- miRetrieve_alldiseases %>% 
    mutate(Topic = toupper(Topic)) %>% 
    filter(Topic == toupper(all_combis$diseases[mydisease])) %>% 
    pull(TargetName)
  no_mirnas <- 10
  topNmiRetrieve <- results_logmedians %>% 
    # do not filter
    # filter(miR %in% disease_miRetrieve_mirnas) %>% 
    arrange(padj, desc(auc)) %>% 
    slice(1:no_mirnas)
  
  # long format
  rainbow_plot_dat <- data01 %>% 
    select(disease, all_of(topNmiRetrieve$miR)) %>% 
    pivot_longer(cols=-disease, names_to = "feature", values_to = "value") 
  
  # Calculate maximum values for each feature
  df_max <- rainbow_plot_dat %>% 
    group_by(feature, disease) %>% 
    summarise(max_value = max(value, na.rm = TRUE)) %>% 
    arrange(desc(max_value)) %>% 
    ungroup() %>% 
    select(-disease) %>% 
    distinct(feature, .keep_all = TRUE) %>% 
    # label
    left_join(topNmiRetrieve %>% select(miR, sign_indicator, padj), by=c("feature"="miR"))  %>% 
    mutate(feature = gsub("_", "-", feature))
  
  df_max <- df_max %>% 
    mutate(sign_indicator_asterisks = case_when(
      padj < 0.001 ~ "***",
      padj>= 0.001 & padj < 0.01 ~ "**",
      padj>= 0.01 & padj < 0.05 ~ "*",
      .default = "n.s."
    ))
  
  rainbow_plot_dat <- rainbow_plot_dat %>% 
    mutate(feature = gsub("_", "-", feature))
  
  # Order the levels of the "feature" variable based on the median values
  feature_order <- rainbow_plot_dat %>%
    group_by(feature) %>%
    summarise(median_value = median(value)) %>%
    arrange(desc(median_value)) %>%
    pull(feature)
  rainbow_plot_dat$feature <- factor(rainbow_plot_dat$feature, levels = feature_order)
  rainbow_plot_dat$disease <- factor(rainbow_plot_dat$disease, levels = c(all_combis$diseases[mydisease], "control"))
  
  # grouped box plots
  p2 <- ggplot(rainbow_plot_dat, aes(x=feature, y=value, fill = disease)) +
    geom_boxplot(position=position_dodge(width=0.8), outlier.shape = NA, alpha=0.7)+    
    geom_jitter(aes(color=disease),position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=0.1, alpha=0.2) +
    # prevent geom_text() from searching for the disease variable, you can override the fill aesthetic by setting fill = NULL within the aes() function 
    geom_text(data = df_max, aes(y = max_value, x=feature, fill = NULL), label = df_max$sign_indicator_asterisks, vjust = -0.5, size=2)+  
    labs(x = "", y = expression(paste("log"[2], " expression")), title = "") +
    scale_y_continuous(breaks = seq(0, ceiling(max(rainbow_plot_dat$value)), by = 2), #limits = c(0, 1), expand = c(0, 0)
    )+
    # theme_classic()+
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none") +
    scale_color_manual(values = thematic::okabe_ito(6)) +
    #scale_color_brewer(palette = "Set1") + # Choose a color palette
    #scale_fill_brewer(palette = "Set1", labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none")+
    my_base_theme()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p2)
  
  p2_listnew[[mydisease]] <- p2
  if (SAVE.files ==TRUE) {
    fig2_disease_vogel2013_path <- glue("{data_path_bestageing2022}/output/plots/fig02vogel2013/001d_{toupper(all_combis$diseases[mydisease])}_topdysregulatedmiRetrievemiR.svg")
    ggsave(filename = fig2_disease_vogel2013_path, plot = p2, 
           width = 8, height = 8, 
           units = "in"  # default
    )
  }
  #combine fig01 and fig02 for disease
  #plot_0102_combined <- ggarrange(fig1vogel2013b, p2 + labs(title = "D"), nrow = 2, ncol = 1, heights = c(1.5, 1))
  #plot_0102_combined
}

fig2vogel2013 <- ggarrange(plotlist = p2_listnew, ncol = 2, nrow = 2, 
                           labels = c("A", "B", "C", "D"))

# To print the plot
print(fig2vogel2013)
if (SAVE.files ==TRUE) {
  fig2vogel2013_path <- glue("{data_path_bestageing2022}/output/plots/fig02vogel2013/001d_arranged_topdysregulatedmiRetrievemiR.svg")
  ggsave(filename = fig2vogel2013_path, plot = fig2vogel2013, 
         width = 12, height = 12, 
         units = "in"  # default
  )
}



# DE volcanos ----------------------------------------------------------------
# moved here from "main.Rmd"
df_de_mirna_topic <- tibble(
  miRNA = character(0),
  Topic = character(0)
)

for(mydisease in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])) {
  
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001d_{all_combis$diseases[mydisease]}_data01.rds")
  
  if(!file.exists(path2dataprocessed)) {
    next
  }
  
  data01 <- readRDS(file = path2dataprocessed)
  # cave_dot here!!
  de.results <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/001d_de_results_batch_corrected.rds"))
  results_logmedians <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/001d_results_logmedians_batch_corrected.rds"))
  

  ### FDR calculation
  # OVERWRITE in script 001c and 001d
  de.results$padj <- p.adjust(de.results$pval.t.test, method = "holm", n = length(de.results$pval.t.test))  # "BH" p-val-inflation (t-test due to shape of volcano for EDA)
  de.results$padj.glm <- p.adjust(de.results$pval.glm, method = "holm", n = length(de.results$pval.t.test))
  de.results$padj.glm_pca <- p.adjust(de.results$pval.glm_pca, method = "holm", n = length(de.results$pval.t.test))
  
  ## manual calc 
  # n.comparisons <- length(de.results$pval.t.test)
  # inv.rank<-n.comparisons-rank(de.results$pval.t.test)
  # adjusted.results.t.test<-de.results$pval.t.test*(n.comparisons)/(n.comparisons-inv.rank+1)
  
  de.results$diffexpressed <- "NO"
  # if log2Foldchange > 1.0 and pvalue < 0.05, set as "UP" 
  de.results$diffexpressed[de.results$log2FoldChange > 1.0 & 
                             de.results$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de.results$diffexpressed[de.results$log2FoldChange < -1.0 & 
                             de.results$padj < 0.05] <- "DOWN"
  
  log2folds_thresholds <- c(1, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
  for (i in 1:length(log2folds_thresholds)) {
    if ( sum(de.results$diffexpressed != "NO") < 15 ) {
      de.results$diffexpressed[de.results$log2FoldChange > log2folds_thresholds[i] & 
                                 de.results$padj < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      de.results$diffexpressed[de.results$log2FoldChange < - log2folds_thresholds[i] & 
                                 de.results$padj < 0.05] <- "DOWN"
      log2folds_thresh <- log2folds_thresholds[i]
    } else{ 
      log2folds_thresh <- log2folds_thresholds[i]
      break
    }
  }
  
  ## label miRNAs 
  # Step 1: select miRNAs from miRetrieve and disease
  # UPDATE NEW 2023-09, not only top50 per disease but all above miRetrieve biomarker threshold
  disease_miRetrieve_mirnas <- miRetrieve_alldiseases %>% 
    mutate(Topic = toupper(Topic)) %>% 
    filter(Topic == toupper(all_combis$diseases[mydisease])) %>% 
    pull(TargetName)
  
  # init
  de.results$delabel <- NA
  de.results$delabel[de.results$diffexpressed != "NO"] <- 
    de.results$miRNA[de.results$diffexpressed != "NO"]
  # ONLY label miRNAs that are also known from LITERATURE 
  ## old literature miRNAs 
  # de.results$delabel <- ifelse(de.results$delabel %in% researchMiRNAAccession$miRNAName_v21,  de.results$delabel, NA)
  # NEW from miRetrieve
  de.results$delabel <- ifelse(de.results$delabel %in% disease_miRetrieve_mirnas,  de.results$delabel, NA)
  
  ## color all literature miRNAs ----------------------------------------
  #de.results$research_mirna <- ifelse(de.results$miRNA %in% researchMiRNAAccession$miRNAName_v21,  "Literature miRNA", "Not found in Literature")
  de.results$research_mirna <- ifelse(de.results$miRNA %in% disease_miRetrieve_mirnas,  "Literature miRNA", "Not found in miRetrieve Top Hits")
  
  ## # label miRNA 106b-5p (special interest)
  ## freys.mirnas <- c("hsa_mir_106b_3p", "hsa_mir_106b_5p")
  ## de.results$freysmirna <- ifelse(de.results$miRNA %in% freys.mirnas, de.results$miRNA, NA)
  ## de.results$freysmirna.col <- ifelse(de.results$miRNA %in% freys.mirnas, de.results$miRNA, "not selected")
  
  ## PLOT ##
  # custom volcano -----------------------------------------------------------
  subtitle.custom <- paste0("log2foldchange = \u00b1 ", log2folds_thresh, " for significant miRNAs (p_adj < 0.05)")
  
  de.results %>% 
    mutate(delabel = gsub("_", "-", delabel)) %>% 
    ggplot(aes(x=log2FoldChange, y=-log10(pval.t.test), col=diffexpressed, label=delabel)) + 
    geom_point(alpha=ifelse(de.results$diffexpressed != "NO", 1, 0.2), size=1, shape=16) + 
    geom_vline(xintercept=c(-log2folds_thresh, log2folds_thresh), alpha=0.3) +
    geom_hline(
      yintercept=-log10(0.05/length(pval.t.test)),  # bonferroni only correction
      alpha=0.3
    ) + 
    labs(
      title = NULL, # paste0("Volcano plot of DE miRNAs in ", disease_vector[mydisease]),
      subtitle = NULL, #subtitle.custom
    ) +
    xlab("log2 fold change") + 
    ylab("-log10 (p-value)") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))+
    #ggthemes::scale_color_few(name="Differentially expressed")+
    ggrepel::geom_text_repel(show.legend = FALSE, size=5)+
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    scale_color_manual(values = thematic::okabe_ito(6), name="Differentially expressed")+
    my_base_theme() -> volcano_1
  print(volcano_1)
  
  if (SAVE.files ==TRUE ){
    filename.volcano <- glue("{data_path_bestageing2022}/output/plots/de_analysis/{all_combis$diseases[mydisease]}/001d_de_volcano.svg")
    ggsave(filename = filename.volcano, plot = volcano_1, 
           width = 9, height = 5, 
           units = "in"  # default
    )
  }
  
  # Color all literature miRNAs
  de.results %>% 
    mutate(delabel = gsub("_", "-", delabel)) %>% 
    ggplot(aes(x=log2FoldChange, y=-log10(pval.t.test), col=research_mirna, label=delabel)) + 
    geom_point(alpha=ifelse(de.results$research_mirna == "Literature miRNA", 1, 0.2), size=1, shape=16) + 
    geom_vline(xintercept=c(-log2folds_thresh, log2folds_thresh), alpha=0.3) +
    geom_hline(
      yintercept=-log10(0.05/length(pval.t.test)),  # bonferroni only correction
      alpha=0.3
    ) + 
    labs(
      title = NULL, # paste0("Volcano plot of DE miRNAs in ", disease_vector[mydisease]),
      subtitle = NULL, # subtitle.custom
    ) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))+
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    scale_color_manual(values = thematic::okabe_ito(6), name=NULL)+
    ggrepel::geom_text_repel(show.legend = FALSE, size=5) -> volcano_2
  print(volcano_2)
  
  if (SAVE.files ==TRUE ){
    filename.volcano2 <-glue("{data_path_bestageing2022}/output/plots/de_analysis/{all_combis$diseases[mydisease]}/001d_de_volcano_LITERATmiRNAs_ANNOT.svg")
    ggsave(filename = filename.volcano2, plot = volcano_2, 
           width = 9, height = 5, 
           units = "in"  # default
    )
  }
  
  # append list for venn diagram (intersection of de mirnas in different diseases)
  df_de_mirna_topic_tmp <- de.results %>% 
    filter(padj.glm < 0.05) %>% 
    select(miRNA) %>% 
    mutate(Topic = all_combis$diseases[mydisease])
  df_de_mirna_topic <- df_de_mirna_topic %>% 
    bind_rows(df_de_mirna_topic_tmp)
}

# Venn diagramm 4 miRetrieve -------------------------------------------------------------


# Venn plot DE-mirnas shared ----------------------------------------------

list_venn_de_mirnas <- list(
  ACS = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "acs"],
  CAD = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "cad"],
  DCM = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "dcm"]
  #HFrEF = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "hfref"]
)

# Create the Venn diagram
venn_plot_de_mirnas <- ggvenn(list_venn_de_mirnas) +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6)) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
venn_plot_de_mirnas

if (SAVE.files ==TRUE ){
  filename.venn_de_mirnas_shared <-glue("{data_path_bestageing2022}/output/plots/venn_dia/001d_venn_de_mirnas_shared.svg")
  ggsave(filename = filename.venn_de_mirnas_shared, plot = venn_plot_de_mirnas, 
         width = 8, height = 8, 
         units = "in"  # default
  )
}

