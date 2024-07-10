

require(dplyr)
require(tidyr)
require(stringr)
require(tibble)
require(forcats)
require(janitor)
require(readxl)
require(glue)
require(RColorBrewer)
require(ggthemes)
require(limma)
require(ggpubr)
require(tidytext)  # for ordering facets https://juliasilge.github.io/tidytext/reference/reorder_within.html 
require(patchwork)
require(broom)
require(tidymodels)
require(finetune)

# Define library and data paths based on system
#

# Get system name
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

# fns -------------------------------------

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))

convert_mir_name <- function(name) {
  name <- gsub("hsa_", "", name)  # Remove 'hsa-' prefix
  name <- gsub("mir", "miR", name)  # Replace 'mir' with 'miR'
  name <- gsub("_", "-", name)  # Replace underscores with hyphens
  return(name)
}
convert_mir_name_V <- Vectorize(convert_mir_name)

remove_digits <- function(vec) {
  # Use regex to replace patterns like '-DIGITS' with an empty string
  return(sub("-\\d+$", "", vec))
}

# functions ---------------------------------------------------------------
mode_function <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

# set analysis parameters -------------------
filterDetMatrix <- TRUE
runLimmaRemoveBatchEffect <- TRUE
SAVE.files <-  TRUE


# load data ---------------------------------
echo_abfrage_ba_join_nearest_exam <- readRDS(file = "./data/abfrage2024/echo_abfrage_ba_join_nearest_exam.rds")

lab2023 <- read_excel(path = glue("{data_path_bestageing2022}/data/elham2023/2023-11-15-elham_metadat_combined_LAB2023_b.xlsx")) %>% 
  select(best_ageing_code, HSTNT, HSTNTHP, NTBNP, INR, BILI, KREA, LEUKO, CRP, HB, CHOL) %>% 
  mutate(
    HSTNT = case_when(is.na(HSTNT) ~ HSTNTHP,
                      .default = HSTNT)
  ) %>% 
  # drop duplicates
  distinct(best_ageing_code, .keep_all = TRUE) %>% 
  # convert to numeric and change "," -> "." for decimals
  mutate(across(!c(best_ageing_code), function(x) as.numeric(gsub(",", ".", x)) )) %>% 
  select(-HSTNTHP) # is second TNT, first TNT measured in HSTNT, kinetics


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
# detection matrix
det_mat_all_mirnas <- read.table(glue("{data_path_bestageing2022}/data/kahraman2023/det_mat_all_mirnas.txt")) %>% 
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

# also cleaned diagnoses data for table 01 to get unique identifiers (e.g. could have diagnoses CAD and HFREF.. )
clean_all_meta_2023_table01 <- readRDS(file = glue("{data_path_bestageing2022}/data/disease_identifier_table01/clean_all_meta_2023_table01.rds"))


# preprocess miRNA data -----------------------------------------------------

clean_all_meta_2023_table01 <- readRDS(file = glue("{data_path_bestageing2022}/data/disease_identifier_table01/clean_all_meta_2023_table01.rds"))

analysis_dat_ukhd <- all_mirnas %>%
  inner_join(clean_all_meta_2023_table01 %>% select(patID, disease, sex, age), by = c("pat_id" = "patID")) %>% 
  relocate(disease, sex, age, .after = pat_id) %>% 
  filter(
    grepl("UKL-HD", pat_id)
  )

analysis_dat_ukhd$disease %>% table()

## FILTERING DET MATRIX ----------------------------------------------------
if (filterDetMatrix == TRUE){
  filter_dat <- analysis_dat_ukhd %>% select(-c(age,sex))
  
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
  print(density_plot_filter_eda_stratified)
  
  # density_ridges_plot_filter_eda_stratified <- ggplot(data_long, aes(x=expression, y=disease, fill=disease)) + 
  #   geom_density_ridges(alpha=0.5) +  # , jittered_points = TRUE, position = position_points_jitter(width = 0.05)) + # too much data involved if points plotted
  #   labs(x="Log2 Expression", y="Density") +  # title="Density Ridge plot of log2 expressions", 
  #   theme_minimal(base_size = 16, base_family = 'Arial')+
  #   scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
  #   theme(axis.text.y = element_blank())
  # #print(density_ridges_plot_filter_eda_stratified)
  
  
  threshold <- 0.7  # not too conservative
  
  # Filter 
  lowExprAll <- colMeans(det_mat_all_mirnas_tmp[, -(1:2)]) <= threshold
  sum(!lowExprAll)
  
  message(glue("After filtering by detection matrix we receieved {sum(!lowExprAll)} miRNA features which we will use for our analysis."))
  
  # Filter for disease samples
  #lowExprDisease <- colMeans(det_mat_all_mirnas_tmp[det_mat_all_mirnas_tmp$disease == disease, -(1:2)]) <= threshold
  
  # Identify miRNAs that are lowly expressed in both conditions
  #toFilter <- lowExprControl & lowExprDisease
  
  # Filter out those miRNAs from the dataset
  filteredData <- filter_dat[, c(TRUE, TRUE, !lowExprAll)] ## TRUE TRUE corresponds to keep pat_id and disease
  
  data01 <- analysis_dat_ukhd %>% 
    select(1:4) %>% 
    dplyr::bind_cols(filteredData[, -(1:2)])
}


## limmma BatchEffect ----------------------------------------

if(runLimmaRemoveBatchEffect == TRUE) {
  batch_center <- remove_digits(data01$pat_id)
  table(batch_center, data01$disease)
  exprs_data <- as.matrix(data01[,-c(1:4)])
  
  # add metadata
  exprs_metadat <- data01 %>% 
    # batch 1
    mutate(batch_center = batch_center) %>% 
    # batch 2
    left_join(hbdx_metadat %>% select(customer_id, slide_id, sample_qc, array_qc),
              by = c("pat_id"="customer_id"))
  batch_slide_id <- exprs_metadat %>% pull(slide_id)
  
  # pheno data
  pheno_data <- data01 %>% 
    select(disease, age, sex) %>% 
    as.data.frame()
  rownames(pheno_data) <- data01$pat_id
  pheno_data$age[is.na(pheno_data$age)] <- mean(pheno_data$age, na.rm = TRUE)
  pheno_data$sex[is.na(pheno_data$sex)] <- mode_function(pheno_data$sex)
  
  
  ### PCA after batch removal -------------------------------------------------
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
    
    #palette <- brewer.pal(12, "Set3")
    palette <- gdocs_pal()(6)
    # takes some time
    exprs_metadat %>% 
      mutate(disease = toupper(disease) ) %>% 
      ggplot(aes(PC1, PC2, color = disease)) +
      geom_point(alpha=0.7, aes(shape=batch_center)) +
      #coord_fixed(ratio=1)+
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC2") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL) -> plot_PC1_PC2
    # plot_PC1_PC2
    
    exprs_metadat %>% 
      mutate(disease = toupper(disease) ) %>% 
      ggplot(aes(PC1, PC3, color = disease)) +
      geom_point(alpha=0.7, aes(shape=batch_center)) +
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC3") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL,) -> plot_PC1_PC3
    #plot_PC1_PC3
    
    if (SAVE.files ==TRUE) {
      pca12_beforebatch <- glue("{data_path_bestageing2022}/output/analysis2024/preprocessing/001_pca12_before_batch.svg")
      ggsave(filename = pca12_beforebatch, plot = plot_PC1_PC2, 
             width = 8, height = 6, 
             units = "in"  # default
      )
      pca13_beforebatch <- glue("{data_path_bestageing2022}/output/analysis2024/preprocessing/001_pca13_before_batch.svg")
      ggsave(filename = pca13_beforebatch, plot = plot_PC1_PC3, 
             width = 8, height = 6, 
             units = "in"  # default
      )
    }
  }
  
  ## run LIMMA -------------------------------------------
  # run limma batch remove https://evayiwenwang.github.io/Managing_batch_effects/adjust.html#accounting-for-batch-effects
  # modlimma <- model.matrix(~1, data=pheno_data) 
  modlimma <- model.matrix( ~ 1, data=pheno_data)  # do not use disease here this time
  adjustedMatrix <- removeBatchEffect(t(exprs_data), batch=batch_slide_id, design = modlimma)  # collinearity issues with center batch:  one-to-one mapping with any condition
  limma_edata <-  adjustedMatrix %>% t() %>%  as.data.frame() %>% as_tibble() %>% 
    mutate(pat_id = data01$pat_id) %>% 
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
    theme(axis.text.x = element_blank())
  # Plot batch-corrected data
  p2 <- ggplot(df_corrected, aes(x = Sample, y = Expression)) +
    geom_boxplot(aes(group = Sample), alpha=0.5) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    labs(title = NULL, x=NULL) +
    ylim(y_limits) +
    scale_color_manual(values = palette, name=NULL) +
    theme(axis.text.x = element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank()  # remove y labels from 2nd plot
    )
  
  # Combine plots
  combined_plot <- p1 | p2
  if (SAVE.files ==TRUE) {
    boxplot_batch_correct <- glue("{data_path_bestageing2022}/output/analysis2024/preprocessing/001_boxplot_before_after_batch.svg")
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
             disease = data01$disease)
    
    # takes some time
    ggplot(exprs_metadat, aes(PC1, PC2, color = disease)) +
      geom_point(alpha=0.9, aes(shape=batch_center)) +
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC2") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL) -> plot_PC1_PC2_corrected
    # plot_PC1_PC2_corrected
    
    ggplot(exprs_metadat, aes(PC1, PC3, color = disease)) +
      geom_point(alpha=0.9, aes(shape=batch_center)) +
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC3") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL)  -> plot_PC1_PC3_corrected
    # plot_PC1_PC3_corrected
    
    if (SAVE.files ==TRUE) {
      pca12_afterbatch <- glue("{data_path_bestageing2022}/output/analysis2024/preprocessing/001_pca12_after_batch.svg")
      ggsave(filename = pca12_afterbatch, plot = plot_PC1_PC2_corrected, 
             width = 8, height = 6, 
             units = "in"  # default
      )
      pca13_afterbatch <- glue("{data_path_bestageing2022}/output/analysis2024/preprocessing/001_pca13_after_batch.svg")
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

## STORE PROCESSED DATA -------------------------------------------------------------------

data01$age[is.na(data01$age)] <- mean(data01$age, na.rm = TRUE)
data01$sex[is.na(data01$sex)] <- mode_function(data01$sex)

if (SAVE.files ==TRUE) {
  path2dataprocessed <- glue("{data_path_bestageing2022}/output/analysis2024/data/001_data01_2024_features.rds")
  saveRDS(object = data01, file = path2dataprocessed)
}

###
# ANALYSIS ----------------------------------------------------------
###

echo_dat <- echo_abfrage_ba_join_nearest_exam %>% 
  mutate(max_wall_thickness = pmax(Septum, Hinterwand, na.rm=TRUE) ) %>% 
  mutate(SV_Funktion_MAPSE = as.numeric(gsub(",", ".", SV_Funktion_MAPSE)) ) %>%
  select(
    best_ageing_code.y, diff, 
    Aortenwurzel, Linker_Vorhof, max_wall_thickness, 
    LV_EDD, LV_EF, RV, sPA, SV_Funktion_EE_Ruhe, SV_Funktion_MAPSE
  )

analysis_dat_ukhd2024_merged <-  data01 %>% 
  left_join(echo_dat, by=c("pat_id"="best_ageing_code.y")) %>% 
  left_join(lab2023, by=c("pat_id" = "best_ageing_code")) %>% 
  relocate(diff:CHOL, .after=age) |> 
  # filter echo/ lab within half a year
  filter(diff < 365/2)

colnames(analysis_dat_ukhd2024_merged)
dim(analysis_dat_ukhd2024_merged)

# Define the range of features
echo_traits <- names(analysis_dat_ukhd2024_merged)[6:14]  # Echo traits
lab_traits <- names(analysis_dat_ukhd2024_merged)[15:23]  # Lab traits
miRNAs <- names(analysis_dat_ukhd2024_merged)[24:518]  # miRNA features

# SANITY CHECK ECHO OUTLIERS !!!
analysis_dat_ukhd2024_merged <- analysis_dat_ukhd2024_merged %>%
  mutate(
    Aortenwurzel = ifelse(Aortenwurzel < 20 | Aortenwurzel > 55, NA, Aortenwurzel),
    Linker_Vorhof = ifelse(Linker_Vorhof < 25 | Linker_Vorhof > 60, NA, Linker_Vorhof),
    max_wall_thickness = ifelse(max_wall_thickness < 5 | max_wall_thickness > 30, NA, max_wall_thickness),
    LV_EDD = ifelse(LV_EDD < 30 | LV_EDD > 90, NA, LV_EDD),
    LV_EF = ifelse(LV_EF < 10 | LV_EF > 75, NA, LV_EF),
    RV = ifelse(RV < 20 | RV > 50, NA, RV),
    sPA = ifelse(sPA < 15 | sPA > 95, NA, sPA),
    SV_Funktion_EE_Ruhe = ifelse(SV_Funktion_EE_Ruhe < 3 | SV_Funktion_EE_Ruhe > 35, NA, SV_Funktion_EE_Ruhe),
    SV_Funktion_MAPSE = ifelse(SV_Funktion_MAPSE < 0.5 | SV_Funktion_MAPSE > 3, NA, SV_Funktion_MAPSE)
  )


# Define the number of top features to display
m <- 20

# CORRELATION ANALYSIS --------------------------------------------
calculate_correlations <- function(traits, miRNAs, data) {
  cor_results <- data.frame(Trait = character(), miRNA = character(), Correlation = numeric(), 
                            P_Value = numeric(), Conf_Lower = numeric(), Conf_Upper = numeric(), 
                            Sample_Size = integer(), stringsAsFactors = FALSE)
  
  for (trait in traits) {
    for (miRNA in miRNAs) {
      subset_data <- data %>% select(all_of(trait), all_of(miRNA)) %>% drop_na()
      sample_size <- nrow(subset_data)
      
      if (nrow(subset_data) > 2) {
        cor_test <- cor.test(subset_data[[trait]], subset_data[[miRNA]], method = "pearson")
        cor_results <- rbind(cor_results, data.frame(
          Trait = trait, 
          miRNA = miRNA, 
          Correlation = cor_test$estimate, 
          P_Value = cor_test$p.value,
          Conf_Lower = cor_test$conf.int[1],
          Conf_Upper = cor_test$conf.int[2],
          Sample_Size = sample_size
        ))
      }
    }
  }
  
  cor_results <- cor_results %>%
    group_by(Trait) %>%
    mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH")) %>%
    ungroup() %>%
    mutate(Significance = case_when(
      Adjusted_P_Value < 0.001 ~ "***",
      Adjusted_P_Value < 0.01 ~ "**",
      Adjusted_P_Value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  return(cor_results)
}

# Calculate correlations for echo and lab traits
echo_cor_results <- calculate_correlations(echo_traits, miRNAs, analysis_dat_ukhd2024_merged)
lab_cor_results <- calculate_correlations(lab_traits, miRNAs, analysis_dat_ukhd2024_merged)

# Filter the top significant features for visualization
top_echo_features <- echo_cor_results %>%
  group_by(Trait) %>%
  arrange(Adjusted_P_Value) %>%
  slice_head(n = m) %>%
  ungroup()

top_lab_features <- lab_cor_results %>%
  group_by(Trait) %>%
  arrange(Adjusted_P_Value) %>%
  slice_head(n = m) %>%
  ungroup()

# Display the top features for echo traits
print(top_echo_features)

# Display the top features for lab traits
print(top_lab_features)



### Dotplot with errorbars | Correlations -------------------------------------

palette <- ggthemes::gdocs_pal()(10)

trait_name_mapping <- c(
  "Aortenwurzel" = "Aortenwurzel",
  "Linker_Vorhof" = "Linker Vorhof",
  "max_wall_thickness" = "Max Wall Thickness",
  "LV_EDD" = "LV EDD",
  "LV_EF" = "LV EF",
  "RV" = "RV",
  "sPA" = "sPA",
  "SV_Funktion_EE_Ruhe" = "E/E'-Ratio",
  "SV_Funktion_MAPSE" = "MAPSE"
)

rename_traits <- function(trait) {
  return(trait_name_mapping[trait])
}
rename_traits_V <- Vectorize(rename_traits)

# Function to prepare data for ordered plotting
prepare_plot_data <- function(cor_results) {
  cor_results %>%
    mutate(miRNA = convert_mir_name_V(miRNA)) %>%
    group_by(Trait) %>%
    arrange(Correlation) %>%
    mutate(miRNA = factor(miRNA, levels = unique(miRNA))) %>%
    ungroup() 
}

top_echo_features_ordered <- prepare_plot_data(top_echo_features)
top_lab_features_ordered <- prepare_plot_data(top_lab_features)

# Calculate the sample size for each trait
sample_sizes_echo <- top_echo_features_ordered %>%
  group_by(Trait) %>%
  summarise(Sample_Size = first(Sample_Size)) %>%
  mutate(New_Trait = paste(rename_traits(Trait), "\nn =", Sample_Size))

sample_sizes_lab <- top_lab_features_ordered %>%
  group_by(Trait) %>%
  summarise(Sample_Size = first(Sample_Size)) %>%
  mutate(New_Trait = paste(Trait, "\nn =", Sample_Size))


# Merge the modified trait names back into the main data
top_echo_features_ordered <- top_echo_features_ordered %>%
  left_join(sample_sizes_echo %>% select(Trait, New_Trait), by = "Trait") 

top_lab_features_ordered <- top_lab_features_ordered %>%
  left_join(sample_sizes_lab %>% select(Trait, New_Trait), by = "Trait") 

# Visualize the correlations for echo traits
echo_plot <- top_echo_features_ordered %>% 
  ggplot(aes(reorder_within(x = miRNA, by=Correlation, within=New_Trait), y = Correlation, ymin = Conf_Lower, ymax = Conf_Upper)) +
  geom_point(size = 3, aes(color = New_Trait)) +
  geom_errorbar(width = 0.2, aes(color = New_Trait)) +
  geom_text(aes(label = Significance), vjust = -3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_reordered() +
  facet_wrap(~ New_Trait, scales = "free_x") +
  labs(
    #title = "Top miRNA-Echo Trait Correlations", 
    x = "miRNA", 
    y = "Correlation Coefficient"
  )+
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_color_manual(values = palette, name=NULL) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 
print(echo_plot)

# Visualize the correlations for lab traits
lab_plot <- top_lab_features_ordered %>% 
  ggplot(aes(reorder_within(x = miRNA, by=Correlation, within=New_Trait), y = Correlation, ymin = Conf_Lower, ymax = Conf_Upper)) +
  geom_point(size = 3, aes(color = New_Trait)) +
  geom_errorbar(width = 0.2, aes(color = New_Trait)) +
  geom_text(aes(label = Significance), vjust = -3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_reordered() +
  facet_wrap(~ New_Trait, scales = "free_x") +
  labs(
    #title = "Top miRNA-Lab Trait Correlations", 
    x = "miRNA", 
    y = "Correlation Coefficient"
  ) +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_color_manual(values = palette, name=NULL) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 
print(lab_plot)

echo_cor_feature_plot <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_echo_cor_feature.svg")
ggsave(filename = echo_cor_feature_plot, plot = echo_plot, 
       width = 18, height = 15, 
       units = "in"  # default
)
lab_cor_feature_plot <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_lab_cor_feature.svg")
ggsave(filename = lab_cor_feature_plot, plot = lab_plot, 
       width = 18, height = 15, 
       units = "in"  # default
)

###
## Volcano plot adaption ----------------------
###

### a) for all traits combined ----------------

# Function to create a volcano plot with two colors based on significance level
create_volcano_plot <- function(cor_results, title) {
  # Calculate the significance level for the plot
  cor_results <- cor_results %>%
    mutate(Significance_Level = -log10(Adjusted_P_Value),
           # Mark significant based on Adjusted P-Value
           Is_Significant = ifelse(Adjusted_P_Value < 0.05, "Significant", "Not Significant"))
  
  # Create and return the volcano plot
  plot <- ggplot(cor_results, aes(x = Correlation, y = Significance_Level)) +
    geom_point(aes(color = Is_Significant), alpha = 0.7) +
    scale_color_manual(values = c("Significant" = "#dc3912", "Not Significant" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    labs(x = "Correlation Coefficient", y = "-log10 Adjusted P-Value", title = title) +
    theme_minimal(base_size = 16, base_family = 'Arial') +
    theme(legend.title = element_blank(), legend.position = "none")
  
  return(plot)
}

# Create and print volcano plot for echo traits
echo_plot <- create_volcano_plot(echo_cor_results, "") # "Volcano Plot of miRNA Correlations for Echo Traits")
print(echo_plot)

volcano_cor_echo <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_volcano_cor_echo.svg")
ggsave(filename = volcano_cor_echo, plot = echo_plot, 
       width = 5, height = 4, 
       units = "in"  # default
)

# Create and print volcano plot for lab traits
lab_plot <- create_volcano_plot(lab_cor_results, "") # "Volcano Plot of miRNA Correlations for Lab Traits")
print(lab_plot)

volcano_cor_lab <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_volcano_cor_lab.svg")
ggsave(filename = volcano_cor_lab, plot = lab_plot, 
       width = 5, height = 4, 
       units = "in"  # default
)


###  b) facet wrap traits ----------------

create_combined_volcano_plot_facet_wrap <- function(data, title) {
  # Calculate the significance level for the plot
  data <- data %>%
    mutate(Significance_Level = -log10(Adjusted_P_Value),
           Is_Significant = ifelse(Adjusted_P_Value < 0.05, "Significant", "Not Significant"))
  
  # Create the combined volcano plot using facet_wrap
  plot <- ggplot(data, aes(x = Correlation, y = Significance_Level)) +
    geom_point(aes(color = Is_Significant), alpha = 0.7) +
    scale_color_manual(values = c("Significant" = "#dc3912", "Not Significant" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    labs(x = "Correlation Coefficient", y = "-log10 Adjusted P-Value", title = title) +
    theme_minimal(base_size = 16, base_family = 'Arial') +
    theme(legend.title = element_blank(), legend.position = "none") +
    facet_wrap(~ Trait, scales = "free_y")  # Use facet_wrap to create individual plots for each trait
  
  return(plot)
}

# Create and print the combined volcano plot
volcano_facet_wrap_echo <- create_combined_volcano_plot_facet_wrap(echo_cor_results, "")
print(volcano_facet_wrap_echo)

volcano_cor_echo_facet <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_volcano_facet_cor_echo.svg")
ggsave(filename = volcano_cor_echo_facet, plot = volcano_facet_wrap_echo, 
       width = 8, height = 8, 
       units = "in"  # default
)


volcano_facet_wrap_lab <- create_combined_volcano_plot_facet_wrap(lab_cor_results, "")
print(volcano_facet_wrap_lab)

volcano_cor_lab_facet <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_volcano_facet_cor_lab.svg")
ggsave(filename = volcano_cor_lab_facet, plot = volcano_facet_wrap_lab, 
       width = 8, height = 8, 
       units = "in"  # default
)


### Heatmap miRNA ~ traits ---------------------------------------------------

create_heatmap <- function(cor_results) {
  # Pivot to a wide format suitable for heatmaps
  heatmap_data <- cor_results %>%
    dplyr::select(miRNA, Trait, Correlation) %>%
    reshape2::acast(miRNA ~ Trait, value.var = "Correlation")
  
  heatmap_data <- as.data.frame(heatmap_data)
  rownames(heatmap_data) <- heatmap_data$miRNA
  heatmap_data$miRNA <- NULL
  
  pheatmap::pheatmap(heatmap_data,
                     color = colorRampPalette(c("#3366cc", "white", "#dc3912"))(100),
                     show_rownames = FALSE,  # Do not show row names (miRNA)
                     show_colnames = TRUE,  # Show column names (Traits)
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean", 
                     fontsize = 16,
                     angle_col = 45)
}


heatmap_cor_coefficient_traits <- create_heatmap(combined_cor_results)
heatmap_cor_coefficient_traits

heatmap_cor_traits <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_heatmap_cor_features.svg")
ggsave(filename = heatmap_cor_traits, plot = heatmap_cor_coefficient_traits, 
       width = 8, height = 8, 
       units = "in"  # default
)



###
# Univariable REGRESSION Analysis-------------------------------------------------------------------------
### 

# Function to calculate regression coefficients and prepare results
calculate_regressions <- function(traits, miRNAs, data) {
  reg_results <- data.frame(Trait = character(), miRNA = character(), 
                            Regression_Coefficient = numeric(), 
                            P_Value = numeric(), Conf_Lower = numeric(), 
                            Conf_Upper = numeric(), stringsAsFactors = FALSE)
  
  for (trait in traits) {
    for (miRNA in miRNAs) {
      # Ensure sex and age are also selected for the subset data
      subset_data <- data %>% select(all_of(c(trait, miRNA, "sex", "age"))) %>% drop_na()
      
      if (nrow(subset_data) > 2) {
        # Fit linear model: trait ~ miRNA + sex + age
        model <- lm(formula = paste(trait, "~", miRNA, "+ sex + age"), data = subset_data)
        # Use broom::tidy with conf.int = TRUE to get confidence intervals
        tidy_model <- broom::tidy(model, conf.int = TRUE)
        
        # Filter to get the coefficient for miRNA (exclude intercept and other covariates)
        coef_info <- tidy_model %>% filter(term == miRNA)
        
        if (nrow(coef_info) == 1) {  # Ensure there is exactly one row of data for the miRNA term
          reg_results <- rbind(reg_results, data.frame(
            Trait = trait, 
            miRNA = miRNA, 
            Regression_Coefficient = coef_info$estimate,
            P_Value = coef_info$p.value,
            Conf_Lower = coef_info$conf.low,  # Ensure these columns exist
            Conf_Upper = coef_info$conf.high
          ))
        }
      }
    }
  }
  
  reg_results <- reg_results %>%
    group_by(Trait) %>%
    mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH")) %>%
    ungroup() %>%
    group_by(Trait) %>%
    arrange(Adjusted_P_Value) %>%
    ungroup() %>%
    mutate(Significance = case_when(
      Adjusted_P_Value < 0.001 ~ "***",
      Adjusted_P_Value < 0.01 ~ "**",
      Adjusted_P_Value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  return(reg_results)
}

# Calculate
echo_reg_results <- calculate_regressions(echo_traits, miRNAs, analysis_dat_ukhd2024_merged)
lab_reg_results <- calculate_regressions(lab_traits, miRNAs, analysis_dat_ukhd2024_merged)


## Volcano plot facet wrap ----------------------------------------------------

create_volcano_plot <- function(reg_results, title) {
  # Modify the data to include a color categorization for significance
  reg_results <- reg_results %>%
    mutate(Significant = ifelse(Adjusted_P_Value < 0.05, "Significant", "Not Significant"))
  
  plot <- ggplot(reg_results, aes(x = Regression_Coefficient, y = -log10(Adjusted_P_Value))) +
    geom_point(aes(color = Significant), alpha = 0.7) +
    scale_color_manual(values = c("Significant" = "#dc3912", "Not Significant" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    labs(x = "Regression Coefficient", y = "-log10 Adjusted P-Value", title = title) +
    theme_minimal(base_size = 16, base_family = 'Arial') +
    theme(legend.title = element_blank(), legend.position = "none") +
    facet_wrap(~ Trait, scales = "free") #+
    #geom_text_repel(aes(label = ifelse(Adjusted_P_Value < 0.05, miRNA, '')), size = 3)
  
  return(plot)
}

# Example for echo traits
echo_volcano_plot <- create_volcano_plot(echo_reg_results, "")
print(echo_volcano_plot)

# Example for lab traits
lab_volcano_plot <- create_volcano_plot(lab_reg_results, "")
print(lab_volcano_plot)


volcano_regression_echo_facet <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_volcano_facet_regression_echo.svg")
ggsave(filename = volcano_regression_echo_facet, plot = echo_volcano_plot, 
       width = 8, height = 8, 
       units = "in"  # default
)

volcano_regression_lab_facet <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_volcano_facet_regression_lab.svg")
ggsave(filename = volcano_regression_lab_facet, plot = lab_volcano_plot, 
       width = 8, height = 8, 
       units = "in"  # default
)



### DOTPLOT -------------------------------------------------------------

# Function to create a dot plot with min-max scaling on regression coefficients
create_dot_plot <- function(reg_results) {
  # Min-Max Scaling function
  min_max_scale <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  # Apply Min-Max Scaling to the absolute regression coefficients
  reg_results <- reg_results %>%
    group_by(Trait) %>% 
    dplyr::mutate(Scaled_Coefficient = min_max_scale(abs(Regression_Coefficient))) %>% 
    ungroup()
  
  # Create the dot plot
  ggplot(reg_results, aes(x = Trait, y = miRNA, size = Scaled_Coefficient, color = Significance)) +
    geom_point(shape=16, alpha = 0.5) +
    scale_size_continuous(range = c(1, 8)) +
    ggthemes::scale_color_gdocs() +
    # scale_color_manual(values = c("***" = "red", "**" = "blue", "*" = "green", "" = "grey")) +
    labs(
      x = "", # Trait", 
      y = "miRNA", 
      # title = "Dot Plot of miRNA Effects on Traits"
    ) +
    theme_minimal(base_size = 16, base_family = 'Arial') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank()
      )
}

# Example usage with your dataset
dot_plot_regression <- create_dot_plot(combined_reg_results)

dotplot_regression_traits <- glue("{data_path_bestageing2022}/output/analysis2024/correlation_features/001_dotplot_regression_features.svg")
ggsave(filename = dotplot_regression_traits, plot = dot_plot_regression, 
       width = 12, height = 10, 
       units = "in"  # default
)


###
# ML regression analysis -------------------------------------------------------------
###

# load helper fn
source("./scripts/helper/regression_tidy_function.R")

# Define the range of features
echo_traits <- names(analysis_dat_ukhd2024_merged)[6:14]  # Echo traits
lab_traits <- names(analysis_dat_ukhd2024_merged)[15:23]  # Lab traits
miRNAs <- names(analysis_dat_ukhd2024_merged)[24:518]  # miRNA features

modeldata = analysis_dat_ukhd2024_merged
modeldata <- modeldata |> 
  mutate(sex = factor(sex))

run_regression <- TRUE
save_files = TRUE

## start ##
if (run_regression) {
  # Ensure the echo_traits variable is available and correctly formatted
  if (is.null(echo_traits) || !is.character(echo_traits)) {
    stop("echo_traits must be a character vector of trait names.")
  }
  
  # Iterate over each trait in the echo_traits list
  for (trait in echo_traits) {
    # Attempt to run the regression analysis and handle potential errors
    results <- tryCatch({
      # Run the regression analysis for the current trait
      run_regression_analysis(
        modeldata = modeldata, 
        .outcome_var = trait, 
        .predictors = c(miRNAs, "age", "sex"), 
        .id_col = "pat_id", 
        .exclude_cols = c("age", "sex"), 
        split_proportion = 0.75, 
        n_feature_select = 25, 
        no_folds = 5, 
        no_repeats = 3, 
        grid_size = 100
      )
    }, error = function(e) {
      cat("Error processing trait", trait, ": ", e$message, "\n")
      return(NULL)  # Return NULL on error
    })
    
    # Only proceed with saving if race_results is not NULL
    if (!is.null(race_results)) {
      # Save the regression results to an RDS file (includes race_results, split, modeldata, recipe)
      filename <- glue("./output/analysis2024/regression/tuning_results/001_list_tune_race_results_{trait}.rds")
      saveRDS(object = results, file = filename)
    }
  }
}

## Load results and plot --------------------------
data_path <- "./output/analysis2024/regression/tuning_results/"
plots <- list()  # Store plots in a list

for (trait in echo_traits) {
  filename <- glue("{data_path}001_list_tune_race_results_{trait}.rds")
  
  # Check if the file exists before trying to read it
  if (file.exists(filename)) {
    # Load the race results for each trait
    results <- readRDS(filename)
    
    # Check if results or race_results is NULL
    if (!is.null(results) && !is.null(results$race_results)) {
      # Create plots using autoplot from the workflows package
      plot <- autoplot(
        results$race_results,
        rank_metric = "rsq",  # Ranking models by R-squared
        metric = c("rsq", "rmse"),  # Visualizing R-squared and RMSE
        select_best = TRUE  # Display only the best point per workflow
      ) +
        geom_text(aes(y = mean - 0.05, label = wflow_id), angle = 45, hjust = 1, size = 3) +
        labs(title=paste("Model Evaluation for", trait),
             subtitle= "Comparison of R-Squared and RMSE") +
        theme_minimal(base_size = 16, base_family = 'Arial') +
        theme(legend.position = "none")
      
      # Store the plot in a list
      plots[[trait]] <- plot
      
      # Optionally print the plot immediately
      print(plot)
    } else {
      message("No results available for trait: ", trait)
    }
  } else {
    message("File does not exist for trait: ", trait)
  }
}

## evaluate results on test set ------------------------------------

# Main loop, takes some time to bootstrap CI for performance measures
results_list <- list()  # Store results for each trait

for (trait in echo_traits) {
  filename <- glue("{data_path}001_list_tune_race_results_{trait}.rds")
  results <- readRDS(filename)
  if(is.null(results)) next
  
  dat_train <- results$dat_train
  dat_test <- results$dat_test  # Assuming dat_test is stored in the results
  
  # Initialize trait-specific results list
  results_list[[trait]] <- list()
  
  # Extract best model and its performance on test data
  for (model_idx in 1:nrow(results$race_results)) {
    wflow_id <- results$race_results[[1]][model_idx]
    cat("Processing wflow_id:", wflow_id, "for trait:", trait, "\n")
    # fns in "regression_tidy_function.R"
    results_list[[trait]][[wflow_id]] <- process_model(results, trait, model_idx)
  }
}

# results_list

saveRDS(results_list, file = glue("{data_path_bestageing2022}/output/analysis2024/regression/vip_results/vip_results_list.rds"))


### look at results -----------

r_squared_list <- list()
for(trait in names(results_list)) {
  # For each trait, extract the R_squared values from each model
  r_squared_values <- sapply(results_list[[trait]], function(model) model$R_squared)
  r_squared_values_traintest <- sapply(results_list[[trait]], function(model) model$train_test_metrics_grouped$.estimate[7])
  # Store the extracted R_squared values in the list, named by trait
  r_squared_list[[trait]] <- r_squared_values
}
r_squared_list

performance_measures_combined <- lapply(names(results_list), function(trait) {
  lapply(names(results_list[[trait]]), function(model) {
    if ("train_test_metrics_grouped" %in% names(results_list[[trait]][[model]])) {
      df <- results_list[[trait]][[model]]$train_test_metrics_grouped
      df$trait <- trait
      df$model <- model
      return(df)
    } else {
      return(NULL)  # Return NULL if the dataframe does not exist
    }
  })
}) %>% bind_rows()

performance_measures_combined |> 
  filter(.metric=="rsq", model=="normalized_LM") |> 
  print(n=100)
  
## VIP ---------

vip_results_list <- list()  # Store results for each trait
# cave takes hours!!
for (trait in echo_traits) {
  filename <- glue("{data_path}001_list_tune_race_results_{trait}.rds")
  results <- readRDS(filename)
  if(is.null(results)) next
  
  # trait specific lists for all models
  vip_results_list[[trait]] <- list()
  for (model_idx in 1:nrow(results$race_results)) {
    wflow_id <- results$race_results[[1]][model_idx]
    cat("VIP: Processing wflow_id:", wflow_id, "for trait:", trait, "\n")
    # fns in "regression_tidy_function.R"
    vip_results_list[[trait]][[wflow_id]] <- calculate_feature_importance(results, trait, model_idx, n_bootstraps=100)
  }
}

saveRDS(vip_results_list, file = glue("{data_path_bestageing2022}/output/analysis2024/regression/vip_results/vip_results_list.rds"))






