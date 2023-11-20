
### INFO ----------------------------------------------------------------------
# this script is sourced from `scripts/render_param_reports.R`
# selection provided by `all_combis$diseases` and `all_combis$analysis`

# update in 003c_full_analysis script: received detection matrix from hummingbird dx, filtered based on this matrix
# full analysis with feature selection and model tuning

# scripts saves tuned models (not finalized) to: "./output/tuning_results/"

# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
SAVE.files <- TRUE

# Define library and data paths based on system
if (system_name == "MacBook-Pro-CR-2065.local" | grepl("laptop-zim.uni-heidelberg.de", system_name)) {
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
require(readxl, lib.loc = lib_path)
require(glue, lib.loc = lib_path)
require(janitor, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(tidyr, lib.loc = lib_path)
require(tibble, lib.loc = lib_path)
require(stringr, lib.loc = lib_path)
require(purrr, lib.loc = lib_path)
require(ggplot2, lib.loc = lib_path)
require(svglite, lib.loc = lib_path)
require(ggrepel, lib.loc = lib_path)
require(ggthemes, lib.loc = lib_path)
require(skimr, lib.loc = lib_path)
require(tableone, lib.loc = lib_path)
require(pROC, lib.loc = lib_path)
require(dials, lib.loc = lib_path)
require(infer, lib.loc = lib_path)
require(modeldata, lib.loc = lib_path)
require(tidymodels, lib.loc = lib_path)
require(rsample, lib.loc = lib_path)
options(tidymodels.dark = TRUE)
## used within tidymodels
# "kknn", "glmnet", "ranger", "naivebayes", "kernlab", "xgboost", "nnet"
require(kknn, lib.loc = lib_path)
require(glmnet, lib.loc = lib_path)
require(ranger, lib.loc = lib_path)
require(naivebayes, lib.loc = lib_path)
require(kernlab, lib.loc = lib_path)
require(xgboost, lib.loc = lib_path)
require(nnet, lib.loc = lib_path)
require(colino, lib.loc = lib_path)
#library(keras, lib.loc = lib_path)  # ## ERRORS with keras
# Fold2: preprocessor 1/1, model 100/100:
# Error: Python shared library not found, Python bindings not loaded.
# Use reticulate::install_miniconda() if you'd like to install a Miniconda Python environment.
#library(reticulate, lib.loc = lib_path)
#reticulate::install_miniconda(path = "/mnt/users/reich/programs/miniconda/", update = TRUE, force = FALSE)
#keras::install_keras(conda = "/mnt/users/reich/miniconda/bin/conda", version = "default")
# reticulate::conda_list()
# reticulate::use_python("/mnt/users/reich/programs/miniconda/bin/python3", required = TRUE)
# Sys.setenv(RETICULATE_PYTHON = "/mnt/users/reich/programs/miniconda/bin/python3")
require(discrim, lib.loc = lib_path)  # naive bayes
require(finetune, lib.loc = lib_path)
require(vetiver, lib.loc = lib_path)
require(workflowsets, lib.loc = lib_path)
require(baguette, lib.loc = lib_path)
require(rules, lib.loc = lib_path)
tidymodels_prefer()
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::slice)
conflicted::conflicts_prefer(dplyr::filter)


source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))

cores <- parallel::detectCores()
if (!grepl("mingw32", R.Version()$platform)) {
  library(doMC, lib.loc = lib_path)
  registerDoMC(cores = cores/2)
} else {
  library(doParallel, lib.loc = lib_path)
  cl <- makePSOCKcluster(cores/2)
  registerDoParallel(cl)
}

# set analysis variables - alternatively sys args -------------------------------
args <- commandArgs(trailingOnly = TRUE)
miRetrieveBiomarker <- as.logical(args[1])  # if TRUE "selected" analysis will use top 50 biomarkers from miRetrieve research
random_selection <- as.logical(args[2])  # random selection only runs as a "selected" analysis!  # default FALSE
no_folds <- 5
no_repeats <- as.integer(args[3])  # default 10

if(length(args) == 0) {
  no_repeats <- 10
  miRetrieveBiomarker <- TRUE
  random_selection <- FALSE
}


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="full")

# if ( !exists("all_combis")  ) {  # check if variable name exists
#   diseases <- c("dcm", "acs", "cad", "hfref")
#   analysis <- c("selected", "full")
#   all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
#     filter(analysis=="selected")
# }


###
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

## seq metadat 09-2023
# batch info
hbdx_metadat <- clean_names(readxl::read_excel(path = glue('{data_path_bestageing2022}/data/kahraman2023/230907_annotation_chrstoph_reich.xlsx')))  # has also multiclass col + diagnoses
#RIN
rin_mean <- mean(as.numeric(hbdx_metadat$rin), na.rm=TRUE)
rin_sd <- sd(as.numeric(hbdx_metadat$rin), na.rm=TRUE)

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
# START LOOP for specified combinations -------------------------------------
###
for (i in 1:nrow(all_combis)) {
  ## DATA FIRST
  cat("\n\n")
  print(glue("|||------------------------------------------------------------------------------------------------------------|||"))
  print(glue("|||-----------------------Start modelling for disease: ", stringr::str_to_upper(all_combis$diseases[i]), ", and selection of miRNAs: ", stringr::str_to_upper(all_combis$analysis[i])," -----------------------|||") )
  print(glue("Modeling parameters: {no_folds}-fold-cross validation with {no_repeats} repeats"))
  
  if(random_selection==TRUE){
    print(glue("|||-----------------------!!!CAVE!!! used random miRNA selection for analysis-----------------------|||"))
  }
  
  if(miRetrieveBiomarker == TRUE & all_combis$analysis[i] == "selected"){
    print(glue("|||-----------------------Used literature miRNAs from miRetrieve Text Mining-----------------------|||"))
  }
  if(miRetrieveBiomarker == FALSE & all_combis$analysis[i] == "selected"){
    print(glue("|||-----------------------Used literature miRNAs from Reviews--------------------------------------|||"))
  }
  print(glue("|||------------------------------------------------------------------------------------------------------------|||"))
  cat("\n")
  # create parameter specific {disease}_ids
  
  filename <- glue("{data_path_BestAgeing}/data/pheno_{all_combis$diseases[i]}.xlsx")
  disease_vector <- paste0(all_combis$diseases[i], "_ids")
  assign(x = disease_vector, 
         value = read_excel(filename) %>% 
           dplyr::pull(BestAgeingCode))
  
  
  disease_ident_df <- tibble(pat_id = c(eval(as.symbol(disease_vector)),control_ids), 
                             disease = c(rep(all_combis$diseases[i], length(eval(as.symbol(disease_vector)))),
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
  
  
  # load disease specific data from 001c ----------------------------------------------------------
  # specific to this script (preprocessing, filtered miRNAs with Detection Matrix, batch effects removed)
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_{all_combis$diseases[i]}_data01.rds")
  if(!file.exists(path2dataprocessed)) {
    next
  }
  
  data01 <- readRDS(file = path2dataprocessed)
  
  
  ## CHOOSE miRetrieve Biomarkers? -----------
  if (miRetrieveBiomarker==TRUE) {
    path2biomarker <- glue("{data_path_bestageing2022}/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-27-human-disease_biomarker_with_accession.rds")
    path2biomarker_count <- glue("{data_path_bestageing2022}/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-27-human-df_count_both_with_accession.rds")
    human_disease_biomarker <- readRDS(file = path2biomarker)
    human_disease_biomarker_count <- readRDS(file = path2biomarker_count)
    
    # 1) most hits on miRetrieve
    human_disease_biomarker_count_tidy <- 
      human_disease_biomarker_count %>%
      drop_na(Accession) %>%
      arrange(desc(miRetrieve)) %>%
      dplyr::slice(1:30)
    
    # 2) highest BM score <-- used this
    human_disease_biomarker <- human_disease_biomarker %>% 
      drop_na(Accession) %>% 
      arrange(desc(Biomarker_score)) %>% 
      left_join(human_disease_biomarker_count %>% 
                  select(Accession, miRetrieve, PubMed), 
                by=c("Accession"="Accession")) %>% 
      # drop_na(miRetrieve) %>%   # some PMIDs not in most hits ;)
      mutate(miRetrieve = ifelse(is.na(miRetrieve), 0, miRetrieve) ) %>% 
      ## !! GROUP BY distinct miRNAs!!
      group_by(Accession) %>% 
      mutate(max_value = (Biomarker_score + miRetrieve)) %>% 
      arrange(desc(max_value)) %>% 
      slice(1) %>% 
      ungroup() %>% 
      arrange(desc(max_value)) %>% 
      drop_na(TargetName) 
    
    human_disease_biomarker$TargetName <-  make_clean_names(human_disease_biomarker$TargetName) %>% 
      str_replace(pattern = "mi_r", replacement = "mir")
    # sanity check if miRNA names are in sequenced df
    human_disease_biomarker <- human_disease_biomarker %>% 
      filter(TargetName %in% colnames(data01) ) %>% 
      slice(1:50)  # top 50 
    
    # save for documentation
    #saveRDS(object = human_disease_biomarker, 
    #        file = glue("/mnt/users/reich/rockerprojects/bestageing2022/data-literature/miRetrieve/{all_combis$diseases[i]}/{Sys.Date()}-human-disease_biomarker_with_accession_max_value.rds"))
  }
  
  # FULL ANALYSIS HERE
  
  # data01 <- data01 %>% 
  #   select(pat_id, disease, age, sex, human_disease_biomarker$TargetName)
  
  # no.mirnas <- length(human_disease_biomarker$TargetName)
  no.mirnas <- ncol(data01) -4
  text_disease <- stringr::str_to_upper(all_combis$diseases[i])
  
  ###
  # MODEL -------------------------------------------------------------------
  ###
  
  ### Initial Split
  set.seed(123)
  
  modeldat <- data01
  
  # make control first factor for all analyses ;)
  modeldat <- modeldat %>% 
    mutate(disease = factor(disease, levels = c("control", all_combis$diseases[i])))
  
  dat_split <- rsample::initial_split(modeldat, strata = disease)
  dat_train <- training(dat_split)
  dat_test <- testing(dat_split)
  
  folds <- vfold_cv(dat_train, strata = disease, v = no_folds, repeats = no_repeats)
  
  ###
  # recipe -------------------------------------------------
  # A
  
  n_feature_select <- 20  # how many mirna select in training?
  normalized_rec <- 
    recipe(disease ~ ., data = dat_train) %>%
    ### https://recipes.tidymodels.org/articles/Ordering.html
    #  To make sure we don’t get any unexpected results, it’s best to use 
    #  the following ordering of high-level transformations:
    #     Skewness Transformations - step_YeoJohnson()
    #     Centering, Scaling, or Normalization on Numeric Predictors
    #     Dummy Variables for Categorical Data
    update_role(pat_id, new_role="ID") %>% 
    step_zv(all_predictors()) %>%  
    step_impute_mean(all_numeric_predictors()) %>% 
    step_impute_mode(all_nominal_predictors(), -disease) %>% 
    step_corr(all_numeric_predictors(), threshold = 0.9) %>% 
    step_YeoJohnson() %>% 
    step_normalize(all_numeric_predictors()) %>%  
    step_dummy(all_nominal_predictors(),-disease) %>% 
    step_select_forests(all_predictors(), -c(age, sex_Male), outcome = "disease", top_p = n_feature_select)
  # B
  poly_rec <- 
    recipe(disease ~ ., data = dat_train) %>%
    update_role(pat_id, new_role="ID") %>% 
    step_zv(all_predictors()) %>%  
    step_impute_mean(all_numeric_predictors()) %>% 
    step_impute_mode(all_nominal_predictors(), -disease) %>% 
    step_corr(all_numeric_predictors(), threshold = 0.9) %>% 
    step_normalize(all_numeric_predictors()) %>%  
    step_poly(all_numeric_predictors()) %>% 
    step_dummy(all_nominal_predictors(),-disease) %>% 
    step_select_forests(all_predictors(), -c(starts_with("age"), starts_with("sex")), outcome = "disease", top_p = n_feature_select)
  #step_interact( ~all_predictors():all_predictors())
  
  # C
  simple_rec <- 
    recipe(disease ~ ., data = dat_train) %>%
    update_role(pat_id, new_role="ID") %>% 
    # ZERO VARIANCE
    step_zv(all_predictors()) %>%  
    # IMPUTE
    step_impute_mean(all_numeric_predictors()) %>% 
    step_impute_mode(all_nominal_predictors(), -disease) %>% 
    # DECORRELATE
    step_corr(all_numeric_predictors(), threshold = 0.9) %>% 
    step_dummy(all_nominal_predictors(),-disease) %>% 
    step_select_forests(all_predictors(), -c(age, sex_Male), outcome = "disease", top_p = n_feature_select)
  
  # SANITY CHECK! takes a while with a lot of features in full analysis
  # prepped_rec <- prep(normalized_rec, dat_train, strings_as_factors = FALSE)  # https://community.rstudio.com/t/how-to-specify-a-column-to-be-unaffected-in-recipes/23056/6
  # test_baked_train <- bake(prepped_rec, new_data = dat_train)
  
  # prepped_rec <- prep(poly_rec, dat_train, strings_as_factors = FALSE)  # https://community.rstudio.com/t/how-to-specify-a-column-to-be-unaffected-in-recipes/23056/6
  # test_baked_train <- bake(prepped_rec, new_data = dat_train)
  
  
  ###
  # specs-parsnip ----------------------------------------------------------
  #parsnip::set_dependency("kknn", "glmnet", "ranger", "naivebayes") #, "kernlab", "xgboost", "nnet")
  
  nearest_neighbor_kknn_spec <-
    nearest_neighbor(neighbors = tune(), weight_func = tune(), 
                     dist_power = tune()) %>%
    set_engine('kknn') %>%
    set_mode('classification')
  
  logistic_reg_glmnet_spec <-
    logistic_reg(penalty = tune(), mixture = tune()) %>%
    set_mode("classification") %>% 
    set_engine('glmnet')
  
  rand_forest_ranger_spec <-
    rand_forest(mtry = tune(), min_n = tune()) %>%
    set_engine('ranger') %>%
    set_mode('classification')
  
  naive_Bayes_naivebayes_spec <-
    naive_Bayes(smoothness = tune(), Laplace = tune()) %>%
    set_engine('naivebayes') %>% 
    set_mode('classification')
  
  svm_linear_kernlab_spec <-
    svm_linear(cost = tune(), margin = tune()) %>%
    set_engine('kernlab') %>%
    set_mode('classification')
  
  svm_poly_kernlab_spec <-
    svm_poly(cost = tune(), degree = tune(), 
             scale_factor = tune(), margin = tune()) %>%
    set_engine('kernlab') %>%
    set_mode('classification')
  
  svm_rbf_kernlab_spec <-
    svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>%
    set_engine('kernlab') %>%
    set_mode('classification')
  
  boost_tree_xgboost_spec <-  # https://parsnip.tidymodels.org/reference/details_boost_tree_xgboost.html
    boost_tree(trees = tune(), tree_depth = tune(), 
               min_n = tune(),                              ## first three: model complexity
               loss_reduction = tune(), 
               stop_iter = tune(), 
               sample_size = tune(), mtry = tune(),         ## randomness
               learn_rate = tune(),                         ## step size
    ) %>%
    set_engine('xgboost') %>%
    set_mode('classification')
  
  mlp_nnet_spec <-
    mlp(hidden_units = tune(), penalty = tune(), 
        epochs = tune()) %>%
    # x Fold4: preprocessor 1/1, model 21/100: Error in nnet.default(x, y, w, entropy = TRUE, ...): too many (18937) we...
    set_engine('nnet', MaxNWts = 20000) %>%  #`keras`
    set_mode('classification')
  
  ##
  # workflow set --------------------------------------------------------------
  ## A
  normalized <- 
    workflow_set(   # IDs are generated {name of preproc} + {name of model}
      preproc = list(normalized = normalized_rec),
      models = list(
        KNN = nearest_neighbor_kknn_spec, 
        SVM_radial = svm_rbf_kernlab_spec, 
        SVM_poly = svm_poly_kernlab_spec, 
        SVM_linear = svm_linear_kernlab_spec,
        neural_network = mlp_nnet_spec
      )
    )
  
  ## B
  simple <- 
    workflow_set(
      preproc = list(simple = simple_rec),
      models = list(
        #naive_bayes = naive_Bayes_naivebayes_spec,    ## ERROR when tuning "Error in pkgs$pkg[[1]] : subscript out of bounds"
        RF = rand_forest_ranger_spec,
        XGB = boost_tree_xgboost_spec
      )
    )
  
  ## C
  with_features <- 
    workflow_set(
      preproc = list(full_quad = poly_rec),
      models = list(
        logistic_reg = logistic_reg_glmnet_spec,
        KNN = nearest_neighbor_kknn_spec
      )
    )
  
  # normalized %>% extract_workflow(id = "normalized_KNN")
  
  ## these objects are tibbles -> row binding does not affect the state of the sets
  ### result is itself a workflow set
  all_workflows <- 
    bind_rows(simple, normalized, with_features) %>% 
    # make workflow IDs a little more simple:
    mutate(wflow_id = gsub("(simple_)|(normalized_)", "", wflow_id))
  
  ###
  # define-grids --------------------------------------------------------------
  ###
  # need to finalize mtry - data dependent (unknown values - how many predictors?)
  ## data dependent: mtry(), sample_size(), num_terms(), num_comp()
  set.seed(123)
  
  grid_size <- 10 #500
  # CAVE RF and XGB max mtry must be in range of predictors!
  grid_RF <- rand_forest_ranger_spec %>%   # 2 hyperparams
    extract_parameter_set_dials() %>% 
    # data dependent
    update(mtry = mtry(range = c(1, n_feature_select+2)) ) %>%  # mirnas + age + sex
    grid_latin_hypercube(size=grid_size) 
  
  grid_XGB <- boost_tree_xgboost_spec %>%   # 2 hyperparams
    extract_parameter_set_dials() %>% 
    # data dependent
    update(mtry = mtry(range = c(1, n_feature_select+2)) ) %>%  # mirnas + age + sex
    grid_latin_hypercube(size=grid_size) 
  
  
  grid_KNN <- nearest_neighbor_kknn_spec %>%  # 3 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=grid_size)
  
  grid_SVM_radial <- svm_rbf_kernlab_spec %>%   # 3 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=grid_size)
  
  grid_SVM_poly <- svm_poly_kernlab_spec %>%   # 4 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=grid_size)
  
  grid_SVM_linear <- svm_linear_kernlab_spec %>%   # 2 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=grid_size)
  
  grid_neural_network <- mlp_nnet_spec %>%   # 3 hyperparams
    extract_parameter_set_dials() %>% 
    update(epochs = epochs() %>% range_set(c(10, 150))) %>%   # epochs()  Range: [10, 1000] (default)
    grid_latin_hypercube(size=grid_size)
  
  grid_full_quad_logistic_reg <- logistic_reg_glmnet_spec %>%  # 2 hyperparams
    extract_parameter_set_dials() %>% 
    grid_latin_hypercube(size=grid_size)
  
  ## supply grid to workflow options
  # https://github.com/tidymodels/workflowsets/issues/37
  all_workflows <- all_workflows %>% 
    option_add(grid = grid_RF, id = "RF") %>% 
    option_add(grid = grid_XGB, id = "XGB") %>% 
    option_add(grid = grid_KNN, id = "KNN") %>% 
    option_add(grid = grid_SVM_radial, id = "SVM_radial") %>% 
    option_add(grid = grid_SVM_poly, id = "SVM_poly") %>% 
    option_add(grid = grid_SVM_linear, id = "SVM_linear") %>% 
    option_add(grid = grid_neural_network, id = "neural_network") %>% 
    option_add(grid = grid_full_quad_logistic_reg, id = "full_quad_logistic_reg") %>% 
    option_add(grid = grid_KNN, id = "full_quad_KNN")   # same hyperparams
  
  #all_workflows <- all_workflows %>% slice(8:9)
  
  # # debug
  # all_workflows <- all_workflows %>% 
  #   filter(wflow_id=="neural_network")
  
  ###
  # tune-race-anova -----------------------------------------------------------
  ###
  
  set.seed(123)
  race_ctrl <-
    control_race(
      verbose = TRUE,
      verbose_elim = FALSE,
      allow_par = TRUE,
      save_pred = TRUE,
      burn_in = 3,
      num_ties = 10,
      alpha = 0.05,
      randomize = TRUE,
      pkgs = NULL,
      event_level = "second",  # changed 12-07-2023
      parallel_over = "everything",
      save_workflow = TRUE
    )
  
  time1 <- Sys.time()
  race_results <-  all_workflows %>%
    workflow_map(
      "tune_race_anova",
      # options to `tune_race_anova()`
      resamples = folds,
      # grid = 25,  # grid specified to workflow_map above
      metrics = metric_set(roc_auc),
      control = race_ctrl,
      # options to `workflow_map()`
      seed = 20221111,
      verbose = TRUE
    )
  time2 <- Sys.time()
  time.diff.race <- time2-time1
  
  print(glue("|||----------------------- It took {round(time.diff.race,4)} mins to tune all grids with the race approch for {stringr::str_to_upper(all_combis$diseases[i])} and selection of miRNAs {stringr::str_to_upper(all_combis$analysis[i])} -----------------------|||"))
  
  ## SAVE RACE RESULTS -----------------------------------------
  if(SAVE.files == TRUE) {
    filename_tune_race_results <- glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_full_analysis_tune_race_results_repeats_{no_repeats}_folds_{no_folds}_{text_disease}_analysis_randomMIR_{random_selection}.rds")
    saveRDS(object = race_results, file = filename_tune_race_results)
  }
  
  
  num_race_models <- sum(collect_metrics(race_results)$n)
  
  # PLOT results ----------------------------------------------
  race_results %>% 
    rank_results() %>% 
    filter(.metric == "roc_auc") %>% 
    select(model, .config, roc_auc=mean, rank) -> rankings_race
  
  top2models_race <- rankings_race$model[1:2]
  for (j in 1:length(all_workflows$info)) {
    if (all_workflows$info[[j]]$model == top2models_race[1]) {
      topmodel1_race <- all_workflows$wflow_id[i]
    }
    
    if (all_workflows$info[[j]]$model == top2models_race[2]) {
      topmodel2_race <- all_workflows$wflow_id[j]
    }
  }
  
  autoplot(
    race_results,
    rank_metric = "roc_auc",  # <- how to order models
    metric = "roc_auc",       # <- which metric to visualize
    select_best = TRUE     # <- one point per workflow
  ) +
    geom_text(aes(y = mean - 0.05, label = wflow_id), angle = 45, hjust = 1, size =3) +
    lims(y = c(0.6, 1.0)) +
    ylab("ROC-AUC")+
    labs(title=paste0("Ranking models: ", text_disease, ", miRNA_set: ", 
                      str_to_upper(all_combis$analysis[i])),
         subtitle= paste0("Racing approach | n_grids evaluated: ", num_race_models))+
    scale_x_continuous(breaks = seq(1, nrow(race_results), by = 1)) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme() +
    theme(legend.position = "none") -> plot_tune_race_ranking
  
  filename_plot_tune_race_ranking <- glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_{Sys.Date()}_tune_race_ranking_repeats_{no_repeats}_folds_{no_folds}_{text_disease}_analysis_{str_to_upper(all_combis$analysis[i])}_miRetrieve_{miRetrieveBiomarker}_randomMIR_{random_selection}.svg")
  ggsave(filename = filename_plot_tune_race_ranking, plot = plot_tune_race_ranking, 
         width = 14, height = 10, 
         units = "in"  # default
  )
  
  
  print(paste0("|||-----------------------Run finished for disease: ", stringr::str_to_upper(all_combis$diseases[i]), 
               ", and selection of miRNAs: ", stringr::str_to_upper(all_combis$analysis[i]),
               " -----------------------|||") )
} # END LOOP


