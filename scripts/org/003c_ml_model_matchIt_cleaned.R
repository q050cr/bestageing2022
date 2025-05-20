### INFO ----------------------------------------------------------------------
# Machine Learning Model with MatchIt for Balanced Samples
# Author: Christoph Reich
# Date: 2024-04-20
#
# This script builds machine learning models with sample matching
# Saves tuned models to: "./output/tuning_results/"

# Load configuration
source("scripts/config/config.R")

# Get project paths
paths <- get_project_paths()
data_path_bestageing2022 <- paths$project_path
data_path_BestAgeing <- paths$data_path_BestAgeing
lib_path <- paths$lib_path

# Global settings
SAVE.files <- TRUE

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)
miRetrieveBiomarker <- as.logical(args[1]) # if TRUE "selected" analysis will use top 50 biomarkers from miRetrieve research
random_selection <- as.logical(args[2]) # random selection only runs as a "selected" analysis
no_folds <- 5
no_repeats <- as.integer(args[3]) # default 10

# Default values if no args provided
if (length(args) == 0) {
  no_repeats <- 10
  miRetrieveBiomarker <- TRUE
  random_selection <- FALSE
}

# Load required libraries
required_packages <- c(
  # Data manipulation
  "readxl", "glue", "janitor", "dplyr", "tidyr", "tibble", "stringr", "purrr",
  # Visualization
  "ggplot2", "svglite", "ggrepel", "ggthemes",
  # Analysis
  "MatchIt", "skimr", "tableone", "pROC", "dials", "infer", "modeldata",
  "tidymodels", "rsample",
  # Model engines
  "kknn", "glmnet", "ranger", "naivebayes", "kernlab", "xgboost",
  "nnet", "colino", "discrim", "finetune", "vetiver", "workflowsets",
  "baguette", "rules", "doMC"
)

# Load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = lib_path)
    library(pkg, character.only = TRUE)
  }
}

# Configure tidymodels
tidymodels_prefer()
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::slice)
conflicted::conflicts_prefer(dplyr::filter)

# Load helper functions
source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))

# Set up parallel processing
cores <- parallel::detectCores()
if (!grepl("mingw32", R.Version()$platform)) {
  library(doMC, lib.loc = lib_path)
  registerDoMC(cores = cores / 2)
} else {
  library(doParallel, lib.loc = lib_path)
  cl <- makePSOCKcluster(cores / 2)
  registerDoParallel(cl)
}

# Define diseases and analysis types
diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%
  filter(analysis == "selected")

#-----------------------------------------------------------------------------
# Load data
#-----------------------------------------------------------------------------

# Load miRNA data
model_data1 <- clean_names(readRDS(file = glue("{data_path_BestAgeing}/data_new/model_data1.RDS")))
load(file = glue("{data_path_BestAgeing}/data/mirnas.rda")) # "UKL-HD" n=765
load(file = glue("{data_path_BestAgeing}/data/data.rda")) # "UKL-HD" n=731

# Clean miRNA data
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

# Load sequencing metadata
hbdx_metadat <- clean_names(readxl::read_excel(path = glue("{data_path_bestageing2022}/data/kahraman2023/230907_annotation_chrstoph_reich.xlsx")))
rin_mean <- mean(as.numeric(hbdx_metadat$rin), na.rm = TRUE)
rin_sd <- sd(as.numeric(hbdx_metadat$rin), na.rm = TRUE)

# Load detection matrix
det_mat_all_mirnas <- read.table(glue("{data_path_bestageing2022}/data/kahraman2023/det_mat_all_mirnas.txt")) %>%
  t() %>%
  as.data.frame() %>%
  clean_names()
colnames(det_mat_all_mirnas) <- str_replace(string = colnames(det_mat_all_mirnas), pattern = "mi_r", replacement = "mir")
rownames(det_mat_all_mirnas) <- gsub("\\.", "-", rownames(det_mat_all_mirnas))
det_mat_all_mirnas <- det_mat_all_mirnas %>%
  rownames_to_column(var = "pat_id") %>%
  as_tibble()

# Process metadata
unique_slide_ids <- length(unique(hbdx_metadat$slide_id))
arrays_per_slide <- nrow(hbdx_metadat) / unique_slide_ids
duplicate_ids <- hbdx_metadat[duplicated(hbdx_metadat$customer_id), ] %>% pull(customer_id)
hbdx_metadat <- hbdx_metadat %>% distinct(customer_id, .keep_all = TRUE)

# Load research miRNAs
load(glue("{data_path_BestAgeing}/data_research/fromR/researchMiRNAAccession.rda"))
researchMiRNAAccession$miRNAName_v21 <- make_clean_names(researchMiRNAAccession$miRNAName_v21) %>%
  str_replace(pattern = "mi_r", replacement = "mir")
# Fix inconsistent naming
researchMiRNAAccession$miRNAName_v21[researchMiRNAAccession$miRNAName_v21 == "hsa_mir_106a_5p"] <- "hsa_mir_106b_5p"

# Load diagnoses and survival data
load(glue("{data_path_BestAgeing}/data/diagnoses_df.rda"))
survival_dat <- clean_names(readRDS(glue("{data_path_bestageing2022}/data/202211908_XMELD_abfrage_best_ageing.rds")))

# Load metadata
load(glue("{data_path_BestAgeing}/data/clean_all_meta.rda"))
clean_all_meta <- clean_all_meta %>% mutate(age = ifelse(age < 18, NA, age))

# Load control IDs
control_ids <- read_excel(glue("{data_path_BestAgeing}/data/pheno_controls.xlsx")) %>%
  dplyr::pull(BestAgeingCode)
# "UKL-HD-00318" both in Control and CAD dataset
control_ids <- control_ids[control_ids != "UKL-HD-00318"]

#-----------------------------------------------------------------------------
# Main analysis loop
#-----------------------------------------------------------------------------

for (i in 1:nrow(all_combis)) {
  # Display info header
  cat("\n\n")
  print(glue("|||------------------------------------------------------------------------------------------------------------|||"))
  print(glue(
    "|||-----------------------Start modelling for disease: ", stringr::str_to_upper(all_combis$diseases[i]),
    ", and selection of miRNAs: ", stringr::str_to_upper(all_combis$analysis[i]), " -----------------------|||"
  ))
  print(glue("Modeling parameters: {no_folds}-fold-cross validation with {no_repeats} repeats"))

  if (random_selection == TRUE) {
    print(glue("|||-----------------------!!!CAVE!!! used random miRNA selection for analysis-----------------------|||"))
  }

  if (miRetrieveBiomarker == TRUE & all_combis$analysis[i] == "selected") {
    print(glue("|||-----------------------Used literature miRNAs from miRetrieve Text Mining-----------------------|||"))
  }
  if (miRetrieveBiomarker == FALSE & all_combis$analysis[i] == "selected") {
    print(glue("|||-----------------------Used literature miRNAs from Reviews--------------------------------------|||"))
  }
  print(glue("|||------------------------------------------------------------------------------------------------------------|||"))
  cat("\n")

  # Load disease-specific data
  filename <- glue("{data_path_BestAgeing}/data/pheno_{all_combis$diseases[i]}.xlsx")
  disease_vector <- paste0(all_combis$diseases[i], "_ids")
  assign(
    x = disease_vector,
    value = read_excel(filename) %>% dplyr::pull(BestAgeingCode)
  )

  # Create disease identification dataframe
  disease_ident_df <- tibble(
    pat_id = c(eval(as.symbol(disease_vector)), control_ids),
    disease = c(
      rep(all_combis$diseases[i], length(eval(as.symbol(disease_vector)))),
      rep("control", length(control_ids))
    )
  )

  # Load preprocessed data
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_{all_combis$diseases[i]}_data01.rds")
  if (!file.exists(path2dataprocessed)) {
    print(glue("Skipping {all_combis$diseases[i]} - processed data not found."))
    next
  }

  data01 <- readRDS(file = path2dataprocessed)

  # Select miRNAs for analysis
  if (miRetrieveBiomarker == TRUE) {
    # Load miRetrieve biomarkers
    path2biomarker <- glue("{data_path_bestageing2022}/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-27-human-disease_biomarker_with_accession.rds")
    path2biomarker_count <- glue("{data_path_bestageing2022}/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-27-human-df_count_both_with_accession.rds")
    human_disease_biomarker <- readRDS(file = path2biomarker)
    human_disease_biomarker_count <- readRDS(file = path2biomarker_count)

    # Process biomarkers
    human_disease_biomarker <- human_disease_biomarker %>%
      drop_na(Accession) %>%
      arrange(desc(Biomarker_score)) %>%
      left_join(
        human_disease_biomarker_count %>% select(Accession, miRetrieve, PubMed),
        by = c("Accession" = "Accession")
      ) %>%
      mutate(miRetrieve = ifelse(is.na(miRetrieve), 0, miRetrieve)) %>%
      group_by(Accession) %>%
      mutate(max_value = (Biomarker_score + miRetrieve)) %>%
      arrange(desc(max_value)) %>%
      slice(1) %>%
      ungroup() %>%
      arrange(desc(max_value)) %>%
      drop_na(TargetName)

    # Clean miRNA names
    human_disease_biomarker$TargetName <- make_clean_names(human_disease_biomarker$TargetName) %>%
      str_replace(pattern = "mi_r", replacement = "mir")

    # Filter and select top 50
    human_disease_biomarker <- human_disease_biomarker %>%
      filter(TargetName %in% colnames(data01)) %>%
      slice(1:50)
  }

  # Select miRNAs for analysis
  data01 <- data01 %>%
    select(pat_id, disease, age, sex, human_disease_biomarker$TargetName)

  no.mirnas <- length(human_disease_biomarker$TargetName)
  text_disease <- stringr::str_to_upper(all_combis$diseases[i])

  #-----------------------------------------------------------------------------
  # Machine Learning Model Training
  #-----------------------------------------------------------------------------

  # Data preparation
  set.seed(123)
  modeldat <- data01 %>% na.omit()
  modeldat$disease <- factor(modeldat$disease, levels = c("control", all_combis$diseases[i]))

  # Perform matching
  m.out <- matchit(disease ~ age + sex, data = modeldat, method = "nearest")
  summary(m.out)

  # Get matched data
  modeldat <- match.data(m.out) %>% select(-c("distance", "weights", "subclass"))

  # Split data
  dat_split <- rsample::initial_split(modeldat, strata = disease)
  dat_train <- training(dat_split)
  dat_test <- testing(dat_split)

  # Create cross-validation folds
  folds <- vfold_cv(dat_train, strata = disease, v = no_folds, repeats = no_repeats)

  #-----------------------------------------------------------------------------
  # Recipe definitions
  #-----------------------------------------------------------------------------

  # Calculate feature selection parameters
  n_feature_select <- ncol(dat_train) - 10

  # Normalized recipe
  normalized_rec <-
    recipe(disease ~ ., data = dat_train) %>%
    update_role(pat_id, new_role = "ID") %>%
    step_zv(all_predictors()) %>%
    step_impute_mean(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors(), -disease) %>%
    step_corr(all_numeric_predictors(), threshold = 0.9) %>%
    step_YeoJohnson() %>%
    step_normalize(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors(), -disease)

  # Polynomial recipe
  poly_rec <-
    recipe(disease ~ ., data = dat_train) %>%
    update_role(pat_id, new_role = "ID") %>%
    step_zv(all_predictors()) %>%
    step_impute_mean(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors(), -disease) %>%
    step_corr(all_numeric_predictors(), threshold = 0.9) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_poly(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors(), -disease)

  # Simple recipe
  simple_rec <-
    recipe(disease ~ ., data = dat_train) %>%
    update_role(pat_id, new_role = "ID") %>%
    step_zv(all_predictors()) %>%
    step_impute_mean(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors(), -disease) %>%
    step_corr(all_numeric_predictors(), threshold = 0.9) %>%
    step_dummy(all_nominal_predictors(), -disease)

  #-----------------------------------------------------------------------------
  # Model specifications
  #-----------------------------------------------------------------------------

  # K-nearest neighbors
  nearest_neighbor_kknn_spec <-
    nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) %>%
    set_engine("kknn") %>%
    set_mode("classification")

  # Logistic regression
  logistic_reg_glmnet_spec <-
    logistic_reg(penalty = tune(), mixture = tune()) %>%
    set_mode("classification") %>%
    set_engine("glmnet")

  # Random forest
  rand_forest_ranger_spec <-
    rand_forest(mtry = tune(), min_n = tune()) %>%
    set_engine("ranger") %>%
    set_mode("classification")

  # Naive Bayes
  naive_Bayes_naivebayes_spec <-
    naive_Bayes(smoothness = tune(), Laplace = tune()) %>%
    set_engine("naivebayes") %>%
    set_mode("classification")

  # SVM linear
  svm_linear_kernlab_spec <-
    svm_linear(cost = tune(), margin = tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")

  # SVM polynomial
  svm_poly_kernlab_spec <-
    svm_poly(cost = tune(), degree = tune(), scale_factor = tune(), margin = tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")

  # SVM radial
  svm_rbf_kernlab_spec <-
    svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")

  # XGBoost
  boost_tree_xgboost_spec <-
    boost_tree(
      trees = tune(), tree_depth = tune(), min_n = tune(),
      loss_reduction = tune(), stop_iter = tune(),
      sample_size = tune(), mtry = tune(), learn_rate = tune()
    ) %>%
    set_engine("xgboost") %>%
    set_mode("classification")

  # Neural network
  mlp_nnet_spec <-
    mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>%
    set_engine("nnet", MaxNWts = 20000) %>%
    set_mode("classification")

  #-----------------------------------------------------------------------------
  # Workflow sets
  #-----------------------------------------------------------------------------

  # Define workflow sets
  normalized <-
    workflow_set(
      preproc = list(normalized = normalized_rec),
      models = list(
        KNN = nearest_neighbor_kknn_spec,
        SVM_radial = svm_rbf_kernlab_spec,
        SVM_poly = svm_poly_kernlab_spec,
        SVM_linear = svm_linear_kernlab_spec,
        neural_network = mlp_nnet_spec,
        logistic_reg_norm = logistic_reg_glmnet_spec
      )
    )

  simple <-
    workflow_set(
      preproc = list(simple = simple_rec),
      models = list(
        RF = rand_forest_ranger_spec,
        XGB = boost_tree_xgboost_spec,
        logistic_reg_simple = logistic_reg_glmnet_spec
      )
    )

  with_features <-
    workflow_set(
      preproc = list(full_quad = poly_rec),
      models = list(
        logistic_reg = logistic_reg_glmnet_spec,
        KNN = nearest_neighbor_kknn_spec
      )
    )

  # Combine workflow sets
  all_workflows <-
    bind_rows(simple, normalized, with_features) %>%
    mutate(wflow_id = gsub("(simple_)|(normalized_)", "", wflow_id))

  #-----------------------------------------------------------------------------
  # Parameter tuning grids
  #-----------------------------------------------------------------------------

  set.seed(123)
  grid_size <- 500

  # Random Forest grid
  grid_RF <- rand_forest_ranger_spec %>%
    extract_parameter_set_dials() %>%
    update(mtry = mtry(range = c(1, n_feature_select + 2))) %>%
    grid_latin_hypercube(size = grid_size)

  # XGBoost grid
  grid_XGB <- boost_tree_xgboost_spec %>%
    extract_parameter_set_dials() %>%
    update(
      mtry = mtry(range = c(1, n_feature_select)),
      trees = trees(range = c(20, 500)),
      min_n = min_n(range = c(2, 40)),
      stop_iter = stop_iter(range = c(3, 20)),
      learn_rate = learn_rate(range = c(-1.5, -0.5), trans = log10_trans()),
      loss_reduction = loss_reduction(range = c(-3, 1.5), trans = log10_trans()),
      sample_size = sample_prop(range = c(0.5, 1)),
      stop_iter = stop_iter(range = c(3, 30))
    ) %>%
    grid_latin_hypercube(size = grid_size)

  # KNN grid
  grid_KNN <- nearest_neighbor_kknn_spec %>%
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size = grid_size)

  # SVM grids
  grid_SVM_radial <- svm_rbf_kernlab_spec %>%
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size = grid_size)

  grid_SVM_poly <- svm_poly_kernlab_spec %>%
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size = grid_size)

  grid_SVM_linear <- svm_linear_kernlab_spec %>%
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size = grid_size)

  # Neural network grid
  grid_neural_network <- mlp_nnet_spec %>%
    extract_parameter_set_dials() %>%
    update(epochs = epochs() %>% range_set(c(10, 150))) %>%
    grid_latin_hypercube(size = grid_size)

  # Logistic regression grid
  grid_full_quad_logistic_reg <- logistic_reg_glmnet_spec %>%
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size = grid_size)

  # Add grids to workflows
  all_workflows <- all_workflows %>%
    option_add(grid = grid_RF, id = "RF") %>%
    option_add(grid = grid_XGB, id = "XGB") %>%
    option_add(grid = grid_KNN, id = "KNN") %>%
    option_add(grid = grid_SVM_radial, id = "SVM_radial") %>%
    option_add(grid = grid_SVM_poly, id = "SVM_poly") %>%
    option_add(grid = grid_SVM_linear, id = "SVM_linear") %>%
    option_add(grid = grid_neural_network, id = "neural_network") %>%
    option_add(grid = grid_full_quad_logistic_reg, id = "full_quad_logistic_reg") %>%
    option_add(grid = grid_KNN, id = "full_quad_KNN") %>%
    option_add(grid = grid_full_quad_logistic_reg, id = "logistic_reg_simple") %>%
    option_add(grid = grid_full_quad_logistic_reg, id = "logistic_reg_norm")

  #-----------------------------------------------------------------------------
  # Tune models with race ANOVA approach
  #-----------------------------------------------------------------------------

  set.seed(123)
  race_ctrl <-
    control_race(
      verbose = FALSE,
      verbose_elim = FALSE,
      allow_par = TRUE,
      save_pred = TRUE,
      burn_in = 3,
      num_ties = 10,
      alpha = 0.05,
      randomize = TRUE,
      pkgs = NULL,
      event_level = "second",
      parallel_over = "everything",
      save_workflow = TRUE
    )

  # Time the model tuning process
  time1 <- Sys.time()
  race_results <- all_workflows %>%
    workflow_map(
      "tune_race_anova",
      resamples = folds,
      metrics = metric_set(roc_auc),
      control = race_ctrl,
      seed = 20221111,
      verbose = TRUE
    )
  time2 <- Sys.time()
  time.diff.race <- time2 - time1

  print(glue("|||----------------------- It took {round(time.diff.race,4)} mins to tune all grids with the race approch for {stringr::str_to_upper(all_combis$diseases[i])} and selection of miRNAs {stringr::str_to_upper(all_combis$analysis[i])} -----------------------|||"))

  # Save results
  if (SAVE.files == TRUE) {
    filename_tune_race_results <- glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_20240122_tune_race_results_repeats_{no_repeats}_folds_{no_folds}_{text_disease}_analysis_{str_to_upper(all_combis$analysis[i])}_miRetrieve_{miRetrieveBiomarker}_randomMIR_{random_selection}_more_logit_matchIt.rds")
    saveRDS(object = race_results, file = filename_tune_race_results)
  }

  # Count models evaluated
  num_race_models <- sum(collect_metrics(race_results)$n)

  # Rank and plot results
  race_results %>%
    rank_results() %>%
    filter(.metric == "roc_auc") %>%
    select(model, .config, roc_auc = mean, rank) -> rankings_race

  # Get top models
  top2models_race <- rankings_race$model[1:2]
  for (j in 1:length(all_workflows$info)) {
    if (all_workflows$info[[j]]$model == top2models_race[1]) {
      topmodel1_race <- all_workflows$wflow_id[i]
    }

    if (all_workflows$info[[j]]$model == top2models_race[2]) {
      topmodel2_race <- all_workflows$wflow_id[j]
    }
  }

  # Create ranking plot
  autoplot(
    race_results,
    rank_metric = "roc_auc",
    metric = "roc_auc",
    select_best = TRUE
  ) +
    geom_text(aes(y = mean - 0.05, label = wflow_id), angle = 45, hjust = 1, size = 3) +
    lims(y = c(0.6, 1.0)) +
    ylab("ROC-AUC") +
    labs(
      title = paste0("Ranking models: ", text_disease, ", miRNA_set: ", str_to_upper(all_combis$analysis[i])),
      subtitle = paste0("Racing approach | n_grids evaluated: ", num_race_models)
    ) +
    scale_x_continuous(breaks = seq(1, nrow(race_results), by = 1)) +
    theme_minimal(base_size = 16, base_family = "Arial") +
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme() +
    theme(legend.position = "none") -> plot_tune_race_ranking

  # Save plot
  filename_plot_tune_race_ranking <- glue("{data_path_bestageing2022}/output/tuning_results/{all_combis$diseases[i]}/003c_20240122_tune_race_ranking_repeats_{no_repeats}_folds_{no_folds}_{text_disease}_analysis_{str_to_upper(all_combis$analysis[i])}_miRetrieve_{miRetrieveBiomarker}_randomMIR_{random_selection}_more_logit_matchIt.svg")
  ggsave(
    filename = filename_plot_tune_race_ranking, plot = plot_tune_race_ranking,
    width = 14, height = 10, units = "in"
  )

  # Print completion message
  print(paste0(
    "|||-----------------------Run finished for disease and matchIt: ", stringr::str_to_upper(all_combis$diseases[i]),
    ", and selection of miRNAs: ", stringr::str_to_upper(all_combis$analysis[i]),
    " -----------------------|||"
  ))
} # END LOOP
