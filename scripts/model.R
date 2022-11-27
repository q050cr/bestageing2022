
## INFO -----------------------
# this script is sourced from `scripts/render_param_reports.R`

# dependencies ---------------------------------------------------------------
library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(kableExtra)
library(skimr)
library(tableone)
library(pROC)
library(tidymodels)
options(tidymodels.dark = TRUE)
library(discrim)  # naive bayes
library(finetune)
library(vetiver)
library(workflowsets)
library(baguette)
library(rules)
tidymodels_prefer()
conflicted::conflict_prefer("expand", "tidyr")

cores <- parallel::detectCores()
if (!grepl("mingw32", R.Version()$platform)) {
  library(doMC)
  registerDoMC(cores = cores/2)
} else {
  library(doParallel)
  cl <- makePSOCKcluster(cores/2)
  registerDoParallel(cl)
}


tune_grid_eval <- FALSE

###
# load data ---------------------------------------------------------------

# MIRNA DAT
model_data1 <- clean_names(readRDS(file = '../BestAgeing/data_new/model_data1.RDS'))  # has also multiclass col + diagnoses

load(file = "../BestAgeing/data/mirnas.rda")  # "UKL-HD" n=765
load(file = "../BestAgeing/data/data.rda")  # "UKL-HD" n=731
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

# load-mirnas-from_research
# create vector of described mirnas
load("../BestAgeing/data_research/fromR/researchMiRNAAccession.rda")
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
load(file = "../BestAgeing/data/diagnoses_df.rda")

## SURVIVAL DAT
survival_dat <- clean_names(readRDS(file = 'data/202211908_XMELD_abfrage_best_ageing.rds'))# %>% 
# original path "../../XMeldPortal_neu/meldeportal-tools-meldeportalclient-9.3/Rout/202211908_XMELD_abfrage_best_ageing.rds"

## metadata from DB
# https://www.bestageing.org/Pages/Login.aspx?ReturnUrl=%2f&AspxAutoDetectCookieSupport=1
load(file = "../BestAgeing/data/clean_all_meta.rda")  # created in "scripts/_prepare_metadata.R"
clean_all_meta <- clean_all_meta %>% 
  mutate(age = ifelse(age < 18, NA, age))  # wrong age remove
# cath data? "hkdb"

## load all original metadat xlsx files again to make sure that also overlapped 
#patients (e.g. dcm+cad) are in each group
control_ids <- read_excel("../BestAgeing/data/pheno_controls.xlsx") %>% 
  dplyr::pull(BestAgeingCode)


###
# start loop for specified combinations -------------------------------------
###
for (i in 1:nrow(all_combis)) {
  
  filename <- paste0("../BestAgeing/data/pheno_", all_combis$diseases[i], ".xlsx")
  disease_vector <- paste0(all_combis$diseases[i], "_ids")
  assign(x = disease_vector, 
         value = read_excel(filename) %>% 
           dplyr::pull(BestAgeingCode))
  
  disease_ident_df <- tibble(pat_id = c(eval(as.symbol(disease_vector)),control_ids), 
                             disease = c(rep(all_combis$diseases[i], length(eval(as.symbol(disease_vector)))),
                                         rep("control", length(control_ids)))
  )
  
  # create parameter specific {disease}_ids
  filename <- paste0("../BestAgeing/data/pheno_", all_combis$diseases[i], ".xlsx")
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
  
  selected_mirna_dat <- paste0("data01_", length(researchMiRNAAccession$miRNAName_v21), "mirnas")  
  assign(x = selected_mirna_dat, 
         value = data01 %>% 
           select(pat_id, disease, age, sex, researchMiRNAAccession$miRNAName_v21))
  
  if (all_combis$analysis[i] == "selected") {
    no.mirnas <- length(researchMiRNAAccession$miRNAName_v21)
    #text_intro <- paste0(no.mirnas, " miRNAs, that have been investigated in cardiac pathophysiology")
  } else{
    no.mirnas <- ncol(all_mirnas)-1  # minus pat_id
    #text_intro <- paste0(no.mirnas, " miRNAs, that have been investigated in cardiac pathophysiology")
  }
  
  text_disease <- stringr::str_to_upper(all_combis$diseases[i])
  
  ###
  # MODEL -------------------------------------------------------------------
  ###

  ### Initial Split
  set.seed(123)
  
  if (all_combis$analysis[i] == "selected") {
    # only 114 mirnas in analysis
    modeldat <- eval(as.symbol(selected_mirna_dat)) 
  } else{
    # full analysis with all mirnas
    modeldat <- data01
  }
  
  # make control first factor for all analyses ;)
  modeldat <- modeldat %>% 
    mutate(disease = factor(disease, levels = c("control", all_combis$diseases[i])))
  
  dat_split <- rsample::initial_split(modeldat, strata = disease)
  dat_train <- training(dat_split)
  dat_test <- testing(dat_split)
  
  folds <- 
    vfold_cv(dat_train, strata = disease, v = 5)
  
  ###
  # recipe -------------------------------------------------
  # A
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
    step_dummy(all_nominal_predictors(),-disease)
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
    step_dummy(all_nominal_predictors(),-disease) #%>% 
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
    step_dummy(all_nominal_predictors(),-disease)
  
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
    set_engine('nnet') %>%
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
  if (all_combis$analysis[i] == "selected") {
    # only 114 mirnas in analysis
    grid_RF <- rand_forest_ranger_spec %>%   # 2 hyperparams
      extract_parameter_set_dials() %>% 
      # data dependent
      update(mtry = mtry(range = c(1, ncol(dat_train))) ) %>% 
      grid_latin_hypercube(size=50) 
  } else{
    # full analysis with all mirnas: limit `mtry()`
    grid_RF <- rand_forest_ranger_spec %>% 
      extract_parameter_set_dials() %>% 
      # data dependent
      update(mtry = mtry(range = c(1, ncol(dat_train)/2)) ) %>%  # not using all features
      grid_latin_hypercube(size=50) 
  }
  
  if (all_combis$analysis[i] == "selected") {
    # only 114 mirnas in analysis
    grid_XGB <- boost_tree_xgboost_spec %>%   # 2 hyperparams
      extract_parameter_set_dials() %>% 
      # data dependent
      update(mtry = mtry(range = c(1, ncol(dat_train))) ) %>% 
      grid_latin_hypercube(size=50) 
  } else{
    # full analysis with all mirnas: limit `mtry()`
    grid_XGB <- boost_tree_xgboost_spec %>% 
      extract_parameter_set_dials() %>% 
      # data dependent
      update(mtry = mtry(range = c(1, ncol(dat_train)/2)) ) %>%  # not using all features
      grid_latin_hypercube(size=50) 
  }
  
  grid_KNN <- nearest_neighbor_kknn_spec %>%  # 3 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=50)
  
  grid_SVM_radial <- svm_rbf_kernlab_spec %>%   # 3 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=50)
  
  grid_SVM_poly <- svm_poly_kernlab_spec %>%   # 4 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=75)
  
  grid_SVM_linear <- svm_linear_kernlab_spec %>%   # 2 hyperparams
    extract_parameter_set_dials() %>%
    grid_latin_hypercube(size=50)
  
  grid_neural_network <- mlp_nnet_spec %>%   # 3 hyperparams
    extract_parameter_set_dials() %>% 
    grid_latin_hypercube(size=75)
  
  grid_full_quad_logistic_reg <- logistic_reg_glmnet_spec %>%  # 2 hyperparams
    extract_parameter_set_dials() %>% 
    grid_latin_hypercube(size=50)
  
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
  
  ###
  # tune-grid-workflows -----------------------------------------
  ###
  
  if (tune_grid_eval == TRUE) {
    set.seed(123)
    # the workflow_map() function will apply the same function to all of the workflows in the set
    ## the default is fn="tune_grid"
    grid_ctrl <- 
      control_grid(
        save_pred = TRUE,
        pkgs = NULL,
        event_level = "first",  # default
        parallel_over = "everything",
        save_workflow = TRUE
      )
    
    time1 <- Sys.time()
    grid_results <- 
      all_workflows %>% 
      workflow_map(
        # options to `tune_grid()`
        resamples = folds,
        #grid = 25,  # grid specifed above
        metrics = metric_set(roc_auc),
        control = grid_ctrl,
        # options to `workflow_map()`
        seed = 20221111,
        verbose = TRUE
      )
    time2 <- Sys.time()
    time.diff.grid <- time2-time1
    
    ## SAVE
    filename_tune_grid_results <- paste0("output/tuning_results/", Sys.Date(), "_tune_grid_results_",
                                         text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                         ".rds")
    saveRDS(object = grid_results, file = filename_tune_grid_results)
    
    num_grid_models <- nrow(collect_metrics(grid_results, summarize = FALSE))
  
  
    # grid_results
    grid_results %>% 
      rank_results() %>% 
      filter(.metric == "roc_auc") %>% 
      select(model, .config, roc_auc=mean, rank) -> rankings
    
    top2models_grid <- rankings$model[1:2]
    for (i in 1:length(all_workflows$info)) {
      if (all_workflows$info[[i]]$model == top2models_grid[1]) {
        topmodel1_grid <- all_workflows$wflow_id[i]
      }
      
      if (all_workflows$info[[i]]$model == top2models_grid[2]) {
        topmodel2_grid <- all_workflows$wflow_id[i]
      }
    }
    
    # workflow-sets-plot-rank}
    autoplot(
      grid_results,
      rank_metric = "roc_auc",  # <- how to order models
      metric = "roc_auc",       # <- which metric to visualize
      select_best = TRUE     # <- one point per workflow
    ) +
      geom_text(aes(y = mean - 0.05, label = wflow_id), angle = 45, hjust = 1, size =3) +
      lims(y = c(0.6, 1.0)) +
      labs(title=paste0("Ranking models: ", text_disease, ", miRNA_set: ", 
                        str_to_upper(all_combis$analysis[i])),
           subtitle="Grid approach")+
      ylab("ROC-AUC")+
      ggthemes::theme_few()+
      theme(legend.position = "none") -> plot_tune_grid_ranking
    
    filename_plot_tune_grid_ranking <- paste0("output/plots/", Sys.Date(), "_tune_grid_ranking_", 
                                              text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                              ".png")
    ggsave(filename = filename_plot_tune_grid_ranking, plot = plot_tune_grid_ranking, 
           width = 14, height = 10, 
           units = "in"  # default
    )
    
  
    # inspect hyperparameter results for specific model
    autoplot(grid_results, id = topmodel1_grid, metric = "roc_auc")+
      ylab("ROC-AUC")+
      labs(title=paste0("Hyperparameter performance: ", text_disease, ", miRNA_set: ", 
           str_to_upper(all_combis$analysis[i]) ),
           subtitle = stringr::str_to_upper(topmodel1_grid))+
      ggthemes::theme_few() -> plot_tune_grid_hyperpars_topmodel1
    
    filename_plot_tune_grid_hyperpars_topmodel1 <- paste0("output/plots/", Sys.Date(), "_tune_grid_hyperpars_topmodel1_", 
                                              text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                              ".png")
    ggsave(filename = filename_plot_tune_grid_hyperpars_topmodel1, plot = plot_tune_grid_hyperpars_topmodel1, 
           width = 14, height = 10, 
           units = "in"  # default
    )
    
    autoplot(grid_results, id = topmodel2_grid, metric = "roc_auc")+
      ylab("ROC-AUC")+
      labs(title="Hyperparameter Performance - Grid Approach",
           subtitle = stringr::str_to_upper(topmodel2_grid))+
      ggthemes::theme_few() -> plot_tune_grid_hyperpars_topmodel2
    
    filename_plot_tune_grid_hyperpars_topmodel2 <- paste0("output/plots/", Sys.Date(), "_tune_grid_hyperpars_topmodel1_", 
                                                          text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                                          ".png")
    ggsave(filename = filename_plot_tune_grid_hyperpars_topmodel2, plot = plot_tune_grid_hyperpars_topmodel2, 
           width = 14, height = 10, 
           units = "in"  # default
    )
    
  }
  
  ###
  # tune-race-anova -----------------------------------------------------------
  ###
  
  set.seed(123)
  race_ctrl <-
    control_race(
      verbose = TRUE,
      verbose_elim = TRUE,
      allow_par = TRUE,
      save_pred = TRUE,
      burn_in = 3,
      num_ties = 10,
      alpha = 0.05,
      randomize = TRUE,
      pkgs = NULL,
      event_level = "first",
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
  
  ## SAVE
  filename_tune_race_results <- paste0("output/tuning_results/", Sys.Date(), "_tune_race_results_",
                                       text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                       ".rds")
  saveRDS(object = race_results, file = filename_tune_race_results)
  
  num_race_models <- sum(collect_metrics(race_results)$n)
  
  # PLOT results
  race_results %>% 
    rank_results() %>% 
    filter(.metric == "roc_auc") %>% 
    select(model, .config, roc_auc=mean, rank) -> rankings_race
  
  top2models_race <- rankings_race$model[1:2]
  for (i in 1:length(all_workflows$info)) {
    if (all_workflows$info[[i]]$model == top2models_race[1]) {
      topmodel1_race <- all_workflows$wflow_id[i]
    }
    
    if (all_workflows$info[[i]]$model == top2models_race[2]) {
      topmodel2_race <- all_workflows$wflow_id[i]
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
    labs(title="Ranking models",
         subtitle="Racing approach")+
    ggthemes::theme_few()+
    theme(legend.position = "none") -> plot_tune_race_ranking
  
  filename_plot_tune_race_ranking <- paste0("output/plots/", Sys.Date(), "_tune_race_ranking_", 
                                            text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                            ".png")
  ggsave(filename = filename_plot_tune_race_ranking, plot = plot_tune_race_ranking, 
         width = 14, height = 10, 
         units = "in"  # default
  )
  
  
  # inspect hyperparameter results for specific model
  autoplot(race_results, id = topmodel1_race, metric = "roc_auc")+
    ylab("ROC-AUC")+
    labs(title="Hyperparameter Performance - Racing Approach",
         subtitle = stringr::str_to_upper(topmodel1_race))+
    ggthemes::theme_few() -> plot_tune_race_hyperpars_topmodel1
  
  filename_plot_tune_race_hyperpars_topmodel1 <- paste0("output/plots/", Sys.Date(), "_tune_race_hyperpars_topmodel1_", 
                                                        text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                                        ".png")
  ggsave(filename = filename_plot_tune_race_hyperpars_topmodel1, plot = plot_tune_race_hyperpars_topmodel1, 
         width = 14, height = 10, 
         units = "in"  # default
  )
  
  autoplot(race_results, id = topmodel2_race, metric = "roc_auc")+
    ylab("ROC-AUC")+
    labs(title="Hyperparameter Performance - Racing Approach",
         subtitle = stringr::str_to_upper(topmodel2_race))+
    ggthemes::theme_few() -> plot_tune_race_hyperpars_topmodel2
  
  filename_plot_tune_race_hyperpars_topmodel2 <- paste0("output/plots/", Sys.Date(), "_tune_race_hyperpars_topmodel1_", 
                                                        text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                                        ".png")
  ggsave(filename = filename_plot_tune_race_hyperpars_topmodel2, plot = plot_tune_race_hyperpars_topmodel2, 
         width = 14, height = 10, 
         units = "in"  # default
  )
  
  if (tune_grid_eval == TRUE) {
    matched_results <- 
      rank_results(race_results, select_best = TRUE) %>% 
      select(wflow_id, .metric, race = mean, config_race = .config) %>% 
      inner_join(
        rank_results(grid_results, select_best = TRUE) %>% 
          select(wflow_id, .metric, complete = mean, 
                 config_complete = .config, model),
        by = c("wflow_id", ".metric"),
      ) %>%  
      filter(.metric == "roc_auc")
    
    matched_results %>% 
      ggplot(aes(x = complete, y = race)) + 
      geom_abline(lty = 3) + 
      geom_point() + 
      geom_text_repel(aes(label = model)) +
      coord_obs_pred() + 
      labs(x = "Complete Grid ROC-AUC", y = "Racing ROC-AUC") +
      ggthemes::theme_few() -> plot_matched_grid_race
    
    filename_plot_matched_grid_race <- paste0("output/plots/", Sys.Date(), "_matched_grid_race_", 
                                                          text_disease, "_analysis_", str_to_upper(all_combis$analysis[i]),
                                                          ".png")
    ggsave(filename = filename_plot_matched_grid_race, plot = plot_matched_grid_race, 
           width = 14, height = 10, 
           units = "in"  # default
    )
  }

  print(paste0("|||-----------------------Run finished for disease: ", stringr::str_to_upper(all_combis$diseases[i]), 
               ", and selection of miRNAs: ", stringr::str_to_upper(all_combis$analysis[i]),
               " -----------------------|||") )
} # END LOOP
