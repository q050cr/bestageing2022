

# shortened script 001c_model_de_analysis_DET_MATRIX to recalculate CIs

runTests <- TRUE

# univariate tests with auc CI ------------------------------------------------

for (i in 1:nrow(all_combis)) { 
  # reassign disease since only full analysis here
  disease <- all_combis[all_combis[["analysis"]] == "full", ]$diseases[i]
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_{all_combis$diseases[i]}_data01.rds")
  if(!file.exists(path2dataprocessed)) {
    next
  }
  
  data01 <- readRDS(file = path2dataprocessed)
  all_filtered_mirnas <- data01 %>% select(-c(disease, age, sex))  # bring to same str as "all_filtered_mirnas" to work with existing code
  
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
  aucs_lowerci <- rep(NA,ncol(all_filtered_mirnas)-1)
  aucs_upperci <- rep(NA,ncol(all_filtered_mirnas)-1)
  
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
      # Calc AUC CI
      roc_obj <- roc(controls = cont, cases = case)
      auc_conf <- ci(roc_obj)
      
      aucs_lowerci[miRNA] <- auc_conf[1]
      aucs_upperci[miRNA] <- auc_conf[3]
      
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
      
      # status bar
      setTxtProgressBar(pb, miRNA)
    }
    
    ## GATHER Results -----------------------------------------------------------
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
                         aucs_univariate_lowerci = aucs_lowerci,
                         aucs_univariate_upperci = aucs_upperci,
                         aucs_glm = auc_glm,
                         aucs_glm_sva = auc_glm_sva,
                         pval.glm_pca = pval.glm_pca,
                         auc_glm_pca = auc_glm_pca,
                         pval.glm_rob_pca =pval.glm_rob_pca,
                         auc_glm_rob_pca=auc_glm_rob_pca)
    
    
    results_logmedians <- tibble(miR =name.mir, 
                                 auc=aucs, 
                                 aucs_univariate_lowerci = aucs_lowerci,
                                 aucs_univariate_upperci = aucs_upperci,
                                 aucs_glm = auc_glm, aucs_glm_sva = auc_glm_sva,
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
      mutate(sign_indicator = ifelse(padj.glm_pca < 0.05, "p.adj≤0.05", "n.s."),   # before padj chosen
             sign_indicator_sva = ifelse(padj.glm_sva < 0.05, "p.adj≤0.05", "n.s."))
    
    if (SAVE.files ==TRUE) {
      filename.de.tibble <- glue("{data_path_bestageing2022}/output/de_results/{disease}/20240125_001c_de_results_batch_corrected.rds")
      saveRDS(object = de.results, file = filename.de.tibble)
      filename.logmedians <- glue("{data_path_bestageing2022}/output/de_results/{disease}/20240125_001c_results_logmedians_batch_corrected.rds")
      saveRDS(object = results_logmedians, file = filename.logmedians)
    }
  }
}



# MATCHED ANALYSIS --------------------------------------------------------

runTests <- TRUE

for (i in 1:nrow(all_combis)) { 
  # reassign disease since only full analysis here
  disease <- all_combis[all_combis[["analysis"]] == "full", ]$diseases[i]
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_{all_combis$diseases[i]}_data01.rds")
  if(!file.exists(path2dataprocessed)) {
    next
  }
  
  data01 <- readRDS(file = path2dataprocessed)
  
  set.seed(123)
  modeldat <- data01
  
  ## MATCHING -------------------------
  modeldat <- na.omit(modeldat) # matching without missings
  modeldat$disease <- as.factor(modeldat$disease)
  
  modeldat$disease <- factor(modeldat$disease, levels=c("control", all_combis$diseases[i]))
  
  m.out <- matchit(disease ~ age + sex, data = modeldat, method = "nearest")
  summary(m.out)
  
  # plot(m.out, type = "jitter")
  data01 <- match.data(m.out) %>% select(-c("distance",        "weights",         "subclass"))
  
  
  ## continue ----------------------
  all_filtered_mirnas <- data01 %>% select(-c(disease, age, sex))  # bring to same str as "all_filtered_mirnas" to work with existing code
  
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
  aucs_lowerci <- rep(NA,ncol(all_filtered_mirnas)-1)
  aucs_upperci <- rep(NA,ncol(all_filtered_mirnas)-1)
  
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
      # Calc AUC CI
      roc_obj <- roc(controls = cont, cases = case)
      auc_conf <- ci(roc_obj)
      
      aucs_lowerci[miRNA] <- auc_conf[1]
      aucs_upperci[miRNA] <- auc_conf[3]
      
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
      
      # status bar
      setTxtProgressBar(pb, miRNA)
    }
    
    ## GATHER Results -----------------------------------------------------------
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
                         aucs_univariate_lowerci = aucs_lowerci,
                         aucs_univariate_upperci = aucs_upperci,
                         aucs_glm = auc_glm,
                         aucs_glm_sva = auc_glm_sva,
                         pval.glm_pca = pval.glm_pca,
                         auc_glm_pca = auc_glm_pca,
                         pval.glm_rob_pca =pval.glm_rob_pca,
                         auc_glm_rob_pca=auc_glm_rob_pca)
    
    
    results_logmedians <- tibble(miR =name.mir, 
                                 auc=aucs, 
                                 aucs_univariate_lowerci = aucs_lowerci,
                                 aucs_univariate_upperci = aucs_upperci,
                                 aucs_glm = auc_glm, aucs_glm_sva = auc_glm_sva,
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
      mutate(sign_indicator = ifelse(padj.glm_pca < 0.05, "p.adj≤0.05", "n.s."),   # before padj chosen
             sign_indicator_sva = ifelse(padj.glm_sva < 0.05, "p.adj≤0.05", "n.s."))
    
    if (SAVE.files ==TRUE) {
      filename.de.tibble <- glue("{data_path_bestageing2022}/output/de_results/{disease}/20240125_001c_de_results_batch_corrected_matched.rds")
      saveRDS(object = de.results, file = filename.de.tibble)
      filename.logmedians <- glue("{data_path_bestageing2022}/output/de_results/{disease}/20240125_001c_results_logmedians_batch_corrected_matched.rds")
      saveRDS(object = results_logmedians, file = filename.logmedians)
    }
  }
}




